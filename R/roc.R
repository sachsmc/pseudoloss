#' Compute ROC curve based on pseudo observations
#'
#' Given a predictor X and a matrix of pseudo observations, with one
#' column per cause, and only the pseudo observations for a fixed time
#' computes the true and false positive fractions for each unique value
#' of the predictor.
#'
#' @param X the predicted pseudo observations
#' @param pseus A matrix or data frame, pseudo obs of interest in the first column
#' @param ramp A ramp function, first argument is a vector, and second is a scalar
#' @export
#' @return A data frame with columns for the false and true positive fractions

calc_roc <-
    function(X, pseus, ramp = function(X.vect, c.scal) as.numeric(X.vect > c.scal)) {
        ## compute roc curve using pseudo values

        pseu <- pseus[, 1]

        cuts <- sort(unique(X))

        if(length(cuts) == 1) return(data.frame(tpf = c(0, 0), fpf = c(1, 1)))

        cutin <- cuts[-c(1, length(cuts))]
        res.0 <- rbind(c(1, 1), t(sapply(cutin, function(cth) {
            c(sum((1 - rowSums(pseus)) * ramp(X , cth)) / sum(1 - rowSums(pseus)),
              sum(pseu * ramp(X, cth)) / sum(pseu))
        })), c(0, 0))
        res <- as.data.frame(res.0)
        res <- data.frame(lapply(res, function(x)
            pmax(pmin(x, 1), 0)))
        colnames(res) <- c("fpf", "tpf")
        res

    }

#' A faster version of the ROC calculation
#'
#' Does not handle ties

fast_empirical_roc <- function(X, pseus) {

    sumpseus <- rowSums(pseus)
    data.frame(tpf = rev(cumsum(pseus[order(X, decreasing = TRUE), 1]) / sum(pseus[, 1])),
               fpf = rev(cumsum(1 - sumpseus[order(X, decreasing = TRUE)]) / (sum(1 - sumpseus))))

}


#' Compute AUROC using the trapezoidal rule
#'
#' @param fpf Vector of false positive fractions/1 - specificity
#' @param tpf Vector of true positive fractions/sensitivy
#' @export

calc_auc <- function(fpf, tpf) {
    n <- length(fpf)
    stopifnot(n == length(tpf))
    dx <- fpf[-n] - fpf[-1]
    mid.y <-  (tpf[-n] + tpf[-1]) / 2
    sum(dx * mid.y)

}

#' Ramp auc calculation
#'
#' Calculates the approximate pseudo AUC using the specified ramp function
#' to be used as the indicator. Using the indicator function as the ramp will
#' result in the empirical AUC.
#'
#' @param X numeric vector of predictions
#' @param pseus matrix of pseudo observations
#' @param rampfunc Ramp function
#' @export

ramp_auc <- function(X, pseus, rampfunc = stepramp) {

    K1 <- pseus[, 1]
    K2 <- 1 - rowSums(pseus)

    fixc <- sort(unique(X))

    traps <- unlist(lapply(1:(length(fixc) - 1), function(i) {

        sum(outer(K2 * ((rampfunc(X, fixc[i])) - (rampfunc(X, fixc[i + 1]))),
                  K1 * (rampfunc(X, fixc[i]))))

    }))

    sum(traps) / (sum(K1) * sum(K2))



}

#' Gradient of the ramp pseudo AUC
#'
#' Gradient of the AUC using the specified ramp function and the first derivative of the ramp.
#' Returns the gradient with respect to the prediction vector. This is for use with xgboost.
#'
#'
#' @param X numeric vector of predictions
#' @param pseus matrix of pseudo observations
#' @param rampfunc Ramp function
#' @param d.rampfunc first derivative of ramp function
#' @export

auc_gradient <- function(X, pseus, rampfunc = smoothramp, d.rampfunc = d.smoothramp) {  # currently only works with smoothramp

    K1 <- pseus[, 1]
    K2 <- 1 - rowSums(pseus)

    fixc <- sort(unique(X))

    traps <- lapply(1:(length(fixc) - 1), function(i) {

        colSums(outer(K2 * (smoothramp(X, fixc[i]) - smoothramp(X, fixc[i + 1])),
              K1 * d.smoothramp(X, fixc[i]))) +
        rowSums(outer(K2 * (d.smoothramp(X, fixc[i]) - d.smoothramp(X, fixc[i + 1])),
                      K1 * smoothramp(X, fixc[i])))

    })

    numer <- rowSums(do.call(cbind, traps))

    #if(any(is.nan(numer))) stop("Ramp function too large, adjust sigma or recale X to prevent NaNs")

    numer / (sum(K1) * sum(K2))

}

#' Hessian of the ramp pseudo AUC
#'
#' Hessian of the AUC using the specified ramp function and the first and second derivatives
#' of the ramp. Returns the hessian with respect to the prediction vector. This is mainly
#' for use with xgboost.
#'
#'
#' @param X numeric vector of predictions
#' @param pseus matrix of pseudo observations
#' @param rampfunc Ramp function
#' @param d.rampfunc first derivative of ramp function
#' @param d2.rampfunc second derivative of ramp function (only the diagonal terms are needed)
#' @export

auc_hessian <- function(X, pseus, rampfunc = smoothramp, d.rampfunc = d.smoothramp, d2.rampfunc = d2.smoothramp) {

    K1 <- pseus[, 1]
    K2 <- 1 - rowSums(pseus)

    fixc <- sort(unique(X))

    traps <- lapply(1:(length(fixc) - 1), function(i) {

        colSums(outer(K2 * (smoothramp(X, fixc[i]) - smoothramp(X, fixc[i + 1])),
                      K1 * d2.smoothramp(X, fixc[i]))) +
            rowSums(outer(K2 * (d2.smoothramp(X, fixc[i]) - d2.smoothramp(X, fixc[i + 1])),
                          K1 * smoothramp(X, fixc[i])))

    })

    numer <- rowSums(do.call(cbind, traps))

    #if(any(is.nan(numer))) stop("Ramp function too large, adjust sigma or recale X to prevent NaNs")

    numer / (sum(K1) * sum(K2))


}

#' Gradient of AUC with respect to vector of predictors
#'
#' Given a vector of numeric predictors for each subject, this calculates the gradient
#' of the pseudo-AUC with respect to the coefficients for the linear combination of
#' the predictors. This can be used for AUC regression, wherein the predictors are independent
#' biomarkers, and also with SuperLearner, wherein the predictors are the prediction
#' results of a library of prediction algorithms.
#'
#' @param X.wide matrix of predictions
#' @param pseus matrix of pseudo observations
#' @param rampfunc ramp function
#' @param d.rampfunc first derivative of the ramp function with respect to the coefficients
#'
#' @export

multivar_auc_gradient <- function(X.wide, pseus, rampfunc = smoothramp, d.rampfunc = d.smoothramp) { # returns a function of beta

    K1 <- pseus[, 1]
    K2 <- 1 - rowSums(pseus)

    pp <- ncol(X.wide)

    function(beta) {

        if(is.null(ncol(X.wide))) { X <- X.wide * beta } else { X <- X.wide %*% beta }
        fixc <- sort(unique(X))
        g.res <- rep(NA, ncol(X.wide))

        for(k in 1:ncol(X.wide)) {

            traps <- lapply(1:(length(fixc) - 1), function(i) {

                outer(K1 * X.wide[, k] * d.smoothramp(X, fixc[i]),
                      K2 * (smoothramp(X, fixc[i]) - smoothramp(X, fixc[i + 1]))) +
                    outer(K1 * smoothramp(X, fixc[i]), K2 * X.wide[, k] * (d.smoothramp(X, fixc[i]) - d.smoothramp(X, fixc[i + 1])))

            })

        numer <- sum(unlist(traps))

        if(any(is.nan(numer))) stop("Ramp function too large, adjust sigma or recale X to prevent NaNs")

        g.res[k] <- numer / (sum(K1) * sum(K2))

        }

        g.res

    }

}
