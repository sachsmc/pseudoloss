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


fast_empirical_roc <- function(X, pseus) {

    sumpseus <- rowSums(pseus)
    data.frame(tpf = rev(cumsum(pseus[order(X, decreasing = TRUE), 1]) / sum(pseus[, 1])),
               fpf = rev(cumsum(1 - sumpseus[order(X, decreasing = TRUE)]) / (sum(1 - sumpseus))))

}


#' Compute AUROC using the trapezoidal rule
#'
#' @param fpf Vector of false positive fractions/1 - specificity
#' @param tpf Vector of true positive fractions/sensitivy

calc_auc <- function(fpf, tpf) {
    n <- length(fpf)
    stopifnot(n == length(tpf))
    dx <- fpf[-n] - fpf[-1]
    mid.y <-  (tpf[-n] + tpf[-1]) / 2
    sum(dx * mid.y)

}


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

    if(any(is.nan(numer))) stop("Ramp function too large, adjust sigma or recale X to prevent NaNs")

    numer / (sum(K1) * sum(K2))

}
