
#' Pseudo observation AUC objective function
#'
#' Designed for use with xgboost, given a vector of predictions and
#' a training matrix, returns the gradient and hessian so that it can
#' be optimized by xgboost
#'
#' @param preds Vector of predictions
#' @param dtrain An object of class xgb.Dmatrix, with additional attributes otherc and tdex
#' @export

pseudoauc_objective <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    otherc <- attr(dtrain, "otherc")

    pseus <- cbind(labels, otherc)
    if(length(table(preds)) == 1) preds <- preds + runif(length(preds), -.1, .1)

    trans.preds <- preds #1 / (1 + exp(-preds))


    grad <- auc_gradient(trans.preds, pseus)
    hess <- auc_hessian(trans.preds, pseus)

    grad[is.nan(grad)] <- min(grad, na.rm = TRUE) * 1.1
    hess[is.nan(hess)] <- min(hess, na.rm = TRUE) * 1.1

    return(list(grad = -grad, hess = -hess))
}
