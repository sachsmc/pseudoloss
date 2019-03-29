
#' Make pseudo observation evaluation matrix
#'
#' Given a pseudo object, and a time of interest,
#' stacks in columns the pseudo observations of the
#' different causes at that time. Result is suitable for use
#' with the ROC and AUC functions.
#'
#' @param pseudo an object as returned from pseudoci, pseudosurv, or similar
#' @param time index of the time point of interest

make_pseudo_matrix <- function(pseudo, time = 1) {

    do.call(cbind, lapply(pseudo$pseudo, function(x) x[, time]))

}

#' Make pseudo observation predictor data frame
#'
#' Given a pseudo object, stacks in rows the pseudo observations for
#' a fixed cause over the different times at which the pseudo
#' observations were calculated, along with the time at which they
#' were calculated.
#'
#' @param pseudo an object as returned from pseudoci, pseudosurv, or similar
#' @param cause index of the cause of interest

make_pseudo_predictor <- function(pseudo, cause = 1) {

    do.call(rbind, lapply(1:length(pseudo$time), function(x) {

        data.frame(time = pseudo$time[x], pseudo.Y = pseudo$pseudo[[cause]][, x])

        }))

}
