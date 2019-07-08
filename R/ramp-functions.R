
#' Sigmoid ramp function
#'
#'@export
smoothramp <- function(X.vect, c.scal, sigma = .01) {

    1 / (1 + exp(-(X.vect - c.scal) / sigma))

}

#' Step ramp function
#'
#'@export

stepramp <- function(X.vect, c.scal) {

    as.numeric(X.vect > c.scal)

}

#' Slope ramp function
#'
#'@export

sloperamp <- function(X.vect, c.scal, sigma = .05) {

    dif <- -(X.vect - c.scal) / sigma

    ifelse(dif <= -.5, 1, ifelse(dif > .5, 0, -dif + .5))

}


#' First derivative of slope ramp function
#'
#'@export

d.sloperamp <- function(X.vect, c.scal, sigma = .05) {

    dif <- -(X.vect - c.scal) / sigma

    ifelse(abs(dif) > .5, 0, 1 / sigma)

}

#' Second derivative of slope ramp function
#'
#'@export

d2.sloperamp <- function(X.vect, c.scal, sigma = .05) {

    dif <- -(X.vect - c.scal) / sigma

    ifelse(abs(dif) > .5, 0, 1 / sigma)

}
#' First derivative of sidmoid ramp function
#'
#'@export

d.smoothramp <- function(X.vect, c.scal, sigma = .05) {

    exp(-(X.vect - c.scal) / sigma) / (sigma * (1 + exp(-(X.vect - c.scal) / sigma))^2)

}

#' Second derivative of sigmoid ramp function
#'
#'@export

d2.smoothramp <- function(X.vect, c.scal, sigma = .05) {

      ifelse(X.vect > c.scal, -1, 1) * (exp((2 * c.scal + X.vect) / sigma) + exp((c.scal + 2 * X.vect) / sigma)) /
        (sigma^2 * (exp(c.scal / sigma) + exp(X.vect / sigma))^3)

}


