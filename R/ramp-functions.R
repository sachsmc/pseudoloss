

smoothramp <- function(X.vect, c.scal, sigma = .01) {

    1 / (1 + exp(-(X.vect - c.scal) / sigma))

}


stepramp <- function(X.vect, c.scal) {

    as.numeric(X.vect > c.scal)

}


sloperamp <- function(X.vect, c.scal, sigma = .05) {

    dif <- -(X.vect - c.scal) / sigma

    ifelse(dif <= -.5, 1, ifelse(dif > .5, 0, -dif + .5))

}


## gradients

d.sloperamp <- function(X.vect, c.scal, sigma = .05) {

    dif <- -(X.vect - c.scal) / sigma

    ifelse(abs(dif) > .5, 0, 1 / sigma)

}


d2.sloperamp <- function(X.vect, c.scal, sigma = .05) {

    dif <- -(X.vect - c.scal) / sigma

    ifelse(abs(dif) > .5, 0, 1 / sigma)

}

d.smoothramp <- function(X.vect, c.scal, sigma = .01) {

    exp(-(X.vect - c.scal) / sigma) / (sigma * (1 + exp(-(X.vect - c.scal) / sigma))^2)

}


d2.smoothramp <- function(X.vect, c.scal, sigma = .01) {

      ifelse(X.vect > c.scal, -1, 1) * (exp((2 * c.scal + X.vect) / sigma) + exp((c.scal + 2 * X.vect) / sigma)) /
        (sigma^2 * (exp(c.scal / sigma) + exp(X.vect / sigma))^3)

}


