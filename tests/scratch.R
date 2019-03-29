
library(pseudo)
library(survival)
library(survivalROC)

pseus <- make_pseudo_matrix(pseudoci(lung$time, lung$status - 1, tmax = 400))


roc.them <- survivalROC(lung$time, lung$status - 1, marker = lung$age, predict.time = 400,
                        method = "KM")


system.time(roc.fast <- fast_empirical_roc(lung$age, pseus))  # about twice as fast but doesn't handle ties

auc.trap <- auc.ramp <- rep(NA, 1000)
for(i in 1:1000) {
fake.x <- rnorm(nrow(lung))

roc.me <- calc_roc(fake.x, pseus)
auc.trap[i] <- calc_auc(roc.me$fpf, roc.me$tpf)
auc.ramp[i] <- ramp_auc(fake.x, pseus, sloperamp)
}



numeric_grad_auc <- function(X, pseus, delta = .005) {

    out.err <- as.numeric(ramp_auc(X, pseus))

    sapply(1:length(X), function(i) {
        pre.in <- X
        pre.in[i] <- X[i] + delta

        err <- as.numeric(ramp_auc(pre.in, pseus))
        (err - out.err) / delta

    })

}

system.time(numt <- numeric_grad_auc(X, pseus))
system.time(relt <- auc_gradient(X, pseus))
