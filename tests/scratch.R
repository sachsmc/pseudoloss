
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

X <- lung$age / max(lung$age)

system.time(numt <- numeric_grad_auc(X, pseus))
system.time(relt <- auc_gradient(X, pseus))


## multivariate testing

X1 <- lung$age / max(lung$age)
X2 <- pseus + rnorm(nrow(lung))

fn.auc <- function(beta) {

    X <- X1 * beta[1] + X2 * beta[2]
    1 - ramp_auc(X, pseus) + 10 * sum(beta^2)
}


gn.auc <- function(beta) -multivar_auc_gradient(cbind(X1, X2), pseus)(beta) + 20 * beta

system.time(ramp_auc(cbind(X1, X2) %*% optim(par = c(0.05, .05) , fn = fn.auc, gr = gn.auc, method = "BFGS")$par, pseus))
system.time(ramp_auc(cbind(X1, X2) %*% optim(par = c(0.05, .05) , fn = fn.auc, method = "BFGS")$par, pseus))


## xgboost testing
library(xgboost)

noise <- matrix(rnorm(nrow(pseus) * 25), ncol = 25)
dtrain.po <- xgb.DMatrix(as(cbind(X1, X2, lung$sex, noise), "dgCMatrix"), label = c(pseus[, 1]))
attr(dtrain.po, "otherc") <- NULL

preds <- runif(nrow(pseus), .4, .6)
pseudoauc_objective(preds, dtrain.po)

library(xgboost)
param <- list(max_depth=6,eta = .3, nthread = 2, verbose = 1,
              objective = pseudoauc_objective, base_score = .5,
              booster = "gbtree", maximize = FALSE,
              early_stopping_rounds = NULL)

po.test <- xgb.train(param, dtrain.po, nrounds = 50, verbose = 1)

plot(calc_roc(predict(po.test, newdata = dtrain.po), pseus))


## superlearner



