qpoints = 30
cc <- pracma::gaussHermite(qpoints);
venll <- function(par,
                  score = mstbl_std,
                  wts   = wts_tbl,
                  cc= cc) {

    # Setup the parameters
    nr <- nrow(score);

    nc <- ncol(score);

    # Reconstruction of the parameters
    mu    <- par[gsub("\\d","", names(par)) %in% "mu"]         #
    fl    <- par[gsub("\\d","", names(par)) %in% "fl"]         # factor loading
    err   <- par[gsub("\\d","", names(par)) %in% "err"]

    # For rows of the mstbl_std and wts_tbl
    coefs  <- sapply(1:nr, function(idx) {
        # For hessian computation only
        fn <- function(fv) {
            means <- mu + fl * fv
            out1  <- wts_row * (dnorm(score_row, mean=means, sd = err, log=TRUE))
            out2  <- dnorm(fv, mean=0,sd=1,log=TRUE)
            -sum(out1,out2,na.rm=TRUE)
        }
        # Calculation
        wts_row   <- as.matrix(wts[idx,])
        score_row <- as.matrix(score[idx,])
        (coef   <- sqrt(2 / pracma::hessian(fn,0))) # pracma::hessian or numDeriv::hessian
    })
    coefs
}
