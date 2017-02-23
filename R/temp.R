venll7 <- function(par,score,wts,cc,qpoints) {

    # Reconstruction of the parameters
    nr <- nrow(score); nc <- ncol(score)
    mu    <- par[grepl("mu", names(par))]         #
    fl    <- par[grepl("fl", names(par))]         # factor loading
    err   <- par[grepl("err", names(par))]

    # 2nd derivative
    coefs <- sqrt(2/(rowSums((fl^2 * wts/err^2), na.rm = TRUE) + 1))

    # fv matrix
    fv_mtx  <- cc$x %o% coefs

    # 3D array:
    wts_arr   <- aperm(array(wts,  dim=c(nr,nc,qpoints)),c(2,3,1))
    score_arr <- aperm(array(score,dim=c(nr,nc,qpoints)),c(2,3,1))
    means_arr <- array(mu,         dim=c(nc,qpoints,nr)) + fl %o% fv_mtx

    # Weighted log likelyhood
    wll_mtx   <- colSums(wts_arr * denfn2(score_arr, mean=means_arr, sd = err),na.rm=TRUE)

    # Joint probability
    joint_mtx <- wll_mtx +  denfn2(fv_mtx, mean=0,sd=1)

    # Gaussian quadrature integral approximation
    # gqi <- log(sum(exp(joint + log(cc$w) +(cc$x)^2),na.rm=TRUE))
    gqi <- matrixStats::colLogSumExps(joint_mtx + log(cc$w) +(cc$x)^2,na.rm=TRUE)
    -sum(log(coefs)+gqi, na.rm=TRUE)
}
