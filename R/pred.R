#' Group Score Prediction Interface
#'
#' Predict the group scores (random effect latent variabl) from the extimated
#' parameters
#'
#' @param pms A named vector of parameters to be optimized.
#' @param score_tbl The data table.
#' @param wts_tbl The corresponding weigting table.
#' @return A data.frame containing pred and stderr. pred is the predicted group score.
#'   score.
#'
#' @export
#'
predict.relvm <- function(object,newdata,level = 0.95){
    object$tag  <- "predict.relvm"

    if (missing(newdata) || is.null(newdata)) {
        score_tbl <- as.data.frame(object$mstbl_std)
        wts_tbl   <- as.data.frame(object$wtbl)
    }
    measure_idx <- rstarating::measure(score_tbl)
    score_tbl <- as.matrix(score_tbl[measure_idx])
    wts_tbl   <- as.matrix(wts_tbl[paste(measure_idx,"wt",sep="_")])

    object$pred <- pred(score_tbl,wts_tbl,pms = object$par)
    if (exists("provider_id",as.data.frame(object$mstbl_std))) {
        object$pred <- cbind(as.data.frame(object$mstbl_std)["provider_id"],object$pred)
    }

    object
}

#' Group Score Prediction
#'
#' @param score_tbl A matrix of standadized measure scores.
#' @param wts_tbl A matrix of weights.
#' @
#'
pred <- function(score_tbl,wts_tbl,pms){

    len = nrow(score_tbl);
    par = stderr = c()

    out <- as.data.frame(sapply(1:len, function(idx) {
        score_row    <- score_tbl[idx,]
        weight_row   <- wts_tbl[idx,]

        # fit the function
        fv <- 1;
        fit <- optim(par     = fv,     # factor variable / latent variable
                     fn      = pnll,   # prediction function
                     gr      = NULL,
                     method  = "BFGS",
                     control = list(maxit=5200),
                     hessian = TRUE,
                     score_row   = score_row,
                     wts_row     = weight_row,
                     pms     = pms)

        c( fit$par,sqrt(diag(solve(fit$hessian))))
    }))

    rownames(out) <- c("pred","stderr")
    (out <- as.data.frame(t(out)))
}


# Prediction negLogLik Function
pnll <- function(fv, score_row, wts_row, pms) {

    # Reconstruction of the parameters
    mu    <- pms$mu         #
    fl    <- pms$fl         # factor loading
    err   <- pms$err        #

    # negtive log likelyhood
    means <- mu + fl * fv
    out1  <- wts_row * (dnorm(score_row, mean=means, sd = err, log=TRUE))
    out2  <- dnorm(fv, mean=0,sd=1,log=TRUE)

    -sum(out1,out2,na.rm=TRUE)
}
