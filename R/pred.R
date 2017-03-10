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
#' @param pms A list or a data frame containing mu, fl, and err.
#'
pred <- function(score_tbl,wts_tbl,pms){
    nr  <- nrow(score_tbl)
    out <- vapply(1:nr, function(idx) {

        # fit the function
        fv <- 1;
        fit <- optim(par     = fv,     # factor variable / latent variable
                     fn      = pnll_cpp,   # prediction function
                     gr      = NULL,
                     method  = "BFGS",
                     control = list(maxit=5200),
                     hessian = TRUE,
                     score_row = score_tbl[idx,],
                     wts_row = wts_tbl[idx,],
                     err     = pms[["err"]],
                     mu      = pms[["mu"]],
                     fl      = pms[["fl"]])

        c( fit$par,1/sqrt(fit$hessian))
    },FUN.VALUE=numeric(2))

    rownames(out) <- c("pred","stderr")
    (out <- t(out))
}


# Prediction negLogLik Function
pnll <- function(fv, score_row, wts_row, err,mu,fl) {

    # negtive log likelyhood
    # out1  <- wts_row * (dnorm(score_row, mean=means, sd = err, log=TRUE))
    # out2  <- dnorm(fv, mean=0,sd=1,log=TRUE)
    # -sum(out1,out2,na.rm=TRUE)
    # -sum(wts_row*dnorm2(score_row, mean=means, sd = err),dnorm_cpp(fv, mean=0,sd=1),na.rm=TRUE)
    sum(wts_row*(0.79817986835+2*log(err)+((score_row-mu - fl * fv)/err)^2)/2,
        (0.79817986835+(fv)^2)/2,na.rm=TRUE)
    # -sum(wts_row * dnorm(score_row,
    #                      mean=pms[["mu"]] + pms[["fl"]] * fv,
    #                      sd = pms[["err"]], log=TRUE),
    #      dnorm(fv, mean=0,sd=1,log=TRUE),na.rm=TRUE)
}
