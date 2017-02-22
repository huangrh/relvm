#' Group Score Prediction
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

    if (missing(newdata) || is.null(newdata)) {
        score_tbl <- object$score_tbl
        wts_tbl   <- object$wts_tbl
    }

    par_names <- names(object$par);
    theta <- list(mu = object$par[grepl("mu", par_names)],
                  fl = object$par[grepl("fl", par_names)],
                  err= object$par[grepl("err",par_names)])

    object$pred <- pred(score_tbl,wts_tbl,pms=theta)
    object$tag <- "test"
    object
}

#
pred <- function(score_tbl = data_tbl$score,
                 wts_tbl   = data_tbl$wts,
                 pms       = theta){

    len = nrow(score_tbl);
    score_tbl<- as.matrix(score_tbl)
    wts_tbl   <- as.matrix(wts_tbl)
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
                     score   = score_row,
                     wts     = weight_row,
                     pms     = pms)

        c( fit$par,sqrt(diag(solve(fit$hessian))))
    }))

    rownames(out) <- c("pred","stderr")
    (out <- as.data.frame(t(out)))
}


# Prediction Function
pnll <- function(fv,
                 score = score_row,
                 wts   = weight_row,
                 pms   = theta) {

    # Reconstruction of the parameters
    mu    <- pms$mu         #
    fl    <- pms$fl         # factor loading
    err   <- pms$err        #
    score <- as.matrix(score)
    wts   <- as.matrix(wts)
    # negtive log likelyhood
    means <- mu + fl * fv
    out1  <- wts * (dnorm(score, mean=means, sd = err, log=TRUE))
    out2  <- dnorm(fv, mean=0,sd=1,log=TRUE)

    -sum(out1,out2,na.rm=TRUE)
}
