#' Multiple Estimation Of The Random Effect Latent Variable Model Parameters
#'
#' Estimate multiple groups of the random effect latent variable model
#'
#' @param object A mstbl object.
#' @param groups A vector of measure group names. The default is NULL, in which
#'   case a vector of all groups will be generated accordingly.
#' @param fit A list of fitting parameters. \itemize{ \item qpoints: The numbe
#'   of the quadrature points. \item init: Initial values for mu, fl, and err
#'   term in a list. fl is the factor loading. They will be initialized
#'   generally if it is null. The default is a list with for all mu and one for
#'   others. \item predict: The default is TRUE. }
#'
#' @return An list of S3 object of class "relvm" with estimated parametes.
#' @seealso \code{\link{mstbl}}
#' @importFrom pracma hessian
#'
#' @export
#'
relvm <- function(object,groups=NULL, fit=control(qpoints=30,init=NULL,predict=TRUE)) {

    # -------------------------------------------------------
    # Prepare to call relvm
    alldf  <- merge(x=object$mstbl_std, y=object$wtbl, all=TRUE)

    # Check & update "groups"
    mtbl       <- create_measure_tbl(alldf)
    all_groups <- unique(mtbl$group)
    if (is.null(groups)) {
        groups <- all_groups
    } else if (any(groups %in% all_groups)) {
        groups <- groups[groups %in% all_groups]
    } else stop("The group name do not match.")

    # control
    if (is.null(fit)) {
        qpoints = 30; init = NULL; predict =TRUE
    } else {
        qpoints = fit$qpoints
        init    = fit$init
        predict = fit$predict
    }

    # ------------------------------------------------------------------#
    # Call relvm_single
    allout <- sapply(groups, relvm_single, df=alldf, qpoints=qpoints,
                      init = init, predict = predict, simplify = FALSE)

    # ------------------------------------------------------------------#
    # After Relvm:
    # Merge the predicted group score if there is multiple group.
    preds <- alldf[,1,drop=FALSE]
    for (group in allout) {preds <- merge(x=preds,y=group$pred,all=TRUE)}

    # Merge factor loadings and other parametes.
    pars <- data.frame()
    for (group in allout) {pars = rbind(pars,group$par)}

    #output
    allout$groups <- structure(list(preds=preds,pars=pars),class="relvm")
    (allout)
}

# Set the relvm fitting control parametes.
control <- function(qpoints = 30,init=NULL,predict=TRUE){
    {(control <- list(qpoints=qpoints,init=init,predict=predict))}
}


#' Estimation Of The Random Effect Latent Variable Model Parameters
#'
#' Estimate the random effect latent variable model
#'
#' @param mstbl_std The standized measure score table.
#' @param wts_tbl The measure score weight table.
#' @param qpoints The numbe of the quadrature points.
#' @param init Initial values for mu, fl, and err term in a list. fl is the
#'   factor loading. They will be initialized generally if it is null. The
#'   default is a list with for all mu and one for others.
#' @param predict The default is TRUE.
#'
#' @return An object of S3 class "relvm" with estimated parametes.
#'
#' @export
#'
relvm_single <- function(group, df = alldf,
                         qpoints   = 30,
                         init      = NULL,
                         predict   = TRUE) {
    # -------------------------------------------------------#
    # Prepare to fit
    # start of the cycle
    start_time = Sys.time()
    cat("Fitting:",group,"=>")

    # data table & weight table
    subdat    <- sub1group(group,df)
    mstbl_std <- subdat$mstbl_std
    wts_tbl   <- subdat$wtbl

    # Setup and initialize the parameters
    nc <- ncol(mstbl_std);
    if (is.null(init)) {
        init <- unlist(list(mu  = rep(0.5, nc),
                            fl  = rep(0.5, nc),
                            err = rep(0.5, nc)))
    }

    #--------------------------------------------------------#
    # Fit the function
    fit <- optim(par     = init,     # Model parameter
                 fn      = venll,    # Estimation function
                 gr      = NULL,
                 method  = "L-BFGS-B",
                 control = list(maxit=1000),
                 hessian = FALSE,
                 score   = mstbl_std,
                 wts     = wts_tbl,
                 qpoints = qpoints)

    #--------------------------------------------------------#
    # Output the fitting
    #
    # Format fit$par
    theta_names <- names(fit$par);
    fit$par     <- data.frame(name = names(mstbl_std),
                              mu = fit$par[grepl("mu", theta_names)],
                              fl = fit$par[grepl("fl", theta_names)],
                              err= fit$par[grepl("err",theta_names)],
                              row.names=NULL)
    # Prediction
    if (identical(predict, TRUE)) {
        pred_out       <- relvm:::pred(mstbl_std,wts_tbl,pms=fit$par);
        names(pred_out)<- paste(names(pred_out),group,sep="_")
        fit$pred       <- cbind(subdat$pid,pred_out)
    }

    # Add three fields to the output
    init_names  <- names(init);
    fit$init    <- data.frame(name = names(mstbl_std),
                              mu = init[grepl("mu", init_names)],
                              fl = init[grepl("fl", init_names)],
                              err= init[grepl("err",init_names)], row.names=NULL)
    fit$mstbl_std = cbind(subdat$pid,mstbl_std)
    fit$wtbl      = cbind(subdat$pid,wts_tbl)

    # Output
    cat("\tTime:",Sys.time() - start_time,"\n")
    structure(fit,class="relvm")
}

# Estimation function
enll <- function(par,
                 score = mstbl_std,
                 wts   = wts_tbl,
                 qpoints= qpoints) {

    # Setup the parameters
    nr <- nrow(score);
    cc <- pracma::gaussHermite(qpoints);
    nc <- ncol(score);

    # Reconstruction of the parameters
    mu    <- par[gsub("\\d","", names(par)) %in% "mu"]         #
    fl    <- par[gsub("\\d","", names(par)) %in% "fl"]         # factor loading
    err   <- par[gsub("\\d","", names(par)) %in% "err"]

    # For rows of the mstbl_std and wts_tbl
    out  <- sapply(1:nr, function(idx) {

        # Calculation
        wts_row   <- as.matrix(wts[idx,])
        score_row <- as.matrix(score[idx,])

        # For hessian computation only
        fn <- function(fv) {
            means <- mu + fl * fv
            out1  <- wts_row * (dnorm(score_row, mean=means, sd = err, log=TRUE))
            out2  <- dnorm(fv, mean=0,sd=1,log=TRUE)
            -sum(out1,out2,na.rm=TRUE)
        }
        coef   <- sqrt(2) * numDeriv::hessian(fn, x=0)^(-1/2)
        fv     <- coef * cc$x

        # Weighted log likelyhood
        wts_m     <- matrix(wts_row,  nrow=qpoints, ncol=nc,byrow=TRUE)
        score_m   <- matrix(score_row,nrow=qpoints, ncol=nc,byrow=TRUE)
        means_m   <- matrix(mu,       nrow=qpoints, ncol=nc,byrow=TRUE) + fv %o% fl
        err_m     <- matrix(err,      nrow=qpoints, ncol=nc,byrow=TRUE)

        wll   <- rowSums(wts_m * (dnorm(score_m, mean=means_m, sd = err_m, log=TRUE)),na.rm=TRUE)

        # Joint probability
        joint <- wll +  (dnorm(fv, mean=0,sd=1,log=TRUE))

        # Gaussian quadrature integral approximation
        # gqi <- log(sum(exp(joint + log(cc$w) +(cc$x)^2),na.rm=TRUE))
        gqi <- matrixStats::logSumExp(joint + log(cc$w) +(cc$x)^2,na.rm=TRUE)
        gqi <- log(coef) + gqi
    })

    -sum(out)
}

# Vectorized Estimation function
venll <- function(par,
                  score  = mstbl_std,
                  wts    = wts_tbl,
                  qpoints= qpoints) {

    # Setup the parameters
    nr <- nrow(score);
    cc <- pracma::gaussHermite(qpoints);
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

    #
    fv_m  <- cc$x %o% coefs
    wts   <- as.matrix(wts)
    score <- as.matrix(score)

    # Weighted log likelyhood
    wts_a   <- aperm(array(wts,  dim=c(nr,nc,qpoints)),c(2,3,1))
    score_a <- aperm(array(score,dim=c(nr,nc,qpoints)),c(2,3,1))
    means_a <- array(mu,         dim=c(nc,qpoints,nr)) + fl %o% fv_m

    wll_m   <- colSums(wts_a * (dnorm(score_a, mean=means_a, sd = err, log=TRUE)),na.rm=TRUE)

    # Joint probability
    joint_m <- wll_m +  (dnorm(fv_m, mean=0,sd=1,log=TRUE))

    # Gaussian quadrature integral approximation
    # gqi <- log(sum(exp(joint + log(cc$w) +(cc$x)^2),na.rm=TRUE))
    gqi <- matrixStats::colLogSumExps(joint_m + log(cc$w) +(cc$x)^2,na.rm=TRUE)
    -sum(log(coefs)+gqi, na.rm=TRUE)
}



