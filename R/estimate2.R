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
relvm2 <- function(object,groups=NULL, fit=control(qpoints=30,init=NULL,predict=TRUE,adaptive=c("noad","ad"))) {

    # -------------------------------------------------------
    # Merge both tables of the measure score and weights.
    alldf  <- merge(x=object$mstbl_std, y=object$wtbl, all=TRUE)

    # Check & update "groups"
    mtbl       <- create_measure_tbl(alldf)
    all_groups <- unique(mtbl$group)
    if (is.null(groups)) {
        groups <- all_groups
    } else if (any(groups %in% all_groups)) {
        groups <- groups[groups %in% all_groups]
    } else stop("The group name do not match.")

    # Fit control
    if (is.null(fit)) {
        qpoints = 30; init = NULL; predict =TRUE;adaptive="noad"
    } else {
        qpoints = fit[["qpoints"]]
        init    = fit[["init"]]
        predict = fit[["predict"]]
        adaptive= fit[["adaptive"]][1]
    }

    # ------------------------------------------------------------------#
    # Call relvm_single
    allout <- sapply(groups, relvm_single2, df=alldf, qpoints=qpoints,
                      init = init, predict = predict, adaptive=adaptive,simplify = FALSE)

    # ------------------------------------------------------------------#
    # After Relvm:
    # Merge the predicted group score if there is multiple group.
    preds <- alldf[,1,drop=FALSE] # take the column "provider_id"
    for (group in allout) {preds <- merge(x=preds,y=group$pred,all=TRUE)}

    # Merge factor loadings and other parametes.
    pars <- data.frame()
    for (group in allout) {pars = rbind(pars,group$par)}

    # convergence
    convergence<- data.frame(convergence=vapply(allout,function(x) {x$convergence},c(0)))
    value      <- data.frame(value=vapply(allout,  function(x) x$value,c(0)))
    message    <- data.frame(message=vapply(allout,function(x) x$message,"0"),stringsAsFactors = FALSE)
    counts     <- t(as.data.frame(vapply(allout,   function(x) x$counts,c(0L,0L))))

    #output
    allout$groups <- structure(list(preds=preds,pars=pars,
                                    counts= as.matrix(counts),
                                    value = as.matrix(value),
                                    message=as.matrix(message),
                                    convergence = as.matrix(convergence)),class="relvm")
    (allout)
}

# Set the relvm fitting control parametes.
control <- function(qpoints = 30,init=NULL,predict=TRUE,adaptive=c("noad","ad")){
    {(control <- list(qpoints=qpoints,init=init,predict=predict,adaptive=adaptive))}
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
relvm_single2 <- function(group, df, qpoints,init,predict,adaptive) {
    # -------------------------------------------------------#
    # Prepare to fit
    # start of the cycle
    start_time = Sys.time()
    cat(sprintf("Fitting: %-15s =>",group))

    # data table & weight table
    subdat    <- sub1group(group,df)
    mstbl_std <- as.matrix(subdat$mstbl_std)
    wts_tbl   <- as.matrix(subdat$wtbl)

    # Setup and initialize the parameters
    nc <- ncol(mstbl_std);
    if (is.null(init)) {
        init <- unlist(list(mu  = rep(0.5, nc),
                            fl  = rep(0.5, nc),
                            err = rep(0.5, nc)))}
    # cc$x & cc$w
    cc <- pracma::gaussHermite(qpoints);

    #--------------------------------------------------------#
    # Fit the function
    fit <- optim(par     = init,      # Model parameter
                 fn      = venll10a,    # Estimation function
                 gr      = NULL,
                 method  = "L-BFGS-B",
                 control = list(maxit=1000),
                 hessian = FALSE,
                 score   = mstbl_std,
                 wts     = wts_tbl,
                 cc      = cc,
                 qpoints = qpoints,
                 adaptive= adaptive)

    #--------------------------------------------------------#
    # Output the fitting
    #
    # Format fit$par
    theta_names <- names(fit$par);
    fit$par     <- data.frame(name = names(subdat$mstbl_std),
                              fl = fit$par[grepl("fl", theta_names)],
                              mu = fit$par[grepl("mu", theta_names)],
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
    fit$init    <- data.frame(name = names(subdat$mstbl_std),
                              fl = init[grepl("fl", init_names)],
                              mu = init[grepl("mu", init_names)],
                              err= init[grepl("err",init_names)], row.names=NULL)

    fit$mstbl_std = cbind(subdat$pid,mstbl_std)
    fit$wtbl      = cbind(subdat$pid,wts_tbl)

    # Output
    cat(" : ", as.character.Date(Sys.time() - start_time),"\n")
    structure(fit,class="relvm")
}

# Simplified normal density function.
dnorm2 <- function(x,mean=0,sd=1) -(log(2 * pi) +2*log(sd)+((x-mean)/sd)^2)/2

# adaptive estimate function
venll10a <- function(par,score,wts,cc,qpoints,adaptive) {
    # Reconstruction of the parameters
    nr <- nrow(score); nc <- ncol(score)
    mu <- par[grepl("mu", names(par))]         #
    fl <- par[grepl("fl", names(par))]         # factor loading
    err<- par[grepl("err", names(par))]

    # Sigma
    sigma <- switch(adaptive,
                    ad = {
                        fitall <- relvm:::pred(score_tbl=score, wts_tbl=wts, pms=list(mu=mu,fl=fl,err=err));
                        u_hat  <- t(array(fitall[,"pred"],dim=c(nr,qpoints)));
                        fitall[,"stderr"];},
                    noad= {
                        u_hat <- array(0,dim=c(qpoints,nr));
                        sqrt(1/(rowSums((fl^2 * wts/err^2), na.rm = TRUE) + 1));})

    coefs <- 1.41421356237 * sigma # sqrt(2) = 1.4142
    # coefs <- rep(mean(coefs),length(coefs))

    # fv matrix
    fv_mtx  <-  cc$x %o% coefs + u_hat

    # 3D array:
    wts_arr   <- aperm(array(wts,  dim=c(nr,nc,qpoints)),c(2,3,1))
    score_arr <- aperm(array(score,dim=c(nr,nc,qpoints)),c(2,3,1))
    means_arr <- array(mu,         dim=c(nc,qpoints,nr)) + fl %o% fv_mtx

    # Weighted log likelyhood
    wll_mtx   <- colSums(wts_arr * dnorm2(score_arr, mean=means_arr, sd = err),na.rm=TRUE)

    # Joint probability
    joint_mtx <- wll_mtx +  dnorm_cpp(fv_mtx, mean=0,sd=1)

    # Gaussian quadrature integral approximation
    gqi <- matrixStats::colLogSumExps(joint_mtx + log(cc$w) +(cc$x)^2,na.rm=TRUE)
    -sum(log(coefs)+gqi, na.rm=TRUE)
}


# Fast Vectorized Estimation function
venll7 <- function(par,score,wts,cc,qpoints) {

    # Reconstruction of the parameters
    nr    <- nrow(score); nc <- ncol(score)
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
    wll_mtx   <- colSums(wts_arr * dnorm2(score_arr, mean=means_arr, sd = err),na.rm=TRUE)

    # Joint probability
    joint_mtx <- wll_mtx +  dnorm_cpp(fv_mtx, mean=0,sd=1)

    # Gaussian quadrature integral approximation
    gqi <- matrixStats::colLogSumExps(joint_mtx + log(cc$w) +(cc$x)^2,na.rm=TRUE)
    -sum(log(coefs)+gqi, na.rm=TRUE)
}


