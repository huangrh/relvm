# Copyright (C) 2016-2017 Ren-Huai Huang <huangrenhuai@gmail.com>
#
# This file is part of relvm.
#
# relvm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# relvm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with relvm.  If not, see <http://www.gnu.org/licenses/>.



#' Multiple Estimation Of The Random Effect Latent Variable Model Parameters
#'
#' Estimate multiple groups of the random effect latent variable model
#'
#' Ad = adaptive
#'
#' @param object A mstbl object.
#' @param groups A vector of measure group names. The default is NULL, in which
#'   case a vector of all groups will be generated accordingly.
#' @param fit A list of fitting parameters. \itemize{ \item qpoints: The numbe
#'   of the quadrature points. \item inits: Initial values for mu, fl, and err
#'   term in a list. fl is the factor loading. They will be initialized
#'   generally if it is null. The default is a list with for all mu and one for
#'   others. \item predict: The default is TRUE. \item adaptive: noad, no
#'   adaptive or ad, use adaptive. }
#'
#' @return An list of S3 object of class "relvm" with estimated parametes.
#' @seealso \code{\link{mstbl}}
#' @importFrom pracma hessian
#'
#' @export
#'
relvm_adaptive <- function(object,
                     groups=NULL,
                     fit=list(qpoints=15,
                              inits=NULL,
                              predict=TRUE,
                              adaptive=c("ad","noad"))) {

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

    # Fit control default
    fit_default   <- list(qpoints = 15,
                          inits   = NULL,
                          predict = TRUE,
                          adaptive= c("ad","noad"))
    extra_default <- fit_default[!(names(fit_default) %in% names(fit))]
    fit[names(extra_default)] <- extra_default

    #
    qpoints = fit[["qpoints"]]
    inits   = fit[["inits"]]
    predict = fit[["predict"]]
    adaptive= fit[["adaptive"]][1]


    # ------------------------------------------------------------------#
    start_time = Sys.time()
    cat(sprintf("Fitting start at: %-15s by method: %-15s\n",start_time,adaptive))

    # Call relvm_single
    allout <- sapply(groups, relvm_adaptive_single, df=alldf, qpoints=qpoints,
                      inits = inits, predict = predict, adaptive=adaptive,simplify = FALSE)

    cat("\n","Total time: ", as.character.Date(Sys.time() - start_time),"\n")
    # ------------------------------------------------------------------#
    # After Relvm:
    # Merge the predicted group score if there is multiple group.
    preds <- alldf[,1,drop=FALSE] # take the column "ccnid"
    for (group in allout) {preds <- merge(x=preds,y=group$pred,all=TRUE)}
    colnames(preds) <- gsub("pred_","",colnames(preds))

    # Calculate the summary score.
    hospital_score <- rstarating::sum_score(preds)
    hospital_score <- merge.data.frame(x=hospital_score,y=object$report_indicator,
                                       by='ccnid',all.x=TRUE)
    hospital_score <- subset(hospital_score, report_indicator == 1)

    # Merge factor loadings and other parametes.
    pars <- data.frame()
    for (group in allout) {pars = rbind(pars,group$par)}

    # groups: convergence, value, message, and counts
    convergence<- data.frame(convergence=vapply(allout,function(x) {x$convergence},c(0)))
    value      <- data.frame(value=vapply(allout,  function(x) x$value,c(0)))
    message    <- data.frame(message=vapply(allout,function(x) x$message,"0"),stringsAsFactors = FALSE)
    counts     <- t(as.data.frame(vapply(allout,   function(x) x$counts,c(0L,0L))))

    #output
    allout$groups <- structure(list(preds=preds,pars=pars,
                                    summary_score=hospital_score,
                                    counts  = as.matrix(counts),
                                    value   = as.matrix(value),
                                    message = as.matrix(message),
                                    convergence = as.matrix(convergence)),class="relvm")
    (allout)
}


#' Estimation Of The Random Effect Latent Variable Model Parameters
#'
#' Estimate the random effect latent variable model
#'
#' @param mstbl_std The standized measure score table.
#' @param wts_tbl The measure score weight table.
#' @param qpoints The numbe of the quadrature points.
#' @param inits Initial values for mu, fl, and err term in a list. fl is the
#'   factor loading. They will be initialized generally if it is null. The
#'   default is a list with for all mu and one for others.
#' @param predict The default is TRUE.
#'
#' @return An object of S3 class "relvm" with estimated parametes.
#'
relvm_adaptive_single <- function(group, df, qpoints,inits,predict,adaptive) {
    # -------------------------------------------------------#
    # Prepare to fit
    # start of the cycle
    start_time = Sys.time()
    cat(adaptive,"-",qpoints,"qpts",": ")
    cat(sprintf("Fitting: %-15s =>",group))

    # data table & weight table in a measure group.
    subdat    <- relvm:::sub1group(group,df)
    mstbl_std <- as.matrix(subdat$mstbl_std)
    wts_tbl   <- as.matrix(subdat$wtbl)

    # Setup and initialize the parameters
    nc <- ncol(mstbl_std);
    init_default <- unlist(list(mu  = rep(0.5, nc),
                                fl  = rep(0.5, nc),
                                err = rep(0.5, nc)))

    if (is.null(inits)) {
        init <- init_default
    } else {
        inits_group <- inits[match(colnames(mstbl_std),inits$name),]
        init <- unlist(list(mu  = inits_group[,"mu"],
                            fl  = inits_group[,"fl"],
                            err = inits_group[,"err"]))
        init <- ifelse(is.na(init),init_default,init)
    }

    # cc$x & cc$w
    cc <- pracma::gaussHermite(qpoints);

    #--------------------------------------------------------#
    # Fit the function
    fit <- optim(par     = init,      # Model parameter
                 fn      = venll10a2,    # Estimation function
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
        pred_out           <- relvm:::pred(mstbl_std,wts_tbl,pms=fit$par);

        colnames(pred_out) <- paste(colnames(pred_out),group,sep="_")
        pred_out           <- cbind(subdat$ccnid,pred_out)
        fit$pred           <- pred_out[,1:2]

        #
        fit$stderr         <- pred_out[,c(1,3)]
    }

    # Add three fields to the output
    init_names  <- names(init);
    fit$init    <- data.frame(name = names(subdat$mstbl_std),
                              fl = init[grepl("fl", init_names)],
                              mu = init[grepl("mu", init_names)],
                              err= init[grepl("err",init_names)], row.names=NULL)

    fit$mstbl_std = cbind(subdat$ccnid,mstbl_std)
    fit$wtbl      = cbind(subdat$ccnid,wts_tbl)

    # Output
    cat(" : ", as.character.Date(Sys.time() - start_time),"\n")
    structure(fit,class="relvm")
}



# adaptive estimate function
venll10a2 <- function(par,score,wts,cc,qpoints,adaptive) {

    # Simplified normal density function.
    dnorm2 <- function(x,mean=0,sd=1) -(log(2 * pi) +2*log(sd)+((x-mean)/sd)^2)/2

    # Reconstruction of the parameters
    nr <- nrow(score); nc <- ncol(score)
    mu <- par[grepl("mu", names(par))]         #
    fl <- par[grepl("fl", names(par))]         # factor loading
    err<- par[grepl("err", names(par))]
    err <- abs(err)

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
    joint_mtx <- wll_mtx +  dnorm2(fv_mtx, mean=0,sd=1)

    # Gaussian quadrature integral approximation
    gqi <- matrixStats::colLogSumExps(joint_mtx + log(cc$w) +(cc$x)^2,na.rm=TRUE)
    -sum(log(coefs)+gqi, na.rm=TRUE)
}


