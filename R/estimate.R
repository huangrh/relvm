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
#' Estimate multiple measure groups of the random effect latent variable model
#'
#' @param object A mstbl object.
#' @param groups A vector of measure group names. The default is NULL, in which
#'   case a vector of all groups will be generated accordingly.
#' @param fit A list of fitting parameters. \itemize{ \item init: Initial values for mu, fl, and err
#'   term in a list. fl is the factor loading. They will be initialized
#'   generally if it is null. The default is a list with for all mu and one for
#'   others. \item predict: The default is TRUE.}
#'
#' @return An list of S3 object of class "relvm" with estimated parametes.
#' @seealso \code{\link{mstbl}}
#' @importFrom pracma hessian
#'
#' @export
#'
relvm <- function(object,groups=NULL,fit=list(init=NULL)) {

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
    fit.default = list(init=NULL,predict=TRUE)
    fit         = merge_default(x=fit,y=fit.default)

    init    = fit[["init"]]
    predict = fit[["predict"]]

    # ------------------------------------------------------------------#
    # Call relvm_single
    # snowfall::sfInit(parallel=TRUE,cpus=2);snowfall::sfExportAll()
    # snowfall::sfExport(create_measure_tbl)

    allout <- sapply(groups, relvm_single_true, df=alldf,
                     init = init, predict = predict,simplify = FALSE)

    # snowfall::sfRemoveAll()
    # snowfall::sfStop()

    # ------------------------------------------------------------------#
    # After Relvm:
    # Merge the predicted group score if there is multiple group.
    preds <- alldf[,1,drop=FALSE] # take the column "provider_id"
    for (group in allout) {preds <- merge(x=preds,y=group$pred,all=TRUE)}
    colnames(preds) <- gsub("pred_","",colnames(preds))

    # Calculate the summary score.
    hospital_score <- rstarating::sum_score(preds)


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
                                    summary_score=hospital_score,
                                    counts= as.matrix(counts),
                                    value = as.matrix(value),
                                    message=as.matrix(message),
                                    convergence = as.matrix(convergence)),class="relvm")
    (allout)
}


# Default parameter helper:
merge_default = function(x,y) {
    i = is.na(match(names(y), names(x)))
    if (any(i)) {
        iy = names(y)[i];
        x[iy] = y[iy]}
    x}


#' Estimation Of The Random Effect Latent Variable Model Parameters
#'
#' Estimate the random effect latent variable model
#'
#' @param mstbl_std The standized measure score table.
#' @param wts_tbl The measure score weight table.
#' @param init Initial values for mu, fl, and err term in a list. fl is the
#'   factor loading. They will be initialized generally if it is null. The
#'   default is a list with for all mu and one for others.
#' @param predict The default is TRUE.
#'
#' @return An object of S3 class "relvm" with estimated parametes.
#'
relvm_single_true <- function(group, df, init, predict) {
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
        init <- unlist(list(mu  = rep(0, nc),
                            fl  = rep(0.8, nc),
                            err = rep(0.8, nc)))}


    #--------------------------------------------------------#
    # Fit the function
    fit <- optim(par     = init,      # Model parameter
                 fn      = venll18, # venll11m,   # Estimation function
                 gr      = NULL,
                 method  = "L-BFGS-B", # "L-BFGS-B"
                 control = list(maxit=1000), # set factr=1e-8
                 hessian = FALSE,
                 score   = mstbl_std,
                 wts     = wts_tbl)

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
        pred_out          <- relvm:::pred(mstbl_std,wts_tbl,pms=fit$par);

        colnames(pred_out)<- paste(colnames(pred_out),group,sep="_")
        pred_out          <- cbind(subdat$pid,pred_out)
        fit$pred          <- pred_out[,1:2]

        #
        fit$stderr        <- pred_out[,c(1,3)]
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

# Simplified Vectorized Estimation function
venll18 <- function(par,score,wts) {
    nr <- nrow(score); nc <- ncol(score)

    # matrix: score, wts, mu,fl,err
    mu  <- matrix(par[grepl("mu", names(par))], nrow=nr,ncol=nc,byrow=TRUE) #
    fl  <- matrix(par[grepl("fl", names(par))], nrow=nr,ncol=nc,byrow=TRUE) # loading_m
    err <- abs(matrix(par[grepl("err", names(par))],nrow=nr,ncol=nc,byrow=TRUE)) # delta_m

    #
    log_fh <- - 1/2*log(2*pi) + rowSums(-wts/2 * (log(2*pi)+2*log(err)), na.rm=TRUE)
    ah     <- - 1/2 * (1+rowSums(wts * fl^2 / err^2,na.rm=TRUE))
    bh     <- rowSums(wts/err^2 * (score-mu) * fl,  na.rm=TRUE)
    ch     <- rowSums(-wts / (2 * err^2) * (score-mu)^2,na.rm=TRUE)

    #
    -sum(log_fh + ch -bh^2/(4*ah) + log(-pi/ah)/2, na.rm=TRUE)
}


