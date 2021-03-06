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
#' @param fit A list of fitting parameters. \itemize{ \item init: Initial values
#'   for mu, fl, and err term in a data frame. fl is the factor loading. They
#'   will be initialized generally if it is null. The default is zero for mu and
#'   one for the rest. \item predict: The default is TRUE, which calculate the
#'   group score for each hospital. \item use_wt: A vector of strings, eg.
#'   c("step1", "step2"), corresponding to step 1 (optimize the model
#'   parameters) and step 2 (predict the group score) in the latent variable
#'   model, respectively. For example, fit = list(use_wt = "step2") tells the program to
#'   ignore the weight in step 1 and only use the weight in step 2. }
#' @param file File name in a rds format to cache the results on hard drive.
#'
#' @return An list of S3 object of class "relvms" with estimated parametes.
#'
#' @seealso \code{\link{mstbl}}
#'
#' @examples
#' # require(rstarating)
#' # To fit just one group: 'outcome_mort'
#' fit <- relvm(mstbl(cms_star_rating_input_2017dec),groups="outcome_mort")
#'
#' # use weight in step 2 only (by default, weigt are used in step 1 and step 2)
#' fit <- relvm(mstbl(cms_star_rating_input_2017dec),groups="outcome_mort", fit=list(use_wt="step2"))
#'
#' # Use weight in step 1 only (for test purpose only).
#' fit <- relvm(mstbl(cms_star_rating_input_2017dec),groups="outcome_mort", fit=list(use_wt="step1"))
#'
#' @export
#'
relvm <- function(object,groups=NULL,fit=list(inits=NULL,use_wt=c("step2")), file = NULL) {

    if (!is.null(file) && file.exists(file)) {
        # if the cached file exists
        allout <- readRDS(file)
    } else {

        # -------------------------------------------------------
        # Merge both tables of the measure score and weights.
        alldf  <- merge(x=object$mstbl_std, y=object$wtbl,by = "ccnid",all=TRUE)

        # Check & update "groups"
        mtbl       <- create_measure_tbl(alldf)
        all_groups <- unique(mtbl$group)
        if (is.null(groups)) {
            groups <- all_groups
        } else if (any(groups %in% all_groups)) {
            groups <- groups[groups %in% all_groups]
        } else stop("The group name do not match.")

        # Fit control
        fit_default   = list(inits=NULL,predict=TRUE,use_wt=c("step1","step2"))
        extra_default <- fit_default[!(names(fit_default) %in% names(fit))]
        fit[names(extra_default)] <- extra_default

        inits    = fit[["inits"]]
        predict  = fit[["predict"]]
        use_wt   = fit[["use_wt"]]

        # ------------------------------------------------------------------#
        # Call relvm_single
        start_time = Sys.time()
        cat(sprintf("Fitting start at: %-15s\n",start_time))

        allout <- sapply(groups,
                         relvm_single_true,
                         df       =   alldf,
                         inits    =   inits,
                         predict  =   predict,
                         use_wt   =   use_wt,
                         simplify =   FALSE)

        cat("\n","Total time: ", as.character.Date(Sys.time() - start_time),"\n")

        # ------------------------------------------------------------------#
        # After Relvm:
        # Merge the predicted group score if there is multiple group.
        preds <- alldf[,"ccnid",drop=FALSE] # take the column "ccnid"
        for (group in allout) {preds <- merge(x=preds,y=group$pred,by="ccnid",all=TRUE)}
        colnames(preds) <- gsub("pred_","",colnames(preds))

        # Calculate the summary score.
        hospital_score <- rstarating::sum_score(preds)
        hospital_score <- merge.data.frame(x=hospital_score,y=object$report_indicator,
                                           by='ccnid',all.x=TRUE)
        # hospital_score <- subset(hospital_score, report_indicator == 1)

        # Merge factor loadings and other parametes.
        pars <- data.frame()
        for (group in allout) {pars = rbind(pars,group$par)}

        # convergence
        convergence<- data.frame(convergence=vapply(allout,function(x) {x$convergence},c(0)))
        value      <- data.frame(value=vapply(allout,      function(x) x$value,c(0)))
        message    <- data.frame(message=vapply(allout,    function(x) x$message,"0"),stringsAsFactors = FALSE)
        counts     <- t(as.data.frame(vapply(allout,       function(x) x$counts,c(0L,0L))))

        #output
        allout$groups <- structure(list(
            preds        = preds,
            pars         = pars,
            summary_score= hospital_score,
            counts       = as.matrix(counts),
            value        = as.matrix(value),
            message      = as.matrix(message),
            convergence  = as.matrix(convergence),
            mstbl_std    = object$mstbl_std,
            wtbl         = object$wtbl
        ),class="relvm")

        # save the data
        if (!is.null(file)) {
            # check if the file directory exists. Create one if not.
            file_dir <- dirname(file)
            if (!dir.exists(file_dir)) file_dir_created <- dir.create(file_dir,recursive =TRUE)
            saveRDS(allout,file=file)
        }
    }

    # Output
    structure(allout, class="relvms")
}




#' Estimation Of The Random Effect Latent Variable Model Parameters
#'
#' Estimate the random effect latent variable model
#'
#' @param group A measure group name.
#' @param df The standardized measure score and measure weight table (alldf).
#' @param inits Initial values for mu, fl, and err term in a list. fl is the
#'   factor loading. They will be initialized generally if it is null. The
#'   default is a list with for all mu and one for others.
#' @param predict The default is TRUE.
#' @param use_wt A vector of strings, possible contains "step1" and "step2",
#'   corresponding to step 1 (optimize the model parameters) and step 2 (predict
#'   the group score) in the latent variable model, respectively.
#'
#' @return An object of S3 class "relvm" with estimated parametes.
#'
relvm_single_true <- function(group, df, inits, predict, use_wt) {
    # -------------------------------------------------------#
    # Prepare to fit
    # start of the cycle
    start_time = Sys.time()
    cat(sprintf("Fitting: %-15s =>",group))

    # -------------------------------------------------------#
    # data table & weight table
    subdat    <- relvm:::sub1group(group,df)
    mstbl_std <- as.matrix(subdat$mstbl_std)
    wts_tbl   <- wts_tbl_step1 <- wts_tbl_step2 <- as.matrix(subdat$wtbl)


    # -------------------------------------------------------#
    # Measure weighting in step 1 and step 2
    #
    # if not use wt in step 1, just set the weight to one.
    #
    if (!("step1" %in% use_wt)) {
        wts_tbl_step1 = lapply(as.data.frame(wts_tbl_step1), function(x) {
            ifelse(is.na(x),NA,1)
        }) %>% as.data.frame() %>% as.matrix()
    }

    # if no use use wt in step2
    if (!("step2" %in% use_wt)) {
        wts_tbl_step2 = lapply(as.data.frame(wts_tbl_step2), function(x) {
            ifelse(is.na(x),NA,1)
        }) %>% as.data.frame() %>% as.matrix()
    }

    # -------------------------------------------------------#
    # Setup and initialize the parameters
    nc <- ncol(mstbl_std);
    init_default <- unlist(list(mu  = rep(0, nc),
                                fl  = rep(0.7, nc),
                                err = rep(0.8, nc)))

    if (is.null(inits)) {
        init <- init_default
    } else {
        inits_group <- inits[match(colnames(mstbl_std),inits$name),]
        init <- unlist(list(mu  = inits_group[,"mu"],
                            fl  = inits_group[,"fl"],
                            err = inits_group[,"err"]))
        init <- ifelse(is.na(init),init_default,init) # set as the default if missing.
    }


    #--------------------------------------------------------#
    # Fit the function
    fit <- optim(par     = init,      # Model parameter
                 fn      = venll18,   # venll11m,   # Estimation function
                 gr      = NULL,
                 method  = "L-BFGS-B", # "L-BFGS-B"
                 control = list(maxit=1000), # set factr=1e-8
                 hessian = FALSE,
                 score   = mstbl_std,
                 wts     = wts_tbl_step1)

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
        pred_out          <- relvm:::pred(mstbl_std, wts_tbl_step2, pms=fit$par);

        colnames(pred_out)<- paste(colnames(pred_out),group,sep="_")
        pred_out          <- cbind(subdat$ccnid,pred_out)
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

    fit$mstbl_std = cbind(subdat$ccnid, mstbl_std)
    fit$wtbl_step1= cbind(subdat$ccnid, wts_tbl_step1)
    fit$wtbl_step2= cbind(subdat$ccnid, wts_tbl_step2)
    fit$wtbl      = cbind(subdat$ccnid, wts_tbl)

    # Output
    cat(" : ", as.character.Date(Sys.time() - start_time))
    cat(", ",fit$message,       "\n")

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


