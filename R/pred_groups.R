# Copyright (C) 2016-2018 Ren-Huai Huang <huangrenhuai@gmail.com>
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


#' Multiple Group score
#'
#' Predict the group scores (random effect latent variabl) from the extimated
#' parameters
#'
#' @param x An object, whose class attribute is "relvms", i.e. multiple group of
#'   relvm. see also relvm function.
#' @param newpars The matrix of modified parameters (factor loadings, etc.).
#' @param newdata The new data table. this table include standardized measure
#'   score and measure weight. See example below.
#' @param level The level of confident interval. This parameter may not be used
#'   any more.
#' @return A updated object.
#'
#' @examples
#' # predict with new parameters
#' fit  <- relvm(mstbl(rstarating::cms_star_rating_input_2017dec))
#' pred <- predict(fit,newpars=fit$groups$pars) # update the newpars.
#' star <- rating(x=pred$groups$summary_score,method="rclus2",score_col="sum_score",iter.max=5000)
#'
#' # predict with new data
#' newdata = list(mstbl_std = fit17$groups$mstbl_std, wtbl = fit17$groups$wtbl)
#' pred_newdata <- predict(fit,newdata=newdata) # update the newpars.
#'
#' @export
#'
predict.relvms <- function(x, newpars=NULL, newdata = NULL, level = 0.95) {
    begin_time = Sys.time()
    cat(toString(begin_time),": Runing predict.relvms\n")

    # New parameters (if there is no new parameters, take the pars from the object)
    if (is.null(newpars)) newpars <- x$groups$pars

    # new data
    if (!is.null(newdata)) {
        # merge the measure score table with the weight table, so the ccnid is matched for a hospital.
        if (exists("mstbl_std", newdata)) {
            x$groups$mstbl_std = newdata$mstbl_std
        } else {
            message("Warning: newdata doesn't have mstbl_std table.")
        }

        if (exists("wtbl",      newdata)) {
            x$groups$wtbl      = newdata$wtbl
        } else {
            message("Warning: newdata doesn't have wtbl table.")
        }
    }

    # mtbl
    mtbl = rstarating::create_measure_tbl(x$groups$mstbl_std)

    # merge
    alldf <- merge(x$groups$mstbl_std, x$groups$wtbl, by = "ccnid", all=TRUE)

    # Find the unique measure groups.
    groups        <- unique(rstarating::create_measure_tbl(alldf)[["group"]])

    # ---------------------------------------
    # call function pred_single_measure_group
    allout <- sapply(groups,
                     relvm:::pred_single_measure_group,
                     alldf= alldf,
                     pars = newpars,
                     mtbl = mtbl,
                     simplify = FALSE)


    # ---------------------------------------
    # after pred_single_...
    # 1. Merge the predicted group score if there is multiple group.
    preds <- alldf[,"ccnid",drop=FALSE] # take the column "ccnid"
    for (group in allout) {preds <- merge(x=preds,y=group$pred,by="ccnid",all=TRUE)}
    colnames(preds) <- gsub("pred_","",colnames(preds))


    # 2. Calculate the summary score.
    hospital_score <- rstarating::sum_score(preds)
    hospital_score <- merge.data.frame(x=hospital_score,
                                       y=x$groups$summary_score[c("ccnid","report_indicator")],
                                       by='ccnid',all.x=TRUE)
    # hospital_score <- subset(hospital_score, report_indicator == 1)

    # ---------------------------------------
    # Output
    x$groups <- structure(list(
        preds        = preds,
        pars         = newpars,
        summary_score= hospital_score,
        mstbl_std    = x$groups$mstbl_std,
        wtbl         = x$groups$wtbl,
        predict_tag  = "predict.relvms"
    ), class="relvm")

    cat("Total time: ", as.character.Date(Sys.time() - begin_time),"\n")

    #
    (x)
}



#' Group score prediction from a single group
#'
#' Predict/calculate a group score from a single group.
#'
#' @param group The measure group specified.
#' @param alldf  A data frame merged from both standardized measure scores and
#'   measure weights.
#' @param pars The parameters in a data frame, including factor loading, err and
#'   offset, etc.
#' @return relvm class object, which includes the predicted group scores and the
#'   corresponding standard deviations.
#'
#'
pred_single_measure_group <- function(group, alldf, pars, mtbl) {

    # ------------------
    # data table & weight table
    subdat    <- relvm:::sub1group(group, alldf)
    mstbl_std <- as.matrix(subdat$mstbl_std)
    wts_tbl   <- as.matrix(subdat$wtbl)



    # Extract the pms (factor loadings, stderr of the factor loading, etc)
    pms <- subset(pars, name %in% mtbl[mtbl$group %in% group,"name"])

    # Run the pred
    pred_out <- cbind(subdat$ccnid, relvm:::pred(mstbl_std,wts_tbl,pms))
    colnames(pred_out)[-1] <- paste(colnames(pred_out)[-1],group,sep="_")

    # Output
    structure(list(
        pred   = pred_out[,1:2],   # predicted group score
        stderr = pred_out[,c(1,3)] # predicted standard deviation
    ),class="relvm")

}
