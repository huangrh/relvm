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
                     fn      = pnll,   # prediction function
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
    sum(wts_row*(0.79817986835 + log(err^2) + ((score_row-mu - fl * fv)/err)^2)/2,
        (0.79817986835+(fv)^2)/2,na.rm=TRUE)
}
