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


#' S3 method of print for relvm class
#'
#'
#'
#' @export
#'
print.relvm <- function(object,...) {
    if (exists("par",object) | exists("pars",object)) {
        cat("\n..$par:\n");        print.data.frame(object$par)}
    if (exists("value",object)) {
        cat("\n..$value:\n");      print.default(object$value)}
    if (exists("counts",object)) {
        cat("\n..$counts:\n");     print.default(object$counts)}
    if (exists("convergence",object)){
        cat("\n..$convergence:\n");print.default(object$convergence)}
    if (exists("message",object)) {
        cat("\n..$message:\n");    print.default(object$message)}
    cat("\n..$field:\n",       names(object),"\n")
}



#' Pipe
#'
#' Pipe function from magrittr package.
#'
`%>%` <- magrittr:::pipe()



