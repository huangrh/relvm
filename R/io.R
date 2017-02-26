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


#' Output Directory
#'
#' Check the input string as a directory. Create a directory if it doesn't exists.
#'
#' @param x A full path directory string.
#'
#' @return A directory string.
#'
#' @export
out_dir <- function(x) {
    #x <- gsub("\\","/",x)
    if (!dir.exists(x)) dir.create(x,recursive=TRUE)
    x
}

#' Pipe
#'
#' Pipe function from magrittr package.
#'
#'
`%>%` <- magrittr:::pipe()
