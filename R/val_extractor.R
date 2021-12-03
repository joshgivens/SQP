#function to extract value at a point x from deriv
#' Title
#'
#' @param f_deriv
#' @param x_args
#'
#' @return
#' @export
#'
#' @examples
val_extractor<-function(f_deriv,x_args){
  do.call(f_deriv,as.list(x_args))[1]
}



