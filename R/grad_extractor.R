#function to extract gradient at a point x from deriv
#' Title
#'
#' @param f_deriv
#' @param x_args
#'
#' @return
#' @export
#'
#' @examples
grad_extractor<-function(f_deriv,x_args){
  attr(do.call(f_deriv,as.list(x_args)),"gradient")[1,]
}
