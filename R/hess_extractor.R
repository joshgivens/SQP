#function to extract hessian at a point x from deriv
#' Title
#'
#' @param f_deriv
#' @param x_args
#'
#' @return
#' @export
#'
#' @examples
hess_extractor<-function(f_deriv,x_args){
  attr(do.call(f_deriv,as.list(x_args)),"hessian")[1,,]
}
