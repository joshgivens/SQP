#Extract derivates (and hessians) from list of expressions
#' Title
#'
#' @param expression_list
#' @param xdim
#' @param hessian
#'
#' @return
#' @export
#'
#' @examples
deriv_extractor<-function(expression_list,xdim=2,hessian=TRUE){
  temp_deriv<-function(expr,xdim,hessian){deriv(expr,namevec=paste0("x",1:xdim),
                                                function.arg=paste0("x",1:xdim),
                                                hessian=hessian)}
  return(lapply(as.list(expression_list),temp_deriv,xdim=xdim,hessian=hessian))
}
