#' Title
#'
#' @param init
#' @param f_exp
#' @param h_exp
#' @param xdim
#' @param meq
#' @param method
#' @param max_iter
#'
#' @return
#' @export
#'
#'
#' @examples
SQP_auto<-function(init,f_exp,h_exp,xdim=2,meq=0,method="standard",max_iter=10000){
  if (method=="standard"){
    hessian=TRUE
  }
  else if (method=="BFGS"){
    hessian=FALSE
  }

  #Get derivatives
  f_full<-deriv_extractor(f_exp,xdim = xdim)[[1]]
  h_full<-deriv_extractor(h_exp,xdim = xdim)

  #Get required functions from derivatives

  #Get values
  f_val<-function(x){val_extractor(f_full,x_args=x)}
  h_val<-function(x)unlist(lapply(h_full,val_extractor,x_args=x))

  #Get gradients
  grad_f<-function(x){
    grad_extractor(f_full,x_args=x)
  }
  grad_h<-function(x) {
    do.call(rbind,lapply(h_full,grad_extractor,x_args=x))
  }

  if (method=="standard"){
    #Get hessian
    hess_L<-function(x,mu,lambda) {
      hess_f<-hess_extractor(f_full,x)
      hess_h_list<-lapply(h_full,hess_extractor,x_args=x)
      hess_h_sum=Reduce(`+`,mapply(`*`,c(mu,lambda),hess_h_list,SIMPLIFY=FALSE))

      return(hess_f+hess_h_sum)
    }
    #Run SQP
    out<-SQP(init,f = f_val,grad_f = grad_f,h= h_val,grad_h = grad_h ,hess_L = hess_L,meq=meq,method=method)
  }
  else if (method=="BFGS"){
    #Run SQP
    out<-SQP(init,f = f_val,grad_f = grad_f,h= h_val,grad_h = grad_h ,meq=meq,method=method,max_iter=max_iter)
  }


}

