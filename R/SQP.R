
#' Title
#'
#' @param init
#' @param f
#' @param grad_f
#' @param hess_L
#' @param h
#' @param grad_h
#' @param meq
#' @param method
#' @param max_iter
#'
#' @return
#' @export
#' @import quadprog
#'
#' @examples
SQP<-function(init,f,grad_f,hess_L=NULL,h,grad_h,meq=0,method="standard",max_iter=10000){
  # init : list of init params, 1st element x_0, 2nd element mu_0, 3rd element, lambda_0)
  # f: function to be minismed
  # grad_f: gradient function of f
  # h: constraint functions (function should return vector output for each constraint)
  #         equality constraints followed by inequality constraints)
  #grad_h: matrix representing gradient of each constrain function.
  #           Each row is a gradient function
  #meq: number of constrains which are equality
  nconstr<-length(h(init[[1]]))
  meq_vec<-1:meq
  #Sort names of initialial values by case
  if(meq==0){
    #names(init)<-c("x","lambda")
    meq_vec=NULL
    mineq_vec<-1:nconstr
  }
  else if(meq==nconstr){
    #names(init)<-c("x","mu")
    meq_vec=1:meq
    mineq_vec=NULL
  }
  else{
    #names(init)<-c("x","mu","lambda")
    meq_vec=1:meq
    mineq_vec=(meq+1):nconstr
  }

  #Inititalise before loop
  x_k<-init$x
  mu_k<-init$mu
  lambda_k<-init$lambda
  if(method=="BFGS"){
    B_k=init$B
  }

  x_mat<-matrix(rep(NA,max_iter*length(x_k)),ncol=max_iter)
  mu_mat<-matrix(rep(NA,max_iter*length(mu_k)),ncol=max_iter)
  lambda_mat<-matrix(rep(NA,max_iter*length(lambda_k)),ncol=max_iter)
  f_vec<-rep(NA,max_iter)

  x_mat[,1]<-x_k
  mu_mat[,1]<-mu_k
  lambda_mat[,1]<-lambda_k
  f_vec[1]<-f(x_k)

  for (k in 2:max_iter){
    ##Set uo parameters for quadratic programme
    grad_fk=grad_f(x_k)
    grad_hk=grad_h(x_k)
    hk=h(x_k)

    #Calculate BFGS
    if(method=="standard"){
      hess_Lk=hess_L(x_k,mu_k,lambda_k)
    }
    else{
      hess_Lk=B_k
    }

    print(k)
    #if(k==37){browser()}
    #Construct QP
    QP_k<-solve.QP(Dmat=hess_Lk, dvec=-grad_fk, Amat=-grad_hk,bvec = hk, meq = meq)

    P_k<-QP_k$solution
    u_k<-(QP_k$Lagrangian)[meq_vec]
    l_k<-(QP_k$Lagrangian)[mineq_vec]

    #Do Line search
    x_k<-x_k+P_k
    mu_k<-u_k
    lambda_k<-l_k

    #Calculate BFGS if used
    if (method=="BFGS"){
      #browser()
      s_k<-x_k-x_mat[,k-1]
      y_k<-grad_f(x_k)-grad_f(x_mat[,k-1])+t(grad_h(x_k)-grad_h(x_mat[,k-1]))%*%c(mu_k,lambda_k)

      #Damping coefficient
      if(t(s_k)%*%y_k>=0.2*t(s_k)%*%B_k%*%s_k){
        theta_k=1
      }
      else{
        theta_k=as.vector((0.8*t(s_k)%*%B_k%*%s_k)/(t(s_k)%*%B_k%*%s_k-t(s_k)%*%y_k))
      }
      r_k=theta_k*y_k+(1-theta_k)*B_k%*%s_k

      #Actual BFGS update
      B_k<-B_k-
        (B_k%*%s_k%*%t(s_k)%*%B_k)/(as.vector(t(s_k)%*%B_k%*%s_k))+
        (r_k%*%t(r_k))/(as.vector(t(s_k)%*%r_k))
    }

    x_mat[,k]<-x_k
    mu_mat[,k]<-mu_k
    lambda_mat[,k]<-lambda_k
    f_vec[k]<-f(x_k)

    #Break point
    if (sum(abs(x_mat[,k]-x_mat[,k-1]))<1e-16){
      break
    }
  }
  return(list(x_mat=x_mat[,1:k],lambda_mat=lambda_mat[,1:k],mu_mat=mu_mat[,1:k],f_vec=f_vec[1:k],iter=k))

}
