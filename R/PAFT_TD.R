library(nloptr)
library(survival)
library(dplyr)
library(nnls) # updated weight
#------------ Function ------------#

#' Parametric AFT model with time dependent covariates
#'
#' @param dat Input data
#' @param X List of vectors indicating the parameters of any covariate that the user needs to introduce in this model
#' @param initi The initial values of betas
#' @param beta List of vectors indicating the effect of the corresponding covariate
#' @param gab The initial support(=SE)
#' @param dist Distribution assumed for log(T_0)
#' @param tol Allowable margin of error
#' @param maxiter The maximum number of iteration
#' @param algorithm Algorithm used for iteration
#'
#' @return $Estimated Betas $itration $status $message
#' @export
#'
#' @examples PAFT_TD(dat=dat1,X=c("Z1","X"),dist="log_normal") # True data: log-normal with b0=1, b1=1, b2=1)
#' @examples PAFT_TD(dat=dat4,X=c("Z1","X"),dist="log_normal") # True data: log-t(3) with b0=1, b1=1, b2=1
#'
PAFT_TD <- function(dat,X,initi=F,beta=NA,gab=NA,dist,tol=1.0e-5,maxiter=2000,algorithm="NLOPT_LN_COBYLA"){

  ## inital values
  if(initi==T){
    beta <- beta
    gab <- gab
  }else{
    dat1 <- as.data.frame(dat%>%group_by(ID)%>%arrange(stop)%>%slice(n()))
    X.mat <- as.matrix(dat1[,X])
    PAFT <- survreg(Surv(dat1[,"stop"],dat1[,"delta"]) ~ X.mat,dist="lognormal")
    beta <- PAFT$coefficients # initial betas
    gab <- PAFT$scale  # initial support set( = SE)
  }

  if(dist=="log_normal"){
    ## optimization function
    PAFT_opt <- function(theta){

      ## 1. Calculation of the observation time
      obs_time <- as.data.frame(dat%>%group_by(ID)%>%slice(n()))[,c("ID","stop")]
      names(obs_time) <- c("ID","obs_time")
      dat <- merge(dat,obs_time,by="ID",all.x=T)
      dat <- dat[order(dat$ID,dat$start),] # ordering


      ## 3. Likelihood function of the AFT with time-depnedent

      ## Term1 : Calculation of the Capital psi
      # : integrated from 0 to yi(observed time) of exp(-z(u)*beat^T)


      ## 1.
      # AA <- -theta[2]*dat[,X[1]]-theta[3]*dat[,X[2]]

      ## 2.
      # A1 <- -theta[2]*dat[,X[1]]
      # A2 <- -theta[3]*dat[,X[2]]
      # AA <- A1+A2

      ## 3.
      AA <- 0
      for(i in 1:length(X)){
        AA <- -theta[i+1]*dat[,X[i]]+AA
      }

      dat$psi <- (dat$stop-dat$start)*exp(AA)
      # dat$psi <- (dat$stop-dat$start)*exp(rowSums(-theta[2:(1+length(X))] *dat[,X]))
      final <- aggregate(dat[,c("psi")],list(dat[,"ID"]),FUN=sum)
      names(final) <- c("ID","psi")
      final <- final[order(final$ID),] # ordering


      ## Term2 : Calculation of psi : exp(-beta*z(yi))
      last_dat <- as.data.frame(dat%>%group_by(ID)%>%slice(n()))[,c("ID",X)]
      BB <- 0
      for(i in 1:length(X)){
        BB <- -theta[i+1]*last_dat[,X[i]]+BB
      }
      final$psi_d <- exp(BB)
      # final$psi_d <- exp(-theta[2]*last_dat[,"Z1"]-theta[3]*last_dat[,"X"])
      # final$psi_d <- exp(rowSums(-theta[2:(1+length(X))] *last_dat[,X]))
      final <- final[order(final$ID),] # ordering

      ## Joint delta
      delta <- as.data.frame(dat%>%group_by(ID)%>%slice(n()))[,c("ID","delta")]
      final <- merge(final,delta,by="ID",all.x=T)
      final <- final[order(final$ID),] # ordering

      ## likelihood function ##
      sigma <- theta[length(theta)]
      logPsi <- log(final$psi)
      dnorm_2 <- dnorm((logPsi-theta[1])/sigma)*(final$psi_d/(sigma*final$psi))
      S <- final$delta*dnorm_2 + (1-final$delta)*(1 - pnorm((logPsi-theta[1])/sigma))
      S <- ifelse(S<=1.0e-20,1.0e-20,S)

      lik <- sum(log(S))
      lik_inv <- -lik
      return(lik_inv)
    }

    # initial values
    # b <- c(1e-20,rep(1e-20,length(X)),1e-20)
    # b <- c(1,rep(1,length(X)),1)
    b <- as.numeric(c(beta,gab))

    # bounded values
    bound <- c(max(abs(b))*10,rep(max(abs(b))*10,length(X)),max(abs(b))*10)
    # bound <- c(Inf,rep(Inf,length(X)),Inf)


    # ----  nloptr function ----#
    res2 <- nloptr(x0=b,eval_f=PAFT_opt,lb = -bound,ub = bound,
                   opts = list("algorithm"=algorithm,"xtol_rel"=tol,maxeval=maxiter))
    # Results
    Betas_results <- res2$solution[1:(1+length(X)+1)] # Estimated betas
    names(Betas_results) <- c(paste("b",c(0:length(X)),sep=""),"sigma")
    itr <- res2$iterations # number of iterations of the algorithm
    status <- res2$status # integer value with the status of the optimization
    message <- res2$message # more informative message with the status of the optimization
    # --------------------------#

    # ----  optimum function ----#
    # res2 <- optim(b,PAFT_opt,lower=-bound,upper=bound,method = "L-BFGS-B")

    # # Results
    # Betas_results <- res2$par[1:(1+length(X)+1)] # Estimated betas
    # names(Betas_results) <- c(paste("b",c(0:length(X)),sep=""),"sigma")
    # message <- res2$message
    # convergence <- res2$convergence
    # count <- res2$counts

    # A <- list(EstBetas=Betas_results,count=count,convergence=convergence,
    #           message=message)
    A <- list(EstBetas=Betas_results,itr=itr,status=status,
              message=message)
  }else{
    A <- NULL
  }

  return(A)
}
