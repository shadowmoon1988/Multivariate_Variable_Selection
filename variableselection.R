require(huge)
require(extlasso)
require(MRCE)
require(foreach)
require(doParallel)
setwd("~/Brant Ying/Document of NUS/NUS/Research/Sparse")
source("function.R")
source("Extended Lasso.R")
source("Seq_Lasso.R")



#Main Part
## Generation
n=500
q=10
p=200

set.seed(100)

nsim=10
#result=NULL
#for( i in 1:nsim){
Simulation = function() {

    E0=huge.generator(n,q,graph="random",verbose = FALSE) # define own generator
    
    E1=npn_gaussian.generator(E0$data,mu_v=rep(0,q),sigma_v=rep(1,q))
    #E1=npn_power.generator(E0$data,mu_v=rep(0,q),sigma_v=rep(1,q))
    
    prob=floor(4*n^0.16)/p
    B = matrix(abs(rnorm(p*q,0,0.1/qnorm(0.875)))+4*n^(-0.15),nrow=p,ncol=q)*
        matrix(sample(c(-1,1),p*q,replace=TRUE,prob=c(0.4,0.6)),nrow=p,ncol=q)*
        matrix(sample(0:1,p*q,replace=TRUE,prob=c(1-prob,prob)),nrow=p,ncol=q)
    
    X=matrix(rnorm(n*p),nrow=n,ncol=p) #design different X structure
    Y=X%*%B+E1
    
    ##Estimation
    ###seperate lasso

    B_sepLasso=matrix(NA,nrow=p,ncol=q) #B_hat=solve(t(X)%*%X)%*%(t(X)%*%Y)
    for (i in 1:q){
        cv_lambda=cv.extlasso(X,Y[,i],family="normal",plot=FALSE)$lambda # CV get the best lambda
        lasso_reg=extlasso(X,Y[,i],family="normal",intercept=FALSE) 
        B_sepLasso[,i]=predict(lasso_reg,mode="lambda",at=cv_lambda)[-1]
    }
    
    # MRCE
    #lam1.vec=rev(10^seq(from=-2, to=0, by=0.5))
    #lam2.vec=rev(10^seq(from=-2, to=0, by=0.5))
    #B_MRCE1=mrce(Y=Y, X=X, lam1.vec=lam1.vec, lam2.vec=lam2.vec, method="cv")$Bhat
    #B_MRCE2=mrce(Y=Y, X=X, lam1=10^(-1.5), lam2=10^(-1), method="single")$Bhat
    
    # Extenting Lasso
    #B_exlasso=ExLasso(X,Y,t=1)
    
    # Our Method
    #B_seq=Seq_Lasso(X,Y,B_sepLasso)

    B_seq1=Seq_Lasso1(X,Y,B_sepLasso)
    
    #result = rbind(result,evaluation(B_seq,B))
    return(evaluation(B_seq1,B))
    # Evaluation
#     result = rbind(result,
#                     c(evaluation(B_seq,B),
#                     evaluation(B_sepLasso,B),
#                     evaluation(B_MRCE2,B),
#                     evaluation(B_exlasso,B))
#                    )
}


n=100
q=10
p=100

set.seed(100)


no_cores = detectCores() -1
cl<-makeCluster(no_cores)
registerDoParallel(cl)
foreach(i=1:20, .combine=cbind, .export = c("n","p","q"), .packages=c('huge','extlasso','MRCE')) %dopar% Simulation()


stopImplicitCluster()



n=250

# write the parallel
# write the detail of computation.
for ( i in 1:nsim){print(Simulation())}