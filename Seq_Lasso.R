source("Coefficient_Estimation.R")
source("Precision_estimation.R")
source("SLASSO.R")
Seq_Lasso1 = function (X, Y, B_int,nonparanormal = TRUE){
    
    B_old = B_int# initial
    # SLasso with SPSS
    itr_count=0
    B_new=matrix(NA,nrow=p,ncol=q)
    trans_e=matrix(0,n,q)
    repeat{
        e=Y-X%*%B_old
        
        if (nonparanormal){
            for (i in 1:q){
                trans_e[,i]=mest(e[,i],e[,i])
            }
        }
        Omega=PME(trans_e)$Omega
        
        ### maximize the likelihood respect to b
        ### Sequential method
        for (i in 1:q){
            Z=trans_e[,-i]
            cond_mean=Z %*% (-Omega[-i,i]/Omega[i,i])
            B_new[,i]=SLasso1(X,Y[,i],e[,i],cond_mean)
        }
        
        itr_count=itr_count+1
        print(itr_count)
        if ( stoprule(B_new,B_old) || itr_count>50 ){ 
            return(B_old)
        } else {
            B_old=B_new
        }
    }
}