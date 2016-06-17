PME = function(mat){
    #SPSS
    sd=apply(mat,2,sd)
    ys=mat
    Theta=matrix(FALSE,q,q)
    S=Theta
    last_EBIC=Inf
    
    r = function(j,k){
    a1= sum((mat[,k]*ys[,j])^2)/(sd[j]*sum(mat[,k]^2))
    a2= sum((mat[,j]*ys[,k])^2)/(sd[k]*sum(mat[,j]^2))
    return(a1+a2)
    }
    
    EBIC = function(S){
        sum = 0
        for ( i in 1:q){
            sum=sum+log(sum(ys[,i]^2)/n/sd[i])    
        }
        sum=sum*n+sum(S)*log(n)+2*(1-log(n)/2.1/log(q))*log(choose(q*(q-1)/2,sum(S)/2))
        return(sum)
    }
    
    H = function(S,j){
        if ( sum(S[,j]) >0 ){
            Z = as.matrix(mat[,S[,j]])
            return(Z %*% solve(t(Z)%*%Z)%*% t(Z))
        } else {
            return(0)    
        }
    }
    
    repeat{
        temp = -Inf 
        pair = NULL 
        for(k in 2:q){
            for (j in 1:(k-1)){
                new = r(j,k)
                if (new>temp & Theta[j,k]==0 ){
                    temp = new
                    pair = c(j,k)
                }
            }    
        }
        
        S[pair[1],pair[2]]=TRUE
        S[pair[2],pair[1]]=TRUE
        new_EBIC=EBIC(S)
        if (new_EBIC < last_EBIC){
            last_EBIC=new_EBIC
            for (i in 1:q){
                ys[,i]= (diag(1,n)-H(S,i))%*%mat[,i]
            }
            Theta=S
        } else {
            break
        }
        
    }
    
 
    Omega=matrix(0,q,q)
    for ( i in 1:q){
        sigma = mean(((diag(1,n)-H(Theta,i))%*%mat[,i])^2)
        Omega[i,i] = 1/sigma 
        if ( sum(Theta[,i]>0)){
            Z = mat[,Theta[,i]]
            beta =  solve(t(Z)%*%Z)%*%t(Z)%*%mat[,i]
            Omega[Theta[,i],i]=-beta*sigma
        } 
    }
    
    Omega = (Omega+t(Omega))/2
    return(list(Omega=Omega,Theta=Theta))
}