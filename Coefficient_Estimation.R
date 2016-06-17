SLasso = function(X,y,e,cm){
        Fun = function(t,cond){ #1st negative deriavetive of obj function
            cond_X=X[,cond]
            dim(cond_X)=c(n,length(cond))
            dim(t)=c(length(cond),1)
            res=y-cond_X%*%t
    
            sumvec=rep(0,length(cond))
            for (i in 1:n){
                sumvec = sumvec + ( mest(res[i],e) - cm[i] )*
                            f_mest(res[i],e) * cond_X[i,]
            }
            return(sumvec)
        }
        
        obj=function(t,cond=1:p){ # tell S
            var=rep(0,p)
            var[cond]=t
            object=sum((mest(y-X%*%var,e)-cm)^2)
            return(object)
        }
        
        OP = function(cond){
            search = rep(0,length(cond))
            alpha = 0.5
            itr=0
            repeat{
                new_search = search + alpha * Fun(search,cond)
                itr =itr+1
                if (  itr > 10  ) {
                    break
                } else {
                    search=new_search
                }
                
            }
            return(search)
        }
        
    SA=1:p
    S=NULL
    beta=rep(0,p)
    EBIC=Inf
    repeat{  
        ###add new variable depend on R()  
        Z=mest(y-X%*%beta,e)
        R=-Inf
        for (i in SA){
            temp=t(X[,i])%*%(Z-cm)
            if(temp > R){
                R=temp
                temp_s=i
            }
        }
        
        
        S=c(S,temp_s)
        SA=SA[SA!=temp_s]
        ##Estimate beta based on recent variable selection
        beta[S]=OP(S)
        
        new_EBIC=n*log(obj(beta)/n)+length(S)*log(n)+2*(1-log(n)/(2.1*log(p)))*log(choose(p,length(S)))
        
        if(new_EBIC>EBIC | length(S)==p){
            break
        } else {
            EBIC=new_EBIC
        }
    }
    return(beta)
}