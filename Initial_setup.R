## Gaussian Distribution Generator with Given Graph/Covariance

if ( Precision_Structur = "In-Build") {
    E0=huge.generator(n,q,graph="random",vis=TRUE,verbose = FALSE)$data # In-build function generator
} else {
    if ( Precision_Structure = "AR1") { 
        Omega = diag(1,q)+rbind(rep(0,q),cbind(diag(0.45,q-1),rep(0,q-1)))+cbind(rep(0,q),rbind(diag(0.4,q-1),rep(0,q-1)))
    } else if (Precision_Structure = "AR2"){
        Omega = diag(1,q)+rbind(rep(0,q),cbind(diag(0.3,q-1),rep(0,q-1)))+cbind(rep(0,q),rbind(diag(0.4,q-1),rep(0,q-1)))+
                rbind(matrix(0,2,q),cbind(diag(0.4,q-2),matrix(0,q-2,2)))+cbind(matrix(0,2,q),rbind(diag(0.4,q-2),matrix(0,q-2,2)))
    } else if (Precision_Structure = "Circle"){    
        Omega =  diag(1,q)+rbind(rep(0,q),cbind(diag(0.45,q-1),rep(0,q-1)))+cbind(rep(0,q),rbind(diag(0.4,q-1),rep(0,q-1)))
        Omega[1,q] = 0.1
        Omega[q,1] = 0.1
    } else if (Precision_Structure = "Cluster"){    
        block = matrix(0.5,5,5)+diag(0.5,5)
        Omega =  diag(block,q/5)
    } else if (Precision_Structure = "ER"){    
        Omega =  
    } else if (Precision_Structure = "Tridiag"){    
        Omega =  
    } else if (Precision_Structure = "Hub"){    
        Omega =  
    } else if (Precision_Structure = "RP"){    
        Omega =  
    } else if (Precision_Structure = "BA"){    
        Omega =  
    }
    E0 = mvrnorm(1,rep(0,q),solve(Omega))
}    
    
    
    
## Transform to NPN Family
E1=npn_gaussian.generator(E0,mu_v=rep(0,q),sigma_v=rep(1,q))
E1=npn_power.generator(E0,mu_v=rep(0,q),sigma_v=rep(1,q))


## Coefficient Matrix Generator 
prob=floor(4*n^0.16)/p
B = matrix(abs(rnorm(p*q,0,0.1/qnorm(0.875)))+4*n^(-0.15),nrow=p,ncol=q)*
    matrix(sample(c(-1,1),p*q,replace=TRUE,prob=c(0.4,0.6)),nrow=p,ncol=q)*
    matrix(sample(0:1,p*q,replace=TRUE,prob=c(1-prob,prob)),nrow=p,ncol=q)


## Feature Generator with Given Structure
X=matrix(rnorm(n*p),nrow=n,ncol=p) #design different X structure


## Calculate the Response
Y=X%*%B+E1