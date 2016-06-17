ExLasso = function(X,Y,t=1){
    X_tilde= diag(1,q)%x%X
    Q=t(X_tilde) %*% X_tilde
    y_tilde = NULL
    for (i in 1:q){
        y_tilde = c(y_tilde,Y[,i])
    }
    c=-t(X_tilde) %*% y_tilde

    p=dim(X)[2]
    q=dim(Y)[2]

    b=rep(0,p*q)
    z=rep(0,p)
    lambda_u=rep(100,p*q)
    lambda_l=rep(100,p*q)
    tau=100
    su=rep(100,p*q)
    sl=rep(100,p*q)
    kesi=100
    start=c(b,z,lambda_u,lambda_l,tau,su,sl,kesi)
    
    repeat{
        MAT=as(rbind(
          cbind(Q,matrix(0,p*q,p),diag(1,p*q),-diag(1,p*q),matrix(0,p*q,2*p*q+2)),
          cbind(matrix(0,p,p*q+p),-(matrix(1,1,q) %x% diag(1,p)),-(matrix(1,1,q) %x% diag(1,p)),matrix(1,p,1),matrix(0,p,2*p*q+1)),
          cbind(diag(-1,p*q),matrix(1,q,1)%x% diag(1,p),matrix(0,p*q,2*p*q+1),diag(-1,p*q),matrix(0,p*q,p*q+1)),
          cbind(diag(1,p*q),matrix(1,q,1)%x% diag(1,p),matrix(0,p*q,3*p*q+1),diag(-1,p*q),matrix(0,p*q,1)),
          cbind(matrix(0,1,p*q),matrix(-1,1,p),matrix(0,1,4*p*q+1),matrix(-1,1,1)),
          cbind(matrix(0,p*q,p*q+p),diag(su),matrix(0,p*q,p*q+1),diag(lambda_u),matrix(0,p*q,p*q+1)),
          cbind(matrix(0,p*q,2*p*q+p),diag(sl),matrix(0,p*q,p*q+1),diag(lambda_l),matrix(0,p*q,1)),
          cbind(matrix(0,1,3*p*q+p),matrix(kesi,1,1),matrix(0,1,2*p*q),matrix(tau,1,1))
        ),"sparseMatrix")
       

        mu=(sum(lambda_l*sl)+sum(lambda_u*su)+tau*kesi)/(2*p*q+1)
        if (mu < 10^(-6)) {
          break
        }
        
        residual_b= Q %*% b + c +lambda_u -lambda_l 
        residual_z= -(matrix(1,1,q) %x% diag(1,p)) %*% lambda_u -(matrix(1,1,q) %x% diag(1,p)) %*% lambda_l + matrix(tau,p,1)
        residual_u= matrix(1,q,1) %x% z -b -su
        residual_l= matrix(1,q,1) %x% z +b -sl
        residual_t= t-sum(z)-kesi
        
        # Caculate affine scaling direction
        residual_hi=lambda_u*su
        residual_lo=lambda_l*sl
        residual_tk=tau*kesi
        
        residual=c(residual_b,residual_z,residual_u,residual_l,residual_t,residual_hi,residual_lo,residual_tk)
        direction = as.matrix(solve(MAT, -residual))  
        
        alpha_aff = min(ifelse(direction<0,-start/direction,1)[-(1:(p*q+p))])
        temp_aff=start+alpha_aff*direction
        
        #Use the affine scaling direction
        mu_aff=(sum(temp_aff[(p*q+p+1):(2*p*q+p)]*temp_aff[(3*p*q+p+2):(4*p*q+p+1)])+
                sum(temp_aff[(2*p*q+p+1):(3*p*q+p)]*temp_aff[(4*p*q+p+2):(5*p*q+p+1)])+
                temp_aff[3*p*q+p+1]*temp_aff[5*p*q+p+2])/(2*p*q+1)   
        sigma= (mu_aff/mu)^3
        
        hat_rhi=direction[(p*q+p+1):(2*p*q+p)]*direction[(3*p*q+p+2):(4*p*q+p+1)]
        hat_rlo=direction[(2*p*q+p+1):(3*p*q+p)]*direction[(4*p*q+p+2):(5*p*q+p+1)]
        hat_rtk=direction[3*p*q+p+1]*direction[5*p*q+p+2]
        
        residual_hi=residual_hi-sigma * mu * matrix(1,p*q,1) + hat_rhi
        residual_lo=residual_lo-sigma * mu * matrix(1,p*q,1) + hat_rlo
        residual_tk=tau*kesi-sigma*mu+hat_rtk
        
        residual=c(residual_b,residual_z,residual_u,residual_l,residual_t,residual_hi,residual_lo,residual_tk)
        direction = as.matrix(solve(MAT, -residual))   
        
        alpha = min(ifelse(direction<0,-start/direction,1)[-(1:(p*q+p))])
        temp=start+alpha*direction
        start=temp
        #print(min(temp[-(1:(p*q+p))]))
        b=temp[1:(p*q)]
        z=temp[(p*q+1):(p*q+p)]
        lambda_u=temp[(p*q+p+1):(2*p*q+p)]
        lambda_l=temp[(2*p*q+p+1):(3*p*q+p)]
        tau=temp[3*p*q+p+1]
        su=temp[(3*p*q+p+2):(4*p*q+p+1)]
        sl=temp[(4*p*q+p+2):(5*p*q+p+1)]
        kesi=temp[5*p*q+p+2]
    }
    return(matrix(ifelse(b>t*10^(-4),b,0),p,q))
}