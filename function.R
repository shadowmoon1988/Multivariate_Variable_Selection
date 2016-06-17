#function part
##Via Gaussian transform to generate Nonparanormal
npn_gaussian.generator=function(data,mu=0.05,sigma=0.4,mu_v,sigma_v){
    ###inner use funciton for npn gaussian transform
    f1 = function(t,m,s,m0,s0){pnorm((t-m0)/s0)*dnorm((t-m)/s)}  
    f2 = function(t,m,s,m0,s0){(pnorm((t-m0)/s0)-integrate(f1,-Inf,Inf,m=m,s=s,m0=m0,s0=s0)[[1]])^2*dnorm((t-m)/s)}
    
    k=dim(data)[2]
    output=matrix(0,dim(data)[1],k)
    for (i in 1:k){
        numerator= pnorm((data[,i]-mu)/sigma)-  integrate(f1,-Inf,Inf,m=mu_v[i],s=sigma_v[i],m0=mu,s0=sigma)[[1]]
        dominator= sqrt(integrate(f2,-Inf,Inf,m=mu_v[i],s=sigma_v[i],m0=mu,s0=sigma)[[1]])
        output[,i]=numerator/dominator
    }
    return(output)
}

##Via Power transform to generate Nonparanormal
npn_power.generator = function(data,alpha=3,mu_v,sigma_v){
    ###inner use funciton for npn power transform
    f1 = function(t,m,s){abs(t-m)^(2*alpha)*dnorm((t-m)/s)}
    
    k=dim(data)[2]
    output=matrix(0,dim(data)[1],k)
    for (i in 1:k){
        numerator=sign(data[,i]-mu_v[i]) * abs(data[,i]-mu_v[i])^alpha
        dominator=sqrt(integrate(f1,-Inf,Inf,m=mu_v[i],s=sigma_v[i])[[1]])
        output[,i]=numerator/dominator
    }
    return(output)
}

## Marginal Estimation
# mest=function(input,data,thershold=1/(2*n),smooth=FALSE){  
#     temp_fun = function(var){
#         temp=sum(ifelse(data<=var,1,0))/n
#         if (smooth & var>=min(data) & var<max(data) ){
#             lower=max(data[data<=var])
#             upper=min(data[data>var])
#             temp=temp+(var-lower)/(n*(upper-lower))
#         }
#         if (temp<thershold) {temp=thershold}
#         if (temp >1-thershold){temp=1-thershold}
#         return(qnorm(temp))
#     }
#     
#     return(sapply(input,temp_fun))
# }

mest=function(input,data,bandwidth=n^(-0.2)){  
    temp_fun = function(var){
        temp=sum(pnorm((var-data)/bandwidth))/n
        thershold=1/2/n
        if (temp<thershold) {temp=thershold}
        if (temp >1-thershold){temp=1-thershold}
        return(qnorm(temp))
    }
    
    return(sapply(input,temp_fun))
}

## 1st Derivative of Marginal Estimation
f_mest = function(var,data,bandwidth=n^(-0.2)){
    output=sum(dnorm((var-data)/bandwidth))/(n*bandwidth) /dnorm(mest(var,data))
    return(output)
}    
## 2nd Derivative of Marginal Estimation
s_mest = function(var,data,bandwidth=n^(-0.2)){
    temp=(var-data)/bandwidth
    output = (sum(dnorm(temp))*mest(var,data)*f_mest(var,data) -sum(dnorm(temp)*temp/bandwidth))/
        dnorm(mest(var,data))
    return(output)
}    

## Evaluation 
evaluation = function(M1,M2){
    nrow=dim(M1)[1]
    ncol=dim(M1)[2]
    FP=0
    FN=0
    r1=0
    r2=0
    
    for ( i in 1:nrow){
        for ( j in 1:ncol){
            if (M2[i,j] == 0){
                if (M1[i,j]!=0) {
                    FP=FP+1
                    r1=r1+1
                }
            }
            else {
                r2=r2+1
                if (M1[i,j]==0) {
                    FN=FN+1
                } else {
                    r1=r1+1
                }
            }
        }
    }
    output=c(FP/r1,1-FN/r2)
    names(output)=c("FDR","PDR")
    return(output)
}


stoprule = function (M1, M2, threshold=c(0.95,0.95)){
    B1 = ifelse( M1 !=0 , 1, 0)
    B2 = ifelse( M2 !=0 , 1, 0)
    r = sum (B1 & B2) /c(sum (B1), sum (B2))
    result = all(r > threshold)
    return(result)
}