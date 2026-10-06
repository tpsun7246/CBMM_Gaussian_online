##install package
library(coda)
library(psych)
library(quantmod)
library(truncnorm)
library(dplyr)
##copula function
copula=function(a,b,c){
  a[a<=0.0001]=0.0001
  a[a>=0.9999]=0.9999
  b[b<=0.0001]=0.0001
  b[b>=0.9999]=0.9999
  d=log(1+c)+(-1-c)*log(a*b)+((-1/c)-2)*log(a^(-c)+b^(-c)-1)  
  return(d)
}

##calculate the predictive probability
mcc=function(dd,prm.m,prm.s,prs.a,prs.b){
  a=rnorm(n = 10000,mean = prm.m,sd = prm.s)
  b=rgamma(n = 10000,shape = prs.a,rate = prs.b)
  w=dnorm(x = dd,mean = a,sd = sqrt(b))
  p=mean(w[is.finite(w)])
  return(p)
}

posterio_pred.m.s <- function(n,data,alp){
  p=c()
  mm=c()
  ms=c()
  ##compute initial value(integral)
  p[1]=lapply(c(1:100), function(x){
    mcc(dd = data[n],prm.m=prm.m,prm.s=prm.s,prs.a=prs.a,prs.b=prs.b)
  }) %>% unlist() %>% median()
  ##create two vectors for accept rate
  
  for (k in 1:(n-1)) {
    if(k==1){
      ##create four vectors for accept rate
      lam.m=c()
      lam.s=c()
      lam.m[1]=40
      lam.s[1]=20
      
      
      for (j in 1:nMCMC) {
        pp=rnorm(n = 1,mean =lam.m[j] ,sd=pp.sd.m)
        a1=dnorm(x = data[n-1],mean = pp,sd = lam.s[j],log = T)+dnorm(x = pp,mean = prm.m,sd = prm.s,log = T)-dnorm(x = data[n-1],mean = lam.m[j],sd = lam.s[j],log = T)-dnorm(x = lam.m[j],mean = prm.m,sd = prm.s,log = T)
        u=runif(1)
        if(log(u)<min(a1,0)){
          lam.m[j+1]=pp
          
          
        }else{
          lam.m[j+1]=lam.m[j]
          
        }
        p1=rtruncnorm(n = 1,a = 0,b = Inf,mean = lam.s[j],sd = pp.sd.s)
        
        a2=dnorm(x = data[n-1],mean = lam.m[j+1],sd = sqrt(p1),log = T)+dgamma(x = p1,shape = prs.a,rate = prs.b,log = T)+log(dtruncnorm(x = lam.s[j],a = 0,b = Inf,mean = p1,sd = pp.sd.s))-dnorm(x = data[n-1],mean = lam.m[j+1],sd = sqrt(lam.s[j]),log = T)-dgamma(x = lam.s[j],shape = prs.a,rate = prs.b,log = T)-log(dtruncnorm(x = p1,a = 0,b = Inf,mean = lam.s[j],sd = pp.sd.s))
        u=runif(1)
        if(log(u)<min(a2,0)){
          lam.s[j+1]=p1
          
        }else{
          lam.s[j+1]=lam.s[j]
          
        }
      }
    }else{
      lam.m=c()
      lam.s=c()
      lam.m[1]=40
      lam.s[1]=20
      for (j in 1:nMCMC) {
        pp=rnorm(n = 1,mean =lam.m[j] ,sd=pp.sd.m)
        a1=sum(dnorm(data[(n-1):(n-k)],pp,sqrt(lam.s[j]),log = T))+sum(copula(pnorm(data[(n-1):(n-k+1)],pp,sqrt(lam.s[j])),pnorm(data[(n-2):(n-k)],pp,sqrt(lam.s[j])),alp))+dnorm(x = pp,mean = prm.m,sd = prm.s,log = T)
        a2=sum(dnorm(data[(n-1):(n-k)],lam.m[j],sqrt(lam.s[j]),log = T))+sum(copula(pnorm(data[(n-1):(n-k+1)],lam.m[j],sqrt(lam.s[j])),pnorm(data[(n-2):(n-k)],lam.m[j],sqrt(lam.s[j])),alp))+dnorm(x = lam.m[j],mean = prm.m,sd = prm.s,log = T)
        u=runif(1)
        if(log(u)<min(a1-a2,0)){
          lam.m[j+1]=pp
          
          
        }else{
          lam.m[j+1]=lam.m[j]
          
          
        }
        p1=rtruncnorm(n = 1,a = 0,b = Inf,mean = lam.s[j],sd = pp.sd.s)
        
        a3=sum(dnorm(data[(n-1):(n-k)],lam.m[j+1],sqrt(p1),log = T))+sum(copula(pnorm(data[(n-1):(n-k+1)],lam.m[j+1],sqrt(p1)),pnorm(data[(n-2):(n-k)],lam.m[j+1],sqrt(p1)),alp))+dgamma(x = p1,shape = prs.a,rate = prs.b,log = T)+log(dtruncnorm(x = lam.s[j],a = 0,b = Inf,mean = p1,sd = pp.sd.s))
        a4=sum(dnorm(data[(n-1):(n-k)],lam.m[j+1],sqrt(lam.s[j]),log = T))+sum(copula(pnorm(data[(n-1):(n-k+1)],lam.m[j+1],sqrt(lam.s[j])),pnorm(data[(n-2):(n-k)],lam.m[j+1],sqrt(lam.s[j])),alp))+dgamma(x = lam.s[j],shape = prs.a,rate = prs.b,log = T)+log(dtruncnorm(x = p1,a = 0,b = Inf,mean = lam.s[j],sd = pp.sd.s))
        u=runif(1)
        if(log(u)<min(a3-a4,0)){
          lam.s[j+1]=p1
          
          
        }else{
          lam.s[j+1]=lam.s[j]
          
          
        }
      }
    }
    u2=dnorm(data[n],lam.m[-c(1:burn.in)],sqrt(lam.s[-c(1:burn.in)]))*exp(copula(a=pnorm(data[n],lam.m[-c(1:burn.in)],sqrt(lam.s[-c(1:burn.in)])),b=pnorm(data[n-1],lam.m[-c(1:burn.in)],sqrt(lam.s[-c(1:burn.in)])),alp))
    p[k+1]=mean(u2)
    mm[k]=pnorm(q =abs(unname(geweke.diag(x = lam.m[-c(1:burn.in)])$z)),lower.tail = F)*2
    ms[k]=pnorm(q =abs(unname(geweke.diag(x = lam.s[-c(1:burn.in)])$z)),lower.tail = F)*2
  }
  return(list(p=p,pos.mean=mm,pos.var=ms))
}

getSymbols("^GSPC",from="2007-07-13",to="2008-07-15")
a.log.return=diff(log(GSPC$GSPC.Close))[-1]
length(a.log.return)
i.a=index(a.log.return)
d1=as.numeric(a.log.return$GSPC.Close)*100
getSymbols("^GSPC",from="2008-07-16",to="2008-11-15")
aa.log.return=diff(log(GSPC$GSPC.Close))[-1]
dd=as.numeric(aa.log.return$GSPC.Close)*100
ii.a=index(aa.log.return)
##algorithm
cpu.cores=30    
H=1/1260
ell=10
##prior and proposal setting
pp.sd.m=sqrt(10)
pp.sd.s=sqrt(10)
prm.m=0
prm.s=sqrt(10)
prs.a=1
prs.b=10^-2
nMCMC = 15000
burn.in=3000
##estimate alpha
likelihood=function(x){
  dd=d1
  sum(dnorm(dd[-length(dd)],mean = mean(dd[-length(dd)]),sd=sd(dd[-length(dd)]),log = T))+sum(copula(pnorm(dd[-length(dd)],mean(dd[-length(dd)]),sd(dd[-length(dd)])),pnorm(dd[-1],mean(dd[-1]),sd(dd[-1])),x))
}
alpha=optimize(f = likelihood,interval = c(0,20),maximum = T)$maximum
iteration=length(dd)
f=function(j){
  a=posterio_pred.m.s(data=dd[1:j],n=j,alp = alpha)
  return(a)
}
##parallel computation
library(doSNOW)
library(snow)
library(parallel)
cl <- makeCluster(cpu.cores)  
registerDoSNOW(cl)   
pb <- txtProgressBar(max = iteration, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
t1<-Sys.time()
pis <- foreach(i = 2:iteration, .packages = c("snow","coda","dplyr","truncnorm"), .options.snow = opts) %dopar%{
  s <- f(i)
}   
close(pb)
stopCluster(cl)   
t2<-Sys.time()
pp=t2-t1
pp
##compute the proportion of MCMC coveragence p-value
mm=matrix(0,1,2)
colnames(mm)=c("mean","var")
pos.mean=lapply(pis, function(x){
  x[[2]]
}) %>% unlist()
mm[,1]=mean(pos.mean>0.05)
pos.var=lapply(pis, function(x){
  x[[3]]
}) %>% unlist()
mm[,2]=mean(pos.var>0.05)
mm
write.csv(x = mm,file = "converge.csv",row.names = F)
###lag algorithm
data=dd
R=list(ell+1)
for (i in 1:(ell+1)){
  R[[i]] = matrix(0, length(data), length(data))
  R[[i]][1,1] = 1
}
for (t in 2:length(data)) {
  R[[1]][t,2:t]= R[[1]][t-1,1:(t-1)]*pis[[t-1]]$p[-1]*(1-H)
  R[[1]][t,1] = sum(R[[1]][t-1,1:(t-1)])*pis[[t-1]]$p[1]*H
  R[[1]][t,1:t]=R[[1]][t,1:t]/sum(R[[1]][t,1:t])
  for (j in 1:ell){
    # Calculate run length
    if (t>=j+1){ # LEXO-1 begins from t=2, LEXO-2 begins from t=3,...
      R[[j+1]][t-j,1:(t-j)]=R[[j]][t-j+1,2:(t-j+1)]+R[[j]][t-j+1,1]*R[[1]][t-j,1:(t-j)]
      
    }
  }
}
##copula runlength
co=NULL
for (i in 1:(ell+1)) {
  co=cbind(co,apply(R[[i]],1,which.max))
}
rownames(co)=as.character(ii.a)
write.csv(x = co,file = "GSPC.co.csv")

cha=which(diff(as.numeric(co[,1]))<0)+1
cha.ii=ii.a[cha]
plot(1:nrow(co), co[,1],xlab="Date",ylab="Run Length",type = "l",lwd=2,xaxt="n",main="The run length plot of S&P 500")
axis(side = 1,at = seq(1,length(ii.a),by=20),labels = ii.a[seq(1,length(ii.a),by=20)])
abline(v=cha,col="red",lwd=2,lty=2)
legend("topleft",legend = cha.ii,col="red",lwd=2,lty=2,bty="n",title = "changepoint position",cex=1.1)



plot(x=1:length(ii.a),y=dd,type="l",ylab="log return*100",xaxt="n",xlab="Date",main="The log return of S&P 500 in testing period",lwd=2)
axis(side = 1,at = seq(1,length(ii.a),by=20),labels = ii.a[seq(1,length(ii.a),by=20)])
abline(v=which(ii.a %in% cha.ii),col=c("red","blue","purple3"),lwd=2,lty=2)
legend("topleft",legend = "changepoint",col="red",lwd=2,lty=2,bty="n",cex=1.2)


