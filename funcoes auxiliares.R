###Periodograma###
periodograma=function(DATA,grafico="TRUE"){
  require(TSA)
  hh=periodogram(DATA, plot=FALSE);
  zz=1/hh$freq
  zz=cbind(zz,hh$spec)
  zz=zz[order(zz[,1]),]
  colnames(zz)=c("periodo","spec")
  if(grafico=="TRUE") plot(zz,type="l",ylab="espectro")
  return(zz)
}

###teste de fisher###
Fisher.test=function(P,alpha=0.05){
  ordenado=P[order(P[,2],decreasing="TRUE"),] 
  spec=ordenado[,2]
  g=max(spec)/(sum(spec))
  n=nrow(P)
  periodo=ordenado[1,1]
  zalpha=1-(alpha/n)^(1/(n-1))
  valorp=n*(1-g)^(n-1)
  a=cbind(g,zalpha,valorp,periodo)
  rownames(a)=""
  return(a)
}

###gma###
#y= time series, k= entre 8 a 12
gma=function(y,k){
  y1=matrix(y,nrow=k)
  media=apply(y1,2,mean)
  maximo=apply(y1,2,max)
  minimo=apply(y1,2,min)
  amplitude=maximo-minimo
  plot(media,amplitude,xlab="M?dias", ylab="Amplitudes")
  abline(lm(amplitude ~ media))
  a=cor.test(media,amplitude)
  a1=cbind(round(a$estimate,4),round(a$p.value,6))
  colnames(a1)=c("correlacao","valor-p")
  return(a1)
}

cs.test = function (x){
  method = "Cox-Stuart test for trend analysis"
  leng = length(x)
  apross = round(leng) %% 2
  if (apross == 1) {
    delete = (length(x)+1)/2
    x = x[ -delete ] 
  }
  half = length(x)/2
  x1 = x[1:half]
  x2 = x[(half+1):(length(x))]
  difference = x1-x2
  signs = sign(difference)
  signcorr = signs[signs != 0]
  pos = signs[signs>0]
  neg = signs[signs<0]
  if (length(pos) < length(neg)) {
    prop = pbinom(length(pos), length(signcorr), 0.5)
    names(prop) = "Increasing trend, p-value"
    rval <- list(method = method, statistic = prop)
    class(rval) = "htest"
    return(rval)
  }
  else {
    prop = pbinom(length(neg), length(signcorr), 0.5)
    names(prop) = "Decreasing trend, p-value"
    rval <- list(method = method, statistic = prop)
    class(rval) = "htest"
    return(rval)
  }
}


###Diagnostico do Modelo###
diagTAR=function(model){
  opar = par(mfrow = c(3, 1), mar = c(3, 4, 3, 2) + 0.1, oma = c(1, 
                                                                 0, 2, 0))
  e=model$residuals
  e=e/sqrt(var(e))
  normalidade=ks.test("pnorm",mean(e),sd(e))
  normalidade$p.value
  plot(e,ylab="Residuos Padronizados")
  mtext(paste("TAR model with ML = ",ML," and MH = ",MH,sep=""), outer = TRUE, cex = 1)
  n=length(e)
  lag.max=10*log10(n)
  acf=acf(e,lag.max=lag.max,main="")
  lbv=rep(0,lag.max)
  for (i in 1:lag.max){
    lbv[i]=LB.test(model,i)$p.value
  }
  plot(lbv,x = 1:lag.max,ylim = c(0, 1), pch = 21,ylab="p-value",xlab="lag") 
  abline(h = 0.05, lty = 2, col = 2)
  return(normalidade)
}

###Grafico de Previs?o
predictTAR=function(model){
  yprev=model$fitted.values
  y1=model$str$x
  n=length(y1)
  t=seq(1,n)
  plot(t,y1,ylab="Serie",type="l")
  n1=length(yprev)
  yprev=c(y1[1:(n-n1)],yprev)
  lines(yprev,col=2)
}

##setar - fun??o para estimar o modelo TAR auto threshold
##ML - vetor contendo quais sao os lag do regime 1
##MH - vetor contendo quais sao os lag do regime 2
##th - valor do threshold - se th n?o for informado ser? estimado um valor de threshold
## thDelay-  o valor do retador na varialvel threshold

select=function(ML=1:31,MH=1:31,d=1,y,x,criterio="AIC"){
  lags=list(NULL)
  criterios=(NULL)
  ww=1
  n1=length(ML)
  n2=length(ML)
  ML1=NULL
  for(i in 1:n1){
    m1=combn(ML,i)
    m1=t(m1)
    nw1=nrow(m1)
    for(k in 1:nw1){
      for(j in 1:n2){
        m2=combn(MH,j)
        m2=t(m2)  
        nw2=nrow(m2)
        for(l in 1:nw2){
          model=setar(y,ML=m1[k,], MH=m2[l,],thVar=x,thDelay=d)
          pp=summary(model)$coef
          k1=length(which(pp[,4]<0.05))
          gc(reset=TRUE)
          if(model$k==(k1+1)){
            regimes=list("ML"=m1[k,],"MH"=m2[l,])
            lags[[ww]]=regimes
            criterios=rbind(criterios,c(AIC(model),MAPE(model)))
            ww=ww+1
          }
        }
      }
    }
  }
  colnames(criterios)=c("AIC","MAPE")
  return(list(LAGS=lags,criterios=criterios))
}

gdados=function(dados,p){
  n=length(dados)
  n1=n-p
  entrada=matrix(0,n1,p)
  saida=dados[(p+1):n]
  for(i in 1:p){
    nw=n1+i-1
    h=c(i:nw)
    entrada[,i]=dados[h]
  }
  return(list(entrada=entrada,saida=saida))
}
