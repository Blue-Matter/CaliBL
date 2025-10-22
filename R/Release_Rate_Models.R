# Release rate models


#CRvec = seq(0,3.2,length.out=20)
#L = CRvec-(max(CRvec)/2) * 0.1
#VM = exp(L)/(1+exp(L))


# RR predictors ---------------------------------------------------------------------------

RRpois = function(BL,CRvec,n=1E4,CRlim=30,V=0){
  if(length(BL)==1)BL = rep(BL,length(CRvec))
  RRs = rep(NA,length(CRvec))
  crs = 0:CRlim
  for(i in 1:length(CRvec)){
    ENC = dpois(crs,CRvec[i])*crs         # expected number caught
    ERV = ENC * V                         # expected voluntary release
    ERB = dpois(crs,CRvec[i])*(crs-BL[i]) * (1-V)
    RRs[i] = (sum(ERV)+sum(ERB[crs>=BL[i]]))/sum(ENC)
  }
  RRs
}

RRlnV = function(BL,CRvec,CV,V=0,n=1E4,CRlim=c(1E-5,30)){
  if(length(BL)==1)BL = rep(BL,length(CRvec))
  RRs = rep(NA,length(CRvec))
  crs = seq(CRlim[1],CRlim[2],length.out=n)
  for(i in 1:length(CRvec)){
    vals0 = dlnorm(crs,log(CRvec[i]),CV)*crs
    valsV = vals0 * (1-V)
    RRs[i] = 1 - sum(valsV[crs<BL[i]])/sum(vals0)
    ENC = dlnorm(crs,log(CRvec[i]),CV)*crs   # expected number caught
    ERV = ENC * V                            # expected voluntary release
    ERB = dlnorm(crs,log(CRvec[i]),CV)*(crs-BL[i]) * (1-V)
    RRs[i] = (sum(ERV)+sum(ERB[crs>=BL[i]]))/sum(ENC)
  }
  RRs
}

RRln2V = function(BL,CRvec,CV,V1=0,V2=0,n=1E4,CRlim=c(1E-5,30)){
  if(length(BL)==1)BL = rep(BL,length(CRvec))
  RRs = rep(NA,length(CRvec))
  crs = seq(CRlim[1],CRlim[2],length.out=n)
  Vs = rep(V1,n)
  Vs[crs >1]=V2
  for(i in 1:length(CRvec)){
    vals0 = dlnorm(crs,log(CRvec[i]),CV)*crs
    valsV = vals0 * (1-Vs)
    RRs[i] = 1 - sum(valsV[crs<BL[i]])/sum(vals0)
    ENC = dlnorm(crs,log(CRvec[i]),CV)*crs   # expected number caught
    ERV = ENC * Vs                            # expected voluntary release
    ERB = dlnorm(crs,log(CRvec[i]),CV)*(crs-BL[i]) * (1-Vs)
    RRs[i] = (sum(ERV)+sum(ERB[crs>=BL[i]]))/sum(ENC)
  }
  RRs
}


RRlnS = function(BL,CRvec,CV,alpha=0,beta=-4,n=1E4,CRlim=c(1E-5,30)){
  if(length(BL)==1)BL = rep(BL,length(CRvec))
  RRs = rep(NA,length(CRvec))
  crs = seq(CRlim[1],CRlim[2],length.out=n)

  Vs = ilogit(beta+alpha*crs)
  for(i in 1:length(CRvec)){
    vals0 = dlnorm(crs,log(CRvec[i]),CV)*crs
    valsV = vals0 * (1-Vs)
    RRs[i] = 1 - sum(valsV[crs<BL[i]])/sum(vals0)
    ENC = dlnorm(crs,log(CRvec[i]),CV)*crs   # expected number caught
    ERV = ENC * Vs                            # expected voluntary release
    ERB = dlnorm(crs,log(CRvec[i]),CV)*(crs-BL[i]) * (1-Vs)
    RRs[i] = (sum(ERV)+sum(ERB[crs>=BL[i]]))/sum(ENC)
  }
  RRs
}


# Fitting functions

logit = function(p)log(p/(1-p))
ilogit = function(x)exp(x)/(1+exp(x))
cleandat = function(dat){
  if(ncol(dat)==2){# if no time varying bag limit (column 3)
    dat1 = dat[dat$RR > 0 & dat$RR<1 & dat$CR!='Inf',]
  }else{
    cond = dat$RR > 0 & dat$RR<1 & dat$CR!='Inf' & !is.na(dat$BG)
    cond[is.na(cond)] = FALSE
    dat1 = dat[cond,]
  }
  dat1
}

minmax = function(vec,minv = 0.1, maxv = 1){
  vec[vec<minv] = minv
  vec[vec>maxv] = maxv
  vec
}

plotfit = function(dat,CRpred,RRpred, point=F){
  plot(dat$CR,dat$RR,pch=19,col="#99999990",ylim=c(0,0.8),xlim=c(0,quantile(dat$CR,0.98))); grid()
  if(!point){
    lines(CRpred,RRpred,col="#ff000090",lwd=2)
  }else{
    points(CRpred,RRpred,col="#ff000090",pch=19)
  }
}

fit_pois = function(dat, BG, OE = 0.2, plot=F){
  dat=cleandat(dat)
  RR = dat$RR
  CR = dat$CR
  OE = minmax(dat$varp)
  logitV=logit(lm(RR~CR,dat=dat)$coefficients[1])

  pois_int = function(logitV,dat, BG,OE){
    pred= pois_RR(dat$CR,BagLim=BG,ilogit(logitV))
    -sum(dnorm(log(dat$RR),log(pred),OE),log=T)
  }

  fit = optimize(pois_int,logit(c(0.0001,0.9999)),dat=dat, BG=BG, OE = OE, hessian=T)
  CRpred = seq(min(dat$CR),max(dat$CR),length.out=100)
  RRpred = pois_RR(CRpred, BagLim=BG, ilogit(fit$minimum))
  legtext1 = paste("BL =",BG); if(length(BG==1)){legtext1=""}
  if(plot){
    plotfit(dat,CRpred,RRpred)
    legend('bottomright',legend = c(legtext1,
                                    paste("V pois =",round(ilogit(fit$minimum),3)),
                                    paste("O.F. =", round(fit$objective,2))),bty='n',text.col="red")
  }
  list(fit=fit, CRpred = CRpred, RRpred=RRpred, CRobs = dat$CR, RRobs = dat$RR, BG = BG)
}


# RRlnV = function(BL,CRvec,CV,V=0,n=1E4,CRlim=c(1E-5,30)){

LN_pred= function(BG,par,dat,type="LNV"){

  if(type == "LNV"){
    pred = RRlnV(BG, dat$CR, CV = exp(par[1]), V=ilogit(par[2]))
  }else if(type == "LN2V"){
    pred = RRln2V(BG,dat$CR, CV = exp(par[3]), V1=ilogit(par[1]),V2=ilogit(par[2]))
  }else if(type == "LNVS"){
    pred = RRlnS(BG, dat$CR, CV = exp(par[2]), alpha=exp(par[1]))
  }else{
    pred = RRlnS(BG, dat$CR, CV = exp(par[2]), alpha=exp(par[1]),beta = par[3])
  }
  pred
}

fit_LNV = function(dat, BG, OE = 0.4, plot=F){
  dat=cleandat(dat)
  tv = "BG" %in% names(dat)
  if(tv) BG = dat$BG
  RR = dat$RR
  CR = dat$CR
  if(tv) OE = minmax(sqrt(dat$varp))

  LNV_int = function(par, dat, BG, OE){
    pred = RRlnV(BG, dat$CR, CV = exp(par[1]), V=ilogit(par[2]))
    -sum(dnorm(logit(dat$RR),logit(pred),OE),log=T)
  }

  fit = optim(c(log(0.5),logit(0.2)), LNV_int, method = "L-BFGS-B",
              lower = c(log(0.1),logit(0.001)),
              upper = c(log(3),logit(0.999)),dat=dat, BG=BG, OE=OE, hessian=T)

  if(!tv){
    CRpred = seq(min(dat$CR),max(dat$CR),length.out=100)
    RRpred = RRlnV(BL=BG, CRpred, CV = exp(fit$par[1]), V=ilogit(fit$par[2]))
  }
  if(tv){
    CRpred = dat$CR
    RRpred = sapply(1:length(BG),function(X){RRlnV(BL = BG[X],CRpred[X],CV = exp(fit$par[1]), V=ilogit(fit$par[2]))})
  }

  legtext1 = paste("BL =",BG); if(length(BG==1)){legtext1=""}
  if(plot){
    plotfit(dat,CRpred,RRpred,point=tv)
    legend('bottomright',legend = c(legtext1,
                                    paste("V0 =",round(ilogit(fit$par[2]),3)),
                                    paste("CV =",round(exp(fit$par[1]),3)),
                                    paste("O.F. =", round(fit$value,2))),
           bty='n',text.col="red")
  }
  list(fit=fit, CRpred = CRpred, RRpred=RRpred, CRobs = dat$CR, RRobs = dat$RR, BG = BG)
}

#RRln2V = function(BL,CRvec,CV,V1=0,V2=0,n=1E4,CRlim=c(1E-5,30)){
fit_LN2V = function(dat, BG, OE = 0.4, plot=F){

  dat=cleandat(dat)
  tv = "BG" %in% names(dat)
  if(tv)OE = minmax(sqrt(dat$varp))
  if(tv) BG = dat$BG
  #initfit = fit_lnorm1(dat,BG,plot=F)$fit

  LN2V_int = function(par, dat, BG, OE){
    pred = RRln2V(BG,dat$CR, CV = exp(par[3]), V1=ilogit(par[1]),V2=ilogit(par[2]))
    -sum(dnorm(logit(dat$RR),logit(pred),OE),log=T)
  }

  fit = optim(c(log(0.1),-2,log(0.55)),LN2V_int, method = "L-BFGS-B",
              lower = c(logit(0.001),logit(0.001),log(0.25)),
              upper = c(logit(0.5),logit(0.8),log(3)),dat=dat, BG=BG,OE=OE,hessian=T)


  if(!tv){
    CRpred = seq(min(dat$CR),max(dat$CR),length.out=100)
    RRpred = RRln2V(BG, CRpred, CV = exp(fit$par[3]), V1=ilogit(fit$par[1]),V2=ilogit(fit$par[2]))
  }
  if(tv){
    CRpred = dat$CR
    RRpred = sapply(1:length(BG),function(X){RRln2V(BG[X],CRpred[X],CV = exp(fit$par[3]), V1=ilogit(fit$par[1]),V2=ilogit(fit$par[2]))})
  }

  legtext1 = paste("BL =",BG); if(length(BG==1)){legtext1=""}

  if(plot){
    plotfit(dat,CRpred,RRpred,point=tv)
    legend('bottomright',legend = c(legtext1,
                                    paste("V1 =",round(ilogit(fit$par[1]),3)),
                                    paste("V2 =",round(ilogit(fit$par[2]),3)),
                                    paste("CV =",round(exp(fit$par[3]),3)),
                                    paste("O.F. =",round(fit$value,2))),
           bty='n',text.col="blue")
  }

  list(fit=fit, CRpred = CRpred, RRpred=RRpred, CRobs = dat$CR, RRobs = dat$RR, BG = BG)
}



# RRlnS = function(BL,CRvec,CV,alpha=0,beta=-4,n=1E4,CRlim=c(1E-5,30)){
fit_LNVS = function(dat, BG, OE = 0.4, plot=F){

  dat=cleandat(dat)

  tv = "BG" %in% names(dat)
  if(tv) BG = dat$BG
  if(tv) OE = minmax(sqrt(dat$varp))
  LNVS_int = function(par, dat, BG, OE){
    pred = RRlnS(BG, dat$CR, CV = exp(par[2]), alpha=exp(par[1]))
    -sum(dnorm(logit(dat$RR),logit(pred),OE),log=T)
  }

  fit = optim(par=c(-2,log(0.5)),LNVS_int, method = "L-BFGS-B",
              lower = c(-5,log(0.25)),
              upper = c(3,log(4)),dat=dat, BG=BG, OE=OE,hessian=T)

  if(!tv){
    CRpred = seq(min(dat$CR),max(dat$CR),length.out=100)
    RRpred = RRlnS(BG, CRpred,CV = exp(fit$par[2]), alpha=exp(fit$par[1]))
  }
  if(tv){
    CRpred = dat$CR
    RRpred = sapply(1:length(BG),function(X){RRlnS(BG[X],CRpred[X],CV = exp(fit$par[2]), alpha=exp(fit$par[1]))})
  }

  legtext1 = paste("BL =",BG); if(length(BG==1)){legtext1=""}
  if(plot){
    plotfit(dat,CRpred,RRpred,point=tv)
    legend('bottomright',legend = c(legtext1,
                                    paste("alpha =",round(exp(fit$par[1]),3)),
                                    paste("CV =",round(exp(fit$par[2]),3)),
                                    paste("O.F. =",round(fit$value,2))),
           bty='n',text.col="red")
  }
  list(fit=fit, CRpred = CRpred, RRpred=RRpred, CRobs = dat$CR, RRobs = dat$RR, BG = BG)
}



fit_LNVSI = function(dat, BG, plot=F, OE=0.4){

  dat=cleandat(dat)
  tv = "BG" %in% names(dat)
  if(tv) BG = dat$BG
  if(tv) OE = minmax(sqrt(dat$varp))
  LNVSI_int = function(par, dat, BG, OE){
    pred = RRlnS(BG, dat$CR, CV = exp(par[2]), alpha=exp(par[1]),beta = par[3])
    -sum(dnorm(logit(dat$RR),logit(pred),OE),log=T)
  }

  fit = optim(par=c(-2,log(0.5),-4),LNVSI_int, method = "L-BFGS-B",
              lower = c(-5,log(0.25),-10),
              upper = c(3,log(4),1),dat=dat, BG=BG, OE=OE,hessian=T)


  if(!tv){
    CRpred = seq(min(dat$CR),max(dat$CR),length.out=100)
    RRpred = RRlnS(BG, CRpred,CV = exp(fit$par[2]), alpha=exp(fit$par[1]),beta=fit$par[3])
  }
  if(tv){
    CRpred = dat$CR
    RRpred = sapply(1:length(BG),function(X){RRlnS(BG[X],CRpred[X],CV = exp(fit$par[2]), alpha=exp(fit$par[1]),beta=fit$par[3])})
  }

  legtext1 = paste("BL =",BG); if(length(BG==1)){legtext1=""}
  if(plot){
    plotfit(dat,CRpred,RRpred,point=tv)
    legend('bottomright',legend = c(legtext1,
                                    paste("alpha =",round(exp(fit$par[1]),3)),
                                    paste("beta =",round(fit$par[3],3)),
                                    paste("CV =",round(exp(fit$par[2]),3)),
                                    paste("O.F. =",round(fit$value,2))),
           bty='n',text.col="blue")
  }
  list(fit=fit, CRpred = CRpred, RRpred=RRpred, CRobs = dat$CR, RRobs = dat$RR, BG = BG)
}




