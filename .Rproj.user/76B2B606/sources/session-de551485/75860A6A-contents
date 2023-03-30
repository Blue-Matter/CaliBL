
# Source

# Fitting lognormal distributions to observed data by half year --------------------

dlnorm_int= function(pars,mids,prb,cpue,mode='opt'){
  pred = dlnorm(mids,pars[1],exp(pars[2]))
  if(mode=='opt')return(sum((prb - pred)^2))
  if(mode=='plot'){
    plot(mids,prb,ylim=c(0,max(pred,prb)),xlab="Catch rate",ylab="Relative Frequency")
    grid()
    lines(mids,pred,col="red")
    legend('topright',legend=c(paste("Log Mean =",round(pars[1],3)),
                               paste("SD =", round(exp(pars[2]),3))),bty='n',text.col="red")
    #legend('top',legend=c(paste("Mean =",round(mean(log(cpue)),3)),
    #                     paste("SD =", round(sd(log(cpue)),3))),bty='n')
  }
}


dofits = function(cpue,cutoff=15, lab=""){
  
  breaks=0:60
  xseq = 0:cutoff
  frq = hist(cpue,breaks,plot=F)
  keep = frq$mids < cutoff
  mids = frq$mids[keep]
  prb = frq$counts[keep]/sum(frq$counts[keep])
  inits=c(log(1.2),log(1))
  fit = optim(par=inits,dlnorm_int,method = 'L-BFGS-B',lower=c(log(0.3),log(0.2)),upper=c(log(4),log(2)),hessian=T, mids=mids,prb=prb,cpue=cpue)
  dlnorm_int(fit$par,mids=mids,prb=prb,cpue=cpue,mode="plot")
  mtext(lab,line=0.5,cex=0.9)
  exp(fit$par)
  
}

# RR predictors ---------------------------------------------------------------------------

pois_RR = function(CRvec,BagLim,V_pois){
  predTheta_pois = ppois(BagLim, CRvec, lower.tail=F) 
  1-((1-V_pois)*(1-predTheta_pois))
}

lnorm1_RR = function(CRvec, BagLim, CV_ln, V_ln){
  predTheta_ln = plnorm(BagLim,log(CRvec),sd=CV_ln,lower.tail=F)
  1-((1-V_ln)*(1-predTheta_ln))
}

lnorm2_RR = function(CRvec, BagLim, CV_Ln2, V1_ln2, V2_ln2){
  predTheta_ln2 = plnorm(BagLim,log(CRvec),sd=CV_ln2,lower.tail=F)
  predOmega_ln2 = plnorm(1,log(CRvec),sd=CV_ln2,lower.tail=T)
  1-((1-V1_ln2) * predOmega_ln2 + (1-V2_ln2) * (1-predTheta_ln2-predOmega_ln2) )
}


#CRvec = seq(0,3.2,length.out=20)
#L = CRvec-(max(CRvec)/2) * 0.1
#VM = exp(L)/(1+exp(L))

lnormVM1_RR = function(CRvec, BagLim, CV_VM1, Vslp1){
  predTheta_VM1 = plnorm(BagLim,log(CRvec),sd=CV_VM1,lower.tail=F)
  VM1 = CRvec*Vslp1 
  1-((1-VM1)*(1-predTheta_VM1))
}

lnormVM2_RR = function(CRvec, BagLim, CV_VM2, Vint, Vslp2){
  predTheta_VM2 = plnorm(BagLim,log(CRvec),sd=CV_VM2,lower.tail=F)
  VM2 = Vint + CRvec*Vslp2 
  1-((1-VM2)*(1-predTheta_VM2))
}


# source_URL("https://raw.github.com/christophergandrud/christophergandrud.github.com/master/SourceCode/CarsScatterExample.R")

# source_URL("https://raw.github.com/blue-matter/blue-matter.github.com/master/CDFW_Bag_Limits/Code/Source.R")

RDS_from_web <- function(myurl) {
  
  tempFile_location<- tempfile()
  download.file(myurl, tempFile_location)
  b <- readRDS(tempFile_location)
  file.remove(tempFile_location)
  b
}


sourcefunc= function(){
  install.packages('downloader')
  library(downloader)
  source_url("https://raw.githubusercontent.com/Blue-Matter/CDFW_Bag_Limits_IO/master/docs/Code/Source.R",prompt=F)
  source_url("https://raw.githubusercontent.com/Blue-Matter/CDFW_Bag_Limits_IO/master/docs/Code/MPs.R",prompt=F)
  # load(url("https://github.com/Blue-Matter/CDFW_Bag_Limits_IO/raw/master/Blue-Matter/CDFW_Bag_Limits_IO/ddocs/OMs/DemoOM"))
  temp = RDS_from_web("https://github.com/Blue-Matter/CDFW_Bag_Limits_IO/blob/2b4464fcfd946dc0f7d5f4b6d13aae97c847eb85/docs/Data/DemoOM.RDATA")
}


# Fitting 2-parameter bag limit model -----------------------------------------------------------------










# Fitting 3-parameter bag limit model


print('CDFW Bag Limit Source Code Loaded')



