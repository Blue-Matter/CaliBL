




# Source
Proj_Comp_2 = function(MSEobj,MSEobj2,ccols = c("#ff000070","#0000ff70","#99999970"),
                     lcols = c("red","blue","black"),refyr = 2020, MPs = c(1,1),newplot=T,
                     legcode = NA, dolabs=T, labs = NA){

  if(newplot)par(mfrow=c(2,1), mai=c(0.4,0.5,0.02,0.02),omi=c(0.4,0.2,0.05,0.05))

  docloud = function(xs,dat,lcol="black",ccol="#99999990"){
    qs=apply(dat,2,quantile,p=c(0.05,0.5,0.95))
    polygon(c(xs,rev(xs)),c(qs[1,],rev(qs[3,])),border=F,col=ccol)
    lines(xs,qs[2,],col=lcol)
  }

  yind = 2:MSEobj@proyears
  # SSB
  ylim = c(0,range(apply(MSEobj@SB_SBMSY,2:3,quantile,p=c(0.05,0.95)))[2])
  yrs = refyr + 1:MSEobj@proyears
  plot(range(yrs),ylim,col="white",xlab="",ylab=""); grid()
  docloud(xs = yrs[yind], dat=MSEobj@SB_SBMSY[,MPs[1],yind],ccol=ccols[1],lcol=lcols[1])
  if(is.function(legcode))legcode()
  docloud(xs = yrs[yind], dat=MSEobj2@SB_SBMSY[,MPs[2],yind],ccol=ccols[2],lcol=lcols[2])
  if(dolabs) mtext("SSB / SSBMSY",2,line=2.2)
  if(!is.na(labs[1]))mtext(labs[1],adj=0.03,line=-1,cex=0.85)

  # Yield
  ylim = c(0,range(apply(MSEobj@Catch,2:3,quantile,p=c(0.05,0.95)))[2])
  yrs = refyr + 1:MSEobj@proyears
  plot(range(yrs),ylim,col="white",xlab="",ylab=""); grid()
  docloud(xs = yrs[yind], dat=MSEobj@Catch[,MPs[1],yind],ccol=ccols[1],lcol=lcols[1])
  docloud(xs = yrs[yind], dat=MSEobj2@Catch[,MPs[2],yind],ccol=ccols[2],lcol=lcols[2])
  if(dolabs) mtext("Yield (t)",2,line=2.2)
  if(!is.na(labs[1]))mtext(labs[2],adj=0.03,line=-1,cex=0.85)


  if(newplot)mtext("Projection year",1,outer=T,line=0.3)

}


# MSE plotting
Proj_Comp = function(MSEobj,ccols = c("#ff000070","#0000ff70","#99999970"),
                     lcols = c("red","blue","black"),refyr = 2020, newplot=T,
                     legcode = NA, legcode2 = NA, dolabs=T, labs = NA){

  if(newplot)par(mfrow=c(2,1), mai=c(0.4,0.5,0.02,0.02),omi=c(0.4,0.2,0.05,0.05))

  docloud = function(xs,dat,lcol="black",ccol="#99999990"){
    qs=apply(dat,2,quantile,p=c(0.05,0.5,0.95))
    polygon(c(xs,rev(xs)),c(qs[1,],rev(qs[3,])),border=F,col=ccol)
    lines(xs,qs[2,],col=lcol)
  }

  # SSB

  SSB =MSEobj@SSB
  SSBMSY = MSEobj@RefPoint$SSBMSY[,,MSEobj@nyears]
  Bmet = SSB/array(SSBMSY,dim(SSB))
  ylim = c(0,range(apply(Bmet,2:3,quantile,p=c(0.05,0.95)))[2])
  yrs = refyr + 1:MSEobj@proyears
  yind = 2:MSEobj@proyears
  plot(range(yrs),ylim,col="white",xlab="",ylab=""); grid()
  for(i in 1:MSEobj@nMPs)docloud(xs = yrs[yind], dat=Bmet[,i,yind],ccol=ccols[i],lcol=lcols[i])
  if(dolabs)mtext("SSB / SSBMSY",2,line=2.2)
  if(!is.na(labs[1]))mtext(labs[1],adj=0.03,line=-1,cex=0.85)

  if(is.function(legcode))legcode()

  # Yield
  ylim = c(0,range(apply(MSEobj@Catch,2:3,quantile,p=c(0.05,0.95)))[2])
  yrs = refyr + 1:MSEobj@proyears
  plot(range(yrs),ylim,col="white",xlab="",ylab=""); grid()
  for(i in 1:MSEobj@nMPs)docloud(xs = yrs[yind], dat=MSEobj@Catch[,i,yind],ccol=ccols[i],lcol=lcols[i])
  if(dolabs) mtext("Yield",2,line=2.2)
  if(!is.na(labs[1]))mtext(labs[2],adj=0.03,line=-1,cex=0.85)
  if(is.function(legcode2))legcode2()

  if(newplot)mtext("Projection year",1,outer=T,line=0.3)

}


# comparative yield plot
YBTplot = function(MSEs,yrs=31:40,nowyr=2020,textadj = 0.04,textcex = 0.9,labs=T,laby=NA,title="",labinv=NA, xadj=1.3,yadj=0.2, plab ="", doaxislabs = F, newplot=T){

  if(newplot)par(mai=c(0.8,0.8,0.05,0.05))
  tcols = c("darkgrey","red","green","blue")
  cols = c("#999999","#ff000095","#00ff0095","#0000ff99")
  Y30 = lapply(MSEs,function(x)apply(x@Catch[,,yrs],2,mean))
  SSB = lapply(MSEs,function(x)apply(x@SB_SBMSY[,,yrs],2,mean))
  nams = lapply(MSEs,function(x)x@MPs)
  xlim = range(SSB)*c(0.99,1.01)
  ylim = range(Y30)*c(0.99,1.01)
  xrng = range(SSB)[2]-range(SSB)[1]
  yrng = range(Y30)[2]-range(Y30)[1]
  plot(SSB[[1]],Y30[[1]],col="white",ylim=ylim,xlim=xlim,xlab="",ylab="")
  if(is.na(labinv[1]))labinv=rep(1,length(SSB))
  grid()
  for(i in 1:length(SSB)){
    points(SSB[[i]],Y30[[i]],col=cols[i],cex=0.95,pch=19)
    lines(SSB[[i]],Y30[[i]],col=cols[i])
    textx = SSB[[i]]+xadj*labinv[i]*textadj*xrng
    texty = Y30[[i]]+yadj*labinv[i]*textadj*yrng
    if(labs)text(textx,texty,nams[[i]],cex=textcex,col=tcols[i])
  }
  if(!is.na(laby[1]))legend('topright',legend=laby,text.col=tcols,bty='n',title=title,title.col="black")
  mtext(plab,adj=0.03,line=-1,cex=0.85)
  if(doaxislabs) mtext(paste0("Expected yield (proj. yrs. ",min(yrs+nowyr),"-",max(yrs+nowyr),",)"),2,line=2.1)
  if(doaxislabs) mtext(paste0("Expected SSB/SSBMSY (proj. yrs. ",min(yrs+nowyr),"-",max(yrs+nowyr),")"),1,line=2.2)
}

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


