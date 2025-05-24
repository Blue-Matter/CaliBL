

# Fitted Bag limit MPs

# library(mvtnorm)


# the lognormal bag limit model with slope voluntary release
RRlnS = function(BL,CRvec,CV,alpha=0,beta=-4,n=1E4,CRlim=c(1E-5,30)){
  logit = function(p)log(p/(1-p))
  ilogit = function(x)exp(x)/(1+exp(x))

  RRs = rep(NA,length(CRvec))
  crs = seq(CRlim[1],CRlim[2],length.out=n)

  Vs = ilogit(beta+alpha*crs)
  for(i in 1:length(CRvec)){
    vals0 = dlnorm(crs,log(CRvec[i]),CV)*crs
    valsV = vals0 * (1-Vs)
    RRs[i] = 1 - sum(valsV[crs<BL])/sum(vals0)
    ENC = dlnorm(crs,log(CRvec[i]),CV)*crs   # expected number caught
    ERV = ENC * Vs                            # expected voluntary release
    ERB = dlnorm(crs,log(CRvec[i]),CV)*(crs-BL) * (1-Vs)
    RRs[i] = (sum(ERV)+sum(ERB[crs>=BL]))/sum(ENC)
  }
  RRs
}

explore = function(){
  fits = readRDS("Data/Processed_2023/Northern_fits.rda")
  LNSind = 1
  solve(fits$LB_PC[[LNSind]]$fit$hessian) # parameter 1 is alpha, parameter 2 is CV
  solve(fits$CRFS_PC[[LNSind]]$fit$hessian)
  solve(fits$CRFS_PR[[LNSind]]$fit$hessian)

  fits$LB_PC[[LNSind]]$fit$par
  fits$CRFS_PC[[LNSind]]$fit$par
  fits$CRFS_PR[[LNSind]]$fit$par

  LB_PC = readRDS("Data/Processed_2023/CPFV_logbook_summary.rda")
  CRFS_PC = readRDS("Data/Processed_2023/CFRS_PC_summary.rda")
  CRFS_PR = readRDS("Data/Processed_2023/CFRS_PR_summary.rda")

  yrange = 2016:2020
  CRobs = (CRFS_PR$Ret +CRFS_PR$Rel)/CRFS_PR$CNTRBTRS

  CRref = mean(CRobs[CRFS_PR$Year %in% yrange & CRFS_PR$MA == "N"],na.rm=T)
  Data = Hist@Data
}

BL  = function(x, Data, reps=1, BLim=3, RecPerc=0.6,
               Fdisc = 0.05, varfac = 1, maxVar = Inf,
               datname="CRFS_PR", yrange = 2015:2019,
               CRref = 2.393813, FMSYfrac = 1, MinSz = 55,
               seed=T){

  mus = list()
  mus [["LB_PC"]] = c(-0.7649721, -1.4862156)
  mus [["CRFS_PC"]] = c(-1.068589, -1.601482)
  mus [["CRFS_PR"]] = c(-0.001212253, -4.171962267)

  covars = list()  # parameter 1 is alpha, parameter 2 is CV
  covars[["LB_PC"]] =  matrix(c(0.03093144, -0.01888812,  -0.01888812,  0.02054848),byrow=T,nrow=2)
  covars[["CRFS_PC"]] =  matrix(c(0.12265755, -0.06664325, -0.06664325,  0.06775760),byrow=T,nrow=2)
  covars[["CRFS_PR"]] =  matrix(c( 0.00795155, -0.2714206, -0.27142060, 10.1322595),byrow=T,nrow=2)

  mu = mus[match(datname, names(mus))][[1]]
  covar = covars[match(datname, names(covars))][[1]]
  diag(covar) = diag(covar)*varfac
  diag(covar)[diag(covar)>maxVar] = maxVar
  if(seed)set.seed(x) # select the same model parameters for each simulation
  pars = as.vector(mvtnorm::rmvnorm(1,mu,covar))
  qs = CRref / mean(Data@VInd[x,Data@Year%in%yrange])
  thisyr = dim(Data@VInd)[2]
  CRvec = Data@VInd[x,thisyr]*qs
  RR = RRlnV(BL = BLim,CRvec=CRvec, CV = exp(pars[1]),V=ilogit(pars[2]))
  #lines(CRvec,RR)

  FMSYeffort = FMSYref(x,Data,reps=1)@Effort

  Rec=new('Rec')
  Rec@Effort = FMSYeffort * FMSYfrac
  Rec@DR = RR*RecPerc
  Rec@L5 = MinSz * 0.975
  Rec@LFS = MinSz * 1.025
  Rec@Fdisc = Fdisc
  Rec

}
class(BL) = "MP"

BL_1 = BL_2 = BL_3 = BL_4 = BL_5 = BL_6 = BL
formals(BL_1)$BLim = 1
formals(BL_2)$BLim = 2
formals(BL_3)$BLim = 3
formals(BL_4)$BLim = 4
formals(BL_5)$BLim = 5
formals(BL_6)$BLim = 6
class(BL_1) = class(BL_2) = class(BL_3) = class(BL_4) = class(BL_5) = class(BL_6) = "MP"

E_7 = E_8 = E_9 = E_10 = E_11 = E_12 = BL
formals(E_7)$FMSYfrac = 0.7
formals(E_8)$FMSYfrac = 0.8
formals(E_9)$FMSYfrac = 0.9
formals(E_10)$FMSYfrac = 1.0
formals(E_11)$FMSYfrac = 1.1
formals(E_12)$FMSYfrac = 1.2
class(E_7)= class(E_8)= class(E_9)= class(E_10)= class(E_11)= class(E_12)= "MP"


EU_7 = EU_8 = EU_9 = EU_10 = EU_11 = EU_12 = BL_6
formals(EU_7)$FMSYfrac = 0.7
formals(EU_8)$FMSYfrac = 0.8
formals(EU_9)$FMSYfrac = 0.9
formals(EU_10)$FMSYfrac = 1.0
formals(EU_11)$FMSYfrac = 1.1
formals(EU_12)$FMSYfrac = 1.2
class(EU_7)= class(EU_8)= class(EU_9)= class(EU_10)= class(EU_11)= class(EU_12)= "MP"


EL_7 = EL_8 = EL_9 = EL_10 = EL_11 = EL_12 = EL_13 = EL_14 = BL_1
formals(EL_7)$FMSYfrac = 0.7
formals(EL_8)$FMSYfrac = 0.8
formals(EL_9)$FMSYfrac = 0.9
formals(EL_10)$FMSYfrac = 1.0
formals(EL_11)$FMSYfrac = 1.1
formals(EL_12)$FMSYfrac = 1.2
formals(EL_13)$FMSYfrac = 1.3
formals(EL_14)$FMSYfrac = 1.4
class(EL_7)= class(EL_8)= class(EL_9)= class(EL_10)= class(EL_11)= class(EL_12)=class(EL_13)=class(EL_14)= "MP"


S_45 = S_50 = S_55 = S_60 = S_65 = BL
formals(S_45)$MinSz = 45
formals(S_50)$MinSz = 50
formals(S_55)$MinSz = 55
formals(S_60)$MinSz = 60
formals(S_65)$MinSz = 65
class(S_45) =  class(S_50) = class(S_55) = class(S_60) = class(S_65) = "MP"


SU_45 = SU_50 = SU_55 = SU_60 = SU_65 = BL_4
formals(SU_45)$MinSz = 45
formals(SU_50)$MinSz = 50
formals(SU_55)$MinSz = 55
formals(SU_60)$MinSz = 60
formals(SU_65)$MinSz = 65
class(SU_45) =  class(SU_50) = class(SU_55) = class(SU_60) = class(SU_65) = "MP"


SL_45 = SL_50 = SL_55 = SL_60 = SL_65 = BL_2
formals(SL_45)$MinSz = 45
formals(SL_50)$MinSz = 50
formals(SL_55)$MinSz = 55
formals(SL_60)$MinSz = 60
formals(SL_65)$MinSz = 65
class(SL_45) =  class(SL_50) = class(SL_55) = class(SL_60) = class(SL_65) = "MP"





C_PR_1 = BL_1; C_PR_2 = BL_2; C_PR_3 = BL_3; C_PR_4 = BL_4; C_PR_5 = BL_5; C_PR_6 = BL_6

C_PC_1 = BL_1; formals(C_PC_1)$datname = "CRFS_PC"; class(C_PC_1) = "MP"
C_PC_2 = BL_2; formals(C_PC_2)$datname = "CRFS_PC"; class(C_PC_2) = "MP"
C_PC_3 = BL_3; formals(C_PC_3)$datname = "CRFS_PC"; class(C_PC_3) = "MP"
C_PC_4 = BL_4; formals(C_PC_4)$datname = "CRFS_PC"; class(C_PC_4) = "MP"
C_PC_5 = BL_5; formals(C_PC_5)$datname = "CRFS_PC"; class(C_PC_5) = "MP"
C_PC_6 = BL_6; formals(C_PC_6)$datname = "CRFS_PC"; class(C_PC_6) = "MP"

L_PC_1 = BL_1; formals(L_PC_1)$datname = "LB_PC"; class(L_PC_1) = "MP"
L_PC_2 = BL_2; formals(L_PC_2)$datname = "LB_PC"; class(L_PC_2) = "MP"
L_PC_3 = BL_3; formals(L_PC_3)$datname = "LB_PC"; class(L_PC_3) = "MP"
L_PC_4 = BL_4; formals(L_PC_4)$datname = "LB_PC"; class(L_PC_4) = "MP"
L_PC_5 = BL_5; formals(L_PC_5)$datname = "LB_PC"; class(L_PC_5) = "MP"
L_PC_6 = BL_6; formals(L_PC_6)$datname = "LB_PC"; class(L_PC_6) = "MP"



# Demonstration Bag limit MPs

curE_BG  = function(x, Data, reps=1, BL=3, CV=0.4, Vint=0.35, Vslp=0.075, RecPerc=0.4, Fdisc = 0.5){

  ly = dim(Data@VInd)[2]
  CR = Data@VInd[x,ly]

  DR_rec = lnormVM2_RR(CR, BagLim = BL, CV_VM2 = CV, Vint = Vint, Vslp2 = Vslp)

  Rec=new('Rec')
  Rec@Effort = 1
  Rec@DR = DR_rec*0.4
  Rec@Fdisc = Fdisc

  Rec

}


BL1_D2 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=1,Fdisc=0.2)
BL2_D2 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=2,Fdisc=0.2)
BL3_D2 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=3,Fdisc=0.2)
BL4_D2 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=4,Fdisc=0.2)
BL5_D2 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=5,Fdisc=0.2)
BL6_D2 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=6,Fdisc=0.2)


BL1_D5 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=1,Fdisc=0.5)
BL2_D5 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=2,Fdisc=0.5)
BL3_D5 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=3,Fdisc=0.5)
BL4_D5 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=4,Fdisc=0.5)
BL5_D5 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=5,Fdisc=0.5)
BL6_D5 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=6,Fdisc=0.5)


BL1_D8 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=1,Fdisc=0.8)
BL2_D8 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=2,Fdisc=0.8)
BL3_D8 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=3,Fdisc=0.8)
BL4_D8 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=4,Fdisc=0.8)
BL5_D8 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=5,Fdisc=0.8)
BL6_D8 = function(x,Data,reps=1)curE_BG(x=x,Data=Data,reps=reps,BL=6,Fdisc=0.8)

class(BL1_D2) = class(BL2_D2) = class(BL3_D2) = class(BL4_D2) = class(BL5_D2) = class(BL6_D2) = 'MP'
class(BL1_D5) = class(BL2_D5) = class(BL3_D5) = class(BL4_D5) = class(BL5_D5) = class(BL6_D5) = 'MP'
class(BL1_D8) = class(BL2_D8) = class(BL3_D8) = class(BL4_D8) = class(BL5_D8) = class(BL6_D8) = 'MP'

BLMPs = paste0("BL",rep(1:6,3),"_D",rep(c(2,5,8),each=6))

print("CDFW Bag Limit MPs Loaded")

# disused

disused_BLfunc = function(x, Data, reps=1, BLim=3, RecPerc=0.6, # this was when LNVS was used
               Fdisc = 0.05, varfac = 0.5, maxVar = 0.2,
               datname="CRFS_PR", yrange = 2015:2019,
               CRref = 2.393813, FMSYfrac = 1, MinSz = 55 ){

  mus = list()
  mus [["LB_PC"]] = c(-0.5,0.2897611)
  mus [["CRFS_PC"]] = c(-1.36329543,  0.05143813)
  mus [["CRFS_PR"]] = c(-0.3221240, -0.1627627)

  covars = list()  # parameter 1 is alpha, parameter 2 is CV
  covars[["LB_PC"]] =  matrix(c(151.66685218, -0.021578167,  -0.021578167,  0.001323111),byrow=T,nrow=2)
  covars[["CRFS_PC"]] =  matrix(c(0.14015458, -0.021578167, -0.021578167, 0.008575245),byrow=T,nrow=2)
  covars[["CRFS_PR"]] =  matrix(c(0.01188297, -0.007439349, -0.007439349, 0.005708415),byrow=T,nrow=2)

  mu = mus[match(datname, names(mus))][[1]]
  covar = covars[match(datname, names(covars))][[1]]
  diag(covar) = diag(covar)*varfac
  diag(covar)[diag(covar)>maxVar] = maxVar
  pars = as.vector(mvtnorm::rmvnorm(1,mu,covar))

  qs = CRref / mean(Data@VInd[x,Data@Year%in%yrange])
  thisyr = dim(Data@VInd)[2]
  CRvec = Data@VInd[x,thisyr]*qs

  RR = RRlnS(BL = BLim,CRvec=CRvec,CV = exp(pars[2]),alpha=exp(pars[1]))

  FMSYeffort = FMSYref(x,Data,reps=1)@Effort

  Rec=new('Rec')
  Rec@Effort = FMSYeffort * FMSYfrac
  Rec@DR = RR*RecPerc
  Rec@L5 = MinSz * 0.975
  Rec@LFS = MinSz * 1.025
  Rec@Fdisc = Fdisc
  Rec

}
class(BL) = "MP"
