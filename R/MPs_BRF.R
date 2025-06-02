

# Fitted Bag limit MPs

# library(mvtnorm)


explore = function(){
  setwd("C:/Users/tcarruth/Documents/GitHub/CDFW_Bag_Limits")
  setwd("C:/Github/CDFW_Bag_Limits")
  fits = readRDS("Data/Processed_BRF_2024/Northern_fits.rda")
  LNSind = 3

  fits$LB_PC[[LNSind]]$fit$par
  fits$CRFS_PC[[LNSind]]$fit$par
  fits$CRFS_PR[[LNSind]]$fit$par

  solve(fits$LB_PC[[LNSind]]$fit$hessian) # parameter 1 is alpha, parameter 2 is CV
  solve(fits$CRFS_PC[[LNSind]]$fit$hessian)
  solve(fits$CRFS_PR[[LNSind]]$fit$hessian)

  LB_PC = readRDS("Data/Processed_2023/CPFV_logbook_summary.rda")
  CRFS_PC = readRDS("Data/Processed_2023/CFRS_PC_summary.rda")
  CRFS_PR = readRDS("Data/Processed_2023/CFRS_PR_summary.rda")

  yrange = 2016:2020
  CRobs = (CRFS_PR$Ret +CRFS_PR$Rel)/CRFS_PR$CNTRBTRS

  CRref = mean(CRobs[CRFS_PR$Year %in% yrange & CRFS_PR$MA == "N"],na.rm=T)
  Data = Hist@Data
}

BL_BRF = function(x, Data, reps=1, BLim=3, RecPerc=0.78,  # recent stock synthesis assessment has it mostly rec catch (2019)
               Fdisc = 0.2, varfac = 0.5, maxVar = 0.2,
               datname="CRFS_PC", yrange = 2015:2019,
               CRref = 2.393813, FMSYfrac = 1, MinSz = 55 ){

  mus = list()
  mus [["LB_PC"]] = c(-5, -0.598)
  mus [["CRFS_PC"]] = c(-3.188, -0.041)
  mus [["CRFS_PR"]] = c(-0.462, 0.227)

  covars = list()  # parameter 1 is alpha, parameter 2 is CV
  covars[["LB_PC"]] =  matrix(c(12.124, -0.0387,  -0.0387,  0.00214),byrow=T,nrow=2)
  covars[["CRFS_PC"]] =  matrix(c(34.63, -1.605, -1.605, 0.1153),byrow=T,nrow=2)
  covars[["CRFS_PR"]] =  matrix(c(0.0446, -0.02847, -0.02847, 0.01935),byrow=T,nrow=2)

  mu = mus[match(datname, names(mus))][[1]]
  covar = covars[match(datname, names(covars))][[1]]
  diag(covar) = diag(covar)*varfac
  diag(covar)[diag(covar)>maxVar] = maxVar
  set.seed(x)
  pars = as.vector(mvtnorm::rmvnorm(1,mu,covar))

  qs = CRref / mean(Data@VInd[x,Data@Year%in%yrange])
  thisyr = dim(Data@VInd)[2]
  CRvec = Data@VInd[x,thisyr]*qs

  RR = RRlnS(BL = BLim,CRvec=CRvec,CV = exp(pars[2]),alpha=exp(pars[1]))
  FMSYeffort = FMSYref(x,Data,reps=1)@Effort

  Rec=new('Rec')
  Rec@Effort = FMSYeffort * FMSYfrac
  Rec@DR = RR*RecPerc
  #Rec@L5 = MinSz * 0.975
  #Rec@LFS = MinSz * 1.025
  Rec@Fdisc = Fdisc
  Rec

}
class(BL_BRF) = "MP"

BL_BRF_1 = BL_BRF_2 = BL_BRF_3 = BL_BRF_4 = BL_BRF_5 = BL_BRF_6 = BL_BRF_7 = BL_BRF_8 = BL_BRF_9 = BL_BRF_10 = BL_BRF
formals(BL_BRF_1)$BLim = 1
formals(BL_BRF_2)$BLim = 2
formals(BL_BRF_3)$BLim = 3
formals(BL_BRF_4)$BLim = 4
formals(BL_BRF_5)$BLim = 5
formals(BL_BRF_6)$BLim = 6
formals(BL_BRF_7)$BLim = 7
formals(BL_BRF_8)$BLim = 8
formals(BL_BRF_9)$BLim = 9
formals(BL_BRF_10)$BLim = 10
class(BL_BRF_1) = class(BL_BRF_2) = class(BL_BRF_3) = class(BL_BRF_4) = class(BL_BRF_5) =
  class(BL_BRF_6) = class(BL_BRF_7) = class(BL_BRF_8) = class(BL_BRF_9) = class(BL_BRF_10) = "MP"

PRM_BRF_05 = PRM_BRF_10 = PRM_BRF_15 = PRM_BRF_20 = PRM_BRF_25 = PRM_BRF_30 = PRM_BRF_35 = PRM_BRF_40 = BL_BRF
formals(PRM_BRF_05)$Fdisc = 0.05
formals(PRM_BRF_10)$Fdisc = 0.10
formals(PRM_BRF_15)$Fdisc = 0.15
formals(PRM_BRF_20)$Fdisc = 0.20
formals(PRM_BRF_25)$Fdisc = 0.25
formals(PRM_BRF_30)$Fdisc = 0.30
formals(PRM_BRF_35)$Fdisc = 0.35
formals(PRM_BRF_40)$Fdisc = 0.40
class(PRM_BRF_05) = class(PRM_BRF_10) = class(PRM_BRF_15) = class(PRM_BRF_20) = class(PRM_BRF_25) =
  class(PRM_BRF_30) = class(PRM_BRF_35) = class(PRM_BRF_40) = "MP"

PRML_BRF_05 = PRML_BRF_10 = PRML_BRF_15 = PRML_BRF_20 = PRML_BRF_25 = PRML_BRF_30 = PRML_BRF_35 = PRML_BRF_40 = BL_BRF_1
formals(PRML_BRF_05)$Fdisc = 0.05
formals(PRML_BRF_10)$Fdisc = 0.10
formals(PRML_BRF_15)$Fdisc = 0.15
formals(PRML_BRF_20)$Fdisc = 0.20
formals(PRML_BRF_25)$Fdisc = 0.25
formals(PRML_BRF_30)$Fdisc = 0.30
formals(PRML_BRF_35)$Fdisc = 0.35
formals(PRML_BRF_40)$Fdisc = 0.40
class(PRML_BRF_05) = class(PRML_BRF_10) = class(PRML_BRF_15) = class(PRML_BRF_20) = class(PRML_BRF_25) =
  class(PRML_BRF_30) = class(PRML_BRF_35) = class(PRML_BRF_40) = "MP"


PRMU_BRF_05 = PRMU_BRF_10 = PRMU_BRF_15 = PRMU_BRF_20 = PRMU_BRF_25 = PRMU_BRF_30 = PRMU_BRF_35 = PRMU_BRF_40 = BL_BRF_6
formals(PRMU_BRF_05)$Fdisc = 0.05
formals(PRMU_BRF_10)$Fdisc = 0.10
formals(PRMU_BRF_15)$Fdisc = 0.15
formals(PRMU_BRF_20)$Fdisc = 0.20
formals(PRMU_BRF_25)$Fdisc = 0.25
formals(PRMU_BRF_30)$Fdisc = 0.30
formals(PRMU_BRF_35)$Fdisc = 0.35
formals(PRMU_BRF_40)$Fdisc = 0.40
class(PRMU_BRF_05) = class(PRMU_BRF_10) = class(PRMU_BRF_15) = class(PRMU_BRF_20) = class(PRMU_BRF_25) =
  class(PRMU_BRF_30) = class(PRMU_BRF_35) = class(PRMU_BRF_40) = "MP"



E_BRF_7 = E_BRF_8 = E_BRF_9 = E_BRF_10 = E_BRF_11 = E_BRF_12 = BL_BRF
formals(E_BRF_7)$FMSYfrac = 0.7
formals(E_BRF_8)$FMSYfrac = 0.8
formals(E_BRF_9)$FMSYfrac = 0.9
formals(E_BRF_10)$FMSYfrac = 1.0
formals(E_BRF_11)$FMSYfrac = 1.1
formals(E_BRF_12)$FMSYfrac = 1.2
class(E_BRF_7)= class(E_BRF_8)= class(E_BRF_9)= class(E_BRF_10)= class(E_BRF_11)= class(E_BRF_12)= "MP"


EU_BRF_7 = EU_BRF_8 = EU_BRF_9 = EU_BRF_10 = EU_BRF_11 = EU_BRF_12 = BL_BRF_6
formals(EU_BRF_7)$FMSYfrac = 0.7
formals(EU_BRF_8)$FMSYfrac = 0.8
formals(EU_BRF_9)$FMSYfrac = 0.9
formals(EU_BRF_10)$FMSYfrac = 1.0
formals(EU_BRF_11)$FMSYfrac = 1.1
formals(EU_BRF_12)$FMSYfrac = 1.2
class(EU_BRF_7)= class(EU_BRF_8)= class(EU_BRF_9)= class(EU_BRF_10)= class(EU_BRF_11)= class(EU_BRF_12)= "MP"


EL_BRF_7 = EL_BRF_8 = EL_BRF_9 = EL_BRF_10 = EL_BRF_11 = EL_BRF_12 = EL_BRF_13 = EL_BRF_14 = BL_BRF_1
formals(EL_BRF_7)$FMSYfrac = 0.7
formals(EL_BRF_8)$FMSYfrac = 0.8
formals(EL_BRF_9)$FMSYfrac = 0.9
formals(EL_BRF_10)$FMSYfrac = 1.0
formals(EL_BRF_11)$FMSYfrac = 1.1
formals(EL_BRF_12)$FMSYfrac = 1.2
formals(EL_BRF_13)$FMSYfrac = 1.3
formals(EL_BRF_14)$FMSYfrac = 1.4
class(EL_BRF_7)= class(EL_BRF_8)= class(EL_BRF_9)= class(EL_BRF_10)= class(EL_BRF_11)= class(EL_BRF_12)=class(EL_BRF_13)=class(EL_BRF_14)= "MP"






C_BRF_PR_1 = BL_BRF_1; C_BRF_PR_2 = BL_BRF_2; C_BRF_PR_3 = BL_BRF_3; C_BRF_PR_4 = BL_BRF_4; C_BRF_PR_5 = BL_BRF_5; C_BRF_PR_6 = BL_BRF_6

C_BRF_PR_1 = BL_BRF_1; formals(C_BRF_PR_1)$datname = "CRFS_PR"; class(C_BRF_PR_1) = "MP"
C_BRF_PR_2 = BL_BRF_2; formals(C_BRF_PR_2)$datname = "CRFS_PR"; class(C_BRF_PR_2) = "MP"
C_BRF_PR_3 = BL_BRF_3; formals(C_BRF_PR_3)$datname = "CRFS_PR"; class(C_BRF_PR_3) = "MP"
C_BRF_PR_4 = BL_BRF_4; formals(C_BRF_PR_4)$datname = "CRFS_PR"; class(C_BRF_PR_4) = "MP"
C_BRF_PR_5 = BL_BRF_5; formals(C_BRF_PR_5)$datname = "CRFS_PR"; class(C_BRF_PR_5) = "MP"
C_BRF_PR_6 = BL_BRF_6; formals(C_BRF_PR_6)$datname = "CRFS_PR"; class(C_BRF_PR_6) = "MP"

L_BRF_PC_1 = BL_BRF_1; formals(L_BRF_PC_1)$datname = "LB_PC"; class(L_BRF_PC_1) = "MP"
L_BRF_PC_2 = BL_BRF_2; formals(L_BRF_PC_2)$datname = "LB_PC"; class(L_BRF_PC_2) = "MP"
L_BRF_PC_3 = BL_BRF_3; formals(L_BRF_PC_3)$datname = "LB_PC"; class(L_BRF_PC_3) = "MP"
L_BRF_PC_4 = BL_BRF_4; formals(L_BRF_PC_4)$datname = "LB_PC"; class(L_BRF_PC_4) = "MP"
L_BRF_PC_5 = BL_BRF_5; formals(L_BRF_PC_5)$datname = "LB_PC"; class(L_BRF_PC_5) = "MP"
L_BRF_PC_6 = BL_BRF_6; formals(L_BRF_PC_6)$datname = "LB_PC"; class(L_BRF_PC_6) = "MP"


print("CDFW Bag Limit MPs for Black Rockfish Loaded")

