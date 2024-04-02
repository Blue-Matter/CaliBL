

# Fitted Bag limit MPs

# library(mvtnorm)


explore = function(){
  setwd("C:/Users/tcarruth/Documents/GitHub/CDFW_Bag_Limits")
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

BL_BRF = function(x, Data, reps=1, BLim=3, RecPerc=0.6,
               Fdisc = 0.05, varfac = 0.5, maxVar = 0.2,
               datname="CRFS_PC", yrange = 2015:2019,
               CRref = 2.393813, FMSYfrac = 1, MinSz = 55 ){

  mus = list()
  mus [["LB_PC"]] = c(-4.4706441, -0.3949828)
  mus [["CRFS_PC"]] = c(-2.0104686, -0.2249778)
  mus [["CRFS_PR"]] = c(-1.0301354, 0.1988783)

  covars = list()  # parameter 1 is alpha, parameter 2 is CV
  covars[["LB_PC"]] =  matrix(c(37.759409, -0.16618101,  -0.166181,  0.00177351),byrow=T,nrow=2)
  covars[["CRFS_PC"]] =  matrix(c(0.21294789, -0.03922442, -0.03922442, 0.01206332),byrow=T,nrow=2)
  covars[["CRFS_PR"]] =  matrix(c(0.09824655, -0.04240281, -0.04240281, 0.02105550),byrow=T,nrow=2)

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
class(BL_BRF) = "MP"

BL_BRF_1 = BL_BRF_2 = BL_BRF_3 = BL_BRF_4 = BL_BRF_5 = BL_BRF_6 = BL_BRF
formals(BL_BRF_1)$BLim = 1
formals(BL_BRF_2)$BLim = 2
formals(BL_BRF_3)$BLim = 3
formals(BL_BRF_4)$BLim = 4
formals(BL_BRF_5)$BLim = 5
formals(BL_BRF_6)$BLim = 6
class(BL_BRF_1) = class(BL_BRF_2) = class(BL_BRF_3) = class(BL_BRF_4) = class(BL_BRF_5) = class(BL_BRF_6) = "MP"


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

