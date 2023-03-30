# Bag limit MPs


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

BGMPs = paste0("BL",rep(1:6,3),"_D",rep(c(2,5,8),each=6))

print("CDFW Bag Limit MPs Loaded")

