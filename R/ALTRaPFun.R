StanRes<-function(PGfile,UpperT,LowerT,DirT,Bac,DH=NULL){
#set.seed(101) #optional from shiny checkbox
#data output in these arrays - LR, second and third arrays are Pr(sec) and Pr(direct) transfer probabilities respectively
LRmod<-matrix(0,6,4)
Hpmod<-matrix(0,6,4)
Hdmod<-matrix(0,6,4)
LRmod2<-matrix(0,6,4)
Hpmod2<-matrix(0,6,4)
Hdmod2<-matrix(0,6,4)
LRmodT<-matrix(0,6,4)
HpmodT<-matrix(0,6,4)
HdmodT<-matrix(0,6,4)
LRmod2T<-matrix(0,6,4)
Hpmod2T<-matrix(0,6,4)
Hdmod2T<-matrix(0,6,4)
SecSum<-matrix(0,6,4)
LRmod3<-matrix(0,6,4)
Hpmod3<-matrix(0,6,4)
Hdmod3<-matrix(0,6,4)
LRtot<-array(dim=c(6,4,4000))
LRtog<-array(dim=c(6,4,4000))
Hp<-matrix(0,6,4)
Hd<-matrix(0,6,4)
d<-1
Offmean<-c(1:100)
#####################################################################################
#To simulate an example where a couple have cohabited and there is a domestic violence.The only evidence is DNA. Suspect denies assault.
# Can the DNA evidence be useful.
#Hp is direct contact and secondary transfer
#Hd is secondary transfer only
########## INPUTS description with suggested example#####################
#n=4 ### n= number of contacts for a given hour***************CHANGE Manually
##### Decide the number of hours before analysis NOTE ALL TIMES ARE RELATIVE TO TIME OF COLLECTION OF SAMPLE
#We need to know the time since contact Where time=6 means somewhere between 6-7pm and time=7 = between7-8pm
#UpperT=6
#We need to know the time of first contact eg time of shower
#LowerT=18 # ie the victim had a shower 12 hours before contact which removed his/her DNA
#DirT is the time since the direct transfer ie time between the event and collection of samples
#Now we estimate probabilities of Transfer between UpperT and LowerT hours representing the trannsfer interval
#Bac=.2 #Pr Background DNA
####################################################################################
LowerT=LowerT-1 #Adjust LowerT to 1 hr blocks
if(LowerT<UpperT){LowerT=UpperT}
#set plot space
Time=PGfile$Time
PGfile$binLR1m<-mapply(binx,PGfile$Log10LR,1)
log10bin=PGfile$binLR1m
plot(Time,log10bin,xlab="Time (hours since contact)",ylab="Pr(LR>x)",ylim=c(0,.8))

#Extract coefficients for given LR>x
Gtrx<-c(1,2,3,4,6,8)###labels used in file "DirectTransfer.csv"
PlotCol<-c("black","green","red", "blue", "black","green")
Plotlty<-c(2,3,1,3,4,6)
legend("topright",inset=.05, legend=c("LR>10","LR>100", "LR>1000",">10,000","LR>1m","LR>10e8"),
       col=c("black","green","red", "blue", "black","green"), lty=c(2,3,1,3,4,6),lwd=2, cex=1)

for (pa in 1:length(Gtrx)){
##########FUNCTION binx###############
binx <-function(x,Gtrx){if(x< Gtrx){y<-0} else {y<-1}
  return(y)}
#########################################
PGfile$binLR1m<-mapply(binx,PGfile$Log10LR,Gtrx[pa])

log10bin=PGfile$binLR1m#set variables/change according to the input file
#########ADJUST time interval to correspond to direct T
Time=PGfile$Time
dat=as.data.frame(cbind(Time,log10bin))
#Priors can be specified as follows - not used
##priorI=normal(location = PriorIntercept, autoscale = TRUE)##If needed
##priorN=normal(location = PriorTime, autoscale = TRUE)##if needed
Model=stan_glm(log10bin~Time, family=binomial(link="logit"),dat)#logistic regression create model
newdataX<-data.frame(Time=seq(0,24,by=0.2))#generate Time values to test
yweight<-predict(Model,newdata=newdataX,type="resp")#create model
mydat<-as.data.frame(cbind(newdata=newdataX,yweight))
##PLOT
lines(mydat$Time,mydat$yweight,col=PlotCol[pa], lty=Plotlty[pa],lwd=3)
#ADD LINES

###COEFFICIENT SIMULATIONS
###launch_shinystan(Model)#Launch this to test parameter errors
##This creates random parameters in var sims
sims <- as.matrix(Model)#generate 4000 random parameters Intercept and Time
colnames(sims)<-c("I","Tc")#change column names and make dataframe
sims<-as.data.frame(sims)#sims contains 2 coefs
#Testdat<-Coefs[which(Coefs$Threshold==Gtrx[pa]),] #Set the time offset according to LR size
#Offset<-Testdat$Offset
##OFFSET CALC
#Zap<-OffsetData[which(OffsetData$equiv==Gtrx[pa]),]## This command is redundant if exact method is used
#Offset<-OffsetCalc(OffsetData,sims,Zap$Pr)## This command is redundant if exact method is used
Offset<-OffsetOyv(sims,Gtrx[pa])##Use this command for Exact method
Offmean[d]<-mean(Offset)
d<-d+1
# All data are in rows in array sims2
###########################Program
for (n in 1:4){ #cycle through the number of contacts per hour
#UpperT=12#point to start 0
#LowerT=24 #no of hours after first contact eg time since shower 12 for time zero
Tarray<-mapply(PrFun,UpperT,LowerT,sims$I,sims$Tc,Offset)#Calc secondary T module

Tarray<-1-(1-Tarray)^n #adjust for number of contacts per hour
####Now calculate the Pr for BN for secondary transfer - need to take account of all sims
#Check to make sure >1hour sec contact otherwise Poisson routine crashes
ifelse(UpperT!=LowerT,
SecT<-apply(Tarray,2,PoisBin),
SecT<-Tarray)
##########################
#######################################################
#Calculate Pr for secondary and direct transfer Hp
#UpperT is time since sample taken and time of assault
direct<-mapply(Logistic,sims$I,sims$Tc,DirT)#Set to time zero time direct where DirT= time since the offence - time samples collected
AsT<-direct
##########This is used if we need to combine sec and dir transfer for agreed events
#Now calculate the Pr Hd/Hp direct + Pr secondary using Poisson Binomial
#Now add Pr of direct Hd/Hp transfer here
if(!is.null(DH)){
directDH<-mapply(Logistic,sims$I,sims$Tc,DH)#Set to time zero time direct where DH=time since samples collected
#AsTDH<-directDH
TarrayDH<-rbind(SecT,directDH)
#And calculate Hp which includes secondary transfer here
SecT<-apply(TarrayDH,2,PoisBin)
}
#######################
##BN CALC using R HUGIN is 3 lines below - redundant here
#BayesHd<-mapply(TablesBN,0,AsT,SecT,Bac)
#BayesHp<-mapply(TablesBN,1,AsT,SecT,Bac)
#LR<-(BayesHp[1,]+BayesHp[2,])/(BayesHd[1,]+BayesHd[2,])
####################################################
#BN CALC USING FORMULAE -fast method
ResLR<-BNformula(SecT,AsT,Bac)
#Conservative method Calculate a quantile eg 0.05
#Calc POI only
LRq<-quantile(ResLR$LRPOI,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
Hpq<-quantile(AsT,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))#Dir T
Hdq<-quantile(SecT,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))#Sec T
Hpaq<-quantile(ResLR$NumPOI,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))#Hp of the LR
Hdaq<-quantile(ResLR$DenPOI,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))#Hd of the LR
LRmod[pa,n]<-LRq[4]#Median LR
Hpmod[pa,n]<-Hpq[4]
Hdmod[pa,n]<-Hdq[4]
LRmod2[pa,n]<-LRq[2]#Lower 5 percentile LR
Hpmod2[pa,n]<-Hpq[2]# the pr direct T
Hdmod2[pa,n]<-Hdq[2]#the pr secondary T
LRmod3[pa,n]<-LRq[6]#Upper 5 percentile LR
Hpmod3[pa,n]<-Hpq[6]# the pr direct T
Hdmod3[pa,n]<-Hdq[6]#the pr secondary T
LRtot[pa,n,]<-ResLR$LRPOI #POI only
LRtog[pa,n,]<-ResLR$LRTog #Unknown +POI
Hp[pa,n]<-Hpaq[4]#Median numerator for LR
Hd[pa,n]<-Hdaq[4]#Median denom for LR
######################
#Calc POI +U/B###
LRq<-quantile(ResLR$LRTog,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
#Hpq<-quantile(AsT,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
#Hdq<-quantile(SecT,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
LRmodT[pa,n]<-LRq[4]#Median LR
HpmodT[pa,n]<-Hpq[4] #Note this is not the numerator it is Pr direct transfer
HdmodT[pa,n]<-Hdq[4]  #Note this is not the denom it is Pr sec transfer
LRmod2T[pa,n]<-LRq[2]#5 percentile
Hpmod2T[pa,n]<-Hpq[2]
Hdmod2T[pa,n]<-Hdq[2]
SecSum[pa,n]<-median(SecT)
}
}
Res<-(list(LRmod=LRmod,Hpmod=Hpmod,Hdmod=Hdmod,LRmod2=LRmod2,Hpmod2=Hpmod2,Hdmod2=Hdmod2,Offmean=Offmean,Model=Model,SecT=SecT,AsT=AsT,sims=sims,
LRmodT=LRmodT,HpmodT=HpmodT,HdmodT=HdmodT,LRmod2T=LRmod2T,Hpmod2T=Hpmod2T,Hdmod2T=Hdmod2T,ResLR=ResLR$LRPOI,ResLRT=ResLR$LRTog,SecSum=SecSum,Hp=Hp,Hd=Hd,LRmod3=LRmod3,Hpmod3=Hpmod3,Hdmod3=Hdmod3,LRtot=LRtot,LRtog=LRtog))
return(Res)
}
#####################################
###############FUNCTIONS###############
binx <-function(x,Gtrx){if(x< Gtrx){y<-0} else {y<-1}
  return(y)}


######FUNCTION Calculate logistic regression probability FUNCTION
Logistic<-function(I,Tc,H) {
PrH<-1/(1+exp(-(I+Tc*H)))
}

########FUNCTION To calculate Binomials
##direct is the Pr of direct transfer - set to zero under Hd usually
#An array is made from upperT to lowerT ie time period since collection of samples
#For sec T the offset is applied to H and y is estimated via logistic regression
PrFun<-function(UpperT,LowerT,I,Tc,Offset){
xarray<-array(UpperT:LowerT)
z=LowerT-UpperT+1
Tarray<-matrix(0,NROW(I),z)
for(i in 1:z){
H=Offset + xarray[i]
Tarray[,i]<-Logistic(I,Tc,H) #Tarray has Pr for one contact per hour between UpperT and LowerT
}#end for loop
return(Tarray)
}







###CALCULATE OFFSET
#CoEf<-coef(Model)#interceptand time respectively
##FUNCTION OFFSET
OffsetCalc<-function(OffsetData,sims,SecT){
H<-0
xmin<-array(1:nrow(sims))
for (i in 1:nrow(sims)){
I<-sims[i,1]
Tc<-sims[i,2]
xmin[i]<-optimize(ResOff,c(1, 100),Tc,I,SecT)$minimum
}
Offset=as.numeric(xmin)
Offset=Offset-0.5
return(Offset)
}

####FUNCTION To calculate the secondary transfer OFFSET Minimise function
###########
ResOff<-function(H,Tc,I,SecT){
AsT<-Logistic(I,Tc,H)
LR<-AsT/SecT
xmin=(LR-1)^2
return(xmin)
}

##FUNCTION OYVINDS EXACT METHOD TO CALC OFFSET using Pareto
OffsetOyv<-function(sims,SecT){
  k<-0.666666667
  alpha<-4.432757
  beta<-5.779410
  xmin<-array(1:nrow(sims))
  for (i in 1:nrow(sims)){
    b0<-sims[i,1]
    b1<-sims[i,2]
    Yx = (((beta+SecT)/beta)^alpha)/(1-k) - 1
    xmin[i]= -(b0 + log(Yx))/b1
  }
  Offset=as.numeric(xmin)
  Offset<-Offset-0.5
  return(Offset)
}

######FUNCTION Calculate logistic regression probability FUNCTION
Logistic<-function(I,Tc,H) {
PrH<-1/(1+exp(-(I+Tc*H)))
}

##Function LR formulae substitute for BN
BNformula<-function(SecT,AsT,Bac){
U<-(AsT*Bac)+(AsT*(1-Bac))+(Bac*(1-AsT))
NUM<-(SecT*(1-AsT))+AsT
#POI ONLY -no U or background
NumPOI<-NUM*(1-Bac)
DenPOI<-SecT*(1-U)
LRPOI<-NumPOI/DenPOI
#POI and Unknown DNA present
NumTog<-NUM*Bac
DenTog<-SecT*U
LRTog<-NumTog/DenTog
ResLR<-(list(NumPOI=NumPOI, DenPOI=DenPOI, LRPOI=LRPOI, NumTog=NumTog, DenTog=DenTog,LRTog=LRTog))
return(ResLR)
}

##########FUNCTION to calc Poisson binomial
library(ggplot2)
PoisBin<-function(Tarray){
#Now combine the number of hours using Poisson Binomial
library(poibin)
kk=0 #Note this is Pr of zero successes
#simulate multiple contacts at 6-12 hours using 10^3 data
Bi=ppoibin(kk=kk, pp=Tarray, method = "DFT-CF",wts=NULL)
NoH=1-Bi #This is the answer # because we want Pr of 1 or more successes
return(NoH)
}

##VIOLIN PLOT FUNCTION
###VIOLIN PLOT POI ONLY
##Important note - to analyse LRs from POI only input
##violin(Results$LRtot)
##To analyse POI + U input:
##violin(Results$LRtog)
violin<-function(Results,n){
  LRT<-array(dim=c(6,4,4000))
  LRT<-(Results) #POI only
  #LRT<-(Results$LRtog)
  a=LRT[,n,] #extract n=1
  a=t(a)#transpose array
  colnames(a)<-c("x=1","x=2","x=3","x=4","x=6","x=8")
  a=as.data.frame(a)
  a=log10(a)
  newcol<-array(1:24000,dim=c(24000,1))
  newcol[1:4000]="x=1"
  newcol[4001:8000]="x=2"
  newcol[8001:12000]="x=3"
  newcol[12001:16000]="x=4"
  newcol[16001:20000]="x=6"
  newcol[20001:24000]="x=8"
  #newcol<-t(newcol)
  apfile<-(a$"x=1")
  apfile<-Append(apfile,a$"x=2",rows=TRUE)
  apfile<-Append(apfile,a$"x=3",rows=TRUE)
  apfile<-Append(apfile,a$"x=4",rows=TRUE)
  apfile<-Append(apfile,a$"x=5",rows=TRUE)
  apfile<-Append(apfile,a$"x=6",rows=TRUE)
  apfile<-Append(apfile,a$"x=8",rows=TRUE)
  apfile<-cbind(apfile,newcol)
  apfile=as.data.frame(apfile)
  apfile$V2<-as.factor(apfile$V2)
  apfile$apfile<-as.numeric(apfile$apfile)
  #############VIOLIN PLOT
  p<-ggplot(apfile, aes(x=V2,y=apfile)) + stat_summary(
    fun.min = function(x) { quantile(x,0.025) },
    fun.max = function(x) { quantile(x,0.975) },
    fun = median, geom="crossbar", width=0.4, fill="red")
  p<-p+stat_summary(
    fun.min = function(x) { quantile(x,0.05) },
    fun.max = function(x) { quantile(x,0.95) },
    fun = median, geom="crossbar", width=0.4, fill="blue")
  p<-p+geom_violin()
  p<-p+geom_boxplot(width=0.2,fill="green")
  p<-p+coord_flip()
  p<-p+ labs(y="log10 likelihood ratio",x="Logistic regression decision threshold")
  p<-p+ theme(text = element_text(size = 15))
  p<-p + scale_y_continuous(breaks=seq(0,10,0.5))
  p<-p + ggtitle(paste("Quantile plots of log10LRs; n=",n)) + theme(plot.title = element_text(hjust = 0.5))
  p
}
  ###################################################
  ##############FUNCTION VIOLIN FOR SENSITIVITY PLOT
  violinSens<-function(Results,n,Gtrx,emp.data){
    if (Gtrx<8){
      emp.data<-format(round(emp.data,digits=2))
    }else{
      emp.data<-format(round(emp.data,digits=3))
    }
    LRT<-array(dim=c(6,4,4000))
    LRT<-(Results) #POI only
    #LRT<-(Results$LRtog)
    a=LRT[,n,] #extract n=1
    a=t(a)#transpose array
    colnames(a)<-c("x=1","x=2","x=3","x=4","x=5","x=6")
    a=as.data.frame(a)
    a=log10(a)
    newcol<-array(1:24000,dim=c(24000,1))
    newcol[1:4000]="x=1"
    newcol[4001:8000]="x=2"
    newcol[8001:12000]="x=3"
    newcol[12001:16000]="x=4"
    newcol[16001:20000]="x=5"
    newcol[20001:24000]="x=6"
    #newcol<-t(newcol)
    apfile<-(a$"x=1")
    apfile<-Append(apfile,a$"x=2",rows=TRUE)
    apfile<-Append(apfile,a$"x=3",rows=TRUE)
    apfile<-Append(apfile,a$"x=4",rows=TRUE)
    apfile<-Append(apfile,a$"x=5",rows=TRUE)
    apfile<-Append(apfile,a$"x=6",rows=TRUE)
    apfile<-cbind(apfile,newcol)
    apfile=as.data.frame(apfile)
    apfile$V2<-as.factor(apfile$V2)
    apfile$apfile<-as.numeric(apfile$apfile)

    apfile$V2<-mapvalues(apfile$V2,from = c("x=1","x=2","x=3","x=4","x=5","x=6"), to= c(emp.data[Gtrx,2],emp.data[Gtrx,3],emp.data[Gtrx,4],emp.data[Gtrx,5],emp.data[Gtrx,6],emp.data[Gtrx,7]))
    p=violinplot(apfile,n,Gtrx)
    p
  }



#############FUNCTION VIOLIN PLOT
violinplot<-function(apfile,n,x){
  p<-ggplot(apfile, aes(x=V2,y=apfile)) + stat_summary(
    fun.min = function(x) { quantile(x,0.025) },
    fun.max = function(x) { quantile(x,0.975) },
    fun = median, geom="crossbar", width=0.4, fill="red")
  p<-p+stat_summary(
    fun.min = function(x) { quantile(x,0.05) },
    fun.max = function(x) { quantile(x,0.95) },
    fun = median, geom="crossbar", width=0.4, fill="blue")
  p<-p+geom_violin()
  p<-p+geom_boxplot(width=0.2,fill="green")
  p<-p+coord_flip()
  p<-p+ labs(y="log10 likelihood ratio",x="Pr secondary transfer")
  p<-p+ theme(text = element_text(size = 15))
  p<-p + scale_y_continuous(breaks=seq(0,10,0.5))
  #p<-p + ggtitle("Density plots: T=6, S=18,n=1") + theme(plot.title = element_text(hjust = 0.5))
  #p<-p + ggtitle(paste("Quantile plots of log10LRs; n=",n)) + theme(plot.title = element_text(hjust = 0.5))
  p<-p + ggtitle(paste("Sensitivity plot, x=",x,",n=",n)) + theme(plot.title = element_text(hjust = 0.5))

  #p<-p+scale_y_continuous("percentiles",sec_axis(datax))

  p
}

##FUNCTION SENSITIVITY ANALYSIS

StanSens<-function(PGfile,UpperT,LowerT,DirT,Bac,Gtrx,DH=NULL){
  LowerT=LowerT-1 #Adjust LowerT to 1 hr blocks
  if(LowerT<UpperT){LowerT=UpperT}
  #set.seed(101)
  #data output in these arrays - LR, second and third arrays are Pr(sec) and Pr(direct) transfer probabilities respectively
  LRmod<-matrix(0,6,4)
  Hpmod<-matrix(0,6,4)
  Hdmod<-matrix(0,6,4)
  LRmod2<-matrix(0,6,4)
  Hpmod2<-matrix(0,6,4)
  Hdmod2<-matrix(0,6,4)
  LRmodT<-matrix(0,6,4)
  HpmodT<-matrix(0,6,4)
  HdmodT<-matrix(0,6,4)
  LRmod2T<-matrix(0,6,4)
  Hpmod2T<-matrix(0,6,4)
  Hdmod2T<-matrix(0,6,4)
  SecSum<-matrix(0,6,4)
  LRmod3<-matrix(0,6,4)
  Hpmod3<-matrix(0,6,4)
  Hdmod3<-matrix(0,6,4)
  LRtot<-array(dim=c(6,4,4000))
  LRtog<-array(dim=c(6,4,4000))
  Hp<-matrix(0,6,4)
  Hd<-matrix(0,6,4)
  d<-1
  Offmean<-c(1:100)
  #####################################################################################
  #To simulate an example where a couple have cohabited and there is a domestic violence.The only evidence is DNA. Suspect denies assault.
  # Can the DNA evidence be useful.
  #Hp is direct contact and secondary transfer
  #Hd is secondary transfer only
  ########## INPUTS description with suggested example#####################
  #n=4 ### n= number of contacts for a given hour***************CHANGE Manually
  ##### Decide the number of hours before analysis NOTE ALL TIMES ARE RELATIVE TO TIME OF COLLECTION OF SAMPLE
  #We need to know the time since contact Where time=6 means somewhere between 6-7pm and time=7 = between7-8pm
  #UpperT=6
  #We need to know the time of first contact eg time of shower
  #LowerT=18 # ie the victim had a shower 12 hours before contact which removed his/her DNA
  #DirT is the time since the direct transfer ie time between the event and collection of samples
  #Now we estimate probabilities of Transfer between UpperT and LowerT hours representing the trannsfer interval
  #Bac=.2 #Pr Background DNA
  ####################################################################################
  ##input sensitivity data
  ################
  emp.data<-data.frame(
    emp_id=c(1,2,3,4,5,6,7,8),##these are values of x from 1000 bootstraps of the data for each value. See table 7.
    "1"=c(0.1643053,0.08928389,0.05223608,0.03238231,0.02103199,0.01419361,0.009890758,0.007082565),#median
    "2"=c(0.2044461,0.11160486,0.0652951,0.04047788,0.02628999,0.01774201,0.012598177,0.009738527),#75 percentile
    "3"=c(0.2303873,0.13392583,0.08471798,0.05828902,0.04144534,0.0303635,0.022477853,0.017503828),#90 percentile
    "4"=c(0.2492831,0.15230387,0.09956092,0.06636477,0.04461198,0.031711,0.025494703,0.021201174),#95 percentile
    "5"=c(0.2785421,0.1702874,0.10799499,0.07480975,0.05699712,0.04453348,0.033199656,0.025755442),#97.5 percerntile
    "6"=c(0.3045894,0.1918559,0.12430095,0.0852547,0.06232639,0.04613692,0.036475513,0.030509873),#99 percentile
    stringsAsFactors = FALSE
  )
  #OffsetData from pareto distribution below

  OffsetData<-data.frame(
  "equiv"=c(1,2,3,4,5,6,7,8),
  "Pr"=c(0.16266221,0.088391058,0.051713726,0.032058488,0.020821677,0.014051672,0.009791853,0.007011741),
  stringsAsFactors = FALSE
  )


  #######################
  #for (pa in 1:length(Gtrx)){
  ##########FUNCTION binx###############
  binx <-function(x,Gtrx){if(x< Gtrx){y<-0} else {y<-1}
    return(y)}
  #########################################
  PGfile$binLR1m<-mapply(binx,PGfile$Log10LR,Gtrx)

  log10bin=PGfile$binLR1m#set variables/change according to the input file
  #########ADJUST time interval to correspond to direct T
  Time=PGfile$Time
  dat=as.data.frame(cbind(Time,log10bin))
  Model=stan_glm(log10bin~Time, family=binomial(link="logit"),dat)#logistic regression create model
  #XX=summary(Model,digits=5)#check model
  newdataX<-data.frame(Time=seq(0,24,by=0.2))#generate Time values to test
  yweight<-predict(Model,newdata=newdataX,type="resp")#create model
  mydat<-as.data.frame(cbind(newdata=newdataX,yweight))

  ###Note to do COEFFICIENT SIMULATIONS (switched off here)
  #launch_shinystan(Results$Model)#Launch this to test parameter errors
  ##This creates random parameters in var sims

  sims <- as.matrix(Model)#generate 4000 random parameters Intercept and Time
  colnames(sims)<-c("I","Tc")#change column names and make dataframe
  sims<-as.data.frame(sims)#sims contains 2 coefs

  ##OFFSET CALC
  Zap<-emp.data[which(emp.data$emp_id==Gtrx),]### This command is redundant if Oyv method is used
  Zap<-as.numeric(Zap)
  for (pa in 1:6){
    Offset<-OffsetCalc(emp.data,sims,Zap[pa+1])## Or use this command for Optim method
    #Offset<-OffsetOyv(sims,Gtrx)##use for Oyv method
    Offmean[d]<-mean(Offset)
    d<-d+1
    # All data are in rows in array sims
    ###########################Program
    for (n in 1:4){ #cycle through the number of contacts per hour
      #UpperT=12#point to start 0
      #LowerT=24 #no of hours after first contact eg time since shower 12 for time zero
      Tarray<-mapply(PrFun,UpperT,LowerT,sims$I,sims$Tc,Offset)#Calc secondary T module

      Tarray<-1-(1-Tarray)^n #adjust for number of contacts per hour
      ####Now calculate the Pr for BN for secondary transfer - need to take account of all sims
      #Check to make sure >1hour sec contact otherwise Poisson routine crashes
      ifelse(UpperT!=LowerT,
             SecT<-apply(Tarray,2,PoisBin),
             SecT<-Tarray)
      ##########################
      #######################################################
      #Calculate Pr for secondary and direct transfer Hp
      #UpperT is time since sample taken and time of assault
      direct<-mapply(Logistic,sims$I,sims$Tc,DirT)#Set to time zero time direct where DirT= time since the offence - time samples collected
      AsT<-direct
      ##########This is used if we need to combine sec and dir transfer for agreed events
      #Now calculate the Pr Hd/Hp direct + Pr secondary using Poisson Binomial
      #Now add Pr of direct Hd/Hp transfer here
      if(!is.null(DH)){
        directDH<-mapply(Logistic,sims$I,sims$Tc,DH)#Set to time zero time direct where DH=time since samples collected
        #AsTDH<-directDH
        TarrayDH<-rbind(SecT,directDH)
        #And calculate Hp which includes secondary transfer here
        SecT<-apply(TarrayDH,2,PoisBin)
      }
      ##########################################################
      ##BN CALC using R HUGIN is 3 lines below Not used
      #BayesHd<-mapply(TablesBN,0,AsT,SecT,Bac)
      #BayesHp<-mapply(TablesBN,1,AsT,SecT,Bac)
      #LR<-(BayesHp[1,]+BayesHp[2,])/(BayesHd[1,]+BayesHd[2,])
      ##########################################################
      #BN CALC USING FORMULAE -fast method
      ResLR<-BNformula(SecT,AsT,Bac)
      #Conservative method Calculate a quantile eg 0.05
      #Calc POI only
      LRqI<-quantile(ResLR$LRPOI,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
      Hpq<-quantile(AsT,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))#Dir T
      Hdq<-quantile(SecT,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))#Sec T
      Hpaq<-quantile(ResLR$NumPOI,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))#Hp of the LR
      Hdaq<-quantile(ResLR$DenPOI,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))#Hd of the LR
      LRmod[pa,n]<-LRqI[4]#Median
      Hpmod[pa,n]<-Hpq[4]
      Hdmod[pa,n]<-Hdq[4]
      LRmod2[pa,n]<-LRqI[2]#Lower 5 percentile
      Hpmod2[pa,n]<-Hpq[2]
      Hdmod2[pa,n]<-Hdq[2]
      LRmod3[pa,n]<-LRqI[6]#Upper 5 percentile
      Hpmod3[pa,n]<-Hpq[6]
      Hdmod3[pa,n]<-Hdq[6]
      LRtot[pa,n,]<-ResLR$LRPOI #POI only
      LRtog[pa,n,]<-ResLR$LRTog #Unknown +POI
      Hp[pa,n]<-Hpaq[4]#Median numerator for LR
      Hd[pa,n]<-Hdaq[4]#Median den for LR
      #Calc POI +U/B###
      LRqT<-quantile(ResLR$LRTog,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
      #Hpq<-quantile(AsT,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
      #Hdq<-quantile(SecT,probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
      LRmodT[pa,n]<-LRqT[4]#Median LR
      HpmodT[pa,n]<-Hpq[4] #Note this is not the numerator it is Pr direct transfer
      HdmodT[pa,n]<-Hdq[4]  #Note this is not the denom it is Pr sec transfer
      LRmod2T[pa,n]<-LRqT[2]#5 percentile
      Hpmod2T[pa,n]<-Hpq[2]
      Hdmod2T[pa,n]<-Hdq[2]
      SecSum[pa,n]<-median(SecT)
    }
  }
  Res<-(list(LRmod=LRmod,Hpmod=Hpmod,Hdmod=Hdmod,LRmod2=LRmod2,Hpmod2=Hpmod2,Hdmod2=Hdmod2,Offmean=Offmean,Model=Model,SecT=SecT,AsT=AsT,sims=sims,
             LRmodT=LRmodT,HpmodT=HpmodT,HdmodT=HdmodT,LRmod2T=LRmod2T,Hpmod2T=Hpmod2T,Hdmod2T=Hdmod2T,ResLR=ResLR$LRPOI,ResLRT=ResLR$LRTog,SecSum=SecSum,
             Hp=Hp,Hd=Hd,LRmod3=LRmod3,Hpmod3=Hpmod3,Hdmod3=Hdmod3,LRtot=LRtot,LRtog=LRtog,Offset=Offset,sims=sims,emp.data=emp.data,Gtrx=Gtrx,LRqI=LRqI,LRqT=LRqT))
  return(Res)
}



