#Load required libraries and files

library(gsl)
library(deSolve)
library(insol)
library(MASS)
library(mvtnorm)
library(VGAM)

dyn.load("SMC_CR14_NAM.so")
dyn.load("FindForcing.so")
dyn.load("insol.so")

# Read data
DataSet<-read.table("~/HuyberRecords/HuyberODP846.txt")
Obs_Length<-308
Obs<-DataSet[Obs_Length:1,5]
Obs_Depths<-DataSet[Obs_Length:1,1]/100
Obs_Depth_Step<-Obs_Depths[1:(Obs_Length-1)]-Obs_Depths[2:Obs_Length]

#Precompute required terms for the astronomical forcing
nap <- 30
nao <- 30
times <- seq(-1000000,0,by=100)/1e4 

Read_Pre<-.Fortran('read_precession', as.integer(nap), as.numeric(rep(0,nap)), as.numeric(rep(0,nap)), as.numeric(rep(0,nap)), DUP=FALSE)
Read_Obl<-.Fortran('read_obliquity', as.integer(nao), as.numeric(rep(0,nao)), as.numeric(rep(0,nao)), as.numeric(rep(0,nao)), DUP=FALSE)
amppre <- Read_Pre[[2]]
omepre <- Read_Pre[[3]]
angpre <- Read_Pre[[4]]
ampobl <- Read_Obl[[2]]
omeobl <- Read_Obl[[3]]
angobl <- Read_Obl[[4]]

GPre<-rep(0,length(times))
GCop<-rep(0,length(times))
GObl<-rep(0,length(times))
for(i in 1:length(times)){
Astro_Out<-.C("FF",as.numeric(amppre),as.numeric(omepre),as.numeric(angpre),as.numeric(ampobl),as.numeric(omeobl),as.numeric(angobl),as.numeric(times[i]),F1=as.numeric(0),F2=as.numeric(0),F3=as.numeric(0))
GPre[i]<-Astro_Out$F1
GCop[i]<-Astro_Out$F2
GObl[i]<-Astro_Out$F3}

# Set age control point
ACP<--780000

#Set prior hyperparameters
x_1_Min<--1.5
x_1_Max<-1.5
x_2_Min<--2.5
x_2_Max<-2.5

D_Min<-3
D_Max<-5
S_Min<-0.5
S_Max<-2

#Select tuning parameters for the particle filters 
Parameter_Particles<-1000
State_Particles<-1000
MCMC_Length<-10

#Array to store parameter particle information
Results<-matrix(rep(0,19*Parameter_Particles),ncol=19)

#Lists to store state values, anscestors and weights
Base<-matrix(rep(0,Parameter_Particles*State_Particles),ncol=State_Particles)
X_1<-list(Base)
for(i in 2:Obs_Length){
Base<-matrix(rep(0,Parameter_Particles*State_Particles),ncol=State_Particles)
X_1[[i]]<-Base}

X_2<-X_1
T<-X_1
w<-X_1
A<-X_1


#Vectors to store effective sample sizes and information about the particle diversity
ESS_Record<-rep(Parameter_Particles,Obs_Length)
Unique_Record_PreMove<-rep(Parameter_Particles,Obs_Length)
Unique_Record_PostMove<-rep(Parameter_Particles,Obs_Length)

#Vectors to store marginal likelihoods
Likelihood_Vector<-rep(0,Parameter_Particles)
Marginal_Likelihood_Vector<-rep(0,Obs_Length)

#Set output folder location
FOL<-99
FileTop<-as.name(paste("../Results/UnforcedCR12CompactedODP846R",FOL,"/",sep=""))
dir.create(paste(FileTop))

#Sample initial values
i<-1
for(j in 1:Parameter_Particles){

Beta_0<-rnorm(1,mean=0.4,sd=0.3)
Beta_1<-rnorm(1,mean=0,sd=0.4)
Beta_2<-rexp(1,rate=1/0.5)
Delta<-rexp(1,rate=1/0.5)
P<-0
C<-0
E<-0
Sigma_1<-rexp(1,rate=1/0.3)
Sigma_2<-rexp(1,rate=1/0.5)
Sigma_Obs<-rexp(1,rate=1/0.1)
D<-runif(1,min=D_Min,max=D_Max)
S<-runif(1,min=S_Min,max=S_Max)
Mu_S<-rgamma(1,shape=180,scale=1/4000000)
SD_S<-rexp(1,rate=500)
Alpha<-rgamma(1,shape=10,scale=2)
p_0<-rbeta(1,shape1=45,shape2=15)
Grad<-rexp(1,rate=4000)

X_1[[i]][j,]<-runif(State_Particles,min=x_1_Min,max=x_1_Max)
X_2[[i]][j,]<-runif(State_Particles,min=x_2_Min,max=x_2_Max)
T[[i]][j,]<-rnorm(State_Particles,mean=ACP,sd=2000)
UD<-Obs_Depths+(Grad/(1-p_0))*Obs_Depths^2
w[[i]][j,]<-dinv.gaussian(-T[[i]][j,],mu=(UD[1]-UD[Obs_Length])/Mu_S,lambda=(UD[1]-UD[Obs_Length])^2/SD_S^2)*dnorm(Obs[1],mean=S*X_1[[i]][j,]+D,sd=Sigma_Obs)

Likelihood_Vector[j]<-(1/State_Particles)*sum(w[[i]][j,])

Results[j,]<-c(Beta_0,Beta_1,Beta_2,Delta,P,C,E,Sigma_1,Sigma_2,Sigma_Obs,D,S,Mu_S,SD_S,Alpha,p_0,Grad,Likelihood_Vector[j],log(Likelihood_Vector[j]))

}

#Normalize weights and evaluate marginal likelihood
Results[,18]<-Results[,18]/sum(Results[,18])

Marginal_Likelihood_Vector[1]<-sum((1/Parameter_Particles)*Likelihood_Vector)

# Write out particle status and marginal likelihood estimate
write.table(Parameter_Particles,file=paste(FileTop,"ESS_Record.txt",sep=""),col.names=FALSE,row.names=FALSE)
write.table(Parameter_Particles,file=paste(FileTop,"Unique_Record_PreMove.txt",sep=""),col.names=FALSE,row.names=FALSE)
write.table(Parameter_Particles,file=paste(FileTop,"Unique_Record_PostMove.txt",sep=""),col.names=FALSE,row.names=FALSE)
write.table(Marginal_Likelihood_Vector[1],file=paste(FileTop,"Marginal_likelihood_Vector.txt",sep=""),col.names=FALSE,row.names=FALSE)

#Arrays to store initial state values
Int_1<-X_1[[1]]
Int_2<-X_2[[1]]
Int_T<-T[[1]]

#Loop through observations
for(i in 2:Obs_Length){

#Calculate ESS
ESS<-1/(sum(Results[,18]^2))
ESS_Record[i]<-ESS

#Write out ESS
write.table(ESS_Record[i],file=paste(FileTop,"ESS_Record.txt",sep=""),col.names=FALSE,row.names=FALSE,append=TRUE)

#Set resampling condition (typically if the ESS falls below N/2)
if(ESS<(Parameter_Particles/2)){

#Calculated weighted means and covariance for the parameters
Mean_Beta_0<-weighted.mean(Results[,1],w=Results[,18])
Mean_Beta_1<-weighted.mean(Results[,2],w=Results[,18])
Mean_Beta_2<-weighted.mean(Results[,3],w=Results[,18])
Mean_Delta<-weighted.mean(Results[,4],w=Results[,18])
#Mean_P<-weighted.mean(Results[,5],w=Results[,18])
#Mean_C<-weighted.mean(Results[,6],w=Results[,18])
#Mean_E<-weighted.mean(Results[,7],w=Results[,18])
Mean_Sigma_1<-weighted.mean(Results[,8],w=Results[,18])
Mean_Sigma_2<-weighted.mean(Results[,9],w=Results[,18])
Mean_Sigma_Obs<-weighted.mean(Results[,10],w=Results[,18])
Mean_D<-weighted.mean(Results[,11],w=Results[,18])
Mean_S<-weighted.mean(Results[,12],w=Results[,18])
Mean_Mu_S<-weighted.mean(Results[,13],w=Results[,18])
Mean_SD_S<-weighted.mean(Results[,14],w=Results[,18])
Mean_Alpha<-weighted.mean(Results[,15],w=Results[,18])
Mean_p_0<-weighted.mean(Results[,16],w=Results[,18])
Mean_Grad<-weighted.mean(Results[,17],w=Results[,18])

Parameter_Mean<-c(Mean_Beta_0,Mean_Beta_1,Mean_Beta_2,Mean_Delta,Mean_Sigma_1,Mean_Sigma_2,Mean_Sigma_Obs,Mean_D,Mean_S,Mean_Mu_S,Mean_SD_S,Mean_Alpha,Mean_p_0,Mean_Grad)
Parameter_Covariance<-cov.wt(Results[,c(1:4,8:17)],wt=Results[,18])$cov

#Resample particles and update relevant arrays
Sample<-sample(1:Parameter_Particles,size=Parameter_Particles,prob=Results[,18],replace=TRUE)
Results<-Results[Sample,]

for(j in 1:(i-1)){
X_1[[j]]<-X_1[[j]][Sample,]
X_2[[j]]<-X_2[[j]][Sample,]
T[[j]]<-T[[j]][Sample,]
w[[j]]<-w[[j]][Sample,]
A[[j]]<-A[[j]][Sample,]
Int_1<-Int_1[Sample,]
Int_2<-Int_2[Sample,]
Int_T<-Int_T[Sample,]}

#Calculated weighted means and variances for the state particles
Mean_1<-weighted.mean(c(Int_1),w=c(w[[i-1]]))
Mean_2<-weighted.mean(c(Int_2),w=c(w[[i-1]]))
Mean_T<-weighted.mean(c(Int_T),w=c(w[[i-1]]))

SD_1<-sqrt(wtd.var(c(Int_1),weights=c(w[[i-1]]),normwt=TRUE))
SD_2<-sqrt(wtd.var(c(Int_2),weights=c(w[[i-1]]),normwt=TRUE))
SD_T<-sqrt(wtd.var(c(Int_T),weights=c(w[[i-1]]),normwt=TRUE))

#Write out number of unique particles from resampling
Unique_Record_PreMove[i]<-length(unique(Results[,1]))

write.table(Unique_Record_PreMove[i],file=paste(FileTop,"Unique_Record_PreMove.txt",sep=""),col.names=FALSE,row.names=FALSE,append=TRUE)

#Perform PMCMC update to diversify particles
for(j in 1:Parameter_Particles){
for(k in 1:MCMC_Length){
Parameter_Current<-Results[j,c(1:4,8:17)]

#Propose parameters
Parameter_Proposed<-rmvnorm(1,mean=Parameter_Mean,sigma=Parameter_Covariance)

#If proposal is outside of the prior support we can reject immediately
if(Parameter_Proposed[8]>D_Min&&Parameter_Proposed[8]<D_Max&&Parameter_Proposed[9]>S_Min&&Parameter_Proposed[9]<S_Max&&Parameter_Proposed[3]>0&&Parameter_Proposed[4]>0&&Parameter_Proposed[5]>0&&Parameter_Proposed[6]>0&&Parameter_Proposed[7]>0&&Parameter_Proposed[10]>0&&Parameter_Proposed[11]>0&&Parameter_Proposed[12]>0&&Parameter_Proposed[13]>0&&Parameter_Proposed[14]>0){

#Initialise new particle filter
X_1_Proposed<-rnorm(State_Particles,mean=Mean_1,sd=SD_1)
X_2_Proposed<-rnorm(State_Particles,mean=Mean_2,sd=SD_2)
T_Proposed<-rnorm(State_Particles,mean=Mean_T,sd=SD_T)
UD<-Obs_Depths+(Parameter_Proposed[14]/(1-Parameter_Proposed[13]))*Obs_Depths^2
w_Proposed<-(dinv.gaussian(-T_Proposed,mu=(UD[1]-UD[Obs_Length])/Parameter_Proposed[10],lambda=(UD[1]-UD[Obs_Length])^2/Parameter_Proposed[11]^2)*dunif(X_1_Proposed,min=x_1_Min,max=x_1_Max)*dunif(X_2_Proposed,min=x_2_Min,max=x_2_Max)*dnorm(Obs[1],mean=Parameter_Proposed[9]*X_1_Proposed+Parameter_Proposed[8],sd=Parameter_Proposed[7])*dnorm(T_Proposed,mean=ACP,sd=2000))/(dnorm(X_1_Proposed,mean=Mean_1,sd=SD_1)*dnorm(X_2_Proposed,mean=Mean_2,sd=SD_2)*dnorm(T_Proposed,mean=Mean_T,sd=SD_T))
A_Proposed<-1:State_Particles

TempMatrixX1<-X_1_Proposed
TempMatrixX2<-X_2_Proposed
TempMatrixT<-T_Proposed
TempMatrixw<-w_Proposed
TempMatrixA<-A_Proposed

Temp_Int_1<-X_1_Proposed
Temp_Int_2<-X_2_Proposed
Temp_Int_T<-T_Proposed

if(sum(w_Proposed)>0){

#Store the (log) marginal likelihood estimate
LML_Track<-log((1/State_Particles)*sum(w_Proposed))

#Loop through necessary observations
if(i>2){
for(LV in 2:(i-1)){

#If the solution explodes (which can happen with some models with poor prior choices, for example) reject the proposal
if(sum(w_Proposed)>0 && LML_Track!="NaN" && min(X_1_Proposed)>-Inf && max(X_1_Proposed)<Inf && min(X_2_Proposed)>-Inf && max(X_2_Proposed)<Inf){

#Draw seed for the particle filter
Random_Seed<-runif(1,min=1,max=2^30)

#Particle filter update
SMC_Out<-.C("SMC_CR12",as.numeric(UD[LV-1]),as.numeric(UD[LV-1]-UD[LV]),as.numeric(0),as.numeric(0),as.integer(State_Particles),as.numeric(Obs[LV]),Log_M_L=as.numeric(0),as.integer(Random_Seed),X_1=as.numeric(X_1_Proposed),X_2=as.numeric(X_2_Proposed),T=as.numeric(T_Proposed),w=as.numeric(w_Proposed),as.numeric(Parameter_Proposed[1]),as.numeric(Parameter_Proposed[2]),as.numeric(-Parameter_Proposed[3]),as.numeric(Parameter_Proposed[4]),as.numeric(0),as.numeric(0),as.numeric(0),as.numeric(Parameter_Proposed[5]),as.numeric(Parameter_Proposed[6]),as.numeric(Parameter_Proposed[7]),as.numeric(Parameter_Proposed[8]),as.numeric(Parameter_Proposed[9]),as.numeric(Parameter_Proposed[10]),as.numeric(Parameter_Proposed[11]),as.numeric(Parameter_Proposed[12]),A=as.integer(A_Proposed),as.numeric(GPre),as.numeric(GCop),as.numeric(GObl))

#Update marginal likelihood
LML_Track<-LML_Track+SMC_Out$Log_M_L

#Update particle state, weights, and anscestors
X_1_Proposed<-SMC_Out$X_1
X_2_Proposed<-SMC_Out$X_2
T_Proposed<-SMC_Out$T
w_Proposed<-SMC_Out$w
A_Proposed<-SMC_Out$A

TempMatrixX1<-c(TempMatrixX1,X_1_Proposed)
TempMatrixX2<-c(TempMatrixX2,X_2_Proposed)
TempMatrixT<-c(TempMatrixT,T_Proposed)
TempMatrixw<-c(TempMatrixw,w_Proposed)
TempMatrixA<-c(TempMatrixA,A_Proposed)

#Update array of initial values
Temp_Int_1<-Temp_Int_1[A_Proposed]
Temp_Int_2<-Temp_Int_2[A_Proposed]
Temp_Int_T<-Temp_Int_T[A_Proposed]

} else{
w_Proposed<-0
LML_Track<--Inf
break}}}

#Determine acceptance of the MCMC proposal
if(LML_Track!="NaN"){
u<-runif(1)
if(log(u)<LML_Track-Results[j,19]+dmvnorm(Parameter_Current,mean=Parameter_Mean,sigma=Parameter_Covariance,log=TRUE)-dmvnorm(Parameter_Proposed,mean=Parameter_Mean,sigma=Parameter_Covariance,log=TRUE)+dnorm(Parameter_Proposed[1],mean=0.4,sd=0.3,log=TRUE)-dnorm(Parameter_Current[1],mean=0.4,sd=0.3,log=TRUE)+dnorm(Parameter_Proposed[2],mean=0,sd=0.4,log=TRUE)-dnorm(Parameter_Current[2],mean=0,sd=0.4,log=TRUE)+dexp(Parameter_Proposed[3],rate=1/0.5,log=TRUE)-dexp(Parameter_Current[3],rate=1/0.5,log=TRUE)+dexp(Parameter_Proposed[4],rate=1/0.5,log=TRUE)-dexp(Parameter_Current[4],rate=1/0.5,log=TRUE)+dexp(Parameter_Proposed[5],rate=1/0.3,log=TRUE)-dexp(Parameter_Current[5],rate=1/0.3,log=TRUE)+dexp(Parameter_Proposed[6],rate=1/0.5,log=TRUE)-dexp(Parameter_Current[6],rate=1/0.5,log=TRUE)+dexp(Parameter_Proposed[7],rate=1/0.1,log=TRUE)-dexp(Parameter_Current[7],rate=1/0.1,log=TRUE)+dgamma(Parameter_Proposed[10],shape=180,scale=1/4000000,log=TRUE)-dgamma(Parameter_Current[10],shape=180,scale=1/4000000,log=TRUE)+dexp(Parameter_Proposed[11],rate=500,log=TRUE)-dexp(Parameter_Current[11],rate=500,log=TRUE)+dgamma(Parameter_Proposed[12],shape=10,scale=2,log=TRUE)-dgamma(Parameter_Current[12],shape=10,scale=2,log=TRUE)+dbeta(Parameter_Proposed[13],shape1=45,shape2=15,log=TRUE)-dbeta(Parameter_Current[13],shape1=45,shape2=15,log=TRUE)+dexp(Parameter_Proposed[14],rate=4000,log=TRUE)-dexp(Parameter_Proposed[14],rate=4000,log=TRUE)) {
Results[j,c(1:4,8:17)]<-Parameter_Proposed
Results[j,19]<-LML_Track

TempMatrixX1<-matrix(TempMatrixX1,ncol=i-1)
TempMatrixX2<-matrix(TempMatrixX2,ncol=i-1)
TempMatrixT<-matrix(TempMatrixT,ncol=i-1)
TempMatrixw<-matrix(TempMatrixw,ncol=i-1)
TempMatrixA<-matrix(TempMatrixA,ncol=i-1)

for(l in 1:(i-1)){
X_1[[l]][j,]<-TempMatrixX1[,l]
X_2[[l]][j,]<-TempMatrixX2[,l]
T[[l]][j,]<-TempMatrixT[,l]
w[[l]][j,]<-TempMatrixw[,l]
A[[l]][j,]<-TempMatrixA[,l]}

Int_1[j,]<-Temp_Int_1
Int_2[j,]<-Temp_Int_2
Int_T[j,]<-Temp_Int_T}}}}}}

Results[,18]<-rep(1/Parameter_Particles,Parameter_Particles)

#Write out number of unique particles from diversification
Unique_Record_PostMove[i]<-length(unique(Results[,1]))

write.table(Unique_Record_PostMove[i],file=paste(FileTop,"Unique_Record_PostMove.txt",sep=""),col.names=FALSE,row.names=FALSE,append=TRUE)}

for(j in 1:Parameter_Particles){

#C code does not accept infinite values, so in the unlikely event that they occur (which can happen with some models with poor prior choices, for example), they are replaced with finite values
X_1[[i-1]][which(X_1[[i-1]]==Inf)]<-1e6
X_1[[i-1]][which(X_1[[i-1]]==-Inf)]<--1e6
X_2[[i-1]][which(X_2[[i-1]]==Inf)]<-1e6
X_2[[i-1]][which(X_2[[i-1]]==-Inf)]<--1e6

if(sum(w[[i-1]][j,])>0){

#Draw seed for the particle filter
Random_Seed<-runif(1,min=1,max=2^30)

#Particle filter update
UD<-Obs_Depths+(Results[j,17]/(1-Results[j,16]))*Obs_Depths^2
SMC_Out<-.C("SMC_CR12",as.numeric(UD[i-1]),as.numeric(UD[i-1]-UD[i]),as.numeric(0),as.numeric(0),as.integer(State_Particles),as.numeric(Obs[i]),Log_M_L=as.numeric(0),as.integer(Random_Seed),X_1=as.numeric(X_1[[i-1]][j,]),X_2=as.numeric(X_2[[i-1]][j,]),T=as.numeric(T[[i-1]][j,]),w=as.numeric(w[[i-1]][j,]),as.numeric(Results[j,1]),as.numeric(Results[j,2]),as.numeric(-Results[j,3]),as.numeric(Results[j,4]),as.numeric(0),as.numeric(0),as.numeric(0),as.numeric(Results[j,8]),as.numeric(Results[j,9]),as.numeric(Results[j,10]),as.numeric(Results[j,11]),as.numeric(Results[j,12]),as.numeric(Results[j,13]),as.numeric(Results[j,14]),as.numeric(Results[j,15]),A=as.integer(A[[i-1]][j,]),as.numeric(GPre),as.numeric(GCop),as.numeric(GObl))

#Update arrays
if(SMC_Out$Log_M_L!="NaN"){
X_1[[i]][j,]<-SMC_Out$X_1
X_2[[i]][j,]<-SMC_Out$X_2
T[[i]][j,]<-SMC_Out$T
w[[i]][j,]<-SMC_Out$w
A[[i]][j,]<-SMC_Out$A

Int_1[j,]<-Int_1[j,SMC_Out$A]
Int_2[j,]<-Int_2[j,SMC_Out$A]
Int_T[j,]<-Int_T[j,SMC_Out$A]

Results[j,19]<-Results[j,19]+SMC_Out$Log_M_L
} else{
w[[i]][j,]<-rep(0,State_Particles)
Results[j,19]<--Inf}
} else{
w[[i]][j,]<-rep(0,State_Particles)
Results[j,19]<--Inf}

Likelihood_Vector[j]<-(1/State_Particles)*sum(w[[i]][j,])}

Results[,18]<-Results[,18]*Likelihood_Vector
Marginal_Likelihood_Vector[i]<-sum(Results[,18])
write.table(Marginal_Likelihood_Vector[i],file=paste(FileTop,"Marginal_likelihood_Vector.txt",sep=""),col.names=FALSE,row.names=FALSE,append=TRUE)
Results[,18]<-Results[,18]/sum(Results[,18])

}

#Write out results
write.table(Results,file=paste(FileTop,"Results.txt",sep=""),col.names=c("Beta_0","Beta_1","Beta_2","Delta","P","C","E","Sigma_1","Sigma_2","Sigma_Obs","Displacement","Scale","Sed_Mean","Sed_SD","Alpha","p_0","Grad","Weight","Marginal Likelihood (Log)"),row.names=FALSE)
Evidence<-prod(Marginal_Likelihood_Vector)
write.table(Evidence,file=paste(FileTop,"Evidence.txt",sep=""),col.names=FALSE,row.names=FALSE)


Sam_T<-matrix(rep(0,Parameter_Particles*Obs_Length),ncol=Obs_Length)
Sam_1<-Sam_T
Sam_2<-Sam_T

for(i in 1:1000){
j<-Obs_Length
if(sum(w[[j]][i,])>0){
Sam_W<-w[[j]][i,]/sum(w[[j]][i,])

Sam<-sample(1:1000,size=1,prob=Sam_W)

for(j in Obs_Length:1){
Sam_T[i,j]<-T[[j]][i,Sam]
Sam_1[i,j]<-X_1[[j]][i,Sam]
Sam_2[i,j]<-X_2[[j]][i,Sam]
Sam<-A[[j]][i,Sam]}
} else{
Sam_T[,j]<-rep(NA,1000)
Sam_1[,j]<-rep(NA,1000)
Sam_2[,j]<-rep(NA,1000)}}

write.table(Sam_T,file=paste(FileTop,"Sam_T.txt",sep=""),row.names=FALSE,col.names=FALSE)
write.table(Sam_1,file=paste(FileTop,"Sam_1.txt",sep=""),row.names=FALSE,col.names=FALSE)
write.table(Sam_2,file=paste(FileTop,"Sam_2.txt",sep=""),row.names=FALSE,col.names=FALSE)


