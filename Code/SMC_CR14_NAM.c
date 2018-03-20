#include <stdio.h>
#include <math.h>
#include <R.h> // R interface
#include <stdlib.h>  // to be able to use the function exit

#include <gsl/gsl_rng.h> // random numbers
#include <gsl/gsl_randist.h>  // random numbers distributions


void SMC_CR12(double *Obs_Depth, double *Delta_S, double *Target_Depth, double *Target_Time, int *Particles, double *Obs, double *LML, int *Random_Seed, double *X_1, double *X_2, double *X_3, double *w, double *Beta_0, double *Beta_1, double *Beta_2, double *Delta, double *Pre, double *Cop, double *Obl, double *Sigma_1, double *Sigma_2, double *Sigma_Obs, double *Dis, double *Sca, double *Mu_S, double *SD_S, double *Alpha, int *Hist, double *GPre, double *GCop, double *GObl){

//Declare variables
int i, j, k, sign1, sign2, u; 
int NAP=30, NAO=30;
double U[Particles[0]], V[Particles[0]], Time[Particles[0]], W[Particles[0]], Delta_T[Particles[0]], End_Time[Particles[0]];
double x, xs, t1, t2, y1, y2, Num, Den, URV, Sum_Weight, CW, F, Int_Time, Int_Step, G_Ran, G_Ran_2, S_Ran;
double a1=0.254829592, a2=-0.284496736, a3=1.421413741, a4=-1.453152027, a5=1.061405429, p=0.3275911, LB, UB;
double New_Mean, Euler_Mean_u, Euler_Mean_v, Euler_SD_u, Euler_SD_v, GW_Mean_u, GW_Mean_v, GW_SD_u, GW_SD_v, Prop_Mean, Prop_SD;
int Fa, Fb;
double Fc, FPre, FCop, FObl;
double ds1=Delta_S[0],ds2=((Obs_Depth[0]-Delta_S[0])-Target_Depth[0]),ds3=(Obs_Depth[0]-Target_Depth[0]);
double M_1=ds1/Mu_S[0],M_2=ds2/Mu_S[0],M_3=ds3/Mu_S[0];
double L_1=(ds1/SD_S[0])*(ds1/SD_S[0]),L_2=(ds2/SD_S[0])*(ds2/SD_S[0]),L_3=(ds3/SD_S[0])*(ds3/SD_S[0]);
double TR, DTR;

//Initialize RNG
const gsl_rng_type * T;
gsl_rng * r;
T = gsl_rng_default;
r = gsl_rng_alloc (T);
gsl_rng_set(r, Random_Seed[0]);

//Normalize weights
Sum_Weight=0;
for(i=0; i<Particles[0]; i++){
  Sum_Weight+=w[i];}
  
for(i=0; i<Particles[0]; i++){
  W[i]=w[i]/Sum_Weight;}

//Resample
for(i=0; i<Particles[0]; i++){
  URV=gsl_rng_uniform (r);
  CW=W[0];
  k=0;
  while(CW<URV){
    k++;
    CW+=W[k];}
  Hist[i]=k+1;
  U[i]=X_1[k];
  V[i]=X_2[k];
  Time[i]=X_3[k];}

//Sample times for each particle
for(i=0; i<Particles[0]; i++){
  W[i]=1;
  
  if(ds2==0){
    Delta_T[i]=Target_Time[0]-Time[i];}
  else{
      G_Ran=gsl_ran_gaussian (r, 1);
      G_Ran_2=G_Ran*G_Ran;
      
      S_Ran=M_1+((M_1*M_1*G_Ran_2)/(2*L_1))-((M_1*sqrt(4*M_1*L_1*G_Ran_2+M_1*M_1*G_Ran_2*G_Ran_2))/(2*L_1));
    
      URV=gsl_rng_uniform (r);
      if(URV<=(M_1/(M_1+S_Ran))){
	Delta_T[i]=S_Ran;}
      else{
	Delta_T[i]=(M_1*M_1)/S_Ran;}
	
      TR=Target_Time[0]-Time[i];
      DTR=TR-Delta_T[i];
      
      if(DTR<=0){
	W[i]=0;}
      else{
	W[i]=W[i]*sqrt((L_2)/(2*PI*DTR*DTR*DTR))*exp(-((L_2*(DTR-M_2)*(DTR-M_2))/(2*M_2*M_2*DTR)));
	W[i]=W[i]/(sqrt((L_3)/(2*PI*TR*TR*TR))*exp(-((L_3*(TR-M_3)*(TR-M_3))/(2*M_3*M_3*TR))));}}
	
  End_Time[i]=Time[i]+Delta_T[i];}
  
//Simulate Climate for each particle
for(i=0; i<Particles[0]; i++){
Int_Time=Time[i];

if(W[i]>0){

  do{

  //Set Euler-Maruyama step size
  Int_Step=100;
  if(End_Time[i]-Int_Time<100){
    Int_Step=End_Time[i]-Int_Time;}

  //Determine orbital parameters using linear interpolation
  Fa=(int)(Int_Time/100);
  Fb=Fa-1;
  Fc=Int_Time-100*Fb;
  if(Fc==0){
    FPre=GPre[10000+Fb];
    FCop=GCop[10000+Fb];
    FObl=GObl[10000+Fb];}
  else if(Fc==100){
    FPre=GPre[10000+Fa];
    FCop=GCop[10000+Fa];
    FObl=GObl[10000+Fa];}
  else{
  FPre=GPre[10000+Fb]+(Fc/100)*(GPre[10000+Fa]-GPre[10000+Fb]);
  FCop=GCop[10000+Fb]+(Fc/100)*(GCop[10000+Fa]-GCop[10000+Fb]);
  FObl=GObl[10000+Fb]+(Fc/100)*(GObl[10000+Fa]-GObl[10000+Fb]);}

  F=Pre[0]*FPre+Cop[0]*FCop+Obl[0]*FObl;
    
  //Calculate target mean and standard deviation
  Euler_Mean_u=U[i]-(Int_Step/10000)*(Beta_0[0]+Beta_1[0]*U[i]+Beta_2[0]*(U[i]-U[i]*U[i]*U[i])+Delta[0]*V[i]+F);
  Euler_Mean_v=V[i]+(Int_Step/10000)*Alpha[0]*Delta[0]*(U[i]+V[i]-(V[i]*V[i]*V[i])/3);
  Euler_SD_u=sqrt(Int_Step/10000)*Sigma_1[0];
  Euler_SD_v=sqrt(Int_Step/10000)*Sigma_2[0];

  //Calculate proposal mean and standard deviation
  GW_Mean_u=Euler_Mean_u+((Sigma_1[0]*Sigma_1[0]*(Int_Step/10000))/(Sigma_1[0]*Sigma_1[0]*((End_Time[i]-Int_Time)/10000)+((Sigma_Obs[0]*Sigma_Obs[0])/(Sca[0]*Sca[0]))))*((1/Sca[0])*(Obs[0]-Dis[0])-(U[i]-((End_Time[i]-Int_Time)/10000)*(Beta_0[0]+Beta_1[0]*U[i]+Beta_2[0]*(U[i]-U[i]*U[i]*U[i])+Delta[0]*V[i]+F)));
  GW_Mean_v=Euler_Mean_v;
  GW_SD_u=Euler_SD_u*sqrt((Sigma_1[0]*Sigma_1[0]*((End_Time[i]-(Int_Time+Int_Step))/10000)+((Sigma_Obs[0]*Sigma_Obs[0])/(Sca[0]*Sca[0])))/(Sigma_1[0]*Sigma_1[0]*((End_Time[i]-Int_Time)/10000)+((Sigma_Obs[0]*Sigma_Obs[0])/(Sca[0]*Sca[0]))));
  GW_SD_v=Euler_SD_v;

  //Propose new state
  U[i]=GW_Mean_u+GW_SD_u*gsl_ran_gaussian (r, 1);
  V[i]=GW_Mean_v+GW_SD_v*gsl_ran_gaussian (r, 1);

  //Update importance weight
  W[i]=W[i]*(GW_SD_u/Euler_SD_u)*exp((1/(2*GW_SD_u*GW_SD_u))*(U[i]-GW_Mean_u)*(U[i]-GW_Mean_u)-(1/(2*Euler_SD_u*Euler_SD_u))*(U[i]-Euler_Mean_u)*(U[i]-Euler_Mean_u));

  //Update current time
  Int_Time=Int_Time+Int_Step;

  } while(Int_Time!=End_Time[i]);

//Assimilate observation
New_Mean=Sca[0]*U[i]+Dis[0];
W[i] = W[i]*(1/sqrt(2*PI*Sigma_Obs[0]*Sigma_Obs[0]))*exp(-(1/(2*Sigma_Obs[0]*Sigma_Obs[0]))*(New_Mean-Obs[0])*(New_Mean-Obs[0]));

}}

//Estimate the log marginal likelihood
Sum_Weight=0;
for(i=0; i<Particles[0]; i++){
  Sum_Weight+=W[i];}
      
LML[0]+=log((1/(double)Particles[0])*Sum_Weight);
      
//Update particle states and weights
for(i=0; i<Particles[0]; i++){
  X_1[i]=U[i];
  X_2[i]=V[i];
  X_3[i]=End_Time[i];
  w[i]=W[i];}
    
//Free RNG memory
gsl_rng_free (r);}
