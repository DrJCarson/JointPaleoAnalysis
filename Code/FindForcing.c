#include <stdio.h>
#include <math.h>
#include <R.h> // R interface
#include <stdlib.h>  // to be able to use the function exit

#include <gsl/gsl_rng.h> // random numbers
#include <gsl/gsl_randist.h>  // random numbers distributions

void FF(double *AMPPRE, double *OMEPRE, double *ANGPRE, double *AMPOBL, double *OMEOBL, double *ANGOBL, double *Time, double *GPRE, double *GCOP, double *GOBL){

  int NAP=30, NAO=30, m;
  double XPRE, XCOP, XOBL, DXDPREDT, DXDOBLDT, SOMEP[NAP], COMEP[NAP], SOMEO[NAO], COMEO[NAO];

    GPRE[0]=0;
    DXDPREDT=0;
    GCOP[0]=0;
    for(m=0;m<NAP;m++){
      SOMEP[m]=sin(OMEPRE[m]*(Time[0]) + ANGPRE[m]);
      COMEP[m]=cos(OMEPRE[m]*(Time[0]) + ANGPRE[m]);
      GPRE[0]=GPRE[0]+AMPPRE[m]*SOMEP[m];
      GCOP[0]=GCOP[0]+AMPPRE[m]*COMEP[m];
      DXDPREDT=DXDPREDT+OMEPRE[m]*AMPPRE[m]*COMEP[m];}
  
    GOBL[0]=0;
    DXDOBLDT=0;
    for(m=0;m<NAO;m++){
      SOMEO[m]=sin(OMEOBL[m]*(Time[0]) + ANGOBL[m]);
      COMEO[m]=cos(OMEOBL[m]*(Time[0]) + ANGOBL[m]);
      GOBL[0]=GOBL[0]+AMPOBL[m]*COMEO[m];
      DXDOBLDT=DXDOBLDT+OMEOBL[m]*AMPOBL[m]*SOMEO[m];}

    
}