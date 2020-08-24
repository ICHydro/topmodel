/* Nash Sutcliff efficiency */

#include "topmodel.h"

double get_Eff(double *Qsim, double *Qobs, int ntimestep)
{
  int i;
  double Eff, numerator, denominator, Qobs_mean;
  int	j = 0;

  Qobs_mean = 0.0;
  numerator = 0.0;
  for(i=0; i<ntimestep; i++){
    if(Qobs[i]>=0) {
      Qobs_mean += Qobs[i];
      numerator += pow(Qobs[i] - Qsim[i], 2.0);
      j++;
    }
  }
  Qobs_mean /= j;

  denominator = 0.0;
  for(i=0; i<ntimestep; i++)
    if(Qobs[i]>=0) {
      denominator += pow(Qobs[i] - Qobs_mean, 2.0);
    }

  if(denominator == 0.0){
    Eff = -999999;
  }else{
    Eff = 1.0 - numerator / denominator;
  }
  return Eff;
}

