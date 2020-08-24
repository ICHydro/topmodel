#include "topmodel.h"

void c_infiltration(double *rain,
                  double *parameters,
                  int *ntimestep,
                  double *result)
{
  int i;
  double dt,CD,K0,m;

  dt = parameters[0];
  CD = parameters[1];
  K0 = parameters[2];
  m  = parameters[3];

  for(i=0; i<*ntimestep; i++)
    result[i] = dt * get_f((i + 1) * dt, rain[i] / dt, CD, K0, m, dt);

  return;

}
