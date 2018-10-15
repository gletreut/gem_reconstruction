//****************************************************************************
// contact_func_gauss.cpp - Contact probability function with a Gaussian form
// factor
// Date: 2017-05-03
// Created by: Guillaume Le Treut
//****************************************************************************

#include <cmath>

#ifndef BIG
#define BIG 1.0e30
#endif

#ifndef SMALL
#define SMALL 1.0e-15
#endif

double
contact_func_f(double gamma, void *params){
  /*
   * return the contact probability as a function of gamma.
   */
  double a, thres, x, x2;

  a = 1.0;
  thres = *((double*)(params));
  x = a / thres;
  x2 = x*x;

  if (gamma == 0.0)
    return 1.0;
  else
    return pow(1.0+gamma*x2,-1.5);
}

double
contact_func_df(double gamma, void *params){
  /*
   * derivative of f.
   */
  double a, thres, x, x2;

  a = 1.0;
  thres = *((double*)(params));
  x = a / thres;
  x2 = x*x;

  if (gamma == 0.0)
    return 1.0;
  else
    return -1.5*x2*pow(1.0+gamma*x2,-2.5);
}

double
contact_func(double gamma, void *params){
  /*
   * return the contact function at gamma minus some input
   * contact probability.
   */
  double* par;
  double thres, c0;

  par = (double*)(params);
  thres = par[0];
  c0 = par[1];

  return contact_func_f(gamma, &thres) - c0;
}

double
contact_func_finv(double c, void *params){
  /*
   * Return the gamma corresponding to the input contact probability c.
   */
  double a, thres, x, x2;

  a = 1.0;
  thres = *((double*)(params));
  x = a / thres;
  x2 = x*x;

  if (c == 0.0)
    return BIG;
  else if (c == 1.0)
    return 0.0;
  else
    return (pow(c,-2.0/3.0)-1.0) / x2;
}

