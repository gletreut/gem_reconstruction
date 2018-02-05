#include <cmath>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#ifndef PI
#define PI 3.1415926535897931
#endif

#ifndef SQRT_TWO
#define SQRT_TWO 1.4142135623730951
#endif

#ifndef ERF_PREF
#define ERF_PREF 7.9788456080286541e-01 //ERF_PREF = sqrt(2.0/PI);
#endif

#ifndef ITERMAX
#define ITERMAX 100
#endif

#ifndef EPSABS
#define EPSABS 0
#endif

#ifndef EPSREL
#define EPSREL 1.0e-3
#endif

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
  x = thres / a * sqrt(3.0 / gamma);
  x2 = x*x;

  if (x < SMALL)
		return 0.0;
	else if ( (gamma == 0.0) ||  (x > BIG))
		return 1.0;
	else
		return gsl_sf_erf(x/SQRT_TWO)-ERF_PREF*x*exp(-x2*0.5);
}

double
contact_func_df(double gamma, void *params){
  /*
   * derivative of f.
   */
  double a, thres, x, x2;

  a = 1.0;
  thres = *((double*)(params));
  x = thres / a * sqrt(3.0 / gamma);
  x2 = x*x;

	if (x < SMALL)
		return 0.0;
	else if ( (gamma == 0.0) ||  (x > BIG))
		return 0.0;
	else
		return ERF_PREF*x2*exp(-x2*0.5);
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
contact_dfunc(double gamma, void *params){
  /*
   * return the contact function at gamma minus some input
   * contact probability.
   */
  double* par;
  double thres, c0;

  par = (double*)(params);
  thres = par[0];
  c0 = par[1];

  return contact_func_df(gamma, &thres);
}

void
contact_funcdfunc(double gamma, void *params, double *func, double *dfunc){
  /*
   * return both the evaluation at point gamma of both:
   *  (i)   func
   *  (ii)  dfunc
   */

  *func = contact_func(gamma,params);
  *dfunc = contact_dfunc(gamma,params);
  return;
}

double
contact_func_finv(double c, void *params){
  /*
   * Return the gamma corresponding to the input contact probability c.
   */
  int status,powMAX;
  double a, thres, gamma, gamma_lo, gamma_hi, fval;
  double *par(0);
  gsl_root_fsolver *s(0);
  gsl_function F;
  const gsl_root_fsolver_type * T = gsl_root_fsolver_brent;

  if (c == 0.0)
    return BIG;

  if (c == 1.0)
    return 0.0;

  a = 1.0;
  thres = *((double*)(params));
  par = new double[2];
  par[0] = thres;
  par[1] = c;
  powMAX=10;

  /* initialize gsl function */
  F.function = &contact_func;
  F.params = par;

  /* initialize root solver */
  s = gsl_root_fsolver_alloc (T);
  gamma_lo = 0.0;
  gamma_hi = thres;
  for (int iter=0;iter<powMAX;iter++){
    fval = contact_func(gamma_hi,par);
    if ( fval < 0.0 )
      break;
    else
      gamma_hi *= 10.0;
  }
  if (fval > 0.0)
    return BIG;

  gsl_root_fsolver_set(s,&F,gamma_lo,gamma_hi);

  /* look for root */
  for (int iter=0;iter<ITERMAX;iter++)
  {
    status = gsl_root_fsolver_iterate (s);
    gamma = gsl_root_fsolver_root (s);
    gamma_lo = gsl_root_fsolver_x_lower (s);
    gamma_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (gamma_lo, gamma_hi, EPSABS, EPSREL);

    if (status == GSL_SUCCESS)
      break;
  }

  /* free memory */
  gsl_root_fsolver_free(s);
  delete par;

  return gamma;
}

