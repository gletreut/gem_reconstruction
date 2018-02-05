#ifndef INVERSE_KIJ_H
#define INVERSE_KIJ_H

/* standard library */
#include <stdio.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <limits>

/* scientific library */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

/* custom files */
#include "sparse_matrix.h"
#include "linalg.h"
#include "utils_gsl.h"

/* workspace */
struct workspace_struct
{
  /*
   * define the problem workspace, that keeps a pointer to each array
   * in the program.
   */

  /* integer */
  /** integer N: size of the problem **/
  int N;
  int seed;

  /* dense matrices */
  /** size NxN **/
  gsl_matrix *trid,*w,*x,*xinv,*workN_1,*workN_2;
  /** size (N+1)x(N+1) **/
  gsl_matrix *k,*cmat,*gammat,*emat,*workN1_1,*workN1_2;
  gsl_matrix *g0,*k1,*g1,*dk;
  gsl_matrix *ac; // active constraints

  /* sparse matrices */
  SparseMat *B;
  SparseMat *BI;
  SparseMat *C;
  SparseMat *CI;

  /* random number generator */
  gsl_rng *rng;

  workspace_struct(int N_);
  ~workspace_struct();
  void init();
  void alloc(int N_);
  void free();
};

typedef struct workspace_struct workspace;

/* model */
struct model_struct
{
  /*
   * interface to the quantities of the problem, using
   * an input workspace with allocated arrays.
   * Note that:
   *  o NOTHING shall be allocated here because everything must goes through
   *  the workspace.
   */

  int N;
  double thres;

  /* attributes */
  /** workspace **/
  workspace *work;

  /** size NxN **/
  gsl_matrix *x,*xinv,*trid,*w;  // named quantities

  /** size (N+1)x(N+1) **/
  gsl_matrix *k,*cmat,*gammat,*emat;

  /** sparse matrices **/
  SparseMat *B;      // w=B k
  SparseMat *BI;     // k=BI w
  SparseMat *C;      // gamma=C xinv
  SparseMat *CI;      // xinv=CI gamma

  /** form factor function **/
  double (*func)(double,void*);
  double (*dfunc)(double,void*);
  double (*funcinv)(double,void*);

  /* constructor */
  model_struct(gsl_matrix *eij_, double thres_, double (*func_)(double,void*), double (*dfunc_)(double,void*), double (*funcinv_)(double,void*), workspace *work_);

  /* built-in functions */
  /** set the model from one input matrix **/
  void set_k(const gsl_matrix *kk);
  void set_x(const gsl_matrix *xx);
  void set_xinv(const gsl_matrix *xxinv);
  void set_gammat(const gsl_matrix *gammat);
  void set_cmat(const gsl_matrix *ccmat);

  /** passage operations from one matrix to another **/
  void k2w(const gsl_matrix *kk, gsl_matrix *ww);
  void w2k(const gsl_matrix *ww, gsl_matrix *kk);
  void w2x(const gsl_matrix *ww, gsl_matrix *xx);
  void x2w(const gsl_matrix *xx, gsl_matrix *ww);
  void x2xinv(const gsl_matrix *xx, gsl_matrix *xxinv);
  void xinv2x(const gsl_matrix *xxinv, gsl_matrix *xx);
  void xinv2gam(const gsl_matrix *xxinv, gsl_matrix *ggammat);
  void gam2xinv(const gsl_matrix *ggammat, gsl_matrix *xxinv);
  void gam2cmat(const gsl_matrix *ggammat, gsl_matrix *ccmat);
  void cmat2gam(const gsl_matrix *ccmat, gsl_matrix *ggammat);

  /** print **/
  void print_infos();

};

typedef struct model_struct model;

/* minimization */
struct minimization_struct
{
  /*
   * define the problem workspace, that keeps a pointer to each array
   * in the program.
   */

  /* attributes */
  /** workspace **/
  workspace *work;

  /** integer **/
  int N,nac;
  int itermax, converged, dump;

  /** double **/
  double step,step_max,tol;
  double gnorm,f0,k0norm,df,dfrel,dknorm;
  double ftol, frtol, ktol;

  /** matrices **/
  gsl_matrix *k1,*g1,*dk;
  gsl_matrix *ac; // active constraints

  /* constructor */
  minimization_struct(workspace *wkp);

  /* functions */
  double (* L) (const gsl_matrix *k, void *params);
  void (* dL) (const gsl_matrix *k, void *params, gsl_matrix *gradient);
  void (* LdL) (const gsl_matrix *k, void *params, double *ll, gsl_matrix *gradient);
  void iterate(gsl_matrix *k, void *params, double *f, gsl_matrix *gradient);
  void test_convergence();
  void project_k(gsl_matrix *k,bool all=false);
  void get_largest_constraint(gsl_matrix *k, size_t *p, size_t *q, gsl_matrix *kwork);
  void min_J(gsl_matrix *k, double *f, void *params);
  int get_nac();
  void add_active_constraint(int p, int q);
  void add_active_constraint_rng();
  void print_infos();
};

typedef struct minimization_struct minimization;

/* projection function */
void project_matrix(model *mod, gsl_matrix *k_init, minimization *min);
void project_matrix(gsl_matrix *emat, gsl_matrix *k_init, double thres, double (*func)(double, void*),  double (*dfunc)(double, void*), double (*funcinv) (double, void *), minimization *min);

/* function to minimize */
double J(const gsl_matrix *kk, void *params);
void dJ(const gsl_matrix *kk, void *params, gsl_matrix *gradient);
void dJ_finite_diff(const gsl_matrix *kk, void *params, gsl_matrix *gradient);
void JdJ(const gsl_matrix *kk, void *params, double *jj, gsl_matrix *gradient);

/* useful functions */
void make_K2W(SparseMat *B, int N);
void make_W2K(SparseMat *BI, int N);
void make_S2G(SparseMat *C, int N);
void make_G2S(SparseMat *CI, int N);
void make_chain_structure(gsl_matrix *m);

#endif
