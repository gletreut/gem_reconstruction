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
   *  o pointers need not be allocated here.
   *  o NOTHING shall be allocated here because everything must goes through
   *  the workspace.
   *  o pointers can be initialized to other values.
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

/* useful functions */
void make_K2W(SparseMat *B, int N);
void make_W2K(SparseMat *BI, int N);
void make_S2G(SparseMat *C, int N);
void make_G2S(SparseMat *CI, int N);
void make_chain_structure(gsl_matrix *m);

#endif
