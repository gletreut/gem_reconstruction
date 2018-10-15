//****************************************************************************
// linalg.h - Various functions to perform linear algebra operations
// Date: 2017-05-03
// Created by: Guillaume Le Treut
//****************************************************************************

#ifndef LINALG_H
#define LINALG_H

/* standard library */
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <limits>

/* scientific library */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#ifdef GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_linalg.h>
#elif defined(LAPACK)
#include <cblas.h>
#include <spblas.h>
#include <lapacke.h>
#elif defined(MKL)
#include "mkl.h"
#endif

/* custom files */
#include "sparse_matrix.h"

/* linear algebra */
double  linalg_dnrm2   (const gsl_vector *v);
double  linalg_dnrm2   (const gsl_matrix *m);
double  linalg_ddot    (const gsl_vector *a, const gsl_vector *b);
double  linalg_ddot    (const gsl_matrix *a, const gsl_matrix *b);
void    linalg_dscal   (const double s, gsl_vector *a);
void    linalg_dscal   (const double s, gsl_matrix *a);
void    linalg_daxpy   (const double alpha, const gsl_vector *a,
                        gsl_vector *b);
void    linalg_daxpy   (const double alpha, const gsl_matrix *a,
                        gsl_matrix *b);
void    linalg_dgemm  (int transA, int transB,
                       const double alpha, const gsl_matrix *A,
                       const gsl_matrix *B, const double beta,
                       gsl_matrix *C);
void    linalg_dgemv  (int transA, const double alpha,
                       const gsl_matrix *A, const gsl_vector *v,
                       const double beta, gsl_vector *w);
void    linalg_spdgemv (const int transa, const double alpha,
                        const  SparseMat *A, const gsl_vector *v,
                        const double beta, gsl_vector *w);
void    linalg_invert (const gsl_matrix *a, gsl_matrix *ainv);

#endif
