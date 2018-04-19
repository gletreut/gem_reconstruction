/* compilation option */
// g++ -O3 -DGSL gem_direct_mapping.cpp -lgsl -lgslcblas -lm
// g++ -O3 -DLAPACK gem_direct_mapping.cpp -lgsl -lsbtk -llapacke -lcblas -llapack -lblas -lgfortran -lm
// g++ -O3 -DMKL gem_direct_mapping.cpp -lgsl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
// * For MKL, the compilation requires to execute before:
//   source /opt/intel/mkl/bin/mklvars.sh intel64

// Compulsory defines
//#define FFACTOR_GAUSS
//#define FFACTOR_THETA

// optional defines
#define VERBOSE
//#define DEBUG

/* standard imports */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <limits>

/* scientific library */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

/* customs imports */
#include "include/sparse_matrix.h"
#include "include/utils_gsl.h"
#include "include/linalg.h"
#include "include/gem_contacts_kij.h"

  /* specify contact functions */
#ifdef FFACTOR_THETA
#include "include/contact_func_theta.cpp"
#elif defined(FFACTOR_GAUSS)
#include "include/contact_func_gauss.cpp"
#endif

using namespace std;

/* functions */
double
dist_to_pos_def(gsl_matrix *mat_)
{
  /*
   * Return the distance to the closest positive definite matrix.
   * The input matrix must be square and symmetrical.
   */
  int N;
  double dist;
  gsl_eigen_symm_workspace *w(0);
  gsl_vector *eval(0);
  gsl_matrix *mat(0);

  /* copy input matrix */
  N = mat_->size1;
  if (N != mat_->size2)
    throw invalid_argument("Input matrix must be square!");
  mat = gsl_matrix_calloc(N,N);
  gsl_matrix_memcpy(mat,mat_);

  /* diagonalize */
  w = gsl_eigen_symm_alloc(N);
  eval = gsl_vector_calloc(N);

  gsl_eigen_symm(mat,eval,w);

  /* compute distance to positive definite matrix ensemble */
  dist = 0.0;
  for (int i=0; i<eval->size; i++){
    double val = gsl_vector_get(eval,i);
    if (val < 0.0)
      dist += val*val;
  }

  /* free memory and return */
  gsl_eigen_symm_free(w);
  gsl_vector_free(eval);
  gsl_matrix_free(mat);
  return sqrt(dist);
}

/* main */
int main(int argc, char *argv []){

//****************************************************************************
//* DECLARATIONS
//****************************************************************************
  int N,npos;
  double a,thres;
  string pathtocmat,pathtokmat,pathtoemat;
  string ffactorcmap,mystr,basename;
  vector<string> parsechain;
  stringstream convert;
  ifstream fin;
  ofstream fcmat, fkmat;

  gsl_matrix *kmat(0),*cmat(0), *emat(0);
  model *mod(0);
  workspace *wkp(0);

//****************************************************************************
//* INITIALIZATIONS
//****************************************************************************

  /* read arguments */
  if ( argc != 6 )
    throw invalid_argument("Syntax: <N> <thres> <emat_in> <cmat_out> <kmat_out>");

  /** import contact map **/{
    convert.clear();
    convert.str(argv[1]);
    convert >> N;

    convert.clear();
    convert.str(argv[2]);
    convert >> thres;

    convert.clear();
    convert.str(argv[3]);
    convert >> pathtoemat;

    convert.clear();
    convert.str(argv[4]);
    convert >> pathtocmat;

    convert.clear();
    convert.str(argv[5]);
    convert >> pathtokmat;
  }

#if defined(DEBUG)
  cout << "ARGUMENTS GIVEN:" << endl;
  printf("%-20s%-d\n","N",N);
  printf("%-20s%-.2f\n","threshold",thres);
  printf("%-20s%-s\n","emat", pathtoemat.c_str());
  printf("%-20s%-s\n","cmat", pathtocmat.c_str());
  printf("%-20s%-s\n","kmat", pathtokmat.c_str());
  printf("\n");
#endif

  /** file outputs **/
  fcmat.open(pathtocmat.c_str());
  fcmat << left << scientific << setprecision(8);
  fkmat.open(pathtokmat.c_str());
  fkmat << left << scientific << setprecision(8);

  /* create matrices */
  emat = gsl_matrix_calloc(N+1,N+1);
  cmat = gsl_matrix_calloc(N+1,N+1);
  kmat = gsl_matrix_calloc(N+1,N+1);

  /* import contact map */
  fin.open(pathtoemat.c_str());
  load_matrix(fin, emat);
  fin.close();

  /* create workspace and model*/
  wkp = new workspace(N);
  mod = new model(emat, thres, contact_func_f, contact_func_df, contact_func_finv, wkp);

  /* Set the model from the input kmat */
  mod->set_cmat(emat);

  gsl_matrix_memcpy(cmat,mod->cmat);
  gsl_matrix_memcpy(kmat,mod->k);

  /* write obtained matrix */
  print_matrix(fcmat,cmat);
  print_matrix(fkmat,kmat);

  /* exit */
  delete mod;
  delete wkp;
  gsl_matrix_free(emat);
  gsl_matrix_free(cmat);
  gsl_matrix_free(kmat);

  return 0;
}

