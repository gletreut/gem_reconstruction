//****************************************************************************
// minimize_thres_norm.cpp - Routine for the reconstruction of a Gaussian Effective Model
// from a given contact probability matrix.
// The arguments are:
//  * the input matrix of contacts.
//  * the largest matrix index N (matrix size is N+1).
//  * the smallest threshold to be used in the GEM mapping.
//  * the largest threshold to be used in the GEM mapping.
//  * the threshold increment.
//  * the global normalization factor to apply to the input matrix of contacts.
//
// Date: 2017-05-03
// Created by: Guillaume Le Treut
//****************************************************************************
/* compilation option */
//#define VERBOSE
//#define DEBUG
//#define FFACTOR_GAUSS

#ifndef CMASK
#define CMASK -1.
#endif

/* standard imports */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cmath>

/* scientific library */
#include <gsl/gsl_matrix.h>

/* customs imports */
#include "include/utils_gsl.h"
#include "include/linalg.h"
#include "include/min_contacts_kij.h"

/* test */
//#include "include/sparse_matrix.cpp"
//#include "include/utils_gsl.cpp"
//#include "include/linalg.cpp"
//#include "include/min_contacts_kij.cpp"
/* test */

/* specify contact functions */
#ifdef FFACTOR_THETA
#include "include/contact_func_theta.cpp"
#elif defined(FFACTOR_GAUSS)
#include "include/contact_func_gauss.cpp"
#endif

using namespace std;

/* add-ons */

int main(int argc, char *argv []){

//****************************************************************************
//* DECLARATIONS
//****************************************************************************
  int N,npos,nthres;
  double a,znorm,thres,thresmin,thresmax,dthres,thresopt,distopt,val;
  string pathtoemat,pathtocmat,pathtokmat,pathtocmat_opt,pathtokmat_opt,pathtosigma,pathtosigma_opt,pathtosigmainv,pathtosigmainv_opt,pathtodist;
  string ffactorcmap,mystr,basename;
  vector<string> parsechain;
  stringstream convert;
  ifstream fin;
  ofstream fdist;

  gsl_matrix *kmat(0),*cmat(0),*kmat_opt(0),*cmat_opt(0), *sigma(0), *sigma_opt(0), *sigmainv(0), *sigmainv_opt(0), *mat1(0), *mat2(0);
  model *mod(0);
  workspace *wkp(0);
  minimization *min(0);

//****************************************************************************
//* INITIALIZATIONS
//****************************************************************************

  /* read arguments */
  if ( argc != 7 )
    throw invalid_argument("Syntax: <cmap> <N> <thresmin> <thresmax> <dthres> <znorm>");

  /** import contact map **/
  {
    convert.clear();
    convert.str(argv[1]);
    convert >> pathtoemat;

    convert.clear();
    convert.str(argv[2]);
    convert >> N;

    convert.clear();
    convert.str(argv[3]);
    convert >> thresmin;

    convert.clear();
    convert.str(argv[4]);
    convert >> thresmax;

    convert.clear();
    convert.str(argv[5]);
    convert >> dthres;

    convert.clear();
    convert.str(argv[6]);
    convert >> znorm;

    if (thresmin > thresmax){
      thres = thresmin;
      thresmin = thresmax;
      thresmax = thres;
    }
  }

#if defined(DEBUG)
  cout << "ARGUMENTS GIVEN:" << endl;
  printf("%-20s%-d\n","N",N);
  printf("%-20s%-.2f\n","znorm",znorm);
  printf("%-20s%-.2f\n","thresmin",thresmin);
  printf("%-20s%-.2f\n","thresmax",thresmax);
  printf("%-20s%-.2f\n","dthres",dthres);
  printf("%-20s%-s\n","pathtoemat",pathtoemat.c_str());
  printf("\n");
#endif

  /* create workspace */
  wkp = new workspace(N);

  /* read input */
  min=new minimization(wkp);
  {
    printf("Enter ftol:\n");
    npos=scanf("%lf", &min->ftol);
    cin.ignore(numeric_limits<std::streamsize>::max(), '\n');
    printf("Enter frtol:\n");
    npos=scanf("%lf", &min->frtol);
    cin.ignore(numeric_limits<std::streamsize>::max(), '\n');
    printf("Enter ktol:\n");
    npos=scanf("%lf", &min->ktol);
    cin.ignore(numeric_limits<std::streamsize>::max(), '\n');
    printf("Enter itermax:\n");
    npos=scanf("%d", &min->itermax);
    cin.ignore(numeric_limits<std::streamsize>::max(), '\n');
    printf("Enter dump:\n");
    npos=scanf("%d", &min->dump);
    cin.ignore(numeric_limits<std::streamsize>::max(), '\n');

#if defined(DEBUG)
    cout << "PARAMETERS READ:" << endl;
    printf("%-20s%-20.6e\n","ftol",min->ftol);
    printf("%-20s%-20.6e\n","frtol",min->frtol);
    printf("%-20s%-20.6e\n","ktol",min->ktol);
    printf("%-20s%-20d\n","itermax",min->itermax);
    printf("%-20s%-20d\n","dump",min->dump);
    printf("\n");
#endif
  }

  /** file outputs **/
  pathtocmat="cmat.dat";
  pathtokmat="kmat.dat";
  pathtosigma="sigma.dat";
  pathtosigmainv="sigmainv.dat";
  pathtocmat_opt="cmat_opt.dat";
  pathtokmat_opt="kmat_opt.dat";
  pathtosigma_opt="sigma_opt.dat";
  pathtosigmainv_opt="sigmainv_opt.dat";

  pathtodist="distances.dat";
  fdist.open(pathtodist.c_str());
  fdist << left << scientific << setprecision(8);

#if defined(DEBUG)
  cout << "OUTPUT:" << endl;
  printf("%-20s%-s\n","cmat", pathtocmat.c_str());
  printf("%-20s%-s\n","kmat", pathtokmat.c_str());
  printf("%-20s%-s\n","sigma", pathtosigma.c_str());
  printf("%-20s%-s\n","sigmainv", pathtosigmainv.c_str());
  printf("%-20s%-s\n","cmat_opt", pathtocmat_opt.c_str());
  printf("%-20s%-s\n","kmat_opt", pathtokmat_opt.c_str());
  printf("%-20s%-s\n","sigma_opt", pathtosigma_opt.c_str());
  printf("%-20s%-s\n","sigmainv_opt", pathtosigmainv_opt.c_str());
  printf("%-20s%-s\n","distances", pathtodist.c_str());
  printf("\n");
#endif

  /* create matrices */
  cmat = gsl_matrix_calloc(N+1,N+1);
  kmat = gsl_matrix_calloc(N+1,N+1);
  sigma = gsl_matrix_calloc(N,N);
  sigmainv = gsl_matrix_calloc(N,N);
  cmat_opt = gsl_matrix_calloc(N+1,N+1);
  kmat_opt = gsl_matrix_calloc(N+1,N+1);
  sigma_opt = gsl_matrix_calloc(N,N);
  sigmainv_opt = gsl_matrix_calloc(N,N);
  mat1 = gsl_matrix_calloc(N+1,N+1);
  mat2 = gsl_matrix_calloc(N,N);

  /* import contact map */
  fin.open(pathtoemat.c_str());
  load_matrix(fin, cmat);
  fin.close();

  /* normalize contact map */
  gsl_matrix_scale(cmat,1.0/znorm);

  /* set to one values to large */
  for (int i=0;i<=N;i++){
    for (int j=i;j<=N;j++){
      val = gsl_matrix_get(cmat,i,j);
      if (val > 1) val=1;
      gsl_matrix_set(cmat,i,j,val);
      gsl_matrix_set(cmat,j,i,val);
    }
  }
  print_matrix(pathtocmat,cmat);

  /* create workspace and model*/
  mod = new model(cmat, thres, contact_func_f, contact_func_df, contact_func_finv, wkp, CMASK);
  printf("%-20s%-lf\n","CMASK",CMASK);
  gsl_matrix_set_all(kmat_opt,0.0);    // initialize at zero
  distopt = 99.9e99;

  for (thres=thresmax; thres > thresmin-1.0e-10 ; thres -= dthres){

    /* set model threshold */
    mod->thres = thres;

    /* update matrices of the model */
    mod->set_cmat(cmat);
    gsl_matrix_memcpy(kmat,mod->k);
    gsl_matrix_memcpy(sigma,mod->xinv);
    gsl_matrix_memcpy(sigmainv,mod->x);

    /* project matrix of the model */
    project_matrix(mod,kmat_opt,min);

    /* get matrices of the projected model */
    mod->set_k(kmat_opt);
    gsl_matrix_memcpy(cmat_opt,mod->cmat);
    gsl_matrix_memcpy(sigma_opt,mod->xinv);
    gsl_matrix_memcpy(sigmainv_opt,mod->x);

    /* compute distances */
    double fc, fc_rel;
    compute_distances(&fc,&fc_rel,mod);

    if (fc < distopt) {
      distopt = fc;
      thresopt = thres;

      /* write matrices of the model */
      print_matrix(pathtokmat,kmat);
      print_matrix(pathtosigma,sigma);
      print_matrix(pathtosigmainv,sigmainv);
      print_matrix(pathtocmat_opt,cmat_opt);
      print_matrix(pathtokmat_opt,kmat_opt);
      print_matrix(pathtosigma_opt,sigma_opt);
      print_matrix(pathtosigmainv_opt,sigmainv_opt);

      /* write threshold */
      FILE *pfout(0);
      pfout = fopen("threshold.dat","w");
      fprintf(pfout,"%-20s%-.6f\n","threshold",thres);
      fclose(pfout);
    }

    /* print distances */
    printf("DISTANCES BETWEEN INPUT AND OUTPUT MATRICES:\n");
    printf("d(c,c_opt)=%.6e\n",fc);
    printf("d_rel(c,c_opt)=%.6e\n",fc_rel);
    //printf("d(k,k_opt)=%.6e\n",fk);
    //printf("d(sigma,sigma_opt)=%.6e\n",fs);
    fdist << setw(20) << thres;
    fdist << setw(20) << fc;
    fdist << setw(20) << fc_rel;
    //fdist << setw(20) << fk;
    //fdist << setw(20) << fs;
    fdist << endl;

    //printf("%-+20.6f%-+20.6e%-+20.6e%-+20.6e%-+20.6e\n", mod->thres, fc, fc_rel, fk, fs);
    printf("%-+20.6f%-+20.6e%-+20.6e\n", mod->thres, fc, fc_rel);
  }

  /* exit */
  delete min;
  delete mod;
  delete wkp;
  gsl_matrix_free(cmat);
  gsl_matrix_free(kmat);
  gsl_matrix_free(sigma);
  gsl_matrix_free(sigmainv);
  gsl_matrix_free(cmat_opt);
  gsl_matrix_free(kmat_opt);
  gsl_matrix_free(sigma_opt);
  gsl_matrix_free(sigmainv_opt);
  gsl_matrix_free(mat1);
  gsl_matrix_free(mat2);

  return 0;
}

