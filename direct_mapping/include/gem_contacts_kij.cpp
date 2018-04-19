/* compile with: g++ main.cpp $(gsl-config --cflags --libs) */

#include "gem_contacts_kij.h"
#include <ctime>

using namespace std;

int dump=1;
const double macheps=numeric_limits<double>::epsilon();
int NOTRANS = 0;
int TRANS   = 1;

/* useful functions */
void
make_K2W(SparseMat *B, int N){
  /*
   * Generate the matrix of passage for the change of variable:
   *  w=B k,
   * where:
   *  k is to be understood as a matrix of size (N+1)^2 x 1,
   *  w is to be understood as a matrix of size N^2 x 1,
   *  B is to be understood as a matrix of size N^2 x (N+1)^2.
   */

  int mu1,mu2;
  double val;

  if ( (B->nrow != N*N) || (B->ncol != (N+1)*(N+1)) )
    throw invalid_argument("B must be of size N^2 x (N+1)^2!");

  B->set_zero();
  /* fill-in B line per line */
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      mu1=i*N+j;


      /* case 1: i == j */
      if (i==j){
        //p=i+1;
        for (int q=0;q<=N;q++){
          if (q==i+1)
            continue;

          mu2=(i+1)*(N+1)+q;
          val=0.5;
          B->set(mu1,mu2,val);
        }
        //q=i+1;
        for (int p=0;p<=N;p++){
          if (p==i+1)
            continue;
          mu2=p*(N+1)+(i+1);
          val=0.5;
          B->set(mu1,mu2,val);
        }
      }

      /* case 2: i <> i */
      else {
        int p,q;

        // subcase 1
        p=i+1;
        q=j+1;
        val=-0.5;
        mu2=p*(N+1)+q;
        B->set(mu1,mu2,val);

        // subcase 2
        p=j+1;
        q=i+1;
        val=-0.5;
        mu2=p*(N+1)+q;
        B->set(mu1,mu2,val);
      }

    }
  }

  /* exit */
  return;
}

void
make_W2K(SparseMat *BI, int N){
  /*
   * Generate the matrix of passage for the change of variable:
   *  k=BI w,
   * where:
   *  k is to be understood as a matrix of size (N+1)^2 x 1,
   *  w is to be understood as a matrix of size N^2 x 1,
   *  BI is to be understood as a matrix of size (N+1)^2 x N^2.
   */

  int mu1,mu2;
  double val;

  if ( (BI->nrow != (N+1)*(N+1)) || (BI->ncol != N*N) )
    throw invalid_argument("BI must be of size (N+1)^2 x N^2!");

  BI->set_zero();
  /* fill-in BI line per line */
  for (int p=0;p<=N;p++){
    for (int q=0;q<=N;q++){
      if (p  == q){
        // val = 0
        continue;
      }

      mu1=p*(N+1)+q;

      /* case 1: p == 0 */
      if (p==0){
        //i=q-1;
        for (int j=0;j<N;j++){
          if (j == q-1)
            continue;
          mu2=(q-1)*N+j;
          val=0.5;
          BI->set(mu1,mu2,val);
        }

        //j=q-1;
        for (int i=0;i<N;i++){
          if (i == q-1)
            continue;
          mu2=i*N+(q-1);
          val=0.5;
          BI->set(mu1,mu2,val);
        }

        //i=j=q-1
        {
          mu2=(q-1)*N+(q-1);
          val=1.0;
          BI->set(mu1,mu2,val);
        }

      }

      /* case 2: q == 0 */
      else if (q==0){
        //i=p-1;
        for (int j=0;j<N;j++){
          if (j == p-1)
            continue;
          mu2=(p-1)*N+j;
          val=0.5;
          BI->set(mu1,mu2,val);
        }

        //j=p-1;
        for (int i=0;i<N;i++){
          if (i == p-1)
            continue;
          mu2=i*N+(p-1);
          val=0.5;
          BI->set(mu1,mu2,val);
        }

        //i=j=p-1
        {
          mu2=(p-1)*N+(p-1);
          val=1.0;
          BI->set(mu1,mu2,val);
        }
      }

      /* case 3: (p <> 0) ^ (q <> 0)*/
      else {
        int i,j;

        // subcase 1
        i=p-1;
        j=q-1;
        val=-0.5;
        mu2=i*N+j;
        BI->set(mu1,mu2,val);

        // subcase 2
        i=q-1;
        j=p-1;
        val=-0.5;
        mu2=i*N+j;
        BI->set(mu1,mu2,val);
      }
    }
  }

  /* exit */
  return;
}

void
make_S2G(SparseMat *C, int N){
  /*
   * Generate the matrix of passage for the change of variable:
   *  \Gamma=C \Sigma,
   * where:
   *  \Sigma is to be understood as a matrix of size N^2 x 1,
   *  \Gamma is to be understood as a matrix of size (N+1)^2 x 1,
   *  C is to be understood as a matrix of size (N+1)^2 x N^2.
   */

  int mu1,mu2;
  int i,j;
  double val;

  if ( (C->nrow != (N+1)*(N+1)) || (C->ncol != N*N) )
    throw invalid_argument("C must be of size (N+1)^2 x N^2!");

  C->set_zero();
  /* fill-in C line per line */
  for (int p=0;p<N+1;p++){
    for (int q=0;q<N+1;q++){
      mu1=p*(N+1)+q;

      /* case 1: p == q */
      if (p==q)
        continue;

      /* case 2: p <> q */
      else {
        // subcase 1
        i=p-1;
        j=q-1;
        if ( !(i<0) && !(j<0) ){
          val=-2.0;
          mu2=i*N+j;
          C->set(mu1,mu2,val);
        }

        // subcase 2
        i=p-1;
        j=p-1;
        if ( !(i<0) && !(j<0) ){
          val=1.0;
          mu2=i*N+j;
          C->set(mu1,mu2,val);
        }

        // subcase 3
        i=q-1;
        j=q-1;
        if ( !(i<0) && !(j<0) ){
          val=1.0;
          mu2=i*N+j;
          C->set(mu1,mu2,val);
        }
      }
    }
  }

  /* exit */
  return;
}

void
make_G2S(SparseMat *CI, int N){
  /*
   * Generate the matrix of passage for the change of variable:
   *  \Sigma=CI \Gamma,
   * where:
   *  \Sigma is to be understood as a matrix of size N^2 x 1,
   *  \Gamma is to be understood as a matrix of size (N+1)^2 x 1,
   *  CI is to be understood as a matrix of size N^2 x (N+1)^2.
   */

  int mu1,mu2;
  int p,q;
  double val;

  if ( (CI->nrow != N*N) || (CI->ncol != (N+1)*(N+1)) )
    throw invalid_argument("CI must be of size N^2 x (N+1)^2!");

  CI->set_zero();
  /* fill-in CI line per line */
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      mu1=i*N+j;

      /* case 1: (p == 0) ^ (q = i+1) */
      p = 0;
      q = i+1;
      mu2 = p*(N+1)+q;
      val = 0.5;
      CI->set(mu1,mu2,val);

      /* case 1: (p == j+1) ^ (q = 0) */
      p = j+1;
      q = 0;
      mu2 = p*(N+1)+q;
      val = 0.5;
      CI->set(mu1,mu2,val);

      /* case 3: (p = i+1) ^ (q = j+1) */
      p = i+1;
      q = j+1;
      if ( p != q){
        mu2 = p*(N+1)+q;
        val = -0.5;
        CI->set(mu1,mu2,val);
      }
    }
  }

  /* exit */
  return;
}

void
make_chain_structure(gsl_matrix *m){
  /* return structure of the bare Gaussian chain */
  int N=m->size1;

  if ( !(N>0) )
    throw invalid_argument("Matrix size must be greater than 0!");

  if ( m->size2 != N)
    throw invalid_argument("Matrix must be square!");

  /* set diagonal elements */
  for (int i=0;i<N-1;i++){
    gsl_matrix_set(m,i,i,2.0);
  }
  gsl_matrix_set(m,N-1,N-1,1.0);

  /* set second diagonal */
  for (int i=1;i<N;i++){
    gsl_matrix_set(m,i,i-1,-1.0);
    gsl_matrix_set(m,i-1,i,-1.0);
  }

  return;
}

/* workspace */
workspace_struct::workspace_struct(int N_) {
  alloc(N_);
  init();
}

workspace_struct::~workspace_struct(){
  free();
}

void
workspace_struct::init(){
  /* random number generator */
  seed=123;
  gsl_rng_set(rng,seed);

  /* Linalg work variables */
  return;
}

void
workspace_struct::alloc(int N_){
  N=N_;
  /* dense matrices */
  /** size NxN **/
  x=gsl_matrix_calloc(N,N);
  xinv=gsl_matrix_calloc(N,N);
  trid=gsl_matrix_calloc(N,N);
  w=gsl_matrix_calloc(N,N);
  workN_1=gsl_matrix_calloc(N,N);
  workN_2=gsl_matrix_calloc(N,N);
  /** size (N+1)x(N+1) **/
  k=gsl_matrix_calloc(N+1,N+1);
  cmat=gsl_matrix_calloc(N+1,N+1);
  gammat=gsl_matrix_calloc(N+1,N+1);
  emat=gsl_matrix_calloc(N+1,N+1);
  workN1_1=gsl_matrix_calloc(N+1,N+1);
  workN1_2=gsl_matrix_calloc(N+1,N+1);
  g0=gsl_matrix_calloc(N+1,N+1);
  k1=gsl_matrix_calloc(N+1,N+1);
  g1=gsl_matrix_calloc(N+1,N+1);
  dk=gsl_matrix_calloc(N+1,N+1);
  ac=gsl_matrix_calloc(N+1,N+1);

  /* sparse matrices */
  B=new SparseMat(N*N,(N+1)*(N+1));
  BI=new SparseMat((N+1)*(N+1),N*N);
  C=new SparseMat((N+1)*(N+1),N*N);
  CI=new SparseMat(N*N,(N+1)*(N+1));

  /* random number generator */
  rng=gsl_rng_alloc (gsl_rng_taus);

  //*/
  return;
}

void
workspace_struct::free(){
  /* dense matrices */
  gsl_matrix_free(x);
  gsl_matrix_free(xinv);
  gsl_matrix_free(trid);
  gsl_matrix_free(w);
  gsl_matrix_free(workN_1);
  gsl_matrix_free(workN_2);
  gsl_matrix_free(k);
  gsl_matrix_free(cmat);
  gsl_matrix_free(emat);
  gsl_matrix_free(gammat);
  gsl_matrix_free(g0);
  gsl_matrix_free(workN1_1);
  gsl_matrix_free(workN1_2);
  gsl_matrix_free(k1);
  gsl_matrix_free(g1);
  gsl_matrix_free(dk);
  gsl_matrix_free(ac);

  /* sparse matrices */
  delete B;
  delete BI;
  delete C;
  delete CI;

  /* random number generator */
  gsl_rng_free(rng);

  return;
}

/* model */
model_struct::model_struct(gsl_matrix *emat_, double thres_, double (*func_)(double,void*), double (*dfunc_)(double,void*), double (*funcinv_)(double,void*), workspace *work_) {
  /* copy threshold */
  thres=thres_;

  /* copy form factor functors */
  func = func_;
  dfunc = dfunc_;
  funcinv = funcinv_;

  /* copy workspace address */
  work=work_;

  /* and copy some references for further referencing */
  /** size of matrix **/
  N=work->N;

  /** useful matrices of the problem **/
  /*** NxN ***/
  x=work->x;
  xinv=work->xinv;
  trid=work->trid;
  w=work->w;

  /*** NxN ***/
  k=work->k;
  cmat=work->cmat;
  gammat=work->gammat;
  emat=work->emat;

  /** matrices for change of variables **/
  B=work->B;
  BI=work->BI;
  C=work->C;
  CI=work->CI;

  /* emat */
  if (emat_->size1 != N+1)
    throw invalid_argument("emat_ must be of size N+1!");

  if (emat_->size2 != N+1)
    throw invalid_argument("emat_ must be a square matrix!");
  gsl_matrix_memcpy(emat,emat_);

  /* change of variable */
  make_K2W(B,N);
  make_W2K(BI,N);
  make_S2G(C,N);
  make_G2S(CI,N);

  /* chain structure */
  make_chain_structure(trid);

}

void
model_struct::set_k(const gsl_matrix *kk){
  /*
   * set workspace state from k input
   */
  if (k != kk)
    gsl_matrix_memcpy(k,kk);    // copying is subdominant compared to inversion

  /* up */
  k2w(k,w);
  w2x(w,x);
  x2xinv(x,xinv);
  xinv2gam(xinv,gammat);
  gam2cmat(gammat,cmat);

  return;
}

void
model_struct::set_x(const gsl_matrix *xx){
  /*
   * set workspace state from x input
   */
  if (x != xx)
    gsl_matrix_memcpy(x,xx);

  /* down */
  x2w(x,w);
  w2k(w,k);

  /* up */
  x2xinv(x,xinv);
  xinv2gam(xinv,gammat);
  gam2cmat(gammat,cmat);

  return;
}

void
model_struct::set_xinv(const gsl_matrix *xxinv){
  /*
   * set workspace state from xxinv input
   */
  if (xinv != xxinv)
    gsl_matrix_memcpy(xinv,xxinv);

  /* down */
  xinv2x(xinv,x);
  x2w(x,w);
  w2k(w,k);

  /* up */
  xinv2gam(xinv,gammat);
  gam2cmat(gammat,cmat);
  return;
}

void
model_struct::set_gammat(const gsl_matrix *ggammat){
  /*
   * set workspace state from gammat input
   */
  if (gammat != ggammat)
    gsl_matrix_memcpy(gammat,ggammat);

  /* down */
  gam2xinv(gammat,xinv);
  xinv2x(xinv,x);
  x2w(x,w);
  w2k(w,k);

  /* up */
  gam2cmat(gammat,cmat);

  return;
}

void
model_struct::set_cmat(const gsl_matrix *ccmat){
  /*
   * set workspace state from xxinv input
   */
  if (cmat != ccmat)
    gsl_matrix_memcpy(cmat,ccmat);

  /* down */
  cmat2gam(cmat,gammat);
  gam2xinv(gammat,xinv);
  xinv2x(xinv,x);
  x2w(x,w);
  w2k(w,k);

  return;
}

void
model_struct::k2w(const gsl_matrix *kk, gsl_matrix *ww){
  /*
   * update w from k in workspace
   */
  gsl_vector_const_view kview = gsl_vector_const_view_array(kk->block->data,kk->block->size);
  gsl_vector_view wview = gsl_vector_view_array(ww->block->data,ww->block->size);

  linalg_spdgemv(NOTRANS, 1.0, B, &kview.vector, 0.0, &wview.vector);

  return;
}

void
model_struct::w2k(const gsl_matrix *ww, gsl_matrix *kk){
  /*
   * update k from w in workspace
   */
  gsl_vector_const_view wview = gsl_vector_const_view_array(ww->block->data,ww->block->size);
  gsl_vector_view kview = gsl_vector_view_array(kk->block->data,kk->block->size);

  linalg_spdgemv(NOTRANS, 1.0, BI, &wview.vector, 0.0, &kview.vector);

  return;
}

void
model_struct::w2x(const gsl_matrix *ww, gsl_matrix *xx){
  /*
   * update x from w in workspace
   */

  gsl_matrix_memcpy(xx,ww);
  linalg_daxpy(1.0,trid,xx);

  return;
}

void
model_struct::x2w(const gsl_matrix *xx, gsl_matrix *ww){
  /*
   * update w from x in workspace
   */

  gsl_matrix_memcpy(ww,xx);
  linalg_daxpy(-1.0,trid,ww);

  return;
}

void
model_struct::x2xinv(const gsl_matrix *xx, gsl_matrix *xxinv){
  /*
   * invert x in workspace
   */
  linalg_invert(xx,xxinv);
  return;
}

void
model_struct::xinv2x(const gsl_matrix *xxinv, gsl_matrix *xx){
  /*
   * invert xinv in workspace
   */
  linalg_invert(xxinv,xx);
  return;
}

void
model_struct::xinv2gam(const gsl_matrix *xxinv, gsl_matrix *ggammat){
  /*
   * compute ggammat from xxinv
   */

  gsl_vector_const_view xiview = gsl_vector_const_view_array(xxinv->block->data,xxinv->block->size);
  gsl_vector_view gamview = gsl_vector_view_array(ggammat->block->data,ggammat->block->size);

  linalg_spdgemv(NOTRANS, 1.0, C, &xiview.vector, 0.0, &gamview.vector);

  return;
}

void
model_struct::gam2xinv(const gsl_matrix *ggammat, gsl_matrix *xxinv){
  /*
   * compute xxinv from ggammat.
   */

  gsl_vector_const_view gamview = gsl_vector_const_view_array(ggammat->block->data,ggammat->block->size);
  gsl_vector_view xiview = gsl_vector_view_array(xxinv->block->data,xxinv->block->size);

  linalg_spdgemv(NOTRANS, 1.0, CI, &gamview.vector, 0.0, &xiview.vector);

  return;
}

void
model_struct::gam2cmat(const gsl_matrix *ggammat, gsl_matrix *ccmat){
  /*
   * compute ccmat from ggammat
   */

  double cij,gij;
  for (int p=0;p<N+1;p++){
    for (int q=0;q<N+1;q++){
      gij = gsl_matrix_get(ggammat,p,q);
      cij = (*func)(gij,&thres);
      gsl_matrix_set(ccmat,p,q,cij);
    }
  }

  return;
}

void
model_struct::cmat2gam(const gsl_matrix *ccmat, gsl_matrix *ggammat){
  /*
   * compute ggammat from ccmat
   */

  double cij,gij,gmax,cmin;

  for (int p=0;p<N+1;p++){
    for (int q=0;q<N+1;q++){
      cij = gsl_matrix_get(ccmat,p,q);
      gij = (*funcinv)(cij,&thres);

      //* Gaussian chain limit
      // The following ensures that contact probabilities cannot
      // be lower than for a Gaussian chain.
      // This is not to be used with the minimization approach because then
      // the computation of the gradient is not accurate.
      gmax = std::abs(p-q);
      cmin = (*func)(gmax,&thres);
      if (cij < cmin)
        gij = gmax;
      // Gaussian chain limit */

      gsl_matrix_set(ggammat,p,q,gij);
    }
  }

  return;
}

void
model_struct::print_infos(){
  gsl_matrix *y(0);

  y=work->workN_1;
  linalg_dgemm(NOTRANS,NOTRANS,1.0,x,xinv,0.0,y);

  printf("C^exp:\n");
  print_matrix(emat);
  printf("C:\n");
  print_matrix(cmat);
  printf("\\Gamma:\n");
  print_matrix(gammat);
  printf("X:\n");
  print_matrix(x);
  printf("X^-1:\n");
  print_matrix(xinv);
  printf("X X^-1:\n");
  print_matrix(y);
  printf("T:\n");
  print_matrix(trid);
  printf("W:\n");
  print_matrix(w);
  printf("K:\n");
  print_matrix(k);
  return;
}

