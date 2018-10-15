//****************************************************************************
// min_contacts_kij.h - Implementation of the definitions in min_contacts_kij.h.
// Date: 2017-05-03
// Created by: Guillaume Le Treut
//****************************************************************************

#include "min_contacts_kij.h"

using namespace std;

//int dump=1;
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
    throw invalid_argument("C must be of size N^2 x (N+1)^2!");

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

/* function to minimize */
double
J(const gsl_matrix *kk, void *params){
  /* function of the quadratic distance */
  double jj,zij,mij,norm;
  int nrow,ncol,N;

  model *mod;
  gsl_matrix *z1(0),*emat(0),*cmat(0),*mask(0);

  /* initializations */
  mod=static_cast<model *> (params);
  emat=mod->emat;
  cmat=mod->cmat;
  z1=mod->work->workN1_1;
  mask = mod->mask;
  norm = linalg_ddot(mask,mask);
  nrow = emat->size1;
  ncol = emat->size2;
  N=mod->N;

  /* update model from k */
  mod->set_k(kk);

  /* compute function using vector views */
  gsl_matrix_memcpy(z1,cmat);   // z1 = cmat
  linalg_daxpy(-1.0,emat,z1);   // z1 <- z1-emat
  for (int i=0;i<N+1;i++){      // remove masked values
    for (int j=0;j<N+1;j++){
      zij = gsl_matrix_get(z1,i,j);
      mij = gsl_matrix_get(mask,i,j);
      gsl_matrix_set(z1,i,j,zij*mij);
    }
  }
  jj = linalg_ddot(z1,z1);
  jj=0.5*jj/norm;

  return jj;
}

void
dJ(const gsl_matrix *kk, void *params, gsl_matrix *gradient){
  /* compute gradient of J */

  int N,nrow,ncol;
  double thres,val,gij,cij,eij,mij,norm;
  model *mod(0);
  gsl_matrix *emat(0),*cmat(0),*gammat(0),*xinv(0),*y1(0),*y2(0),*z1(0),*z2(0),*mask(0);
  SparseMat *B(0),*C(0);
  gsl_vector_view y1view,y2view,z1view,z2view,gview;

  /* initializations */
  mod=static_cast<model *> (params);
  N=mod->N;
  thres=mod->thres;
  emat=mod->emat;
  cmat=mod->cmat;
  gammat=mod->gammat;
  xinv=mod->xinv;
  y1=mod->work->workN_1;
  y2=mod->work->workN_2;
  z1=mod->work->workN1_1;
  z2=mod->work->workN1_2;
  mask = mod->mask;
  norm = linalg_ddot(mask,mask);
  B=mod->B;
  C=mod->C;
  nrow = emat->size1;
  ncol = emat->size2;

  /* update model from k */
  mod->set_k(kk);

  /* compute gradient according to K */
  y1view = gsl_vector_view_array(y1->block->data,y1->block->size);
  y2view = gsl_vector_view_array(y2->block->data,y2->block->size);
  z1view = gsl_vector_view_array(z1->block->data,z1->block->size);
  z2view = gsl_vector_view_array(z2->block->data,z2->block->size);
  gview  = gsl_vector_view_array(gradient->block->data,gradient->block->size);

  /** compute upsilon matrix **/
  gsl_matrix_memcpy(z1,cmat);     // g=y
  linalg_daxpy(-1.0,emat,z1);     // g = cmat - emat
  for (int i=0;i<N+1;i++){
    for (int j=0;j<N+1;j++){
      gij = gsl_matrix_get(gammat,i,j);
      cij = gsl_matrix_get(cmat,i,j);
      eij = gsl_matrix_get(emat,i,j);
      mij = gsl_matrix_get(mask,i,j);
      val = (cij-eij) * mod->dfunc(gij,&thres) * mij;      // z_ij = (cij-eij) x f'(gamma_ij) x mij
      gsl_matrix_set(z2,i,j,val);
    }
  }
  linalg_spdgemv(TRANS, 1.0, C, &z2view.vector, 0.0, &y1view.vector);

  /** compute gradient according to X **/
  linalg_dgemm(NOTRANS,TRANS,1.0,y1,xinv,0.0,y2);
  linalg_dgemm(TRANS,NOTRANS,-1.0,xinv,y2,0.0,y1);  // y =-X^-1 T \Upsilon X^-1 T

  /** compute gradient according to K with W=B K**/
  linalg_spdgemv(TRANS, 1.0, B, &y1view.vector, 0.0, &gview.vector);

  /** divide by N^2 **/
  linalg_dscal(1.0/norm,gradient);

  /* exit */
  return;
}

void
dJ_finite_diff(const gsl_matrix *kk, void *params, gsl_matrix *gradient){
  /* compute gradient of J */
  int n;
  double EPS,h,df;
  gsl_matrix *ktp(0);
  EPS=1.0e-8;

  n=kk->size1;
  gsl_matrix_set_all(gradient,0.0);
  ktp=gsl_matrix_calloc(n,n);
  gsl_matrix_memcpy(ktp,kk);
  gsl_vector_const_view kkview = gsl_vector_const_view_array(kk->block->data,kk->block->size);
  gsl_vector_view ktpview = gsl_vector_view_array(ktp->block->data,ktp->block->size);
  gsl_vector_view gview = gsl_vector_view_array(gradient->block->data,gradient->block->size);

  for (int mu=0;mu<ktpview.vector.size;mu++){
    double temp=gsl_vector_get(&kkview.vector,mu);
    h=EPS*abs(temp);
    if (h==0.0) h=EPS;
    gsl_vector_set(&ktpview.vector,mu,temp+h);
    h=gsl_vector_get(&ktpview.vector,mu)-temp;

    df=(J(ktp,params)-J(kk,params))/h;
    gsl_vector_set(&gview.vector,mu,df);
    gsl_vector_set(&ktpview.vector,mu,temp);
  }

  gsl_matrix_free(ktp);

  return;
}

void
JdJ(const gsl_matrix *kk, void *params, double *jj, gsl_matrix *gradient){
  /*
   * Evaluate both function and gradient
   */

  int N,nrow,ncol;
  double thres,val,gij,cij,eij,zij,mij,norm;
  model *mod(0);
  gsl_matrix *emat(0),*cmat(0),*gammat(0),*xinv(0),*y1(0),*y2(0),*z1(0),*z2(0),*mask(0);
  SparseMat *B(0),*C(0);
  gsl_vector_view y1view,y2view,z1view,z2view,gview;

  /* initializations */
  mod=static_cast<model *> (params);
  N=mod->N;
  thres=mod->thres;
  emat=mod->emat;
  cmat=mod->cmat;
  gammat=mod->gammat;
  xinv=mod->xinv;
  y1=mod->work->workN_1;
  y2=mod->work->workN_2;
  z1=mod->work->workN1_1;
  z2=mod->work->workN1_2;
  B=mod->B;
  C=mod->C;
  nrow = emat->size1;
  ncol = emat->size2;
  mask = mod->mask;
  norm = linalg_ddot(mask,mask);

  /* update model from k */
  mod->set_k(kk);

  /* compute function */
  gsl_matrix_memcpy(z1,cmat);   // y = cmat
  linalg_daxpy(-1.0,emat,z1);   // y <- y-emat
  for (int i=0;i<N+1;i++){      // remove masked values
    for (int j=0;j<N+1;j++){
      zij = gsl_matrix_get(z1,i,j);
      mij = gsl_matrix_get(mask,i,j);
      gsl_matrix_set(z1,i,j,zij*mij);
    }
  }

  *jj = linalg_ddot(z1,z1);
  *jj=0.5*(*jj)/norm;

  /* compute gradient according to K */
  y1view = gsl_vector_view_array(y1->block->data,y1->block->size);
  y2view = gsl_vector_view_array(y2->block->data,y2->block->size);
  z1view = gsl_vector_view_array(z1->block->data,z1->block->size);
  z2view = gsl_vector_view_array(z2->block->data,z2->block->size);
  gview  = gsl_vector_view_array(gradient->block->data,gradient->block->size);

  /** compute upsilon matrix **/
//  gsl_matrix_memcpy(z1,cmat);     // g=y
//  linalg_daxpy(-1.0,emat,z1);     // g = cmat - emat
  for (int i=0;i<N+1;i++){
    for (int j=0;j<N+1;j++){
      gij = gsl_matrix_get(gammat,i,j);
      cij = gsl_matrix_get(cmat,i,j);
      eij = gsl_matrix_get(emat,i,j);
      mij = gsl_matrix_get(mask,i,j);
      val = (cij-eij) * mod->dfunc(gij,&thres) * mij;      // z_ij = (cij-eij) x f'(gamma_ij) x mij
      gsl_matrix_set(z2,i,j,val);
    }
  }
  linalg_spdgemv(TRANS, 1.0, C, &z2view.vector, 0.0, &y1view.vector); // y = C^T . z2

  /** compute gradient according to X **/
  linalg_dgemm(NOTRANS,TRANS,1.0,y1,xinv,0.0,y2);   // y = \Upsilon X^-1 T
  linalg_dgemm(TRANS,NOTRANS,-1.0,xinv,y2,0.0,y1);  // y =-X^-1 T \Upsilon X^-1 T

  /** compute gradient according to K with W=B K**/
  linalg_spdgemv(TRANS, 1.0, B, &y1view.vector, 0.0, &gview.vector);

  /** divide by N^2 **/
  linalg_dscal(1.0/norm,gradient);

  /* exit */
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
  mask=gsl_matrix_calloc(N+1,N+1);
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
  gsl_matrix_free(mask);
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
model_struct::model_struct(gsl_matrix *emat_, double thres_, double (*func_)(double,void*), double (*dfunc_)(double,void*), double (*funcinv_)(double,void*), workspace *work_, double cmask) {
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

  /*** N+1xN+1 ***/
  k=work->k;
  cmat=work->cmat;
  gammat=work->gammat;
  emat=work->emat;
  mask=work->mask;

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

  /* build mask based on emat */
  gsl_matrix_set_all(mask,1.);
  for (size_t i=0;i<mask->size1;i++){
    for (size_t j=0;j<mask->size2;j++){
      double eij = gsl_matrix_get(emat,i,j);
      if (!(eij > cmask)) gsl_matrix_set(mask,i,j,0.);
    }
  }
  /* TEST
  //gsl_matrix_set(mask,0,1,0);
  printf("MASK:\n");
  print_matrix(mask);
  // TEST */

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


      /* Gaussian chain limit
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

void
compute_distances(double *f, double *f_rel, void *params){
  /* compute various distances given the input model */
  double zij,mij,f1,f2,norm;
  int N;

  model *mod;
  /* matrices of size N+1 */
  gsl_matrix *z1(0),*z2(0),*emat(0),*cmat(0),*mask(0);

  /* initializations */
  mod=static_cast<model *> (params);
  emat=mod->emat;
  cmat=mod->cmat;
  z1=mod->work->workN1_1;
  z2=mod->work->workN1_2;
  mask=mod->mask;
  N=mod->N;

  /* normalization - number of non-zero elements in the mask */
  norm = linalg_ddot(mask,mask);
  /* TEST
  cout << "*******************" << endl;
  cout << "norm = " << norm << endl;
  cout << "*******************" << endl;
  // TEST */

  /* compute f */
  gsl_matrix_memcpy(z1,cmat);   // z1 = cmat
  linalg_daxpy(-1.0,emat,z1);   // z1 <- z1-emat
  for (int i=0;i<N+1;i++){      // remove masked values
    for (int j=0;j<N+1;j++){
      mij = gsl_matrix_get(mask,i,j);
      zij = gsl_matrix_get(z1,i,j);
      gsl_matrix_set(z1,i,j,zij*mij);
    }
  }
  *f = linalg_dnrm2(z1) / sqrt(norm);

  /* compute f_rel */
  gsl_matrix_memcpy(z1,cmat);   // z1 = cmat
  gsl_matrix_memcpy(z2,emat);   // z2 = emat
  for (int i=0;i<N+1;i++){      // remove masked values
    for (int j=0;j<N+1;j++){
      // mask
      mij = gsl_matrix_get(mask,i,j);

      // z1
      zij = gsl_matrix_get(z1,i,j);
      gsl_matrix_set(z1,i,j,zij*mij);

      // z2
      zij = gsl_matrix_get(z2,i,j);
      gsl_matrix_set(z2,i,j,zij*mij);
    }
  }
  f1 = linalg_dnrm2(z1) / sqrt(norm);
  f2 = linalg_dnrm2(z2) / sqrt(norm);
  *f_rel = 2*(*f) / (f1 + f2); // 2 ||M1 - M2|| / (||M1|| + ||M2||)

  return;
}

/* minimization */
minimization_struct::minimization_struct(workspace *wkp){
  /* copy workspace address */
  work=wkp;

  /** matrices for the minimization **/
  k1=work->k1;
  g1=work->g1;
  dk=work->dk;
  ac=work->ac;

  /** matrix size **/
  N=work->N;

  /** active constraints **/
  gsl_matrix_set_all(ac,0.0);
  nac=0;

  /** minimization parameters **/
  itermax = 1.0e4;
  step = step_max = 1.0;
  tol = 0.1;
  ftol = macheps;
  frtol = 3.0e-08;
  ktol = macheps;
  converged = 0;
  gnorm = k0norm = f0 = 1.0;
  df = dfrel = dknorm = 0.0;
}

void
minimization_struct::iterate(gsl_matrix *k, void *params, double *f, gsl_matrix *gradient){
  /* take one step in the steepest descent direction */

  /* declarations */
  double f1;
  int failed;

  /* start */
  failed = 0;

  /* test on the gradient */
  gnorm = linalg_dnrm2(gradient);

  if (gnorm == 0.0)
  {
    df=0.0;
    gsl_matrix_set_zero(dk);
    return;
  }

trial:
  /* next point */
  gsl_matrix_set_zero(dk);
  linalg_daxpy(-step/gnorm,gradient,dk); // descent direction
  gsl_matrix_memcpy(k1,k);
  linalg_daxpy(1.0,dk,k1);

  /* project k1 */
  project_k(k1);

  /* evaluate function at new point */
  LdL(k1, params, &f1, g1);
  df=f1-*f;

  /* adjust step */
  if (df > 0.0){
    /* downhill step failed, reduce step-size and try again */
    failed = 1;
    step *= tol;
    goto trial;
  }

  if (failed)
    step *=tol;
  else
    step *= 2.0;

  /* update attributes */
  dfrel = std::abs(df)*2.0/(std::abs(f1)+std::abs(*f));
  *f=f1;
  gsl_matrix_memcpy(dk,k1);
  linalg_daxpy(-1.0,k,dk);
  dknorm=linalg_dnrm2(dk);
  gsl_matrix_memcpy(k,k1);
  gsl_matrix_memcpy(gradient,g1);

  return;
}

void
minimization_struct::test_convergence(){
  /* stop criterion */

  if (std::fabs(df/f0) < ftol){
    printf("condition satisfied: |df/f0|=%.6e\n",std::abs(df/f0));
    converged = 1;
    return;
  }
  if (dfrel < frtol){
    printf("condition satisfied: |df| 2 /(|fn| + |fnn|)=%.6e\n", dfrel);
    converged = 1;
    return;
  }
  if (dknorm/k0norm < ktol){
    printf("condition satisfied: ||dk||/||k0||=%.6e\n",dknorm/k0norm);
    converged = 1;
    return;
  }
  if (gnorm == 0.0){
    printf("condition satisfied: ||gradient||=%.6e\n",gnorm);
    converged = 1;
    return;
  }
  return;
}

void
minimization_struct::project_k(gsl_matrix *k, bool all){
  /*
   * project inequality constraints on R+
   */

  int p,q;
  int cons;
  double val;

  for (p=0;p<=N;p++){
    for (q=0;q<=N;q++){
      cons=int(gsl_matrix_get(ac,p,q));
      if (  (cons == 1) // it the constraint is active
         || (all) ){    // or if all is activated.
        /* diagonal to zero */
        if (p == q)
          val=0.0;

        /* k0q coefficients
        else if ( (p == 0) || (q == 0) )
          val=0.0;
        //*/

        //* leave bond coefficients to zero
        else if ( (p == q+1) || (p == q-1) )
          val=0.0;
        //*/

        /* project other coefficients on R+ */
        else {
          val=gsl_matrix_get(k,p,q);
          val=max(val,0.0);
        }

        /* change value */
        gsl_matrix_set(k,p,q,val);
      }
    }
  }
  return;
}

void
minimization_struct::get_largest_constraint(gsl_matrix *k, size_t *p, size_t *q, gsl_matrix *kwork){

  size_t pmin,pmax,qmin,qmax;
  double kwmax,kwmin;
  gsl_matrix_memcpy(kwork,k);
  project_k(kwork,true);

  linalg_daxpy(-1.0,k,kwork);
  gsl_matrix_minmax_index(kwork,&pmin,&qmin,&pmax,&qmax);
  kwmin=gsl_matrix_get(kwork,pmin,qmin);
  kwmax=gsl_matrix_get(kwork,pmax,qmax);
  kwmin=std::fabs(kwmin);
  kwmax=std::fabs(kwmax);
  if (kwmin > kwmax){
    *p=pmin;
    *q=qmin;
  }
  else {
    *p=pmax;
    *q=qmax;
  }

  return;
}

void
minimization_struct::min_J(gsl_matrix *k, double *f, void *params){
  /* declarations */
  gsl_matrix *gradient(0);

  /* initializations */
  gradient = work->g0;

  /* initialize steepest descent */
  project_k(k);
  LdL(k,params,f,gradient);
  f0=(*f);                               // store initial function value
  k0norm=linalg_dnrm2(k);  // store norm of initial k
  if (k0norm == 0.0) k0norm=1.0;
  gnorm=linalg_dnrm2(gradient);
  step = step_max;
  converged = 0;

  /* minimization */
  for (int iter=1;iter<itermax;iter++){

    /** iterate **/
    iterate(k,params,f,gradient);

    /** test convergence **/
    test_convergence();

#if defined(VERBOSE)
    if (converged == 1)
      printf("\t%-s\n", "Minimum found!");
#endif

#if defined(VERBOSE)
    if ( (iter % dump == 0) || (converged == 1) ){
      printf("%-20d%-+20.8e%-+20.8e%-+20.8e%-+20.8e%-+20.8e%-+20.8e\n", iter, step, *f , gnorm, std::abs(df/f0), dknorm/k0norm,dfrel);
    }
#endif

    if (iter == itermax-1) {
#if defined(VERBOSE)
      printf("\t%-s\n", "Last iteration reached!");
#endif
    }

    if (converged == 1) {
      break;
    }
  }

  return;
}

int
minimization_struct::get_nac(){
  double dnac(0);

  dnac = linalg_ddot(ac,ac);

  return int(dnac);
}

void
minimization_struct::add_active_constraint_rng(){
  int p,q;
  double r;
  gsl_rng *rng(0);

  rng=work->rng;

  if ( nac == (N+1)*(N+1) ){
    printf("Already all constraints are active!\n");
    return;
  }

  /* draw a pair of indices */
  r=gsl_rng_uniform(rng); //U(0,1)
  p=int((N+1)*r);
  r=gsl_rng_uniform(rng); //U(0,1)
  q=int((N+1)*r);

  /* activate constraint */
  gsl_matrix_set(ac,p,q,1.0);
  gsl_matrix_set(ac,q,p,1.0);

  /* update nac */
  nac=get_nac();

  /* exit */
  return;
}

void
minimization_struct::add_active_constraint(int p, int q){

  if ( nac == (N+1)*(N+1) ){
#ifdef DEBUG
    printf("Already all constraints are active!\n");
#endif
    return;
  }

  /* activate constraint */
  gsl_matrix_set(ac,p,q,1);
  gsl_matrix_set(ac,q,p,1);

  /* update nac */
  nac=get_nac();

  /* exit */
  return;
}

void
minimization_struct::print_infos(){

  printf("N=%d  nac=%d\n",N,nac);
  printf("itermax=%d  converged=%d\n",itermax,converged);
  printf("step=%.2lf  step_max=%.2lf  tol=%.2lf\n",step,step_max,tol);
  printf("ftol=%.4e  ktol=%.4e\n",ftol,ktol);

  printf("K1:\n");
  print_matrix(k1);
  printf("G1:\n");
  print_matrix(g1);
  printf("dK:\n");
  print_matrix(dk);
  printf("AC:\n");
  print_matrix(ac);

  return;
}

/* wrapper function */
void
project_matrix(model *mod, gsl_matrix *K, minimization *min){
  int N;
  size_t p,q;
  double f,kmin;
  workspace *wkp(0);
  gsl_matrix *gradient(0);

  /* initializations */
  N=mod->N;
  wkp=mod->work;
  min->LdL=JdJ;
  gradient=wkp->g0;

  /* test on sparse matrix density
  double frac1,frac2;
  frac1 = double(mod->B->nz) / (mod->B->nzmax);
  frac2 = double(mod->B->nz) / (mod->B->nrow * mod->B->ncol);
  printf("nz=%d  nzmax=%d  P=%d  Q=%d  nz/nzmax=%.6e  nz/PQ=%.6e\n",mod->B->nz,mod->B->nzmax,mod->B->nrow,mod->B->ncol,frac1,frac2);
  //*/

  /** initial condition **/
  {
    /* from emat
    gsl_matrix_memcpy(mod->cmat,mod->emat);       // inverse of sigma in X^0
    mod->set_cmat(mod->cmat);
    gsl_matrix_memcpy(K,mod->k);
    //*/
    /* from the tridiagonal matrix T
    gsl_matrix_memcpy(mod->x,mod->trid);           // trid in X^0
    mod->set_x(mod->x);
    gsl_matrix_memcpy(K,mod->k);
    //*/
    /* random noise
    //mod->set_cmat(emat);
    gsl_matrix *R = wkp->workN_1;
    double eps=1.0e0;
    generate_random_matrix2(R, eps, wkp->rng);
    gsl_matrix_add(mod->x,R);
    mod->set_x(mod->x);
    for (int i=0;i<N+1;i++)
      for (int j=0;j<N+1;j++)
        if (gsl_matrix_get(mod->gammat,i,j) < 0.0)
          gsl_matrix_set(mod->gammat,i,j,std::abs(gsl_matrix_get(mod->gammat,i,j)));
    mod->set_gammat(mod->gammat);
    gsl_matrix_memcpy(K,mod->k);
    //*/

//    f=J(K,mod);                                   // function at the initial point
//    dJ(K,mod,gradient);                                  // gradient at the initial point
    JdJ(K,mod,&f,gradient);
#ifdef DEBUG
    mod->print_infos();
#endif
  }

#if defined(DEBUG)
  {
    printf("==============================================================================\n");
    printf("INITIALIZATION\n");
    printf("==============================================================================\n");
    printf("N=%d  N*N=%g\n", N,double(N*N));
    printf("X^0:\n");
    print_matrix(wkp->x);
    printf("f(X)=%.6e\n",f);
    printf("gradf(X):\n");
    print_matrix(gradient);

    /* test gradient expression
    gsl_matrix *gradient_fd = gsl_matrix_calloc(gradient->size1,gradient->size2);
    dJ_finite_diff(K,mod,gradient_fd);
//    printf("gradf(X) (fin. diff):\n");
//    print_matrix(gradient_fd);
    gsl_vector_view rview;
    gsl_vector *rtp = gsl_vector_calloc(N+1);
    for (int i=0;i<N;i++){
      rview = gsl_matrix_row(gradient,i);
      gsl_vector_memcpy(rtp,&rview.vector);
      rview = gsl_matrix_row(gradient_fd,i);
      linalg_daxpy(-1.0,&rview.vector,rtp);
      printf("gradient[%d]-gradient_fd[%d]\n",i,i);
      print_vector(rtp);
      printf("||gradient[%d]-gradient_fd[%d]||=%.6e\n",i,i,linalg_dnrm2(rtp));

    }
    gsl_vector_free(rtp);

    linalg_daxpy(-1.0,gradient,gradient_fd);
    printf("|| gradient - gradient_fd || = %.6e\n",linalg_dnrm2(gradient_fd));
    printf("|| gradient - gradient_fd || / ||gradient|| = %.6e\n",linalg_dnrm2(gradient_fd)/linalg_dnrm2(gradient));
    printf("|| gradient - gradient_fd || / (N+1) = %.6e\n",linalg_dnrm2(gradient_fd)/(N+1));
    // end test */
  }
#endif

  /* run algorithm */

#if defined(VERBOSE)
    {
      printf("==============================================================================\n");
      printf("STEEPEST DESCENT WITH PROJECTION\n");
      printf("==============================================================================\n");
    }
#endif

    /* first minimization */
    gsl_matrix_set_all(min->ac,1); // full constraints from start (the next loop is irrelevant)
    min->min_J(K, &f, mod);

#ifdef DEBUG
    /* result */
    printf("K\n");
    print_matrix(K);
    printf("C:\n");
    print_matrix(mod->cmat);
    printf("E:\n");
    print_matrix(mod->emat);
#endif

#ifdef VERBOSE
    /* print distance of cmat matrix to emat */
    double distance,distance_rel;
    mod->set_k(K);
    compute_distances(&distance,&distance_rel,mod);
    printf("||C - E|| / (#nz el. mask) =%.6e\n",distance);

    /* print constraint history */
//    printf("Active constraints:\n");
//    print_matrix(min->ac);
#endif

  /* exit */
  //delete min;
  return;
}

void
project_matrix(gsl_matrix *emat, gsl_matrix *K, double thres, double (*func)(double, void*),  double (*dfunc)(double, void*), double (*funcinv) (double, void *), minimization *min){
  /*
   */

  /* declarations */
  int N;
  double f,kmin;
  workspace *wkp(0);
  model *mod(0);

  /* initializations */
  N=emat->size1-1;
  wkp=min->work;
  mod=new model(emat,thres,func,dfunc,funcinv,wkp);

  /* project */
  project_matrix(mod,K, min);

  /* exit */
  delete mod;

  return;
}

