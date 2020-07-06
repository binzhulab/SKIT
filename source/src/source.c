/* History: Jul 15 2019 Initial coding
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define DEBUG 0
#define EQUAL_EPS 1e-6
#define ALMOST_ZERO 1e-150
#define SQRT2PI 2.5066282746310002416
#define SQRT2 1.4142135623730951455
#define NRETURN 5
#define EXP0ARG 1e-16

#define MAX( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}

/*
static void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %g ", vec[i]);
  }
  printf("\n \n");
}
static void print_iVec(vec, n, name)
int *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %d ", vec[i]);
  }
  printf("\n \n");
}
*/

/* Function to allocate memory for a double vector */
static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */


static int * iVec_alloc(n, initFlag, initVal)
int n, initFlag, initVal;
{
  int i, *ret, *p;

  ret = (int *) malloc(n*sizeof(int));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: iVec_alloc */

static void getInfo(x, y, N, xne0, yne0, sdx, sdy)
double *x, *y, *sdx, *sdy;
int N, *xne0, *yne0;
{
  int i, *pixi, *piyi, mx=0, my=0;
  double tmp, tmp2, *pxi, *pyi, sumx=0.0, sumx2=0.0, sumy=0.0, sumy2=0.0;

  *sdx = -1.0;
  *sdy = -1.0;

  /* Determine which are 0, and get the standard deviations from non-zero ones */
  for (i=0, pxi=x, pyi=y, pixi=xne0, piyi=yne0; i<N; i++, pxi++, pyi++, pixi++, piyi++) {
    tmp  = *pxi;
    tmp2 = tmp*tmp;
    if (tmp2 > ALMOST_ZERO) {
      *pixi  = 1;
      sumx  += tmp; 
      sumx2 += tmp2;
      mx++;
    } else {
      *pixi = 0;
    }
    tmp  = *pyi;
    tmp2 = tmp*tmp;
    if (tmp2 > ALMOST_ZERO) {
      *piyi = 1;
      sumy  += tmp;
      sumy2 += tmp2;
      my++;
    } else {
      *piyi = 0;
    }
  }

  if (mx > 1) {
    tmp  = sumx/mx;
    tmp2 = (sumx2 - mx*tmp*tmp)/(mx - 1.0);
    *sdx = (tmp2 > 0.0) ? sqrt(tmp2) : 0.0;
  }
  if (my > 1) {
    tmp  = sumy/my;
    tmp2 = (sumy2 - my*tmp*tmp)/(my - 1.0);
    *sdy = (tmp2 > 0.0) ? sqrt(tmp2) : 0.0;
  }

} /* END: getInfo */

static double getDefaultBandwith(sdx, sdy, pown02, def, retSIG)
double sdx, sdy, def, pown02, *retSIG;
{
  double SIG, bw;

  if ((sdx > -0.5) && (sdy > -0.5)) {
    SIG = MAX(sdx, sdy);
  } else {
    SIG = 1.0;
  }
  if (SIG > 0.0) {
    bw = SIG*pown02;
  } else {
    bw = def;
  }
  *retSIG = SIG;

  return(bw);

} /* END: getDefaultBandwith */

static void TS_main(x, y, N, h, xne0, yne0, retT0, retT1, retT2, retT3, retT4)
double *x, *y, h, *retT0, *retT1, *retT2, *retT3, *retT4;
int N, *xne0, *yne0;
{
  int i, j, *pixi, *piyi, *pixj, *piyj, sum_xeq0=0, sum_yeq0=0, sum_xyeq0=0;
  int xine0, yine0, xjne0, yjne0, xijne0, xijeq0, yijne0, yijeq0, xieq0, yieq0;
  double *pxi, *pyi, *pxj, *pyj, tmp, n, n2, n3, n4, T1, T2, T3, T4, sigma, sigma2;
  double twoSigma2, oneOverSqrt2piSigma, xi, xj, yi, yj, dnormxij, dnormyij;
  double T2_dsum1=0.0, T2_dsum2=0.0, T2_sumi=0.0, T3_dsum1=0.0, T3_dsum2=0.0, T3_sumi=0.0;
  double T4_dsum1=0.0, T4_sumi=0.0, xsumj, ysumj;
  
  n                   = (double) N;
  n2                  = n*n;
  n3                  = n*n2;
  n4                  = n*n3;
  sigma               = SQRT2*h;
  sigma2              = sigma*sigma;
  twoSigma2           = 2.0*sigma2;
  oneOverSqrt2piSigma = 1.0/(SQRT2PI*sigma);

  for (i=0, pixi=xne0, piyi=yne0, pxi=x, pyi=y; i<N; i++, pixi++, piyi++, pxi++, pyi++) {
    xine0      = *pixi;
    yine0      = *piyi; 
    xieq0      = !xine0;
    yieq0      = !yine0;
    
    if (xine0) {
      xi    = *pxi;
      xsumj = 0.0;
    } else {
      sum_xeq0++;
    }
    if (yine0) {
      yi    = *pyi;   
      ysumj = 0.0;
    } else {
      sum_yeq0++;
    }
    if (xieq0 && yieq0) {
      sum_xyeq0++;
      continue; 
    }

    for (j=0, pixj=xne0, piyj=yne0, pxj=x, pyj=y; j<N; j++, pixj++, piyj++, pxj++, pyj++) {
      xjne0    = *pixj;
      yjne0    = *piyj; 
      xijne0   = xine0 && xjne0;
      yijne0   = yine0 && yjne0;
      xijeq0   = xieq0 && !xjne0;
      yijeq0   = yieq0 && !yjne0;

      if (xijne0) {
        xj  = *pxj;
        tmp = xi - xj;
        tmp = tmp*tmp;
        if (tmp > EXP0ARG) {
          dnormxij = exp(-tmp/twoSigma2);
        } else {
          dnormxij = 1.0;
        }
        xsumj     += dnormxij;
        T2_dsum1  += dnormxij*yijeq0;
      } else {
        dnormxij = 0.0;
      }
      
      if (yijne0) {
        yj       = *pyj;
        tmp      = yi - yj;
        tmp      = tmp*tmp;
        if (tmp > EXP0ARG) {
          dnormyij = exp(-tmp/twoSigma2);
        } else {
          dnormyij = 1.0;
        }
        ysumj     += dnormyij;
        T3_dsum1  += dnormyij*xijeq0;
        T4_dsum1  += dnormxij*dnormyij;
      } 

    } /* END: loop j */

    if (xine0) {
      xsumj     *= oneOverSqrt2piSigma;
      T2_dsum2  += xsumj;
      T2_sumi   += xsumj*yieq0;
    }
    if (yine0) { 
      ysumj     *= oneOverSqrt2piSigma;
      T3_dsum2  += ysumj;
      T3_sumi   += ysumj*xieq0;
    }  
    if (xine0 && yine0) T4_sumi += xsumj*ysumj;
  }

  T2_dsum1 *= oneOverSqrt2piSigma;
  T3_dsum1 *= oneOverSqrt2piSigma;
  T4_dsum1 *= oneOverSqrt2piSigma*oneOverSqrt2piSigma;

  tmp = sum_xyeq0/n - sum_xeq0*sum_yeq0/n2;
  T1  = tmp*tmp; 
  T2  = T2_dsum1/n2 + sum_yeq0*sum_yeq0*T2_dsum2/n4 - 2.0*sum_yeq0*T2_sumi/n3;
  T2  = MAX(T2, 0.0);
  T3  = T3_dsum1/n2 + sum_xeq0*sum_xeq0*T3_dsum2/n4 - 2.0*sum_xeq0*T3_sumi/n3;
  T3  = MAX(T3, 0.0);
  T4  = T4_dsum1/n2 + T2_dsum2*T3_dsum2/n4 - 2.0*T4_sumi/n3;
  T4  = MAX(T4, 0.0);

  *retT0 = T1 + T2 + T3 + T4;
  *retT1 = T1;
  *retT2 = T2;
  *retT3 = T3;
  *retT4 = T4;

} /* END: TS_main */


void C_TS(x, y, pn, ph, ret)
double *x, *y, *ret, *ph;
int *pn;
{
  int *xne0, *yne0, n;
  double sdx, sdy, T0, T1, T2, T3, T4;
  
  n    = *pn;
  xne0 = iVec_alloc(n, 0, 0);
  yne0 = iVec_alloc(n, 0, 0);

  getInfo(x, y, n, xne0, yne0, &sdx, &sdy);
  TS_main(x, y, n, *ph, xne0, yne0, &T0, &T1, &T2, &T3, &T4);
  ret[0] = T0;
  ret[1] = T1;
  ret[2] = T2;
  ret[3] = T3;
  ret[4] = T4;

  free(xne0);
  free(yne0);

} /* END: C_TS */

static void getBootVec(vec, n, ret)
double *vec, *ret;
int n;
{
  int i, j;
  double dn, *pret;

  /* runif(0.0, i+1.0) gives a number 0 <= x <= i */

  dn = (double) n;
  for (i=0, pret=ret; i<n; i++, pret++) {
    j = floor(runif(0.0, dn));
    *pret = vec[j];
  }

} /* END: getBootVec */

static void TS_boot(x, y, n, h, h_fixed, SIG, nboot, xne0, yne0, print, TSobs, pown02, retp)
double *x, *y, h, *retp, SIG, *TSobs, pown02;
int h_fixed, n, nboot, print, *xne0, *yne0;
{
  int iter, sum0, sum1, sum2, sum3, sum4;
  double *bootx, *booty, bw_boot, sdx, sdy, retSIG, T0, T1, T2, T3, T4, T00, T10, T20, T30, T40;

  bootx = dVec_alloc(n, 0, 0.0);
  booty = dVec_alloc(n, 0, 0.0);

  if (!h_fixed) {
    if (SIG > 0.0) {
      bw_boot = SIG*pown02;
    } else {
      bw_boot = h;
    }
  } else {
    bw_boot = h;
  }

  T00  = TSobs[0];
  T10  = TSobs[1];
  T20  = TSobs[2];
  T30  = TSobs[3];
  T40  = TSobs[4];
  sum0 = 0;
  sum1 = 0;
  sum2 = 0;
  sum3 = 0;
  sum4 = 0;

  for (iter=1; iter<nboot+1; iter++) {
    if (print && !(iter % nboot)) Rprintf("bootstrap sample %d\n", iter);
    getBootVec(x, n, bootx);
    getBootVec(y, n, booty);
    getInfo(bootx, booty, n, xne0, yne0, &sdx, &sdy);

    if (!h_fixed) bw_boot = getDefaultBandwith(sdx, sdy, pown02, h, &retSIG);
    TS_main(bootx, booty, n, bw_boot, xne0, yne0, &T0, &T1, &T2, &T3, &T4);
    
    sum0 += (T0 >= T00);
    sum1 += (T1 >= T10);
    sum2 += (T2 >= T20);
    sum3 += (T3 >= T30);
    sum4 += (T4 >= T40);
  }

  retp[0] = (1.0 + sum0)/(nboot + 1.0);
  retp[1] = (1.0 + sum1)/(nboot + 1.0);
  retp[2] = (1.0 + sum2)/(nboot + 1.0);
  retp[3] = (1.0 + sum3)/(nboot + 1.0);
  retp[4] = (1.0 + sum4)/(nboot + 1.0);

  free(bootx);
  free(booty);

} /* END: TS_boot */

void C_TS_boot(x, y, pn, ph, phfixed, pnboot, pprint, retp, retbw, retobs)
double *x, *y, *ph, *retp, *retbw, *retobs;
int *pn, *phfixed, *pnboot, *pprint;
{
  int *xne0, *yne0, n, h_fixed;
  double sdx, sdy, bw_obs, SIG, h, TSobs[NRETURN], pown02;

  n       = *pn;
  h       = *ph;
  h_fixed = *phfixed;
  pown02  = pow((double) n, -0.2);

  /* For random number generation */
  GetRNGstate();

  xne0 = iVec_alloc(n, 0, 0);
  yne0 = iVec_alloc(n, 0, 0);
  getInfo(x, y, n, xne0, yne0, &sdx, &sdy);
  if (!h_fixed) {
    bw_obs = getDefaultBandwith(sdx, sdy, pown02, h, &SIG);
  } else {
    bw_obs = h;
  }

  *retbw = bw_obs;

  /* Observed results */
  TS_main(x, y, n, bw_obs, xne0, yne0, &TSobs[0], &TSobs[1], &TSobs[2], &TSobs[3], &TSobs[4]);
  retobs[0] = TSobs[0];
  retobs[1] = TSobs[1];
  retobs[2] = TSobs[2];
  retobs[3] = TSobs[3];
  retobs[4] = TSobs[4];

  TS_boot(x, y, n, bw_obs, *phfixed, SIG, *pnboot, xne0, yne0, *pprint, TSobs, pown02, retp);

  free(xne0);
  free(yne0);
  PutRNGstate();  

  return;

} /* END: C_TS_boot */

