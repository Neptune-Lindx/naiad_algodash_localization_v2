#ifndef PARLAMBDA_H
#define PARLAMBDA_H
#include "carvig.h"
#define PAR_MIN_SAT (4)//3个双差+4个系统
extern double FindMinRatio(double cl, int npar);
/* Q: symmetric covariance, column-major n x n (Q[i+j*n]); 1=PD, 0=not */
extern int checkeig(int n, const double *Q);
extern int parlambda(rtk_t *rtk,int n, int m, const double *a, const double *Q,
                     double *F, double *s,int *sat_record);
/* same interface as PAR_resamb_LAMBDA / resamb_LAMBDA2 */
extern int resamb_parlambda(rtk_t *rtk, const nav_t *nav, const obsd_t *obs,
                            double *bias, double *xa,
                            ddsat_t *ddsat, int *namb,
                            ddsat_t *wlsat, int *nw,
                            double *dxwl, double *Qywl, double *wlb,
                            const int *vflg, int nv,
                            int *nny, int *fixwl);
#endif /* PARLAMBDA_H */