#include "carvig.h"


static double prod(const double *val, int n);

static double normcdf(const double x);

extern void mat_block(const double *A, int n, int m, int row_start, int rows,  int col_start, int cols,  double *C);

static void mat_add(const int *int_part, const double *frac_part, double *full, int n);

extern void fimatmul(const char *tr, int n, int k, int m, double alpha,
					 const double *A, const int *B, double beta, double *C);

extern void ifmatmul(const char *tr, int n, int k, int m, double alpha,
				   const int *A, const double *B, double beta, double *C);

extern int *ieye(int n);
extern int *izeros(int n, int m);
static void fatalerr(const char *format, ...);

static void mat_int_frac(const double *full, int *int_part, double *frac_part, int n);
//将Qz拆分成未固定协方差和固定协方差四块矩阵，n表示Qz的行和列，a表示未固定的行和列
static void mat_a_split(double *Qz,double *Qz1,double *Qz2,double *Qz1z2,double *Qz2z1,int n,int a);

static int PAR_ddidx(rtk_t *rtk, int *ix,int *sat_record);

//static int PAR_reset_ssat_fix(rtk_t *rtk, int *index, int unfixed, int nb);
static int parsearch(int n, int m, const double *a, const double *Qa, int *record, double *L, double *D,
				   double *Z, const double *zs, double *zn, double *s, int * delete_index_in_origin, double thresar);

extern int PAR_lambda(rtk_t *rtk,int n, int m, const double *a, const double *Q, double *F,
                  double *s,int *sat_record);

/* same interface as resamb_LAMBDA2 in rtkpos.cc */
extern int PAR_resamb_LAMBDA(rtk_t *rtk, const nav_t *nav, const obsd_t *obs,
                      double *bias, double *xa,
                      ddsat_t *ddsat, int *namb,
                      ddsat_t *wlsat, int *nw,
                      double *dxwl, double *Qywl, double *wlb,
                      const int *vflg, int nv,
                      int *nny, int *fixwl);

static void perm(int n, double *L, double *D, int j, double del, double *Z);

static void gauss(int n, double *L, double *Z, int i, int j);

static int LD(int n, const double *Q, double *L, double *D);

static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s);
