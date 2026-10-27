/*------------------------------------------------------------------------------
* refsat.cc : RTK reference-satellite selection by observation quality
*
* Li et al., Meas. Sci. Technol. 35 (2024) 066309.
* Indicators: PDOP (negative), C/N0 (positive), pseudorange multipath (negative).
* PDOP is rebuilt with the candidate kept in the double difference as the base.
* GPS, QZSS and SBAS share group m=4. BeiDou GEO is never the reference.
*-----------------------------------------------------------------------------*/
#include "carvig.h"

#define REF_PDOP_BAD 1E3
#define REF_RHO      0.5
#define REF_ELMIN    (20.0*D2R)         /* rover elevation required of a reference */
#define REF_DSCORE   (-0.15)            /* keep the reference while its score drop stays within this */

/* m=1:GLO, 2:GAL, 3:BDS, 4:GPS/QZS/SBS, 5:IRN (same groups as ddres) --------*/
static int ref_in_group(int sys, int m)
{
    switch (sys) {
        case SYS_GPS: return m==4;
        case SYS_SBS: return m==4;
        case SYS_GLO: return m==1;
        case SYS_GAL: return m==2;
        case SYS_CMP: return m==3;
        case SYS_QZS: return m==4;
        case SYS_IRN: return m==5;
    }
    return 0;
}
/* BDS-2 C01-C05 and BDS-3 C59-C63 -------------------------------------------*/
static int ref_bds_geo(int sat)
{
    int prn=0;
    if (satsys(sat,&prn)!=SYS_CMP) return 0;
    return prn<=5||prn>=59;
}
static int ref_validobs(int i, int j, int f, int nf, const double *y)
{
    return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
           (f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
}
static int ref_chkfrq(const obsd_t *obs, int nf)
{
    int i;
    for (i=0;i<nf;i++) if (obs->L[i]>0.0) return 1;
    return 0;
}
/* any carrier slip keeps the satellite out of the reference set --------------*/
static int ref_slip_any(const rtk_t *rtk, int sat, int nf)
{
    int q;

    if (sat<=0||sat>MAXSAT) return 1;
    for (q=0;q<nf&&q<NFREQ;q++) if (rtk->ssat[sat-1].slip[q]) return 1;
    return 0;
}
/* dual-frequency pseudorange multipath (m). f1 is the higher carrier. -------*/
static int ref_pr_mp(const obsd_t *obs, const double *lam, int f, int nf,
                     double *mp)
{
    int f1,f2,k,t;
    double a,phi1,phi2;

    if (f<0||f>=nf||lam[f]<=0.0||obs->L[f]==0.0||obs->P[f]==0.0) return 0;
    f1=f;
    f2=-1;
    for (k=0;k<nf;k++) {
        if (k==f||lam[k]<=0.0||obs->L[k]==0.0||obs->P[k]==0.0) continue;
        f2=k;
        break;
    }
    if (f2<0) return 0;
    if (lam[f2]<lam[f1]) { t=f1; f1=f2; f2=t; }
    a=SQR(lam[f2]/lam[f1]);
    if (fabs(a-1.0)<1E-3) return 0;
    phi1=obs->L[f1]*lam[f1];
    phi2=obs->L[f2]*lam[f2];
    if (f==f1) *mp=obs->P[f1]-(1.0+2.0/(a-1.0))*phi1+(2.0/(a-1.0))*phi2;
    else       *mp=obs->P[f2]-(2.0*a/(a-1.0))*phi1+(2.0*a/(a-1.0)-1.0)*phi2;
    return 1;
}
/* |combination - arc mean|. update=0 does not move the mean. ----------------*/
static double ref_mp_epoch(rtk_t *rtk, const nav_t *nav, const obsd_t *obs,
                           int sat, int ff, int nf, int update)
{
    ssat_t *ss=&rtk->ssat[sat-1];
    double raw,dev;

    if (ff<0||ff>=NFREQ) return -1.0;
    if (ss->slip[ff]&1) {
        ss->mp_n[ff]=0;
        ss->mp_arc[ff]=0.0;
    }
    if (!ref_pr_mp(obs,nav->lam[sat-1],ff,nf,&raw)) return -1.0;
    if (update) {
        ss->mp_n[ff]++;
        ss->mp_arc[ff]+=(raw-ss->mp_arc[ff])/ss->mp_n[ff];
    }
    if (ss->mp_n[ff]==0) return -1.0;
    dev=raw-ss->mp_arc[ff];
    return dev<0.0?-dev:dev;
}
/* group PDOP with candidate p kept as the double-difference base -----------*/
static double ref_dd_pdop(const double *azel, const int *iu, const int *idx,
                          int n, int piv)
{
    double H[3*MAXSAT],Q[9],elp,azp,el,az,de,dn,du,s;
    int i,col;

    elp=azel[1+iu[idx[piv]]*2];
    azp=azel[iu[idx[piv]]*2];
    for (i=col=0;i<n;i++) {
        if (i==piv) continue;
        el=azel[1+iu[idx[i]]*2];
        az=azel[iu[idx[i]]*2];
        de=cos(el)*sin(az)-cos(elp)*sin(azp);
        dn=cos(el)*cos(az)-cos(elp)*cos(azp);
        du=sin(el)-sin(elp);
        H[0+3*col]=de;
        H[1+3*col]=dn;
        H[2+3*col]=du;
        col++;
    }
    if (col<3) return REF_PDOP_BAD;
    matmul("NT",3,3,col,1.0,H,H,0.0,Q);
    if (matinv(Q,3)) return REF_PDOP_BAD;
    s=Q[0]+Q[4]+Q[8];
    return s>0.0?sqrt(s):REF_PDOP_BAD;
}
static void ref_norm_col(const double *x, double *r, int n, int positive)
{
    int i;
    double xmin,xmax,den;

    xmin=xmax=x[0];
    for (i=1;i<n;i++) {
        if (x[i]<xmin) xmin=x[i];
        if (x[i]>xmax) xmax=x[i];
    }
    den=xmax-xmin;
    if (den<1E-12) {
        for (i=0;i<n;i++) r[i]=1.0;
        return;
    }
    for (i=0;i<n;i++) r[i]=positive?(x[i]-xmin)/den:(xmax-x[i])/den;
}
/* entropy weight. A flat column has E=1 and drops out. ----------------------*/
static void ref_entropy_w(const double r[3][MAXSAT], int n, double *w)
{
    int j,i;
    double p,sum,e,den,lnn;

    lnn=log((double)n);
    for (j=0;j<3;j++) {
        for (i=0,sum=0.0;i<n;i++) sum+=r[j][i];
        if (sum<1E-15||lnn<1E-15) { w[j]=0.0; continue; }
        for (i=0,e=0.0;i<n;i++) {
            p=r[j][i]/sum;
            if (p<=1E-15) continue;
            e-=p*log(p);
        }
        w[j]=1.0-e/lnn;
        if (w[j]<0.0) w[j]=0.0;
    }
    for (j=0,den=0.0;j<3;j++) den+=w[j];
    if (den<1E-12) {
        w[0]=w[1]=w[2]=1.0/3.0;
        return;
    }
    for (j=0;j<3;j++) w[j]/=den;
}
/* grey relational degree, rho=0.5. Ideal point is the column maximum. -------*/
static void ref_grey(const double r[3][MAXSAT], int n, const double *w,
                     double *score)
{
    int i,j;
    double ideal[3],d,dmin,dmax,eps;

    for (j=0;j<3;j++) {
        ideal[j]=r[j][0];
        for (i=1;i<n;i++) if (r[j][i]>ideal[j]) ideal[j]=r[j][i];
    }
    dmin=1E9;
    dmax=0.0;
    for (i=0;i<n;i++) for (j=0;j<3;j++) {
        d=fabs(ideal[j]-r[j][i]);
        if (d<dmin) dmin=d;
        if (d>dmax) dmax=d;
    }
    for (i=0;i<n;i++) {
        score[i]=0.0;
        for (j=0;j<3;j++) {
            d=fabs(ideal[j]-r[j][i]);
            eps=dmax<1E-15?1.0:(dmin+REF_RHO*dmax)/(d+REF_RHO*dmax);
            score[i]+=w[j]*eps;
        }
    }
}
/* last grey score of each reference slot. Promoted to the baseline when the
* epoch time changes, so repeated ddres calls in one epoch share one baseline. */
static double ref_cur_score[NUMSYS][2*NFREQ];
static double ref_prev_score[NUMSYS][2*NFREQ];
static int ref_cur_ok[NUMSYS][2*NFREQ];
static int ref_prev_ok[NUMSYS][2*NFREQ];
static gtime_t ref_score_time;
static int ref_score_have;

static void ref_score_epoch(gtime_t time)
{
    if (ref_score_have&&fabs(timediff(time,ref_score_time))<1E-6) return;
    if (ref_score_have) {
        memcpy(ref_prev_score,ref_cur_score,sizeof(ref_prev_score));
        memcpy(ref_prev_ok,ref_cur_ok,sizeof(ref_prev_ok));
    }
    memset(ref_cur_ok,0,sizeof(ref_cur_ok));
    ref_score_time=time;
    ref_score_have=1;
}
/* select one reference satellite for system group m and frequency slot f.
* return index into sat[], or -1 when the group has no candidate. */
extern int sel_refsat(rtk_t *rtk, const nav_t *nav, const obsd_t *obs,
                      const int *sat, const int *iu, const int *ir,
                      const double *azel, const double *y, int ns,
                      int m, int f, int nf)
{
    int j,k,ngeom,ncand,ff,s,best,hold,nworse,prev;
    int geom[MAXSAT],cand[MAXSAT],mpok[MAXSAT];
    double x[3][MAXSAT],r[3][MAXSAT],col[MAXSAT],w[3],score[MAXSAT],mp,el;
    char id[8];

    ff=f%nf;
    if (ff<0||ff>=NFREQ||ns<=0) return -1;
    ref_score_epoch(obs[iu[0]].time);

    for (j=ngeom=0;j<ns;j++) {
        s=sat[j];
        if (!ref_in_group(rtk->ssat[s-1].sys,m)) continue;
        if (rtk->ssat[s-1].slip[ff]&1) {
            rtk->ssat[s-1].mp_n[ff]=0;
            rtk->ssat[s-1].mp_arc[ff]=0.0;
        }
        el=azel[1+iu[j]*2];
        if (el<=0.0) continue;
        if (!ref_validobs(iu[j],ir[j],f,nf,y)) continue;
        if (f<nf&&!ref_chkfrq(obs+ir[j],nf)) continue;
        geom[ngeom++]=j;
    }
    /* prefer satellites above 20 deg; fall back when the group has none */
    for (k=0;k<2;k++) {
        for (j=ncand=0;j<ngeom;j++) {
            s=sat[geom[j]];
            if (ref_bds_geo(s)) continue;
            if (ref_slip_any(rtk,s,nf)) continue;
            if (!k&&azel[1+iu[geom[j]]*2]<=REF_ELMIN) continue;
            cand[ncand++]=geom[j];
        }
        if (ncand>0) break;
    }
    if (ncand<=0) return -1;
    if (ncand==1) return cand[0];

    for (j=0;j<ncand;j++) {
        k=cand[j];
        s=sat[k];
        for (best=0;best<ngeom;best++) if (geom[best]==k) break;
        x[0][j]=ref_dd_pdop(azel,iu,geom,ngeom,best);
        x[1][j]=obs[iu[k]].SNR[ff]*0.25;
        mp=ref_mp_epoch(rtk,nav,obs+iu[k],s,ff,nf,f<nf);
        mpok[j]=mp>=0.0;
        x[2][j]=mpok[j]?mp:0.0;
    }
    for (k=-1,j=0;j<ncand;j++) {
        if (!mpok[j]) continue;
        if (k<0||x[2][j]>x[2][k]) k=j;
    }
    for (j=0;j<ncand;j++) if (!mpok[j]) x[2][j]=k<0?0.0:x[2][k];

    for (k=0;k<3;k++) {
        for (j=0;j<ncand;j++) col[j]=x[k][j];
        ref_norm_col(col,col,ncand,k==1);
        for (j=0;j<ncand;j++) r[k][j]=col[j];
    }
    ref_entropy_w(r,ncand,w);
    ref_grey(r,ncand,w,score);

    for (best=0,j=1;j<ncand;j++) if (score[j]>score[best]) best=j;

    /* keep the previous reference while it has no slip and either stays in
     * the top half or its grey score fell by no more than 0.15 */
    prev=(m>=0&&m<NUMSYS&&f>=0&&f<2*NFREQ)?rtk->refsat[m][f]:0;
    if (prev>0&&prev<=MAXSAT&&!ref_slip_any(rtk,prev,nf)) {
        for (hold=-1,j=0;j<ncand;j++) if (sat[cand[j]]==prev) { hold=j; break; }
        if (hold>=0) {
            for (nworse=0,j=0;j<ncand;j++) if (score[j]<=score[hold]) nworse++;
            if (nworse>ncand/2||
                (ref_prev_ok[m][f]&&score[hold]-ref_prev_score[m][f]>=REF_DSCORE))
                best=hold;
        }
    }
    if (m>=0&&m<NUMSYS&&f>=0&&f<2*NFREQ) {
        ref_cur_score[m][f]=score[best];
        ref_cur_ok[m][f]=1;
    }
    s=sat[cand[best]];
    satno2id(s,id);
    trace(4,"refsat: %s m=%d f=%d pdop=%.3f cn0=%.1f mp=%.3f r=%.3f %.3f %.3f w=%.3f %.3f %.3f p=%.3f\n",
          id,m,ff+1,x[0][best],x[1][best],mpok[best]?x[2][best]:-1.0,
          r[0][best],r[1][best],r[2][best],w[0],w[1],w[2],score[best]);
    return cand[best];
}
