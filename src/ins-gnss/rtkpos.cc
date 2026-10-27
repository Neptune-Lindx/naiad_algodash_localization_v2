/*------------------------------------------------------------------------------
* rtkpos.cc : precise positioning
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/12 1.0  new
*           2007/03/13 1.1  add slip detection by LLI flag
*           2007/04/18 1.2  add antenna pcv correction
*                           change rtkpos argin
*           2008/07/18 1.3  refactored
*           2009/01/02 1.4  modify rtk positioning api
*           2009/03/09 1.5  support glonass, gallileo and qzs
*           2009/08/27 1.6  fix bug on numerical exception
*           2009/09/03 1.7  add check of valid satellite number
*                           add check time sync for moving-base
*           2009/11/23 1.8  add api rtkopenstat(),rtkclosestat()
*                           add receiver h/w bias estimation
*                           add solution status output
*           2010/04/04 1.9  support ppp-kinematic and ppp-static modes
*                           support earth tide correction
*                           changed api:
*                               rtkpos()
*           2010/09/07 1.10 add elevation mask to hold ambiguity
*           2012/02/01 1.11 add extended receiver error model
*                           add glonass interchannel bias correction
*                           add slip detectior by L1-L5 gf jump
*                           output snr of rover receiver in residuals
*           2013/03/10 1.12 add otl and pole tides corrections
*           2014/05/26 1.13 support beidou and galileo
*                           add output of gal-gps and bds-gps time offset
*           2014/05/28 1.14 fix bug on memory exception with many sys and freq
*           2014/08/26 1.15 add functino to swap sol-stat file with keywords
*           2014/10/21 1.16 fix bug on beidou amb-res with pos2-bdsarmode=0
*           2014/11/08 1.17 fix bug on ar-degradation by unhealthy satellites
*           2015/03/23 1.18 residuals referenced to reference satellite
*           2015/05/20 1.19 no output solution status file with Q=0
*           2015/07/22 1.20 fix bug on base station position setting
*           2016/07/30 1.21 suppress single solution if !prcopt.outsingle
*                           fix bug on slip detection of backward filter
*           2016/08/20 1.22 fix bug on ddres() function
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include <carvig.h>
#include "PAR_LAMBDA.h"
#include "parlambda.h"
/* constants/macros ----------------------------------------------------------*/
#define VAR_POS     SQR(30.0)  /* initial variance of receiver pos (m^2) */
#define VAR_VEL     SQR(10.0)  /* initial variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)  /* initial variance of receiver acc ((m/ss)^2) */
#define VAR_AMB     SQR(10.0)  /* initial variance of receiver phase bias */
#define VAR_HWBIAS  SQR(1.0)   /* initial variance of h/w bias ((m/MHz)^2) */
#define CHKNUMERIC  1          /* check numeric for given value */
#define VAR_GRA     SQR(0.001) /* initial variance of gradient (m^2) */
#define INIT_ZWD    0.15       /* initial zwd (m) */

#define PRN_HWBIAS  1E-6       /* process noise of h/w bias (m/MHz/sqrt(s)) */
#define GAP_RESION  120        /* gap to reset ionosphere parameters (epochs) */

#define VAR_HOLDAMB 0.001      /* constraint to hold ambiguity (cycle^2) */
#define VAR_WLCONST 0.001      /* constraint to WL ambiguity */
#define THRES_AMB   1.5        /* threshold of constraint WL ambiguity */
#define THRES_HOLDAMB 3.0      /* threshold of hold ambiguity */
#define THRES_MW_JUMP 3.0      /* threshold of MW cycle slip detect*/
#define MINFAILC      30       /* min counts of solution valid fail */
#define MINFIXC       4        /* min fix count for fix abmguity */
#define UPDNEWSAT     1        /* update new satellites after update last epoch satellites */
#define NOINSJACO     0        /* no add jacobians of ins states to ekf filter */
#define DETECT_OUTLIER 0       /* detect outlier for double-differenced measurement data */
#define THRES_L1L2RES 0.03     /* threshold of L1/L2 double difference residual for detecting outliers*/

#define THRES_INHERIT_TIME 1.5 /* threshold of ambiguity inherit */
#define THRES_INHERIT_BIAS 1.0 /* threshold of ambiguity inherit */
#define INHERIT_AMB 1          /* ambiguity inherit option */

#define NX_SD       9          /* SD epoch LS states: pos(3)+clk+sys(5) */
#define MAXITR_SD  10          /* max iteration for SD epoch slip LS */

#define TTOL_MOVEB  (1.0+2*DTTOL)
/* time sync tolerance for moving-baseline (s) */

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)      /* phase bias (s:satno,f:freq) */

#ifdef EXTGSI
extern int resamb_WLNL(rtk_t *rtk, const obsd_t *obs, const int *sat,const int *iu,
                       const int *ir, int ns, const nav_t *nav,const double *azel);
extern int resamb_TCAR(rtk_t *rtk, const obsd_t *obs, const int *sat,const int *iu,
                       const int *ir, int ns, const nav_t *nav,const double *azel);
#else
extern int resamb_WLNL(rtk_t *rtk, const obsd_t *obs, const int *sat,const int *iu,
                       const int *ir, int ns, const nav_t *nav,const double *azel)
{return 0;}
extern int resamb_TCAR(rtk_t *rtk, const obsd_t *obs, const int *sat,const int *iu,
                       const int *ir, int ns, const nav_t *nav,const double *azel)
{return 0;};
#endif

/* global variables ----------------------------------------------------------*/
static int statlevel=0;          /* rtk status output level (0:off) */
static FILE *fp_stat=NULL;       /* rtk status file pointer */
static char file_stat[1024]="";  /* rtk status file original path */
static gtime_t time_stat={0};    /* rtk status file time */


/* open solution status file ---------------------------------------------------
* open solution status file and set output level
* args   : char     *file   I   rtk status file
*          int      level   I   rtk status level (0: off)
* return : status (1:ok,0:error)
* notes  : file can contains time keywords (%Y,%y,%m...) defined in reppath().
*          The time to replace keywords is based on UTC of CPU time.
* output : solution status file record format
*
*   $POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          posx/posy/posz    : position x/y/z ecef (m) float
*          posxf/posyf/poszf : position x/y/z ecef (m) fixed
*
*   $VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          vele/veln/velu    : velocity e/n/u (m/s) float
*          acce/accn/accu    : acceleration e/n/u (m/s^2) float
*          velef/velnf/veluf : velocity e/n/u (m/s) fixed
*          accef/accnf/accuf : acceleration e/n/u (m/s^2) fixed
*
*   $CLK,week,tow,stat,clk1,clk2,clk3,clk4
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          clk1     : receiver clock bias GPS (ns)
*          clk2     : receiver clock bias GLO-GPS (ns)
*          clk3     : receiver clock bias GAL-GPS (ns)
*          clk4     : receiver clock bias BDS-GPS (ns)
*
*   $ION,week,tow,stat,sat,az,el,ion,ion-fixed
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          sat      : satellite id
*          az/el    : azimuth/elevation angle(deg)
*          ion      : vertical ionospheric delay L1 (m) float
*          ion-fixed: vertical ionospheric delay L1 (m) fixed
*
*   $TROP,week,tow,stat,rcv,ztd,ztdf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          rcv      : receiver (1:rover,2:base station)
*          ztd      : zenith total delay (m) float
*          ztdf     : zenith total delay (m) fixed
*
*   $HWBIAS,week,tow,stat,frq,bias,biasf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          frq      : frequency (1:L1,2:L2,...)
*          bias     : h/w bias coefficient (m/MHz) float
*          biasf    : h/w bias coefficient (m/MHz) fixed
*
*   $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc
*          week/tow : gps week no/time of week (s)
*          sat/frq  : satellite id/frequency (1:L1,2:L2,...)
*          az/el    : azimuth/elevation angle (deg)
*          resp     : pseudorange residual (m)
*          resc     : carrier-phase residual (m)
*          vsat     : valid data flag (0:invalid,1:valid)
*          snr      : signal strength (dbHz)
*          fix      : ambiguity flag  (0:no data,1:float,2:fixed,3:hold,4:ppp)
*          slip     : cycle-slip flag (bit1:slip,bit2:parity unknown)
*          lock     : carrier-lock count
*          outc     : data outage count
*          slipc    : cycle-slip count
*          rejc     : data reject (outlier) count
*
*   $POSI,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          posx/posy/posz    : position x/y/z ecef (m) float/fixed
*
*   $VELACCI,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          vele/veln/velu    : velocity e/n/u (m/s) float/fixed
*          acce/accn/accu    : acceleration e/n/u (m/s^2) float/fixed
*
*   $ATT,week,tow,roll,pitch,yaw,roll,rollf,pitchf,yawf
*          week/tow          : gps week no/time of week (s)
*          roll/pitch/yaw    : attitude euler angle between n-frame and b-frame float/fixed (deg)
*
*   $GBIAS,week,tow,gbiasx,gbiasy,gbiasz,gbiasxf,gbiasyf,gbiaszf
*          week/tow          : gps week no/timw of week (s)
*          gyro bias-x/y/z   : gyro bias float/fixed (rad/s)
*
*   $ABIAS,week,tow,abiasx,abiasy,abiasz,abiasxf,abiasyf,abiaszf
*          week/tow          : gps week no/time of week (s)
*          accl bias-x/y/z   : accl bias float/fixed (m/s^2)
*
*   $INSTA,week,tow,stat
*          week/tow    : gps week no/time of week (s)
*          stat        : ins update status
*-----------------------------------------------------------------------------*/
extern int rtkopenstat(const char *file, int level)
{
    gtime_t time=utc2gpst(timeget());
    char path[1024];

    trace(3,"rtkopenstat: file=%s level=%d\n",file,level);

    if (level<=0) return 0;

    reppath(file,path,time,"","");

    if (!(fp_stat=fopen(path,"w"))) {
        trace(1,"rtkopenstat: file open error path=%s\n",path);
        return 0;
    }
    strcpy(file_stat,file);
    time_stat=time;
    statlevel=level;
    return 1;
}
/* close solution status file --------------------------------------------------
* close solution status file
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkclosestat(void)
{
    trace(3,"rtkclosestat:\n");

    if (fp_stat) fclose(fp_stat);
    fp_stat=NULL;
    file_stat[0]='\0';
    statlevel=0;
}
/* write solution status to buffer -------------------------------------------*/
extern int rtkoutstat(rtk_t *rtk, char *buff)
{
    ssat_t *ssat=NULL;
    double tow,pos[3],vel[3],acc[3],vela[3]={0},acca[3]={0},xa[3];
    int i,j,week,est,nfreq,nf=NF(&rtk->opt);
    char id[32],*p=buff;

    if (rtk->sol.stat<=SOLQ_NONE) {
        return 0;
    }
    /* write ppp solution status to buffer */
    if (rtk->opt.mode>=PMODE_PPP_KINEMA) {
        return pppoutstat(rtk,buff);
    }
    est=rtk->opt.mode>=PMODE_DGPS;
    nfreq=est?nf:1;
    tow=time2gpst(rtk->sol.time,&week);

    /* receiver position */
    if (est) {
        for (i=0;i<3;i++) xa[i]=i<rtk->na?rtk->xa[i]:0.0;
        p+=sprintf(p,"$POSI,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                   rtk->sol.stat,rtk->x[0],rtk->x[1],rtk->x[2],xa[0],xa[1],
                   xa[2]);
    }
    else {
        p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                   rtk->sol.stat,rtk->sol.rr[0],rtk->sol.rr[1],rtk->sol.rr[2],
                   0.0,0.0,0.0);
    }
    /* receiver velocity and acceleration */
    if (est&&rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->x+3,vel);
        ecef2enu(pos,rtk->x+6,acc);
        if (rtk->na>=6) ecef2enu(pos,rtk->xa+3,vela);
        if (rtk->na>=9) ecef2enu(pos,rtk->xa+6,acca);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                   week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],acc[0],acc[1],
                   acc[2],vela[0],vela[1],vela[2],acca[0],acca[1],acca[2]);
    }
    else {
        ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->sol.rr+3,vel);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                   week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    }
    /* receiver clocks */
    p+=sprintf(p,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
               week,tow,rtk->sol.stat,1,rtk->sol.dtr[0]*1E9,rtk->sol.dtr[1]*1E9,
               rtk->sol.dtr[2]*1E9,rtk->sol.dtr[3]*1E9);

    /* ionospheric parameters */
    if (est&&rtk->opt.ionoopt==IONOOPT_EST) {
        for (i=0;i<MAXSAT;i++) {
            ssat=rtk->ssat+i;
            if (!ssat->vs) continue;
            satno2id(i+1,id);
            j=II(i+1,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            p+=sprintf(p,"$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,id,ssat->azel[0]*R2D,ssat->azel[1]*R2D,
                       rtk->x[j],xa[0]);
        }
    }
    /* tropospheric parameters */
    if (est&&(rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG)) {
        for (i=0;i<2;i++) {
            j=IT(i,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            p+=sprintf(p,"$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,i+1,rtk->x[j],xa[0]);
        }
    }
    /* receiver h/w bias */
    if (est&&rtk->opt.glomodear==2) {
        for (i=0;i<nfreq;i++) {
            j=IL(i,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            p+=sprintf(p,"$HWBIAS,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,i+1,rtk->x[j],xa[0]);
        }
    }
    return (int)(p-buff);
}
/*write ins solution status to buffer-----------------------------------------*/
extern int insoutstat(rtk_t *rtk, char *buff)
{
    double tow,pos[3],vel[3],acc[3];
    int week; char *p=buff;

    trace(3,"insoutstat:\n");

    if (rtk->sol.stat<=SOLQ_NONE) {
        return 0;
    }
    /* write ppp solution status to buffer */
    if (rtk->opt.mode==PMODE_INS_UPDATE&&rtk->opt.insopt.tc==INSTC_PPK) {
        return pppoutstat(rtk,buff);
    }
    tow=time2gpst(rtk->sol.time,&week);

    /* receiver position */
    p+=sprintf(p,"$POSI,%d,%.3f,%d,%.4f,%.4f,%.4f\n",week,tow,
               rtk->sol.stat,rtk->sol.rr[0],rtk->sol.rr[1],rtk->sol.rr[2]);

    /* receiver velocity/acceleration */
    ecef2pos(rtk->sol.rr,pos);
    ecef2enu(pos,rtk->sol.rr+3,vel);
    ecef2enu(pos,rtk->sol.rr+6,acc);
    p+=sprintf(p,"$VELACCI,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
               rtk->sol.stat,vel[0],vel[1],vel[2],acc[0],acc[1],acc[2]);

    /* imu body attitude */
    p+=sprintf(p,"$ATT,%d,%.3f,%d,%.4f,%.4f,%.4f\n",week,tow,rtk->sol.stat,
               rtk->sol.att[0]*R2D,rtk->sol.att[1]*R2D,rtk->sol.att[2]*R2D);

    /* gyro/accl bias */
    p+=sprintf(p,"$GBIAS,%d,%.3f,%d,%.4f,%.4f,%.4f\n",week,tow,rtk->sol.stat,rtk->sol.bg[0],rtk->sol.bg[1],rtk->sol.bg[2]);
    p+=sprintf(p,"$ABIAS,%d,%.3f,%d,%.4f,%.4f,%.4f\n",week,tow,rtk->sol.stat,rtk->sol.ba[0],rtk->sol.bg[1],rtk->sol.bg[2]);

    /* ins solution status */
    p+=sprintf(p,"$INSTA,%d,%.3f,%d",week,tow,rtk->sol.ista);
    return (int)(p-buff);
}
/* swap solution status file -------------------------------------------------*/
static void swapsolstat(void)
{
    gtime_t time=utc2gpst(timeget());
    char path[1024];

    if ((int)(time2gpst(time     ,NULL)/INT_SWAP_STAT)==
        (int)(time2gpst(time_stat,NULL)/INT_SWAP_STAT)) {
        return;
    }
    time_stat=time;

    if (!reppath(file_stat,path,time,"","")) {
        return;
    }
    if (fp_stat) fclose(fp_stat);

    if (!(fp_stat=fopen(path,"w"))) {
        trace(2,"swapsolstat: file open error path=%s\n",path);
        return;
    }
    trace(3,"swapsolstat: path=%s\n",path);
}
/* output solution status ----------------------------------------------------*/
static void outsolstat(rtk_t *rtk)
{
    ssat_t *ssat;
    double tow;
    char buff[MAXSOLMSG+1],id[32];
    int i,j,n,week,nfreq,nf=NF(&rtk->opt);

    if (statlevel<=0||!fp_stat||!rtk->sol.stat) return;

    trace(3,"outsolstat:\n");

    /* swap solution status file */
    swapsolstat();

    /* write solution status */
    if (rtk->opt.mode< PMODE_INS_UPDATE) n=rtkoutstat(rtk,buff);
    if (rtk->opt.mode>=PMODE_INS_UPDATE) n=insoutstat(rtk,buff);
    buff[n]='\0';

    fputs(buff,fp_stat);

    if ((rtk->sol.stat==SOLQ_NONE&&rtk->ins.stat==INSS_NONE)||statlevel<=1) return;

    tow=time2gpst(rtk->sol.time,&week);
    nfreq=rtk->opt.mode>=PMODE_DGPS?nf:1;

    /* write residuals and status */
    for (i=0;i<MAXSAT;i++) {
        ssat=rtk->ssat+i;
        if (!ssat->vs) continue;
        satno2id(i+1,id);
        for (j=0;j<nfreq;j++) {
            fprintf(fp_stat,"$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d\n",
                    week,tow,id,j+1,ssat->azel[0]*R2D,ssat->azel[1]*R2D,
                    ssat->resp[j],ssat->resc[j],ssat->vsat[j],ssat->snr[j]*0.25,
                    ssat->fix[j],ssat->slip[j]&3,ssat->lock[j],ssat->outc[j],
                    ssat->slipc[j],ssat->rejc[j]);
        }
    }
}
/* one line per system and phase frequency: week tow sys frq sat az el.
 * Written once per epoch. Read by bin/test_skills/scripts/plot_refsat.py. */
static void log_refsat(rtk_t *rtk, const int refsat[NUMSYS][2*NFREQ], gtime_t time)
{
    static FILE *fp=NULL;
    static const char *sysname[]={"GPS","GLO","GAL","BDS","GPS","IRN"};
    char id[8];
    int m,f,sat,week,nf;
    double tow;

    if (!fp) {
        fp=fopen("bin/test_skills/feature/refsat.log","w");
        if (!fp) {
            trace(2,"log_refsat: open bin/test_skills/feature/refsat.log failed\n");
            return;
        }
        fprintf(fp,"# week tow sys frq sat az_deg el_deg\n");
    }
    nf=NF(&rtk->opt);
    tow=time2gpst(time,&week);
    for (m=0;m<6;m++) for (f=0;f<nf&&f<NFREQ;f++) {
        sat=refsat[m][f];
        if (sat<=0||sat>MAXSAT) continue;
        satno2id(sat,id);
        fprintf(fp,"%d %.3f %s %d %s %.1f %.1f\n",
                week,tow,sysname[m],f+1,id,
                rtk->ssat[sat-1].azel[0]*R2D,
                rtk->ssat[sat-1].azel[1]*R2D);
    }
    fflush(fp);
}
/* one line per rover satellite and frequency, once per epoch.
 * dcn0 = |C/N0-C/N0*|, ncn0 = dcn0/STD*. Read by plot_cn0_varerr.py. ------*/
static void log_cn0var(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const double *azel, int ns, gtime_t time)
{
    static FILE *fp=NULL;
    char id[8];
    int i,f,week,sys,nf;
    double tow,el,cn0,dcn0,ncn0;

    if (!fp) {
        fp=fopen("bin/test_skills/feature/cn0_varerr.log","w");
        if (!fp) {
            trace(2,"log_cn0var: open bin/test_skills/feature/cn0_varerr.log failed\n");
            return;
        }
        fprintf(fp,"# week tow sat frq el_deg cn0 dcn0 ncn0\n");
        fprintf(fp,"# dcn0=|cn0-cn0*|  ncn0=|cn0-cn0*|/std*\n");
    }
    nf=NF(&rtk->opt);
    tow=time2gpst(time,&week);
    for (i=0;i<ns;i++) {
        if (sat[i]<=0||sat[i]>MAXSAT) continue;
        if (!(sys=satsys(sat[i],NULL))) continue;
        el=azel[1+iu[i]*2];
        satno2id(sat[i],id);
        for (f=0;f<nf&&f<NFREQ;f++) {
            cn0=obs[iu[i]].SNR[f]*0.25;
            if (!cn0_varerr(&rtk->opt,sys,f,el,cn0,&dcn0,&ncn0)) continue;
            fprintf(fp,"%d %.3f %s %d %.2f %.2f %.3f %.3f\n",
                    week,tow,id,f+1,el*R2D,cn0,dcn0,ncn0);
        }
    }
    fflush(fp);
}
/* save error message --------------------------------------------------------*/
static void errmsg(rtk_t *rtk, const char *format, ...)
{
    char buff[256],tstr[32];
    int n;
    va_list ap;
    time2str(rtk->sol.time,tstr,2);
    n=sprintf(buff,"%s: ",tstr+11);
    va_start(ap,format);
    n+=vsprintf(buff+n,format,ap);
    va_end(ap);
    n=n<MAXERRMSG-rtk->neb?n:MAXERRMSG-rtk->neb;
    memcpy(rtk->errbuf+rtk->neb,buff,n);
    rtk->neb+=n;
    trace(2,"%s",buff);
}
/* single-differenced observable ---------------------------------------------*/
static double sdobs(const obsd_t *obs, int i, int j, int f)
{
    double pi=f<NFREQ?obs[i].L[f]:obs[i].P[f-NFREQ];
    double pj=f<NFREQ?obs[j].L[f]:obs[j].P[f-NFREQ];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
/* check precious epoch solution status---------------------------------------*/
static int chkpsol(const rtk_t* rtk)
{
    static double threshold=0.1;
    int f1=0,f2=0;

    f1=rtk->sol.pstat==SOLQ_FIX||rtk->sol.pstat==SOLQ_FLOAT;
    f2=SQRT(rtk->sol.pqr[0]+rtk->sol.pqr[1]+rtk->sol.pqr[2])<threshold;
    return f1&&f2;
}
/* single-differenced pseudorange observable by INS---------------------------*/
static double sdobspins(const rtk_t* rtk,const obsd_t *obs,const double *rs,
                        int i,int j,int f)
{
    double rr[3],e[3],sdg,sdi,ri,rj;

    sdg=sdobs(obs,i,j,f);
    if (rtk->opt.mode!=PMODE_INS_TGNSS) return sdg;

    insp2antp(&rtk->ins,rr);

    if ((ri=geodist(rs+6*i,rr,e))<=0.0) return sdg;
    if ((rj=geodist(rs+6*j,rtk->rb,e))<=0.0) return sdg;

    sdi=ri-rj;
    if (chkpsol(rtk)) {
        return sdi;
    }
    else return sdg;
}
/* single-differenced geometry-free linear combination of phase --------------*/
static double gfobs_L1L2(const obsd_t *obs, int i, int j, const double *lam)
{
    double pi=sdobs(obs,i,j,0)*lam[0],pj=sdobs(obs,i,j,1)*lam[1];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
static double gfobs_L1L5(const obsd_t *obs, int i, int j, const double *lam)
{
    double pi=sdobs(obs,i,j,0)*lam[0],pj=sdobs(obs,i,j,2)*lam[2];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
/* single-differenced measurement error variance -----------------------------*/
static double varerr(int sat, int sys, double el, double bl, double dt, int f,
                     double cn0, const prcopt_t *opt)
{
    double a,b,c=opt->err[3]*bl/1E4,d=CLIGHT*opt->sclkstab*dt,fact=1.0;
    double elu,sinel;
    int i=sys==SYS_GLO?1:(sys==SYS_GAL?2:0),nf=NF(opt),ff=nf>0?f%nf:0;

    /* eqs (9)(12)(13). θ̄=0 deletes. A lowered virtual elevation stays
     * below the 5° floor; that floor is only for the original elevation. */
    elu=cn0_el(opt,sys,ff,el,cn0);
    if (!(elu>0.0)) return 1E6;
    /* clamp elev: sin(el)~0 makes weight Inf/0 and can singularize estposbias H */
    if (elu<el) sinel=sin(elu);
    else sinel=sin(elu<5.0*D2R?5.0*D2R:elu);

    /* extended error model */
    if (f>=nf&&opt->exterr.ena[0]) { /* code */
        a=opt->exterr.cerr[i][  (f-nf)*2];
        b=opt->exterr.cerr[i][1+(f-nf)*2];
        if (sys==SYS_SBS) {a*=EFACT_SBS; b*=EFACT_SBS;}
        if (sys==SYS_CMP) {a*=EFACT_CMP; b*=EFACT_CMP;}
    }
    else if (f<nf&&opt->exterr.ena[1]) { /* phase */
        a=opt->exterr.perr[i][  f*2];
        b=opt->exterr.perr[i][1+f*2];
        if (sys==SYS_SBS) {a*=EFACT_SBS; b*=EFACT_SBS;}
        if (sys==SYS_CMP) {a*=EFACT_CMP; b*=EFACT_CMP;}
    }
    else { /* normal error model */
        if (f>=nf) fact=opt->eratio[f-nf];
        if (fact<=0.0)  fact=opt->eratio[0];
        /* Qt varerr uses EFACT_GPS for BeiDou. EFACT_CMP=3 shifted the
           held ambiguities onto the wrong integers through the open section. */
        fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
        a=fact*opt->err[1];
        b=fact*opt->err[2];
    }
    return 2.0*(opt->ionoopt==IONOOPT_IFLC?3.0:1.0)*(a*a+b*b/sinel/sinel+c*c)+d*d;
}
/* baseline length -----------------------------------------------------------*/
static double baseline(const double *ru, const double *rb, double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return norm(dr,3);
}
/* check state/covariance numeric (nan/inf) ----------------------------------*/
static int chk_xp(const double *x, const double *P, int nx)
{
    int i;
#if CHKNUMERIC
    if (!x||!P||nx<=0) return 0;
    for (i=0;i<nx;i++) {
        if (isnan(x[i])||isinf(x[i])) return 0;
        if (isnan(P[i+i*nx])||isinf(P[i+i*nx])) return 0;
    }
#endif
    return 1;
}
/* initialize state and covariance -------------------------------------------*/
extern void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi; for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}
/* initialize ins state and covariance----------------------------------------*/
extern void insinitx(insstate_t *ins,double xi,double var,int i)
{
    int j;
    ins->x[i]=xi; for (j=0;j<ins->nx;j++) {
        ins->P[i+j*ins->nx]=ins->P[j+i*ins->nx]=i==j?var:0.0;
    }
}
/* select common satellites between rover and reference station --------------*/
static int selsat(const obsd_t *obs, double *azel, int nu, int nr,
                  const prcopt_t *opt, int *sat, int *iu, int *ir)
{
    int i,j,k=0;

    trace(3,"selsat  : nu=%d nr=%d\n",nu,nr);

    for (i=0,j=nu;i<nu&&j<nu+nr;i++,j++) {
        if      (obs[i].sat<obs[j].sat) j--;
        else if (obs[i].sat>obs[j].sat) i--;
        else if (azel[1+j*2]>=opt->elmin) { /* elevation at base station */
            sat[k]=obs[i].sat; iu[k]=i; ir[k++]=j;
            trace(4,"(%2d) sat=%3d iu=%2d ir=%2d\n",k-1,obs[i].sat,i,j);
        }
    }
    return k;
}
/* temporal update of position/velocity/acceleration -------------------------*/
static void udpos(rtk_t *rtk, double tt)
{
    double *F,*P,*FP,*x,*xp,pos[3],Q[9]={0},Qv[9],var=0.0;
    int i,j,*ix,nx;

    trace(3,"udpos   : tt=%.3f\n",tt);

    /* fixed mode */
    if (rtk->opt.mode==PMODE_FIXED) {
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
            for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        }
    }
    /* static mode */
    if (rtk->opt.mode==PMODE_STATIC) return;

    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        return;
    }
    /* check variance of estimated postion */
    for (i=0;i<3;i++) var+=rtk->P[i+i*rtk->nx];
    var/=3.0;

    if (var>VAR_POS) {
        /* reset position with large variance */
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        trace(2,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }
    /* generate valid state index */
    ix=imat(rtk->nx,1);
    for (i=nx=0;i<rtk->nx;i++) {
        if (rtk->x[i]!=0.0&&rtk->P[i+i*rtk->nx]>0.0) ix[nx++]=i;
    }
    if (nx<9) {
        free(ix);
        return;
    }
    /* state transition of position/velocity/acceleration */
    F=eye(nx); P=mat(nx,nx);
    FP=mat(nx,nx); x=mat(nx,1); xp=mat(nx,1);

    for (i=0;i<6;i++) {
        F[i+(i+3)*nx]=tt;
    }
    for (i=0;i<3;i++) {
        F[i+(i+6)*nx]=SQR(tt)/2.0;
    }
    for (i=0;i<nx;i++) {
        x[i]=rtk->x[ix[i]];
        for (j=0;j<nx;j++) {
            P[i+j*nx]=rtk->P[ix[i]+ix[j]*rtk->nx];
        }
    }
    /* x=F*x, P=F*P*F+Q */
    matmul("NN",nx,1,nx,1.0,F,x,0.0,xp);
    matmul("NN",nx,nx,nx,1.0,F,P,0.0,FP);
    matmul("NT",nx,nx,nx,1.0,FP,F,0.0,P);

    for (i=0;i<nx;i++) {
        rtk->x[ix[i]]=xp[i];
        for (j=0;j<nx;j++) {
            rtk->P[ix[i]+ix[j]*rtk->nx]=P[i+j*nx];
        }
    }
    /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3])*fabs(tt);
    Q[8]=SQR(rtk->opt.prn[4])*fabs(tt);
    ecef2pos(rtk->x,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
            rtk->P[i+6+(j+6)*rtk->nx]+=Qv[i+j*3];
        }
    free(ix); free(F); free(P);
    free(FP); free(x); free(xp);
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udion(rtk_t *rtk, double tt, double bl, const int *sat, int ns)
{
    double el,fact,*x=NULL;
    int i,j,tc=0;
    insstate_t *ins=NULL;

    trace(3,"udion   : tt=%.3f bl=%.0f ns=%d\n",tt,bl,ns);

    ins=&rtk->ins;
    x=(tc=rtk->opt.mode==PMODE_INS_TGNSS?1:0)?rtk->ins.x:rtk->x;

    for (i=1;i<=MAXSAT;i++) {
        j=tc?xiIo(&rtk->opt.insopt,i):II(i,&rtk->opt);
        if (x[j]!=0.0&&
            rtk->ssat[i-1].outc[0]>GAP_RESION&&
            rtk->ssat[i-1].outc[1]>GAP_RESION)
            x[j]=0.0;
    }
    for (i=0;i<ns;i++) {
        j=tc?xiIo(&rtk->opt.insopt,sat[i]):II(sat[i],&rtk->opt);

        if (x[j]==0.0) {
            tc?
            insinitx(&rtk->ins,1E-6,SQR(rtk->opt.std[1]*bl/1E4),j):
            initx(rtk,1E-6,SQR(rtk->opt.std[1]*bl/1E4),j);
        }
        else {
            /* elevation dependent factor of process noise */
            el=rtk->ssat[sat[i]-1].azel[1];
            fact=cos(el);
            tc?
                    ins->P[j+j*ins->nx]+=SQR(rtk->opt.prn[1]*bl/1E4*fact)*fabs(tt):
                    rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1]*bl/1E4*fact)*fabs(tt);
        }
    }
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop(rtk_t *rtk, double tt, double bl)
{
    int i,j,k,tc=0,n;
    double *x=NULL,*P=NULL;
    insstate_t *ins=NULL;

    trace(3,"udtrop  : tt=%.3f\n",tt);

    ins=&rtk->ins;
    tc=rtk->opt.mode==PMODE_INS_TGNSS;
    x=tc?rtk->ins.x:rtk->x;
    P=tc?rtk->ins.P:rtk->P;
    n=tc?rtk->ins.nx:rtk->nx;

    for (i=0;i<2;i++) {
        j=tc?xiTr((&rtk->opt.insopt),i):IT(i,&rtk->opt);
        if (x[j]==0.0) {
            tc?
            insinitx(ins,INIT_ZWD,SQR(rtk->opt.std[2]),j):
            initx(rtk,INIT_ZWD,SQR(rtk->opt.std[2]),j); /* initial zwd */

            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) {
                    tc?
                    insinitx(&rtk->ins,1E-6,VAR_GRA,++j):
                    initx(rtk,1E-6,VAR_GRA,++j);
                }
            }
        }
        else {
            P[j+j*n]+=SQR(rtk->opt.prn[2])*fabs(tt);

            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) {
                    P[++j*(1+n)]+=SQR(rtk->opt.prn[2]*0.3)*fabs(tt);
                }
            }
        }
    }
}
/* temporal update of receiver h/w biases ------------------------------------*/
static void udrcvbias(rtk_t *rtk, double tt)
{
    int i,j,tc,nx;
    double *x,*P;

    trace(3,"udrcvbias: tt=%.3f\n",tt);

    tc=rtk->opt.mode==PMODE_INS_TGNSS?1:0;
    nx=tc?rtk->ins.nx:rtk->nx;
    x=tc?rtk->ins.x:rtk->x;
    P=tc?rtk->ins.P:rtk->P;

    for (i=0;i<NFREQGLO;i++) {
        j=tc?xiLl((&rtk->opt.insopt),i):IL(i,&rtk->opt);

        if (x[j]==0.0) {
            tc?
            insinitx(&rtk->ins,1E-6,VAR_HWBIAS,j):
            initx(rtk,1E-6,VAR_HWBIAS,j);
        }
            /* hold to fixed solution */
        else if (rtk->nfix>=rtk->opt.minfix
                 &&rtk->sol.ratio>rtk->opt.thresar[0]) {
            tc?
            insinitx(&rtk->ins,rtk->ins.xb[j],rtk->ins.Pb[j+j*rtk->ins.nx],j):
            initx(rtk,rtk->xa[j],rtk->Pa[j+j*rtk->na],j);
        }
        else {
            P[j+j*nx]+=SQR(PRN_HWBIAS)*fabs(tt);
        }
    }
}
/* carrier-phase bias correction by fcb --------------------------------------*/
extern void corr_phase_bias_fcb(obsd_t *obs, int n, const nav_t *nav)
{
    int i,j,k;

    for (i=0;i<nav->nf;i++) {
        if (timediff(nav->fcb[i].te,obs[0].time)<-1E-3) continue;
        if (timediff(nav->fcb[i].ts,obs[0].time)> 1E-3) break;
        for (j=0;j<n;j++) {
            for (k=0;k<NFREQ;k++) {
                if (obs[j].L[k]==0.0) continue;
                obs[j].L[k]-=nav->fcb[i].bias[obs[j].sat-1][k];
            }
        }
        return;
    }
}
/* carrier-phase bias correction by ssr --------------------------------------*/
extern void corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav)
{
    double lam;
    int i,j,code;

    for (i=0;i<n;i++) for (j=0;j<NFREQ;j++) {
            if (!(code=obs[i].code[j])) continue;
            if ((lam=nav->lam[obs[i].sat-1][j])==0.0) continue;

            /* correct phase bias (cyc) */
            obs[i].L[j]-=nav->ssr[obs[i].sat-1].pbias[code-1]/lam;
        }
}
/* frequency (Hz) from nav wavelength table ----------------------------------*/
static double obs_freq(const obsd_t *obs, int f, const nav_t *nav)
{
    double lam;
    if (f<0||f>=NFREQ) return 0.0;
    lam=nav->lam[obs->sat-1][f];
    return lam>0.0?CLIGHT/lam:0.0;
}
/* detect receiver clock jump on one sat/freq (phase-code TD residual) -------
 * type: 1=code, 2=phase, 3=both. Does not touch slip flags. */
static int detclkjmp(rtk_t *rtk, const obsd_t *obs, int i, int rcv, int f,
                     const nav_t *nav, double *Si, int *type)
{
    int sat=obs[i].sat,rcv0=obs[i].rcv-1;
    double tt,dph,dpp,freq,L_real,P_real,th;

    *Si=0.0; *type=0;
    if (rcv0<0||rcv0>1) return 0;
    freq=obs_freq(obs+i,f,nav);
    if (freq<=0.0) return 0;
    if (obs[i].L[f]==0.0||obs[i].P[f]==0.0||rtk->ssat[sat-1].pP[rcv0][f]==0.0||
        rtk->ssat[sat-1].ph[rcv0][f]==0.0) return 0;
    if (fabs(tt=timediff(obs[i].time,rtk->ssat[sat-1].pt[rcv0][f]))<DTTOL) return 0;

    L_real=obs[i].L[f]+1.0e-3*(double)rtk->opt.Js_h[rcv-1]*freq;
    P_real=obs[i].P[f]+1.0e-3*(double)rtk->opt.Js_P[rcv-1]*CLIGHT;
    dph=(L_real-rtk->ssat[sat-1].ph[rcv0][f])*CLIGHT/freq;
    dpp=P_real-rtk->ssat[sat-1].pP[rcv0][f];
    *Si=dph-dpp;

    /* |Si| > ~0.5 ms * c  (receiver clock jump candidate) */
    if (fabs(*Si)<=CLIGHT*0.001*0.5-15.0) {
        *Si=0.0; *type=0;
        return 0;
    }
    th=CLIGHT*0.001-fabs(tt)*fabs(obs[i].D[f]);
    if (th<CLIGHT*0.001*0.3) th=CLIGHT*0.001*0.3;

    if (fabs(dpp)>=th&&fabs(dph)<th) *type=1;      /* code jump */
    else if (fabs(dpp)<th&&fabs(dph)>=th) *type=2; /* phase jump */
    else *type=3;                                  /* both / mixed */
    return 1;
}
/* repair receiver clock jump (update Js_h / Js_P) ---------------------------
 * Common-mode only: majority of sats must agree; never marks cycle slips.
 * Wrong Js (unsigned / sparse false alarms) destroys AR — keep call off unless
 * real ms-level clock jumps are confirmed in the dataset. */
static void repair_clkjmp(rtk_t *rtk, const obsd_t *obs, const int *sat,
                          const int *iu, const int *ir, int ns, const nav_t *nav)
{
    int i,j,idx,rcv,typ,tmp,n[2]={0},ntype[2][4]={{0}},Js,need,best;
    double Si,sum[2]={0},ms,frac;

    (void)sat;
    trace(3,"repair_clkjmp: ns=%d\n",ns);

    for (j=0;j<2;j++) {
        rcv=j+1;
        for (i=0;i<ns;i++) {
            idx=(rcv==1)?iu[i]:ir[i];
            tmp=0; Si=0.0;
            if (!detclkjmp(rtk,obs,idx,rcv,0,nav,&Si,&tmp)) continue;
            if (tmp<1||tmp>3) continue;
            sum[j]+=Si;
            n[j]++;
            ntype[j][tmp]++;
        }
        need=ns/2;
        if (need<4) need=4;
        if (n[j]<need) continue;

        best=1;
        for (typ=2;typ<=3;typ++) if (ntype[j][typ]>ntype[j][best]) best=typ;
        if (ntype[j][best]*2<n[j]) continue; /* no majority type */

        ms=sum[j]/n[j]/CLIGHT; /* seconds */
        frac=fabs(ms)*1000.0;
        if (fabs(frac-floor(frac+0.5))>0.25) continue; /* not near-integer ms */
        Js=(int)floor(ms*1000.0+(ms>=0.0?0.5:-0.5)); /* signed */
        if (Js==0) continue;

        if (best==1) rtk->opt.Js_P[j]+=Js;
        else if (best==2) rtk->opt.Js_h[j]+=Js;
        else { rtk->opt.Js_P[j]+=Js; rtk->opt.Js_h[j]+=Js; }

        trace(2,"repair_clkjmp: rcv=%d type=%d Js=%d n=%d/%d -> h=%d P=%d\n",
              rcv,best,Js,n[j],ns,rtk->opt.Js_h[j],rtk->opt.Js_P[j]);
    }
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int i, int rcv)
{
    unsigned int slip,LLI;
    int f,sat=obs[i].sat;

    trace(3,"detslp_ll: i=%d rcv=%d\n",i,rcv);

    for (f=0;f<rtk->opt.nf;f++) {

        if (obs[i].L[f]==0.0||
            fabs(timediff(obs[i].time,rtk->ssat[sat-1].pt[rcv-1][f]))<DTTOL) {
            continue;
        }
        /* restore previous LLI */
        if (rcv==1) LLI=getbitu(&rtk->ssat[sat-1].slip[f],0,2); /* rover */
        else        LLI=getbitu(&rtk->ssat[sat-1].slip[f],2,2); /* base  */

        /* detect slip by cycle slip flag in LLI */
        if (rtk->tt>=0.0) { /* forward */
            if (obs[i].LLI[f]&1) {
                errmsg(rtk,"slip detected forward  (sat=%2d rcv=%d F=%d LLI=%x)\n",
                       sat,rcv,f+1,obs[i].LLI[f]);
            }
            slip=obs[i].LLI[f];
        }
        else { /* backward */
            if (LLI&1) {
                errmsg(rtk,"slip detected backward (sat=%2d rcv=%d F=%d LLI=%x)\n",
                       sat,rcv,f+1,LLI);
            }
            slip=LLI;
        }
        /* detect slip by parity unknown flag transition in LLI */
        if (((LLI&2)&&!(obs[i].LLI[f]&2))||(!(LLI&2)&&(obs[i].LLI[f]&2))) {
            errmsg(rtk,"slip detected half-cyc (sat=%2d rcv=%d F=%d LLI=%x->%x)\n",
                   sat,rcv,f+1,LLI,obs[i].LLI[f]);
            slip|=1;
        }
        /* save current LLI */
        if (rcv==1) setbitu(&rtk->ssat[sat-1].slip[f],0,2,obs[i].LLI[f]);
        else        setbitu(&rtk->ssat[sat-1].slip[f],2,2,obs[i].LLI[f]);

        /* save slip and half-cycle valid flag */
        rtk->ssat[sat-1].slip[f]|=(unsigned char)slip;
        rtk->ssat[sat-1].half[f]=(obs[i].LLI[f]&2)?0:1;
    }
}
/* detect cycle slip by L1-L2 geometry free phase jump -----------------------*/
static void detslp_gf_L1L2(rtk_t *rtk, const obsd_t *obs, int i, int j,
                           const nav_t *nav)
{
    int sat=obs[i].sat;
    double g0,g1;

    trace(3,"detslp_gf_L1L2: i=%d j=%d\n",i,j);

    if (rtk->opt.nf<=1||(g1=gfobs_L1L2(obs,i,j,nav->lam[sat-1]))==0.0) return;

    g0=rtk->ssat[sat-1].gf; rtk->ssat[sat-1].gf=g1;
    if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {

        rtk->ssat[sat-1].slip[0]|=1;
        rtk->ssat[sat-1].slip[1]|=1;
        errmsg(rtk,"slip detected GF2_jump (sat=%2d GF_L1_L2=%.3f %.3f)\n",sat,g0,g1);
    }
}
/* detect cycle slip by L1-L5 geometry free phase jump -----------------------*/
static void detslp_gf_L1L5(rtk_t *rtk, const obsd_t *obs, int i, int j,
                           const nav_t *nav)
{
    int sat=obs[i].sat;
    double g0,g1;

    trace(3,"detslp_gf_L1L5: i=%d j=%d\n",i,j);

    if (rtk->opt.nf<=2||(g1=gfobs_L1L5(obs,i,j,nav->lam[sat-1]))==0.0) return;

    g0=rtk->ssat[sat-1].gf2; rtk->ssat[sat-1].gf2=g1;
    if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {

        rtk->ssat[sat-1].slip[0]|=1;
        rtk->ssat[sat-1].slip[2]|=1;
        errmsg(rtk,"slip detected GF5_jump (sat=%2d GF_L1_L5=%.3f %.3f)\n",sat,g0,g1);
    }
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
static double mwmeas(const obsd_t *obs, int iu, int ir, const nav_t *nav)
{
    const double *lam=nav->lam[obs->sat-1];
    double l1,l2,p1,p2;
    int i=(satsys(obs->sat,NULL)&(SYS_GAL|SYS_SBS))?2:1;

    if (obs[ir].L[0]==0.0||obs[iu].L[0]==0.0) return 0.0;
    if (obs[ir].P[i]==0.0||obs[iu].P[i]==0.0) return 0.0;

    l1=sdobs(obs,iu,ir,0);
    l2=sdobs(obs,iu,ir,i);
    p1=sdobs(obs,iu,ir,0+NFREQ);
    p2=sdobs(obs,iu,ir,i+NFREQ);

    return lam[0]*lam[i]*(l1-l2)/(lam[i]-lam[0])-
           (lam[i]*p1+lam[0]*p2)/(lam[i]+lam[0]);
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
//mw是利用L1和L2频点的伪距和载波观测数据组合进行周跳探测的一种方法，只能探测双频的数据。
//该方法由于包含伪距噪声，故也只能探测大的周跳值。mw方法适用于在宽项固定中进行周跳探测。
static void detslp_mw(rtk_t *rtk, const obsd_t *obs, int iu, int ir,
                      int sat,const nav_t *nav)
{
    double w0,w1;
    int j;

    trace(3,"detslp_mw:\n");

    if ((w1=mwmeas(obs,iu,ir,nav))==0.0) return;

    w0=rtk->ssat[sat-1].mw;
    rtk->ssat[sat-1].mw=w1;

    trace(4,"detslip_mw: sat=%2d mw0=%8.3f mw1=%8.3f\n",sat,w0,w1);

    if (w0!=0.0&&fabs(w1-w0)>THRES_MW_JUMP) {
        errmsg(rtk,"detslip_mw: slip detected sat=%2d mw=%8.3f->%8.3f\n",
               sat,w0,w1);
        for (j=0;j<rtk->opt.nf;j++) rtk->ssat[sat-1].slip[j]|=1;
    }
}
/* SD epoch phase residual for cycle-slip detection
 * rs: sat positions aligned with obs (ecef+vel, 6 per sat), like udbias/satposs
 * return nv = n_use + n_constraint; *n_use = real phase obs count */
static int resphase(int iter, rtk_t *rtk, const obsd_t *obs, int n,
	const int *iu, const int *ir, int rcv, const nav_t *nav, const double *rs,
	const double *x, double *v, double *H, int *vsat, double *resh,
	int *n_use, double *vv)
{
	gtime_t time;
	double r, pr, freq=0.0, rr[3], pos[3], dtr, e[3], pe[3], azel_[2];
	double *H_, *H1, v1, sig, Ls, li=0.0, mi=0.0, ni=0.0, dt=0.0, lam, el;
	const double *rsi;
	int i, j, f, f_use, nv=0, sat, sys, idx, rcv0, ncon=0, mask[NX_SD-3]={0};
	char id[8];

	H_=mat(NX_SD,NX_SD); H1=mat(NX_SD,1);
	if (!H_||!H1) { free(H_); free(H1); return 0; }

	if (rcv==1) { for (i=0;i<3;i++) rr[i]=rtk->sol.rr[i]; }
	else if (rcv==2) { for (i=0;i<3;i++) rr[i]=rtk->rb[i]; }
	else { free(H_); free(H1); return 0; }
	/* static base: rb must be valid ECEF (else LOS geometry is garbage) */
	if (rcv==2&&norm(rr,3)<1.0) {
		trace(2,"resphase rcv=2 skip: base position unset\n");
		free(H_); free(H1);
		return 0;
	}
	dtr=x[3];
	ecef2pos(rr,pos);

	for (i=*n_use=0;i<n&&i<MAXOBS;i++) {
		if (rcv==1) {
			idx=iu[i]; sat=obs[idx].sat; time=obs[idx].time;
			vsat[idx]=0; resh[idx]=0.0;
		}
		else {
			idx=ir[i]; sat=obs[idx].sat; time=obs[idx].time;
			vsat[idx]=0; resh[idx]=0.0;
		}
		if (!(sys=satsys(sat,NULL))) continue;

		/* reject duplicated observation data */
		if (i<n-1&&i<MAXOBS-1) {
			int idx2=(rcv==1)?iu[i+1]:ir[i+1];
			if (sat==obs[idx2].sat) {
				trace(2,"duplicated obs data %s sat=%d\n",time_str(time,3),sat);
				i++; continue;
			}
		}

		rsi=rs+idx*6;
		if (norm(rsi,3)<1.0) continue; /* no sat position -> skip (avoid identical LOS) */
		if ((r=geodist(rsi,rr,e))<=0.0) continue;

		rcv0=obs[idx].rcv-1;
		if (rcv0<0||rcv0>1) continue;
		/* previous-epoch sat ECEF (same as RTKLIB-rtk-update ssat.prs) */
		if (norm(rtk->ssat[sat-1].prs[rcv0],3)<1.0) continue;
		if ((pr=geodist(rtk->ssat[sat-1].prs[rcv0],rr,pe))<=0.0) continue;

		/* elev from THIS receiver geometry (ssat.azel is rover-only / often unset) */
		el=satazel(pos,e,azel_);
		if (iter>0&&el<rtk->opt.elmin) continue;

		lam=nav->lam[sat-1][0];
		if (lam<=0.0) continue;

		Ls=0.0; freq=0.0; f_use=-1;
		for (f=0;f<rtk->opt.nf;f++) {
			if (obs[idx].L[f]==0.0||rtk->ssat[sat-1].ph[rcv0][f]==0.0) continue;
			lam=nav->lam[sat-1][f];
			if (lam<=0.0) continue;
			freq=CLIGHT/lam;
			/* apply clock-jump phase repair (ms -> cycle) like RTKLIB-rtk-update */
			Ls=(obs[idx].L[f]+1.0e-3*(double)rtk->opt.Js_h[rcv0]*freq
			    -rtk->ssat[sat-1].ph[rcv0][f])*lam-(r-pr);
			f_use=f;
			break; /* use first valid frequency */
		}
		if (freq<=0.0||f_use<0) continue;

		li=(rsi[0]-rr[0])/r;
		mi=(rsi[1]-rr[1])/r;
		ni=(rsi[2]-rr[2])/r;
		v1=Ls-(-li*x[0]-mi*x[1]-ni*x[2]+dtr);

		for (j=0;j<NX_SD;j++) H1[j]=0.0;
		H1[0]=-li; H1[1]=-mi; H1[2]=-ni; H1[3]=1.0;

		if (sys==SYS_GLO) { v1-=x[4]; H1[4]=1.0; mask[1]=1; }
		else if (sys==SYS_GAL) { v1-=x[5]; H1[5]=1.0; mask[2]=1; }
		else if (sys==SYS_CMP) { v1-=x[6]; H1[6]=1.0; mask[3]=1; }
		else if (sys==SYS_IRN) { v1-=x[7]; H1[7]=1.0; mask[4]=1; }
		else mask[0]=1; /* GPS, QZSS and SBAS share one clock */

		dt=timediff(obs[idx].time,rtk->ssat[sat-1].pt[rcv0][f_use]);
		sig=varerr(sat,sys,el,0.0,dt,f_use,obs[idx].SNR[f_use]*0.25,&rtk->opt);
		if (sig<=0.0||isnan(sig)||isinf(sig)) sig=1.0;
		matmul("TN",NX_SD,NX_SD,1,1.0,H1,H1,0.0,H_);
		sig=1.0/sig;
		for (j=0;j<NX_SD*NX_SD;j++) H[j]+=H_[j]*sig;
		for (j=0;j<NX_SD;j++) v[j]+=H1[j]*v1*sig;
		(*vv)+=v1*sig*v1;
		vsat[idx]=1; resh[idx]=fabs(v1)*sqrt(sig);
		nv++; (*n_use)++;

		if (iter==0) {
			satno2id(sat,id);
			trace(4,"resphase rcv=%d sat=%s Ls=%7.3f v1=%7.3f resh=%7.3f el=%.1f dt=%.2f\n",
				rcv,id,Ls,v1,resh[idx],el*R2D,dt);
		}
	}
	/* constraint to avoid rank deficiency when a system is absent */
	for (i=0;i<NX_SD-3;i++) {
		if (mask[i]) continue;
		v1=0.0;
		for (j=0;j<NX_SD;j++) H1[j]=(j==i+3)?1.0:0.0;
		sig=0.01;
		matmul("TN",NX_SD,NX_SD,1,1.0,H1,H1,0.0,H_);
		sig=1.0/sig;
		for (j=0;j<NX_SD*NX_SD;j++) H[j]+=H_[j]*sig;
		for (j=0;j<NX_SD;j++) v[j]+=H1[j]*v1*sig;
		(*vv)+=v1*sig*v1;
		nv++; ncon++;
	}
	/* static base: Δxyz≈0 between epochs (mode != moving-base) */
	if (rcv==2&&rtk->opt.mode!=PMODE_MOVEB) {
		for (i=0;i<3;i++) {
			v1=0.0;
			for (j=0;j<NX_SD;j++) H1[j]=(j==i)?1.0:0.0;
			sig=0.01; /* soft constraint ~1 cm */
			matmul("TN",NX_SD,NX_SD,1,1.0,H1,H1,0.0,H_);
			sig=1.0/sig;
			for (j=0;j<NX_SD*NX_SD;j++) H[j]+=H_[j]*sig;
			for (j=0;j<NX_SD;j++) v[j]+=H1[j]*v1*sig;
			(*vv)+=v1*sig*v1;
			nv++; ncon++;
		}
	}
	trace(3,"resphase rcv=%d iter=%d n_use=%d ncon=%d nv=%d vv=%.3f\n",
		rcv,iter,*n_use,ncon,nv,*vv);
	free(H_); free(H1);
	return nv;
}
/* return 1=no slip, 0=slip suspected; uses real obs count n_use for DOF */
static int slip_exists(int rcv, const double vv, int n_use, int nx)
{
	double sigma;
	(void)rcv;
	if (n_use<=nx) {
		trace(3,"slip_exists rcv=%d skip: n_use=%d <= nx=%d\n",rcv,n_use,nx);
		return 0;
	}
	sigma=sqrt(vv/(n_use-nx));
	trace(3,"slip_exists rcv=%d sigma=%.3f thr=3.0 n_use=%d vv=%.3f -> %s\n",
		rcv,sigma,n_use,vv,sigma>3.0?"SLIP":"OK");
	if (sigma>3.0) return 0;
	return 1;
}
/* return: 1=clean, 0=slip, -1=fail (lack obs / singular) */
static int estposbias(rtk_t *rtk, const obsd_t *obs, const int *iu, const int *ir, int rcv,
	const nav_t *nav, const double *rs, int ns, double *resh, int *vsat,
	int *n_use, int *nv)
{
	double *x,*dx,*v,*H,*Q,vv,ndx;
	int i,j,info,stat=-1;

	trace(3,"estposbias rcv=%d ns=%d\n",rcv,ns);

	x=zeros(NX_SD,1); dx=mat(NX_SD,1); v=mat(NX_SD,1); H=mat(NX_SD,NX_SD); Q=mat(NX_SD,NX_SD);
	if (!x||!dx||!v||!H||!Q) {
		free(x); free(dx); free(v); free(H); free(Q);
		return -1;
	}

	for (i=0;i<MAXITR_SD;i++) {
		memset(H,0,sizeof(double)*NX_SD*NX_SD);
		memset(v,0,sizeof(double)*NX_SD);
		vv=0.0;

		*nv=resphase(i,rtk,obs,ns,iu,ir,rcv,nav,rs,x,v,H,vsat,resh,n_use,&vv);

		/* need enough REAL observations (constraints do not count) */
		if (*n_use<=NX_SD) {
			trace(2,"estposbias rcv=%d lack obs n_use=%d nv=%d (need >%d)\n",
				rcv,*n_use,*nv,NX_SD);
			stat=-1;
			break;
		}
		matcpy(Q,H,NX_SD,NX_SD);
		if ((info=matinv(Q,NX_SD))) {
			trace(2,"estposbias rcv=%d matinv error info=%d n_use=%d\n",rcv,info,*n_use);
			stat=-1;
			break;
		}
		matmul("NN",NX_SD,1,NX_SD,1.0,Q,v,0.0,dx);
		for (j=0;j<NX_SD;j++) x[j]+=dx[j];
		ndx=norm(dx,NX_SD);
		trace(4,"estposbias rcv=%d iter=%d |dx|=%.3e\n",rcv,i,ndx);
		if (ndx<1E-4) {
			stat=slip_exists(rcv,vv,*n_use,NX_SD);
			trace(3,"estposbias rcv=%d converged iter=%d result=%s\n",
				rcv,i,stat==1?"clean":"slip");
			free(x); free(dx); free(v); free(H); free(Q);
			return stat;
		}
	}
	if (stat==-1&&i>=MAXITR_SD) {
		trace(2,"estposbias rcv=%d not converged after %d iter |dx|=%.3e\n",
			rcv,MAXITR_SD,norm(dx,NX_SD));
	}
	free(x); free(dx); free(v); free(H); free(Q);
	return -1;
}

static int cmpresh(const void *p1, const void *p2)
{
	double *q1=(double *)p1,*q2=(double *)p2;
	return (int)(*q2)-(int)(*q1); /* q2-q1: 由大到小 */
}

/* slip_sd: flat [MAXSAT*NFREQ], index = (sat-1)*NFREQ+f; may be NULL */
static void detslp_sd_epoch(rtk_t *rtk, const obsd_t *obs, const int *iu, const int *ir,
	int rcv, const nav_t *nav, const double *rs, int ns, uint8_t *slip_sd)
{
	double *resh_e,*resh_e_max;
	int i,j,k=0,m=0,f=0,ret,n_use=0,nv=0,reject_flag=0,nslip=0;
	int *vsat_e,*i_reject,*iu_e,*ir_e;
	char code[20];

	trace(3,"detslp_sd_epoch: %s rcv=%d ns=%d\n",time_str(obs[0].time,0),rcv,ns);
	if (ns<=0||rs==NULL) return;

	resh_e=zeros(MAXOBS*2,1);
	resh_e_max=zeros(MAXOBS,1);
	vsat_e=imat(MAXOBS*2,1);
	i_reject=imat(MAXOBS*2,1);
	iu_e=imat(ns>0?ns:1,1);
	ir_e=imat(ns>0?ns:1,1);
	if (!resh_e||!resh_e_max||!vsat_e||!i_reject||!iu_e||!ir_e) {
		free(resh_e); free(resh_e_max); free(vsat_e); free(i_reject); free(iu_e); free(ir_e);
		return;
	}
	memset(vsat_e,0,sizeof(int)*MAXOBS*2);
	memset(i_reject,0,sizeof(int)*MAXOBS*2);

	ret=estposbias(rtk,obs,iu,ir,rcv,nav,rs,ns,resh_e,vsat_e,&n_use,&nv);
	if (ret==1) {
		trace(3,"detslp_sd_epoch rcv=%d clean n_use=%d nv=%d\n",rcv,n_use,nv);
		free(resh_e); free(resh_e_max); free(vsat_e); free(i_reject); free(iu_e); free(ir_e);
		return;
	}
	if (ret<0) {
		/* lack of obs / singular: NOT a slip — skip reject loop */
		trace(2,"detslp_sd_epoch rcv=%d skip (est fail) n_use=%d nv=%d ns=%d\n",
			rcv,n_use,nv,ns);
		free(resh_e); free(resh_e_max); free(vsat_e); free(i_reject); free(iu_e); free(ir_e);
		return;
	}

	/* ret==0: residual sigma too large — try exclude largest residuals */
	trace(3,"detslp_sd_epoch rcv=%d slip_suspected n_use=%d -> exclude loop\n",rcv,n_use);
	for (i=1;i<ns;i++) {
		if (n_use-i<=NX_SD) {
			trace(3,"detslp_sd_epoch rcv=%d exclude stop: remain n_use~%d\n",rcv,n_use-i);
			break;
		}
		memset(iu_e,0,sizeof(int)*ns);
		memset(ir_e,0,sizeof(int)*ns);
		memset(resh_e_max,0,sizeof(double)*MAXOBS);
		for (j=0;j<ns;j++) {
			if (rcv==1) resh_e_max[j]=resh_e[iu[j]];
			else if (rcv==2) resh_e_max[j]=resh_e[ir[j]];
		}
		qsort(resh_e_max,ns,sizeof(double),cmpresh);
		for (k=0;k<i;k++) for (j=0;j<ns;j++) {
			if (rcv==1&&resh_e[iu[j]]==resh_e_max[k]) {
				i_reject[k]=iu[j];
				continue;
			}
			else if (rcv==2&&resh_e[ir[j]]==resh_e_max[k]) {
				i_reject[k]=ir[j];
				continue;
			}
		}
		for (j=k=0;j<ns;j++) {
			reject_flag=0;
			for (m=0;m<i;m++) {
				if ((rcv==1&&iu[j]==i_reject[m])||
					(rcv==2&&ir[j]==i_reject[m])) reject_flag=1;
			}
			if (reject_flag) continue;
			if (rcv==1) iu_e[k]=iu[j];
			else if (rcv==2) ir_e[k]=ir[j];
			k++;
		}
		ret=estposbias(rtk,obs,iu_e,ir_e,rcv,nav,rs,ns-i,resh_e,vsat_e,&n_use,&nv);
		if (ret!=1) {
			trace(4,"detslp_sd_epoch rcv=%d exclude i=%d still bad ret=%d\n",rcv,i,ret);
			continue;
		}
		for (m=0;m<i;m++) {
			satno2id(obs[i_reject[m]].sat,code);
			trace(2,"detslp_sd_epoch SLIP sat=%s(%d) rcv=%d resh=%.3f nrej=%d\n",
				code,obs[i_reject[m]].sat,rcv,resh_e_max[m],i);
			for (f=0;f<rtk->opt.nf;f++) {
				rtk->ssat[obs[i_reject[m]].sat-1].slip[f]|=1;
				if (slip_sd) slip_sd[(obs[i_reject[m]].sat-1)*NFREQ+f]=1;
			}
			errmsg(rtk,"slip detected sd epoch (sat=%2d %s, rcv=%d)\n",
				obs[i_reject[m]].sat,code,rcv);
			nslip++;
		}
		trace(3,"detslp_sd_epoch rcv=%d fixed after exclude %d sat(s)\n",rcv,i);
		break;
	}
	if (nslip==0) {
		trace(2,"detslp_sd_epoch rcv=%d slip suspected but exclude failed\n",rcv);
	}
	free(resh_e); free(resh_e_max); free(vsat_e); free(i_reject); free(iu_e); free(ir_e);
}

/* test # frequencys for observation data-------------------------------------*/
static int chkfrq(const obsd_t *obs,const prcopt_t *opt)
{
    // /* only check two frequency for phase observation data */
    // int i; for (i=0;i<NF(opt);i++) if (obs->L[i]==0.0) return 0; return 1;
    int i; for (i=0;i<NF(opt);i++) if (obs->L[i]>0.0) return 1; return 0;
}
/* temporal update of phase biases -------------------------------------------*/
//周跳检测与载波相位初始值的估计
static void udbias(rtk_t *rtk, double tt, const obsd_t *obs, const int *sat, const double *rs,
                   const int *iu, const int *ir, int ns, const nav_t *nav)
{
    double cp,pr,cp1,cp2,pr1,pr2,*bias,lami,lam1,lam2,C1,C2;
    double *x,*P,eps=rtk->opt.eps;
    int i,j,f,slip,reset,nf=NF(&rtk->opt),nx,tc=0,ib;

    trace(3,"udbias  : tt=%.3f ns=%d\n",tt,ns);

    tc=rtk->opt.mode==PMODE_INS_TGNSS?1:0;
    tc?x=rtk->ins.x:x=rtk->x;
    tc?P=rtk->ins.P:P=rtk->P;
    tc?nx=rtk->ins.nx:nx=rtk->nx;

    /* clear slip flags, keep the saved base LLI in bits 2-3 (Qt udbias) */
    for (i=0;i<ns;i++) {
        for (f=0;f<rtk->opt.nf;f++) rtk->ssat[sat[i]-1].slip[f]&=0xFC;
    }
#if 0
    /* clock-jump repair OFF (20260729_2205). Qt udbias does not call it. */
    repair_clkjmp(rtk,obs,sat,iu,ir,ns,nav);
#endif

    for (i=0;i<ns;i++) {

        /* detect cycle slip by LLI */
        detslp_ll(rtk,obs,iu[i],1);
        detslp_ll(rtk,obs,ir[i],2);

        /* detect cycle slip by geometry-free phase jump */
        detslp_gf_L1L2(rtk,obs,iu[i],ir[i],nav);
        detslp_gf_L1L5(rtk,obs,iu[i],ir[i],nav);

#if 0
        /* MW and epoch-SD slips are not in Qt udbias (LLI + GF only). */
        detslp_mw(rtk,obs,iu[i],ir[i],sat[i],nav);
#endif

        /* update half-cycle valid flag */
        for (f=0;f<nf;f++) {
            rtk->ssat[sat[i]-1].half[f]=(unsigned char)!((obs[iu[i]].LLI[f]&2)||(obs[ir[i]].LLI[f]&2));
        }
    }

#if 0
    /* epoch SD residual slips are not in Qt udbias */
    detslp_sd_epoch(rtk,obs,iu,ir,1,nav,rs,ns,NULL);
    detslp_sd_epoch(rtk,obs,iu,ir,2,nav,rs,ns,NULL);
#endif

    /* store current sat ECEF -> ssat.prs for next-epoch TD (r-pr), same as RTKLIB-rtk-update */
    for (i=0;i<ns;i++) {
        for (j=0;j<3;j++) {
            if (obs[iu[i]].L[0]!=0.0)
                rtk->ssat[obs[iu[i]].sat-1].prs[obs[iu[i]].rcv-1][j]=rs[iu[i]*6+j];
            if (obs[ir[i]].L[0]!=0.0)
                rtk->ssat[obs[ir[i]].sat-1].prs[obs[ir[i]].rcv-1][j]=rs[ir[i]*6+j];
        }
    }

    for (f=0;f<nf;f++) {
        /* reset phase-bias if instantaneous AR or expire obs outage counter */
        for (i=1;i<=MAXSAT;i++) {

            if (tc&&rtk->ins.ptct.time!=0) {
                rtk->ssat[i-1].outc[f]+=ROUND(timediff(rtk->sol.time,rtk->ins.ptct)/eps);
            }
            else rtk->ssat[i-1].outc[f]+=ROUND(rtk->tt/eps);

            reset=rtk->ssat[i-1].outc[f]>(unsigned int)rtk->opt.maxout;
            ib=tc?xiBs((&rtk->opt.insopt),i,f):IB(i,f,(&rtk->opt));

            if (rtk->opt.modear==ARMODE_INST&&x[ib]!=0.0) {
                tc?
                insinitx(&rtk->ins,0.0,SQR(rtk->opt.std[0]),ib):
                initx(rtk,0.0,SQR(rtk->opt.std[0]),ib);
            }
            else if (reset&&x[ib]!=0.0) {
                tc?
                insinitx(&rtk->ins,0.0,0.0,ib):
                initx(rtk,0.0,SQR(rtk->opt.std[0]),ib);
                trace(3,"obs outage counter overflow (sat=%3d L%d n=%d)\n",i,f+1,rtk->ssat[i-1].outc[f]);
                rtk->ssat[i-1].outc[f]=0;
            }
            if (rtk->opt.modear!=ARMODE_INST&&reset) {
                rtk->ssat[i-1].lock[f]=-rtk->opt.minlock;
            }
        }
        /* reset phase-bias if detecting cycle slip */
        for (i=0;i<ns;i++) {
            j=tc?xiBs((&rtk->opt.insopt),sat[i],f):IB(sat[i],f,&rtk->opt);

            P[j+j*nx]+=rtk->opt.prn[0]*rtk->opt.prn[0]*fabs(tt);
            if (rtk->ssat[sat[i]-1].sfrq[0]) {
                P[j+j*nx]+=SQR(0.5);
            }
            slip=rtk->ssat[sat[i]-1].slip[f];
            if (rtk->opt.ionoopt==IONOOPT_IFLC) slip|=rtk->ssat[sat[i]-1].slip[1];
            if (rtk->opt.modear==ARMODE_INST||!(slip&1)) continue;
            x[j]=0.0;
            rtk->ssat[sat[i]-1].lock[f]=-rtk->opt.minlock;
        }
        bias=zeros(ns,1);

        /* estimate approximate phase-bias by phase - code */
        for (i=0;i<ns;i++) {

            if (rtk->opt.ionoopt!=IONOOPT_IFLC) {
                cp=sdobs(obs,iu[i],ir[i],f); /* cycle */
                pr=sdobs(obs,iu[i],ir[i],f+NFREQ);

                rtk->ssat[sat[i]-1].sdi[f]=sdobspins(rtk,obs,rs,iu[i],ir[i],f+NFREQ);
                rtk->ssat[sat[i]-1].sdp[f]=pr;
                rtk->ssat[sat[i]-1].sdc[f]=cp;
                lami=nav->lam[sat[i]-1][f];
                if (cp==0.0||pr==0.0||lami<=0.0) continue;

                bias[i]=cp-pr/lami;
            }
            else {
                cp1=sdobs(obs,iu[i],ir[i],0);
                cp2=sdobs(obs,iu[i],ir[i],1);
                pr1=sdobs(obs,iu[i],ir[i],NFREQ);
                pr2=sdobs(obs,iu[i],ir[i],NFREQ+1);
                lam1=nav->lam[sat[i]-1][0];
                lam2=nav->lam[sat[i]-1][1];

                if (cp1==0.0||cp2==0.0||
                    pr1==0.0||pr2==0.0||
                    lam1<=0.0||lam2<=0.0) continue;

                C1= SQR(lam2)/(SQR(lam2)-SQR(lam1));
                C2=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
                bias[i]=(C1*lam1*cp1+C2*lam2*cp2)-(C1*pr1+C2*pr2);
            }
        }
        /* correct phase-bias offset so a new satellite stays coherent with
           the ones already in the filter (Qt udbias). Kinematic only. */
        if (!tc) {
            int ncoh=0;
            double offset=0.0;
            for (i=0;i<ns;i++) {
                ib=IB(sat[i],f,&rtk->opt);
                if (bias[i]==0.0||x[ib]==0.0) continue;
                offset+=bias[i]-x[ib];
                ncoh++;
            }
            if (ncoh>0) {
                for (i=1;i<=MAXSAT;i++) {
                    ib=IB(i,f,&rtk->opt);
                    if (x[ib]!=0.0) x[ib]+=offset/ncoh;
                }
            }
        }
        /* set initial states of phase-bias */
        for (i=0;i<ns;i++) {

            tc?ib=xiBs((&rtk->opt.insopt),sat[i],f):ib=IB(sat[i],f,&rtk->opt);

            if (bias[i]==0.0||x[ib]!=0.0) continue;
            tc?
            insinitx(&rtk->ins,bias[i],SQR(rtk->opt.std[0]),ib):
            initx(rtk,bias[i],SQR(rtk->opt.std[0]),ib);
            /* wait minlock epochs (pos2-arlockcnt) before AR */
            rtk->ssat[sat[i]-1].lock[f]=-rtk->opt.minlock;
        }
        free(bias);
    }
}
/* ins/gnss couple state update------------------------------------------------*/
static void insudsta(rtk_t *rtk,double tt)
{
    trace(3,"insudsta : tt=%6.4lf\n",tt);
}
/* reset estimated filter if counts of solution valid fail is many-------------*/
static void udreset(rtk_t *rtk)
{
    int i,j,f;
    if (rtk->opt.mode<=PMODE_FIXED) {
        for (i=0;i<3;i++) initx(rtk,rtk->x[i],VAR_POS,i);
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) initx(rtk,rtk->x[i],VAR_VEL,i);
            for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        }
        for (i=1;i<=MAXSAT;i++) for (f=0;f<NF(&rtk->opt);f++) {
                j=IB(i,f,&rtk->opt);
                initx(rtk,rtk->x[j],VAR_AMB,j);
            }
    }
    else if (rtk->opt.mode>=PMODE_INS_LGNSS) {
        getP0(&rtk->opt.insopt,rtk->ins.P );
        getP0(&rtk->opt.insopt,rtk->ins.Pa);
    }
}
/* temporal update of states --------------------------------------------------*/
static void udstate(rtk_t *rtk, const obsd_t *obs, const int *sat, const double *rs,
                    const int *iu, const int *ir, int ns, const nav_t *nav)
{
    int tc;
    double tt=rtk->tt,bl,dr[3];

    trace(3,"udstate : ns=%d\n",ns);

    /* temporal update of position/velocity/acceleration */
    rtk->opt.mode==PMODE_INS_TGNSS?insudsta(rtk,tt):udpos(rtk,tt);

    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt>=IONOOPT_EST) {
        bl=baseline(rtk->x,rtk->rb,dr);
        udion(rtk,tt,bl,sat,ns);
    }
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt>=TROPOPT_EST) {
        bl=baseline(rtk->x,rtk->rb,dr);
        udtrop(rtk,tt,bl);
    }
    /* temporal update of eceiver h/w bias */
    if (rtk->opt.glomodear==2&&(rtk->opt.navsys&SYS_GLO)) {
        udrcvbias(rtk,tt);
    }
    /* tc-mode */
    tc=rtk->opt.mode==PMODE_INS_TGNSS&&rtk->opt.insopt.tc==INSTC_RTK;

    /* temporal update of phase-bias */
    if ((rtk->opt.mode>PMODE_DGPS&&rtk->opt.mode<PMODE_INS_UPDATE)||tc) {
        udbias(rtk,tt,obs,sat,rs,iu,ir,ns,nav);
    }
#if 0
    /* check counts of solution valid fail */
    if (rtk->failc>=MINFAILC) udreset(rtk);
#endif
}
/* undifferenced phase/code residual for satellite ---------------------------*/
static void zdres_sat(int base, double r, const obsd_t *obs, const nav_t *nav,
                      const double *azel, const double *dant,
                      const prcopt_t *opt, double *y)
{
    const double *lam=nav->lam[obs->sat-1];
    double f1,f2,C1,C2,dant_if;
    int i,nf=NF(opt);

    if (opt->ionoopt==IONOOPT_IFLC) { /* iono-free linear combination */
        if (lam[0]==0.0||lam[1]==0.0) return;

        if (testsnr(base,0,azel[1],obs->SNR[0]*0.25,&opt->snrmask)||
            testsnr(base,1,azel[1],obs->SNR[1]*0.25,&opt->snrmask)) return;

        f1=CLIGHT/lam[0];
        f2=CLIGHT/lam[1];
        C1= SQR(f1)/(SQR(f1)-SQR(f2));
        C2=-SQR(f2)/(SQR(f1)-SQR(f2));
        dant_if=C1*dant[0]+C2*dant[1];

        if (obs->L[0]!=0.0&&obs->L[1]!=0.0) {
            y[0]=C1*obs->L[0]*lam[0]+C2*obs->L[1]*lam[1]-r-dant_if;
        }
        if (obs->P[0]!=0.0&&obs->P[1]!=0.0) {
            y[1]=C1*obs->P[0]+C2*obs->P[1]-r-dant_if;
        }
    }
    else {
        for (i=0;i<nf;i++) {
            if (lam[i]==0.0) continue;

            /* check snr mask */
            if (testsnr(base,i,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) {
                continue;
            }
            /* residuals = observable - pseudorange */
            if (obs->L[i]!=0.0) y[i   ]=obs->L[i]*lam[i]-r-dant[i];
            if (obs->P[i]!=0.0) y[i+nf]=obs->P[i]       -r-dant[i];
        }
    }
}
/* undifferenced phase/code residuals ----------------------------------------*/
static int zdres(int base, const obsd_t *obs, int n, const double *rs,
                 const double *dts, const int *svh, const nav_t *nav,
                 const double *rr, const prcopt_t *opt, int index, double *y,
                 double *e, double *azel)
{
    double r,rr_[3],pos[3],dant[NFREQ]={0},disp[3];
    double zhd,zazel[]={0.0,90.0*D2R};
    int i,nf=NF(opt);

    trace(3,"zdres   : n=%d\n",n);

    for (i=0;i<n*nf*2;i++) y[i]=0.0;

    if (norm(rr,3)<=0.0) return 0; /* no receiver position */

    for (i=0;i<3;i++) rr_[i]=rr[i];

    /* earth tide correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rr_,opt->tidecorr,&nav->erp,
                 opt->odisp[base],disp);
        for (i=0;i<3;i++) rr_[i]+=disp[i];
    }
    ecef2pos(rr_,pos);

    for (i=0;i<n;i++) {
        /* compute geometric-range and azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr_,e+i*3))<=0.0) continue;
        if (satazel(pos,e+i*3,azel+i*2)<opt->elmin) continue;

        /* excluded satellite? */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;

        /* satellite clock-bias */
        r+=-CLIGHT*dts[i*2];

        /* troposphere delay model (hydrostatic) */
        zhd=tropmodel(obs[0].time,pos,zazel,0.0);
        r+=tropmapf(obs[i].time,pos,azel+i*2,NULL)*zhd;

        /* receiver antenna phase center correction */
        antmodel(opt->pcvr+index,opt->antdel[index],azel+i*2,opt->posopt[1],
                 dant);

        /* undifferenced phase/code residual for satellite */
        zdres_sat(base,r,obs+i,nav,azel+i*2,dant,opt,y+i*nf*2);
    }
    trace(4,"rr_=%.3f %.3f %.3f\n",rr_[0],rr_[1],rr_[2]);
    trace(4,"pos=%.9f %.9f %.3f\n",pos[0]*R2D,pos[1]*R2D,pos[2]);
    for (i=0;i<n;i++) {
        trace(4,"sat=%2d %13.3f %13.3f %13.3f %13.10f %6.1f %5.1f\n",
              obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],dts[i*2],azel[i*2]*R2D,
              azel[1+i*2]*R2D);
    }
    trace(3,"y=\n"); tracemat(3,y,nf*2,n,15,6);

    return 1;
}
/* test valid observation data -----------------------------------------------*/
static int validobs(int i, int j, int f, int nf, double *y)
{
    /* if no phase observable, psudorange is also unusable */
    return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
           (f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
}
/* double-differenced measurement error covariance ---------------------------*/
static void ddcov(const int *nb, int n, const double *Ri, const double *Rj,
                  int nv, double *R)
{
    int i,j,k=0,b;

    trace(3,"ddcov   : n=%d\n",n);

    for (i=0;i<nv*nv;i++) R[i]=0.0;
    for (b=0;b<n;k+=nb[b++]) {

        for (i=0;i<nb[b];i++) for (j=0;j<nb[b];j++) {
                R[k+i+(k+j)*nv]=Ri[k+i]+(i==j?Rj[k+i]:0.0);
            }
    }
    trace(5,"R=\n");
    tracemat(5,R,nv,nv,12,6);
}
/* baseline length constraint ------------------------------------------------*/
static int constbl(rtk_t *rtk, const double *x, const double *P, double *v,
                   double *H, double *Ri, double *Rj, int index)
{
    const double thres=0.1; /* threshold for nonliearity (v.2.3.0) */
    double xb[3],b[3],bb,var=0.0,rr[3];
    int i,nx,tc=0,ix;

    trace(3,"constbl : \n");

    tc=rtk->opt.mode==PMODE_INS_TGNSS;

    /* number of estimate states */
    nx=tc?rtk->ins.nx:rtk->nx;
    ix=tc?xiP(&rtk->opt.insopt):0;

    /* rover station position */
    if (tc) {
        insp2antp(&rtk->ins,rr);
    }
    else matcpy(rr,x,1,3);

    /* no constraint */
    if (rtk->opt.baseline[0]<=0.0) return 0;

    /* time-adjusted baseline vector and length */
    for (i=0;i<3;i++) {
#if 0
        xb[i]=rtk->rb[i]+rtk->rb[i+3]*rtk->sol.age;
#else
        xb[i]=rtk->rb[i];
#endif
        b[i]=rr[i]-xb[i];
    }
    bb=norm(b,3);

    /* approximate variance of solution */
    if (P) {
        for (i=0;i<3;i++) var+=P[(ix+i)+(ix+i)*nx];
        var/=3.0;
    }
    /* check nonlinearity */
    if (var>thres*thres*bb*bb) {
        trace(3,"constbl : equation nonlinear (bb=%.3f var=%.3f)\n",bb,var);
        return 0;
    }
    /* constraint to baseline length */
    v[index]=rtk->opt.baseline[0]-bb;
    if (H) {
        for (i=0;i<3;i++) H[i+index*nx]=b[i]/bb;
    }
    Ri[index]=0.0;
    Rj[index]=SQR(rtk->opt.baseline[1]);

    trace(4,"baseline len v=%13.3f R=%8.6f %8.6f\n",v[index],Ri[index],Rj[index]);
    return 1;
}
/* precise tropspheric model -------------------------------------------------*/
static double prectrop(gtime_t time, const double *pos, int r,
                       const double *azel, const prcopt_t *opt, const double *x,
                       double *dtdx)
{
    double m_w=0.0,cotz,grad_n,grad_e;
    int i,flag=0;

    flag=opt->mode==PMODE_INS_TGNSS?1:0;
    flag?i=xiTr((&opt->insopt),r):i=IT(r,opt);

    /* wet mapping function */
    tropmapf(time,pos,azel,&m_w);

    if (opt->tropopt>=TROPOPT_ESTG&&azel[1]>0.0) {

        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[i+1]+grad_e*x[i+2];
        dtdx[1]=grad_n*x[i];
        dtdx[2]=grad_e*x[i];
    }
    else dtdx[1]=dtdx[2]=0.0;
    dtdx[0]=m_w;
    return m_w*x[i];
}
/* glonass inter-channel bias correction -------------------------------------*/
static double gloicbcorr(int sat1, int sat2, const prcopt_t *opt, double lam1,
                         double lam2, int f)
{
    double dfreq;

    if (f>=NFREQGLO||f>=opt->nf||!opt->exterr.ena[2]) return 0.0;
    dfreq=(CLIGHT/lam1-CLIGHT/lam2)/(f==0?DFRQ1_GLO:DFRQ2_GLO);
    return opt->exterr.gloicb[f]*0.01*dfreq; /* (m) */
}
/* test navi system (m=1:glo,2:gal,3:bds,4:gps/qzs/sbs,5:irn) -----------------*/
static int test_sys(int sys, int m)
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
/* get system index (m)-------------------------------------------------------*/
static int sysind(int sys)
{
    switch (sys) {
        case SYS_GPS: return 4;
        case SYS_SBS: return 4;
        case SYS_GLO: return 1;
        case SYS_GAL: return 2;
        case SYS_CMP: return 3;
        case SYS_QZS: return 4;
        case SYS_IRN: return 5;
    }
    return 4;
}
/* jacobian of double-differncee phase/code by ins position-------------------*/
static void jacob_dd_dp(const insstate_t *ins,const double *ei,const double *ej,
                        double *dddp)
{
    dddp[0]=ei[0]-ej[0];
    dddp[1]=ei[1]-ej[1];
    dddp[2]=ei[2]-ej[2];
}
/* jacobian of perturb rotation wrt. perturb euler angles--------------------*/
static void jacobian_prot_pang(const double *Cbe,double *S)
{
    double rpy[3]={0};
    c2rpy(Cbe,rpy);
    S[0]= cos(rpy[1])*cos(rpy[2]); S[3]=sin(rpy[2]); S[6]=0.0;
    S[1]=-cos(rpy[1])*sin(rpy[2]); S[4]=cos(rpy[2]); S[7]=0.0;
    S[2]= sin(rpy[1]);             S[5]=0.0;         S[8]=1.0;
}
/* jacobian of double-difference phase/code by ins attitude error-------------*/
static void jacob_dd_da(const insstate_t *ins,const double *ei,const double *ej,
                        double *ddda)
{
#if NOINSJACO
    return; /* no add jacobians */
#else
    int i;
    double sl[3],T[9];

    matmul("NN",3,1,3,1.0,ins->Cbe,ins->lever,0.0,sl);
    skewsym3(sl,T);

    for (i=0;i<3;i++) sl[i]=-ei[i]+ej[i];
    matmul("NN",1,3,3,1.0,sl,T,0.0,ddda);
#endif
}
/* jacobian of double-differnce phase/code by lever arm-----------------------*/
static void jacob_dd_dl(const insstate_t *ins,const double *ei,const double *ej,
                        double *dddl)
{
#if NOINSJACO
    return; /* no add jacobians */
#else
    double sl[3];
    sl[0]=ei[0]-ej[0];
    sl[1]=ei[1]-ej[1];
    sl[2]=ei[2]-ej[2];
    matmul("NN",1,3,3,1.0,sl,ins->Cbe,0.0,dddl);
#endif
}
/* check whether it is new double difference satellite-----------------------*/
static int chknews(const rtk_t *rtk,int sat1,int sat2,int freq)
{
    int i; for (i=0;i<rtk->ns;i++) {
        if (sat1==rtk->sat[i].sat1&&
            sat2==rtk->sat[i].sat2&&
            freq==rtk->sat[i].f) {
            return 0;
        }
    }
    if (rtk->ns==0) return 0; return 1;
}
/* double-differenced phase/code residuals -----------------------------------*/
static int ddres(rtk_t *rtk, const nav_t *nav, const obsd_t *obs, double dt,const double *x,
                 const double *P, const int *sat, double *y,
                 double *e, double *azel, const int *iu, const int *ir, int ns,
                 double *v, double *H, double *R, int *vflg,const double *re,
                 double *Ri_,double *Rj_,int *nb_,int *b_,
                 int refsat[NUMSYS][2*NFREQ],int usedualfrq)
{
    prcopt_t *opt=&rtk->opt;
    insopt_t *insopt=&opt->insopt;
    double bl,dr[3],posu[3],posr[3],didxi,didxj,*im,*vc,ddi,ddg,factor=1.0;
    double *tropr,*tropu,*dtdxr,*dtdxu,*Ri,*Rj,lami,lamj,fi,fj,df,*Hi=NULL,rr[3];
    double dp[3]={0},da[3]={0},dl[3]={0},S[9],dap[3];
    int i,j,k,m,f,ff,nv=0,nb[NFREQ*4*2+2]={0},b=0,sysi,sysj,nf=NF(opt),tc,nx;
    int ii,ij,flag=0;

    trace(3,"ddres   : dt=%.1f ns=%d\n",dt,ns);
    trace(3,"base position=%.3lf %.3lf %.3lf\n",rtk->rb[0],rtk->rb[1],rtk->rb[2]);

    /* tc=0: common rtk mode
     * tc=1: tightly coupled mode
     * */
    tc=opt->mode==PMODE_INS_TGNSS?1:0;

    /* flag=0: dgps-tightly coupled mode
     * flag=1: rtk-tightly coupled mode
     * */
    flag=opt->mode==PMODE_INS_TGNSS&&insopt->tc==INSTC_RTK;

    tc?nx=rtk->ins.nx:nx=rtk->nx;
    tc?matcpy(rr,re,1,3):matcpy(rr,x,1,3);

    bl=baseline(rr,rtk->rb,dr);
    ecef2pos(rr,posu); ecef2pos(rtk->rb,posr);

    Ri=mat(ns*nf*2+2,1); Rj=mat(ns*nf*2+2,1); im=mat(ns,1);
    tropu=mat(ns,1); tropr=mat(ns,1);
    dtdxu=mat(ns,3); dtdxr=mat(ns,3);
    vc=mat(ns*nf,1);

    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
            rtk->ssat[i].resp[j]=rtk->ssat[i].resc[j]=0.0;
            rtk->ssat[i].news[j]=0;
        }
    /* compute factors of ionospheric and tropospheric delay */
    for (i=0;i<ns;i++) {
        if (opt->ionoopt>=IONOOPT_EST) {
            im[i]=(ionmapf(posu,azel+iu[i]*2)+ionmapf(posr,azel+ir[i]*2))/2.0;
        }
        if (opt->tropopt>=TROPOPT_EST) {
            tropu[i]=prectrop(rtk->sol.time,posu,0,azel+iu[i]*2,opt,x,dtdxu+i*3);
            tropr[i]=prectrop(rtk->sol.time,posr,1,azel+ir[i]*2,opt,x,dtdxr+i*3);
        }
    }
    for (m=0;m<6;m++) { /* m=1:GLO,2:GAL,3:BDS,4:GPS/QZS/SBS,5:IRN */

        /* first frequency only: PDOP, C/N0 and multipath are L1; every carrier slip is checked */
        i=sel_refsat(rtk,nav,obs,sat,iu,ir,azel,y,ns,m,0,nf);

        for (f=(opt->mode>PMODE_DGPS&&opt->mode<PMODE_INS_UPDATE)||flag?0:nf;f<nf*2;f++) {

            if (i<0) continue;
            /* L2/L5 and the code slots use the first-frequency reference */
            if (refsat) {

                if (refsat[m][f]!=sat[i]) {
                    trace(3,"reference satelite change: sat=%3d->%3d  f=%2d\n",refsat[m][f],sat[i],f);
                }
                /* set double difference reference satelite */
                refsat[m][f]=sat[i];
            }
            if (!validobs(iu[i],ir[i],f,nf,y)) continue;
            /* make double difference */
            for (j=0;j<ns;j++) {
                if (i==j) continue;
                sysi=rtk->ssat[sat[i]-1].sys;
                sysj=rtk->ssat[sat[j]-1].sys;
                if (!validobs(iu[j],ir[j],f,nf,y)) continue;
                if (!test_sys(sysj,m)) continue;

                factor=1.0;

                if (f<nf) {
                    /* check all frequency phase-observation data */
                    if (usedualfrq&&!chkfrq(obs+iu[j],opt)) continue;

                    /* Qt does not down-weight a slipped phase. tc still does. */
                    if (tc&&rtk->ssat[sat[j]-1].slip[f]) {
                        factor*=10.0;
                    }
                }
                ff=f%nf;
                lami=nav->lam[sat[i]-1][ff];
                lamj=nav->lam[sat[j]-1][ff];
                if (lami<=0.0||lamj<=0.0) continue;
                if (H) {
                    Hi=H+nv*nx;
                    for (k=0;k<nx;k++) Hi[k]=0.0;
                }
                /* double-differenced residual */
                v[nv]=(y[f+iu[i]*nf*2]-y[f+ir[i]*nf*2])-
                      (y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2]);

                /* partial derivatives by rover position */
                if (H) {
                    if (!tc) for (k=0;k<3;k++) {
                            Hi[k]=-e[k+iu[i]*3]+e[k+iu[j]*3];
                        }
                    else {
                        /* partial derivations by ins position */
                        jacob_dd_dp(&rtk->ins,&e[iu[i]*3],&e[iu[j]*3],dp);
                        Hi[xiP(insopt)+0]=dp[0];
                        Hi[xiP(insopt)+1]=dp[1];
                        Hi[xiP(insopt)+2]=dp[2];
                    }
                }
                /* partial derivation by ins attitude error */
                if (H&&tc) {
                    jacob_dd_da(&rtk->ins,&e[iu[i]*3],&e[iu[j]*3],da);

#if UPD_IN_EULER
                    jacobian_prot_pang(rtk->ins.Cbe,S);
                    matcpy(dap,da,1,3);
                    matmul("NN",1,3,3,1.0,dap,S,0.0,da);
#endif
                    Hi[xiA(insopt)+0]=da[0];
                    Hi[xiA(insopt)+1]=da[1];
                    Hi[xiA(insopt)+2]=da[2];
                }
                /* partial derivation by lever arm */
                if (H&&tc&&xnLa(insopt)) {
                    jacob_dd_dl(&rtk->ins,&e[iu[i]*3],&e[iu[j]*3],dl);
                    Hi[xiLa(insopt)+0]=dl[0];
                    Hi[xiLa(insopt)+1]=dl[1];
                    Hi[xiLa(insopt)+2]=dl[2];
                }
                /* double-differenced ionospheric delay term */
                if (opt->ionoopt==IONOOPT_EST) {
                    ii=tc?xiIo(insopt,sat[i]):II(sat[i],opt);
                    ij=tc?xiIo(insopt,sat[j]):II(sat[j],opt);

                    fi=lami/lam_carr[0]; fj=lamj/lam_carr[0];
                    didxi=(f<nf?-1.0:1.0)*fi*fi*im[i];
                    didxj=(f<nf?-1.0:1.0)*fj*fj*im[j];

                    v[nv]-=didxi*x[ii]-didxj*x[ij];
                    if (H) {
                        Hi[ii]=didxi; Hi[ij]=-didxj;
                    }
                }
                /* double-differenced tropospheric delay term */
                if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                    v[nv]-=(tropu[i]-tropu[j])-(tropr[i]-tropr[j]);
                    for (k=0;k<(opt->tropopt<TROPOPT_ESTG?1:3);k++) {
                        if (!H) continue;
                        ii=tc?xiTr(insopt,0):IT(0,opt);
                        ij=tc?xiTr(insopt,1):IT(1,opt);
                        Hi[ii+k]= (dtdxu[k+i*3]-dtdxu[k+j*3]);
                        Hi[ij+k]=-(dtdxr[k+i*3]-dtdxr[k+j*3]);
                    }
                }
                /* double-differenced phase-bias term */
                if (f<nf) {

                    ii=tc?xiBs(insopt,sat[i],f):IB(sat[i],f,opt);
                    ij=tc?xiBs(insopt,sat[j],f):IB(sat[j],f,opt);

                    if (opt->ionoopt!=IONOOPT_IFLC) {
                        v[nv]-=lami*x[ii]-lamj*x[ij];
                        if (H) Hi[ii]=lami,Hi[ij]=-lamj;
                    }
                    else {
                        v[nv]-=x[ii]-x[ij];
                        if (H) Hi[ii]=1.0,Hi[ij]=-1.0;
                    }
                }
                /* glonass receiver h/w bias term */
                if (rtk->opt.glomodear==2
                    &&sysi==SYS_GLO&&sysj==SYS_GLO
                    &&ff<NFREQGLO) {

                    /* freq-difference (MHz) */
                    df=(CLIGHT/lami-CLIGHT/lamj)/1E6;
                    ii=tc?xiLl(insopt,ff):IL(ff,opt);

                    v[nv]-=df*x[ii];
                    if (H) Hi[ii]=df;
                }
                    /* glonass interchannel bias correction */
                else if (sysi==SYS_GLO&&sysj==SYS_GLO) {
                    v[nv]-=gloicbcorr(sat[i],sat[j],&rtk->opt,lami,lamj,f);
                }
                if (f<nf) rtk->ssat[sat[j]-1].resc[f   ]=v[nv];
                else      rtk->ssat[sat[j]-1].resp[f-nf]=v[nv];

                if (f>=nf&&tc) {
                    ddi=rtk->ssat[sat[i]-1].sdi[f%nf]-rtk->ssat[sat[j]-1].sdi[f%nf];
                    ddg=rtk->ssat[sat[i]-1].sdp[f%nf]-rtk->ssat[sat[j]-1].sdp[f%nf];
                    if (fabs(ddi-ddg)>1.0) {
                        factor*=5.0;
                    }
                }
                /* test innovation */
                if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
                    if (f<nf) {
                        rtk->ssat[sat[i]-1].rejc[f]++;
                        rtk->ssat[sat[j]-1].rejc[f]++;
                    }
                    trace(2,"outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                          sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv]);
                    /* Qt drops the residual. tc keeps it, down-weighted. */
                    if (!tc) continue;
                    factor*=30.0;
                }
                /* single-difference measurement error variances */
                Ri[nv]=varerr(sat[i],sysi,azel[1+iu[i]*2],bl,dt,f,
                             obs[iu[i]].SNR[f%nf]*0.25,opt);
                Rj[nv]=varerr(sat[j],sysj,azel[1+iu[j]*2],bl,dt,f,
                             obs[iu[j]].SNR[f%nf]*0.25,opt);

                /* set valid data flags */
                if (opt->mode>PMODE_DGPS&&opt->mode<PMODE_INS_UPDATE||flag) {
                    if (f<nf) {
                        rtk->ssat[sat[i]-1].vsat [f]=rtk->ssat[sat[j]-1].vsat[f]=1;
                        rtk->ssat[sat[j]-1].index[f]=(unsigned int)nv;
                    }
                }
                else {
                    rtk->ssat[sat[i]-1].vsat[f-nf]=rtk->ssat[sat[j]-1].vsat[f-nf]=1;
                }
                /* code measurement data valid flag */
                if (f>=nf) {
                    rtk->ssat[sat[i]-1].vsatc[f-nf]=rtk->ssat[sat[j]-1].vsatc[f-nf]=1;
                }
                trace(3,"sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n",
                      sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv],Ri[nv],Rj[nv]);
#if UPDNEWSAT
                /* check this phase double difference is new? */
                if (f<nf&&(rtk->refsat[m][f]&&sat[i]==rtk->refsat[m][f])&&chknews(rtk,sat[i],sat[j],f)) {
                    trace(3,"new satellite come in,sat=%2d-%2d\n",sat[i],sat[j]);

                    /* new double difference satellite flag */
                    rtk->ssat[sat[j]-1].news[f]=1;

                    /* disable double difference measurement */
                    Rj[nv]*=10.0;
                }
#endif
                Rj[nv]*=factor;
                vflg[nv++]=(sat[i]<<16)|(sat[j]<<8)|((f<nf?0:1)<<4)|(f%nf);
                nb[b]++;
            }
            b++;
        }
    }
    /* end of system loop */
#if DETECT_OUTLIER
    static const double r0=re_norm(0.95),r1=re_norm(0.99);
    static double s0;

    /* detect outlier by L1/L2 phase double difference residual */
    if (nv&&nf>=2&&opt->mode>PMODE_DGPS) {
        double *vv=mat(ns,1),avv=0.0,v0=0.0;
        int l,*s=imat(ns,1);

        /* iteration for detect outliers */
        for (m=0,flag=1;m<ns;m++,flag=1) {
            for (l=0,i=0;i<ns;i++) {
                if (!rtk->ssat[sat[i]-1].vsat[0]) continue;

                /* exclude reference satellite */
                if (refsat[sysind(rtk->ssat[sat[i]-1].sys)][0]==sat[i]) continue;
                if (refsat[sysind(rtk->ssat[sat[i]-1].sys)][1]==sat[i]) continue;

                /* L1-L2 */
                j=rtk->ssat[sat[i]-1].index[0];
                k=rtk->ssat[sat[i]-1].index[1];

                vv[l]=v[j]-v[k];

                /* index of dd-res */
                s[l++]=i;
            }
            if (l<3) break;

            /* outliers detect */
            for (i=0;i<l;i++) avv+=vv[i]; avv/=l;
            for (i=0;i<l;i++) vv[i]-=avv;

            /* standard deviation */
            matmul("NT",1,1,l,1.0/(l-1),vv,vv,0.0,&v0);

            /* chi-square detect outliers */
            for (i=0;i<l;i++) {
                if (fabs(v0)>=THRES_L1L2RES&&fabs(vv[i])/SQRT(v0)>=r0) {

                    trace(2,"L1/L2 dd-res detect outlier,sat=%2d\n",sat[s[i]]);
                    j=rtk->ssat[sat[s[i]]-1].index[0];
                    k=rtk->ssat[sat[s[i]]-1].index[1];

                    /* disable this satellite */
                    Rj[j]=Rj[k]=SQR(100.0);

                    for (f=0;f<nf;f++) {
                        rtk->ssat[sat[s[i]]-1].vsat[f]=0;
                        rtk->ssat[sat[s[i]]-1].lock[f]++;
                    }
                    /* outlier flag */
                    flag=0;
                }
            }
            if (flag) break;
        }
        free(s); free(vv);
    }
    /* code double-difference measurement outliers detect */
    int sat1,sat2,*ind;
    for (m=0;m<6&&nv;m++) { /* m=1:GLO,2:GAL,3:BDS,4:GPS/QZS/SBS,5:IRN */
        ind=imat(nv,1);
        for (i=0,k=0,s0=0.0;i<nv;i++) {

            if (((vflg[i]>>4)&0xF)==0) continue;

            /* sat no. */
            sat1=(vflg[i]>>16)&0xFF;
            sat2=(vflg[i]>> 8)&0xFF;

            /* test navi system */
            if (!test_sys(rtk->ssat[sat1-1].sys,m)) continue;
            if (!test_sys(rtk->ssat[sat2-1].sys,m)) continue;

            vc[k]=v[i];
            s0+=v[i];
            ind[k++]=i;
        }
        if (k>2) {
            s0/=k;
            for (i=0;i<k;i++) vc[i]-=s0;
            matmul("NT",1,1,k,1.0/(k-1),vc,vc,0.0,&s0);

            /* outliers detect */
            for (i=0;i<k;i++) {
                if (fabs(vc[i])/SQRT(s0)>=r1) {

                    /* disable this satellite */
                    Rj[ind[i]]=SQR(100.0);
                }
                else if (fabs(vc[i])/SQRT(s0)>=r0) {

                    /* degrade this satellite */
                    Rj[ind[i]]=SQR(10.0);
                }
            }
        }
        free(ind);
    }
#endif
    /* baseline length constraint for moving baseline */
    if (opt->mode==PMODE_MOVEB&&constbl(rtk,x,P,v,H,Ri,Rj,nv)) {
        vflg[nv++]=3<<4;
        nb[b++]++;
    }
    /* double-differenced measurement error covariance */
    ddcov(nb,b,Ri,Rj,nv,R);

    /* output used variables if need */
    if (Ri_&&Rj_) {
        matcpy(Ri_,Ri,1,nv);
        matcpy(Rj_,Rj,1,nv);
    }
    if (nb_&&b_) {
        for (*b_=b,i=0;i<b;i++) nb_[i]=nb[i];
    }
    free(Ri); free(Rj); free(im);
    free(tropu); free(tropr);
    free(dtdxu); free(dtdxr);
    free(vc); return nv;
}
/* time-interpolation of residuals (for post-mission) ------------------------*/
static double intpres(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                      rtk_t *rtk, double *y)
{
    static obsd_t obsb[MAXOBS];
    static double yb[MAXOBS*NFREQ*2],rs[MAXOBS*6],dts[MAXOBS*2],var[MAXOBS];
    static double e[MAXOBS*3],azel[MAXOBS*2];
    static int nb=0,svh[MAXOBS*2];
    prcopt_t *opt=&rtk->opt;
    double tt=timediff(time,obs[0].time),ttb,*p,*q;
    int i,j,k,nf=NF(opt);

    trace(3,"intpres : n=%d tt=%.1f\n",n,tt);

    if (nb==0||fabs(tt)<DTTOL) {
        nb=n; for (i=0;i<n;i++) obsb[i]=obs[i];
        return tt;
    }
    ttb=timediff(time,obsb[0].time);
    if (fabs(ttb)>opt->maxtdiff*2.0||ttb==tt) return tt;

    satposs(time,obsb,nb,nav,opt->sateph,rs,dts,var,svh);

    if (!zdres(1,obsb,nb,rs,dts,svh,nav,rtk->rb,opt,1,yb,e,azel)) {
        return tt;
    }
    for (i=0;i<n;i++) {
        for (j=0;j<nb;j++) if (obsb[j].sat==obs[i].sat) break;
        if (j>=nb) continue;
        for (k=0,p=y+i*nf*2,q=yb+j*nf*2;k<nf*2;k++,p++,q++) {
            if (*p==0.0||*q==0.0) *p=0.0;
            else *p=(ttb*(*p)-tt*(*q))/(ttb-tt);
        }
    }
    return fabs(ttb)>fabs(tt)?ttb:tt;
}
/* L1/B1 and L2/B2 ambiguity to WL ambiguity transformation matrix------------*/
static int ddmat_WL(int na,int nb,const ddsat_t *ddsat,double *D,
                    int *index,ddsat_t *wlsat)
{
    int i,j,k,*flag,ny=na+nb;

    trace(3,"ddmat_WL:\n");

    flag=imat(nb,1); for (i=0;i<nb;i++) flag[i]=1;
    for (i=0;i<na;i++) D[i+i*ny]=1.0;

    /* transformation matrix  */
    for (k=0,i=0;i<nb;i++) {

        for (j=0;j<nb;j++) {
            if (flag[j]&&i!=j
                &&ddsat[j].sat1==ddsat[i].sat1
                &&ddsat[j].sat2==ddsat[i].sat2&&ddsat[j].f==1) {

                D[i+(k+na)*ny]= 1.0;
                D[j+(k+na)*ny]=-1.0;
                wlsat[k].sat1=ddsat[j].sat1;
                wlsat[k].sat2=ddsat[j].sat2;
                wlsat[k].f=0;

                index[2*k+0]=i; index[2*k+1]=j;
                flag[i]=flag[j]=0;
                k++;
            }
        }
    }
    trace(3,"WL D=\n");
    tracemat(3,D,ny,k+na,12,5);
    free(flag); return k;
}
/* add double-difference ambiguity-------------------------------------------*/
static int addddamb(amb_t *amb)
{
    ddamb_t *data;
    if (amb->nmax<=amb->nb) {
        if (amb->nmax<=0) amb->nmax=MAXSAT;
        else amb->nmax*=2;

        if (!(data=(ddamb_t *)realloc(amb->amb,sizeof(ddamb_t)*amb->nmax))) {
            trace(2,"add ambiguity: realloc error n=%dx%d\n",sizeof(ddamb_t),amb->nmax);
            free(amb->amb);
            amb->amb=NULL;
            amb->nb=amb->nmax=0;
            return -1;
        }
        amb->amb=data;
    }
    return 1;
}
/* get double-difference ambiguity -------------------------------------------*/
static ddamb_t *getddamb(amb_t *bias,int sat1,int sat2,int f)
{
    int i; for (i=0;i<bias->nb;i++)
        if (sat1==bias->amb[i].sat1&&
            sat2==bias->amb[i].sat2&&f==bias->amb[i].f) return bias->amb+i;
    return NULL;
}
/* extract double-difference ambiguity----------------------------------------*/
static int storeddamb(rtk_t *rtk, const ddsat_t *ddsat, int nb,const double *bias)
{
    ddamb_t *amb=NULL,amb0={0};
    int i,j;

    trace(3,"storeddamb: nb=%d\n",nb);

    /* update double-difference ambiguity list */
    for (j=0,i=0;i<nb;i++) {

        /* no inherit fix */
        if (ddsat[i].flag==1) continue;

        if ((amb=getddamb(&rtk->bias,ddsat[i].sat1,ddsat[i].sat2,ddsat[i].f))==NULL) {
            addddamb(&rtk->bias);

            /* new ambiguity */
            amb=&rtk->bias.amb[rtk->bias.nb++];
            *amb=amb0;
        }
        /* updates ambiguity */
        amb->time =rtk->sol.time;
        amb->ratio=rtk->sol.ratio;

        /* fix counts */
        if (fabs(amb->bias-bias[i])>1E-3) {
            amb->c=0;
        }
        else {
            amb->c++;
        }
        amb->bias=bias[i];

        /* update satellite and frequency no. */
        amb->sat1=ddsat[i].sat1;
        amb->sat2=ddsat[i].sat2;
        amb->f=ddsat[i].f;
    }
    return j;
}
/* single to double-difference transformation matrix (D') --------------------*/
//ddsat:双差卫星的卫星编号，nv卫星总数
static int ddmat(rtk_t *rtk, double *D,ddsat_t *ddsat,const int *vflg,int nv,
                 int flag)
{
    int i,j,k,m,f,tc,sat1,sat2;
    int nb=0,nx=rtk->nx,na=rtk->na;
    static unsigned char fixflag[NFREQ]={0};
    double *x=rtk->x;

    trace(3,"ddmat   :\n");

    tc=rtk->opt.mode==PMODE_INS_TGNSS;

    if (tc) {nx=rtk->ins.nx; na=rtk->ins.nb;}
    if (tc) {
        x=rtk->ins.x;
    }
    //fix清零
    for (i=0;i<MAXSAT;i++) memcpy(rtk->ssat[i].fix,fixflag,sizeof(unsigned char)*NFREQ);
    for (i=0;i<na;i++) D[i+i*nx]=1.0;

    /* for all double difference measurements */
    for (i=0;i<nv;i++) {

        /* check observation type */
        if (((vflg[i]>>4)&0xF)==1) continue;

        sat1=(vflg[i]>>16)&0xFF; /* reference satellite */
        sat2=(vflg[i]>> 8)&0xFF; /* other satellite */
        f=vflg[i]&0xF;           /* freq no. */

        m=sysind(rtk->ssat[sat2-1].sys);

        if ((m==1&&rtk->opt.glomodear==0)||
            (m==3&&rtk->opt.bdsmodear==0)) {
            continue;
        }
        /* reference satellite index */
        tc?k=xiBs(&rtk->opt.insopt,sat1,f):k=IB(sat1,f,&rtk->opt);

        /* other satellite index */
        tc?j=xiBs(&rtk->opt.insopt,sat2,f):j=IB(sat2,f,&rtk->opt);

        /* fix flag */
        rtk->ssat[sat1-1].fix[f]=2;

        /* check double difference ambiguity value */
        if (x[j]==0.0||x[k]==0.0) continue;

        /* slicp cycle */
        if (rtk->ssat[sat2-1].slip[f]) continue;

#if UPDNEWSAT
        /* no fix new double difference ambiguity */
        if (rtk->ssat[sat2-1].news[f]) continue;
#endif
        //卫星valid后并连续锁定opt.minlock次以后双差才生效
        if ((rtk->ssat[sat2-1].lock[f]>=0||flag)&&
            rtk->ssat[sat2-1].vsat[f]&&
            rtk->ssat[sat2-1].azel[1]>=rtk->opt.elmaskar) {

            D[k+(na+nb)*nx]= 1.0;
            D[j+(na+nb)*nx]=-1.0;

            ddsat[nb].f=f;
            ddsat[nb].sat1=sat1;
            ddsat[nb].sat2=sat2;

            rtk->ssat[sat2-1].fix[f]=2; /* fix */
            nb++;
        }
        else rtk->ssat[sat2-1].fix[f]=1; /* no fix */
    }
    trace(5,"D=\n");
    tracemat(5,D,nx,na+nb,2,0);

    trace(3,"fix double difference satelites:\n");
    for (i=0;i<nb;i++) {
        trace(3,"sat=%3d-%3d lock1(%3d-%3d) lock2(%3d-%3d) el(%4.1lf-%4.1lf) slip(%d-%d)\n",
              ddsat[i].sat1,ddsat[i].sat2,
              rtk->ssat[ddsat[i].sat1-1].lock[0],rtk->ssat[ddsat[i].sat1-1].lock[1],
              rtk->ssat[ddsat[i].sat2-1].lock[0],rtk->ssat[ddsat[i].sat2-1].lock[1],
              rtk->ssat[ddsat[i].sat1-1].azel[1]*R2D,
              rtk->ssat[ddsat[i].sat2-1].azel[1]*R2D,
              rtk->ssat[ddsat[i].sat1-1].slip[0],
              rtk->ssat[ddsat[i].sat2-1].slip[0]);
    }
    return nb;
}
/* restore single-differenced ambiguity --------------------------------------*/
static void restamb(rtk_t *rtk, const double *bias, ddsat_t *ddsat,
                    int nb, double *xa)
{
    int i,f;
    int sat1,sat2,ib1,ib2,tc,nx,na;

    trace(3,"restamb :nb=%d\n",nb);

    tc=rtk->opt.mode==PMODE_INS_TGNSS;

    tc?nx=rtk->ins.nx:nx=rtk->nx;
    tc?na=rtk->ins.nb:na=rtk->na;

    /* initial estimated states */
    for (i=0;i<nx;i++) tc?xa[i]=rtk->ins.x [i]:xa[i]=rtk->x [i];
    for (i=0;i<na;i++) tc?xa[i]=rtk->ins.xb[i]:xa[i]=rtk->xa[i];

    for (i=0;i<nb;i++) {

        sat1=ddsat[i].sat1; /* reference satellite */
        sat2=ddsat[i].sat2; /* other satellite */
        f=ddsat[i].f;       /* frq-no. */

        /* state index of ambiguity */
        ib1=tc?xiBs(&rtk->opt.insopt,sat1,f):
            IB(sat1,f,&rtk->opt);
        ib2=tc?xiBs(&rtk->opt.insopt,sat2,f):
            IB(sat2,f,&rtk->opt);

        /* restore single-differenced ambiguity */
        xa[ib2]=xa[ib1]-bias[i];
    }
}
/* restore single-differenced ambiguity --------------------------------------*/
static void restambWL(rtk_t *rtk,const nav_t *nav,const obsd_t *obs,const int *index,int nw,ddsat_t *ddsat,
                      int nb,const double *wlbias,const double *bias,
                      double *xs)
{
    double *y=mat(nb,1);
    double *lc=mat(nw,1);
    int i,sat1,sat2,f,ib1,ib2,tc;

    trace(3,"restambWL: nw=%d\n",nw);

    tc=rtk->opt.mode==PMODE_INS_TGNSS;

    matcpy(y,bias,1,nb);

    //lc组合
    const double *lam=nav->lam[obs->sat-1];
    double f1,f2;
    if (lam[0]==0.0||lam[1]==0.0) return;
    f1=CLIGHT/lam[0];
    f2=CLIGHT/lam[1];

    double C1 = f1*f1/(f1*f1 - f2*f2);
    double C2 = -f1*f2/(f1*f1 - f2*f2);
    double C3 = 1+ f2/f1;
    double C4 = -f2/(f1-f2);
    for (i=0;i<nw;i++) {
        /*lc 无电离层组合 */
        lc[i]=C1*y[index[2*i]]+C2*y[index[2*i+1]];
    }
    for (i=0;i<nw;i++) {
        /* n1浮点解*/
        y[index[2*i]]=C3*lc[i]+C4*wlbias[i];
    }
    for (i=0;i<nw;i++) {
        /* n2浮点解*/
        y[index[2*i+1]]=y[index[2*i]]-wlbias[i];
    }
    /* restore single ambiguity */
    for (i=0;i<nb;i++) {

        sat1=ddsat[i].sat1; /* reference satellite */
        sat2=ddsat[i].sat2; /* other satellite */
        f=ddsat[i].f;       /* frq-no. */

        /* state index of ambiguity */
        ib1=tc?xiBs(&rtk->opt.insopt,sat1,f):IB(sat1,f,&rtk->opt);
        ib2=tc?xiBs(&rtk->opt.insopt,sat2,f):IB(sat2,f,&rtk->opt);

        /* restore single-differenced ambiguity */
        xs[ib2]=xs[ib1]-y[i];
    }
    free(y);
}
/* hold integer ambiguity ----------------------------------------------------*/
static void holdamb(rtk_t *rtk, insstate_t *ins, const double *xa,const ddsat_t *ddsat,
                    int ns)
{
    double *v,*H,*R,*x,*P,*r;
    int i,f,info;
    int nb,nv=0;
    int sat1,sat2,ib1,ib2,tc,nx;
    ddamb_t *pamb=NULL;

    trace(3,"holdamb :ns=%d\n",ns);

    tc=rtk->opt.mode==PMODE_INS_TGNSS;

    tc?nx=rtk->ins.nx:nx=rtk->nx;
    tc?nb=rtk->ins.nx-rtk->ins.nb:nb=rtk->nx-rtk->na;

    tc?x=rtk->ins.x:x=rtk->x;
    tc?P=rtk->ins.P:P=rtk->P;

    v=mat(nb,1); H=zeros(nb,nx);
    r=mat(nb,1);

    for (i=0;i<ns;i++) {
        sat1=ddsat[i].sat1; /* reference satellite */
        sat2=ddsat[i].sat2; /* other satellite */
        f=ddsat[i].f;       /* frq-no. */

        if (ddsat[i].flag) continue;

        /* state index of ambiguity */
        ib1=tc?xiBs(&rtk->opt.insopt,sat1,f):IB(sat1,f,&rtk->opt);
        ib2=tc?xiBs(&rtk->opt.insopt,sat2,f):IB(sat2,f,&rtk->opt);

        /* elevation mask to hold ambiguity (deg) */
        if (rtk->ssat[sat2-1].azel[1]<rtk->opt.elmaskhold) continue;

        rtk->ssat[sat1-1].fix[f]=3; /* hold */
        rtk->ssat[sat2-1].fix[f]=3; /* hold */

        /* constraint to fixed ambiguity */
        v[nv]=(xa[ib1]-xa[ib2])-(x[ib1]-x[ib2]);
#if 0
        /* get double difference ambiguity of precious epoch */
        pamb=getddamb(&rtk->bias,sat1,sat2,f);
        if (!pamb) {
            /* new dd-ambiguity */
            r[nv]=10.0*VAR_HOLDAMB;
        }
        else {
            /* already hold */
            if (fabs(timediff(pamb->time,rtk->sol.time))<1.5&&fabs(pamb->bias-(xa[ib1]-xa[ib2]))<0.5) {
                r[nv]=VAR_HOLDAMB;
            }
            else r[nv]=3.0*VAR_HOLDAMB;
        }
#else
        r[nv]=VAR_HOLDAMB;
#endif
        /* Qt holdamb constrains every DD. The 3-cycle skip is kept
           on the tightly coupled path only. */
        if (tc&&fabs(v[nv])>THRES_HOLDAMB) continue;

        H[ib1+nv*nx]= 1.0;
        H[ib2+nv*nx]=-1.0;
        nv++;
    }
    if (nv>0) {
        double *x_b,*P_b;
        R=zeros(nv,nv);
        x_b=mat(nx,1); P_b=mat(nx,nx);
        matcpy(x_b,x,nx,1); matcpy(P_b,P,nx,nx);
        for (i=0;i<nv;i++) {
            R[i+i*nv]=(rtk->bias.inherit?3.0*r[i]:r[i]);
        }
        /* close-loop states set to zero */
        if (tc) {
            for (i=0;i<xnCl(&rtk->opt.insopt);i++) x[i]=1E-20;
        }
        /* update states with constraints */
        if ((info=filter(x,P,H,v,R,nx,nv))) {
            errmsg(rtk,"filter error (info=%d)\n",info);
        }
#if CHKNUMERIC
        if (!chk_xp(x,P,nx)) {
            fprintf(stderr,"check numeric error: nan or inf\n");
            errmsg(rtk,"check numeric error: nan or inf\n");
            matcpy(x,x_b,nx,1); matcpy(P,P_b,nx,nx);
            info=-1;
        }
#endif
        if (!info&&tc) {
            /* close loop for ins states */
            clp(ins,&rtk->opt.insopt,x);
        }
        free(R); free(x_b); free(P_b);
    }
    free(v); free(H);
    free(r);
}
/* store WL ambiguity---------------------------------------------------------*/
static int storeambwl(rtk_t *rtk, const ddsat_t *ddsat, int nw,const double *bias)
{
    ddamb_t *amb=NULL,amb0={0};
    int i,j;

    trace(3,"storeambwl: nb=%d\n",nw);

    /* update double-difference ambiguity list */
    for (j=0,i=0;i<nw;i++) {

        /* no inherit */
        if (ddsat[i].flag==1) continue;

        if ((amb=getddamb(&rtk->wlbias,ddsat[i].sat1,ddsat[i].sat2,ddsat[i].f))==NULL) {
            addddamb(&rtk->wlbias);

            /* new ambiguity */
            amb=&rtk->wlbias.amb[rtk->wlbias.nb++];
            *amb=amb0;
        }
        /* updates ambiguity */
        amb->time =rtk->sol.time;
        amb->ratio=rtk->sol.wlratio;

        /* fix counts */
        if (fabs(amb->bias-bias[i])>1E-3) {
            amb->c=0;
        }
        else {
            amb->c++;
        }
        amb->bias=bias[i];

        /* update satellite and frequency no. */
        amb->sat1=ddsat[i].sat1;
        amb->sat2=ddsat[i].sat2;
        amb->f=ddsat[i].f;
    }
    return j;
}
/* inherit WL ambiguity继承宽巷模糊度------------------------------------------------------*/
static int inheritambwl(rtk_t *rtk, ddsat_t *wlsat,const double *wl,
                        int nw,double *b,double *s)
{
    ddamb_t *pamb=NULL;
    int i,k=0;

    trace(3,"inheritambwl:\n");
#if INHERIT_AMB
    for (i=0;i<nw;i++) b[i]=wl[i];
    for (i=0;i<nw;i++) {
        wlsat[i].flag=1;
        if ((pamb=getddamb(&rtk->wlbias,wlsat[i].sat1,wlsat[i].sat2,wlsat[i].f))==NULL) continue;
        if (pamb->ratio<rtk->opt.thresar[0]*1.5) continue;

        if (timediff(rtk->sol.time,pamb->time)>THRES_INHERIT_TIME) continue;
        if (fabs(wl[i]-pamb->bias)>THRES_INHERIT_BIAS) continue;
        if (pamb->c<MINFIXC) continue;
        wlsat[i].flag=0;
        b[i]=pamb->bias;
        k++;
    }
    s[1]=s[0]=1.0;
    if (k>=3) {
        s[1]=33.3;
        s[0]=1.00;

        trace(3,"inherit ambiguity=\n");
        tracemat(3,b,1,nw,12,6);
        return 1;
    }
#endif
    return 0;
}
/* resolve WL integer ambiguity by LAMBDA------------------------------------*/
static int resamb_WL(rtk_t *rtk,const nav_t *nav,const obsd_t *obs, double *Qy, double *y, int ny, int *index,const double *D,
                     int *nw, ddsat_t *ddsat, ddsat_t *wlsat,
                     double *dxwl,double *Qywl,double *wlb,int *nny)
{
    int i,j,k,info=0,na,tc,inherit=0;
    double *Qw,*wl,*b,*v,*R,*r,*H,s[2];

    trace(3,"resamb_WL:\n");

    rtk->sol.wlratio=0.0;

    if (*nw<=2) {
        errmsg(rtk,"no valid WL double-difference\n");
        return 0;
    }
    tc=rtk->opt.mode==PMODE_INS_TGNSS;
    tc?na=rtk->ins.nb:na=rtk->na;

    Qw=mat(ny,ny); wl=mat(ny,1); H=zeros(ny,ny); v=mat(ny,1);
    b=mat(ny,2); r=mat(ny,1);
    R=zeros(ny,ny);

    /* covariance matrix of WL ambiguity */
    matmul33("TNN",D,Qy,D,*nw,ny,ny,*nw,Qw);

    trace(3,"Qw=\n");
    tracemat(3,Qw,*nw,*nw,12,6);

    for (i=0;i<*nw;i++) {

        /* WL float ambiguity: N1-N2 */
        wl[na+i]=y[na+index[2*i]]-y[na+index[2*i+1]];
    }
    /* lambda/mlambda integer least-square estimation */
    if (!lambda(*nw,2,wl+na,Qw,b,s)) {

        trace(3,"WL-N(0)=\n"); tracemat(3,wl+na,1,*nw,10,3);
        trace(3,"WL-N(1)=\n"); tracemat(3,b    ,1,*nw,10,3);
        trace(3,"WL-N(2)=\n"); tracemat(3,b+*nw,1,*nw,10,3);

        rtk->sol.wlratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.wlratio>999.9f) {
            rtk->sol.wlratio=999.9f;
        }
        /* validation by popular ratio-test */
        if ((s[0]>0.0&&s[1]/s[0]>=rtk->opt.thresar[0])||(inherit=inheritambwl(rtk,wlsat,wl+na,*nw,b,s))) {

            if (inherit) {
                rtk->sol.wlratio=(float)(s[1]/s[0]);
            }
            /* WL constraint for N1 and N2 ambiguity */
            for (k=0,i=0;i<*nw;i++) {
                v[k]=b[i]-wl[na+i];

                if (fabs(v[k])>THRES_AMB) continue;
                H[index[2*i+0]+na+ny*k]= 1.0;
                H[index[2*i+1]+na+ny*k]=-1.0;
                r[k]=(inherit?SQR(VAR_WLCONST*3.0):SQR(VAR_WLCONST));
                k++;
            }
            if (k&&ny) {
                for (i=0;i<k;i++) R[i+i*k]=r[i];

                /* filter for constraint */
                if (filter(y,Qy,H,v,R,ny,k)) {
                    trace(2,"filter error\n");
                    info=0;
                }
                else {
                    /* output solution */
                    matcpy(dxwl,y , 1,na);
                    matcpy(Qywl,Qy,ny,ny);

                    if (nny) *nny=ny;

                    /* restore single ambiguity 利用宽巷固定解修正单差浮点模糊度*/
                    restambWL(rtk,nav,obs,index,*nw,ddsat,ny-na,b,y+na,dxwl);

                    /* output WL bias */
                    if (wlb) {
                        matcpy(wlb,b,1,*nw);
                    }
                    info=1; /* fix ok */
                }
            }
        }
        else {
            info=0;
        }
    }
    if (info) {
        trace(3,"WL ambiguity fix ok\n");
    }
    else {
        trace(2,"WL ambiguity fix fail\n");
        *nw=0;
    }
    free(Qw); free(wl); free(b);
    free(r ); free(R ); free(H);
    free(v ); 
    return info;
}
/* inherit ambiguity模糊度继承---------------------------------------------------------*/
static int inheritamb(rtk_t *rtk, ddsat_t *ddsat,const double *y,const double *Qy,
                      int na,int nb,double *b,double *s)
{
    ddamb_t *pamb=NULL;
    int i,k=0;

    trace(3,"inheritamb:\n");
#if INHERIT_AMB
    for (i=0;i<nb;i++) b[i]=y[i];
    for (i=0;i<nb;i++) {
        ddsat[i].flag=1;
        if ((pamb=getddamb(&rtk->bias,ddsat[i].sat1,ddsat[i].sat2,ddsat[i].f))==NULL) continue;

        if (timediff(rtk->sol.time,pamb->time)>THRES_INHERIT_TIME) continue;
        if (fabs(y[i]-pamb->bias)>THRES_INHERIT_BIAS) continue;
        if (pamb->c<MINFIXC) continue;
        ddsat[i].flag=0;
        b[i]=pamb->bias;
        k++;
    }
    s[1]=s[0]=1.0;
    if (k>MINFIXC) {
        s[1]=99.9;
        s[0]=1.00;

        trace(3,"inherit ambiguity=\n");
        tracemat(3,b,1,nb,12,6);
        return 1;
    }
#endif
    return 0;
}
/* resolve integer ambiguity by LAMBDA --------------------------------------*/
static int resamb_LAMBDA(rtk_t *rtk,const nav_t *nav,const obsd_t *obs, double *bias, double *xa, ddsat_t *ddsat,int *namb,
                         ddsat_t *wlsat, int *nw,double *dxwl,double *Qywl,double *wlb,
                         const int *vflg,int nv,int fixnolock,
                         int *nny,int *fixwl)
{
    prcopt_t *opt=&rtk->opt;
    insstate_t *ins=&rtk->ins;
    insopt_t *insopt=&rtk->opt.insopt;
    int i,j,ny,nb,info,nx=rtk->nx,na=rtk->na,*index,tc,flag=0,inherit=0;
    double *D,*DP,*y,*Qy,*b,*db,*Qb,*Qab,*QQ,s[2];
    double *DD,*x,*P,*Pa,*xb;

    trace(3,"resamb_LAMBDA :\n");

    if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||//伪距差分和无模糊度固定模式下不解算整周模糊度
        rtk->opt.thresar[0]<1.0) {//ratio值永远>1,所以ratio阈值设置<1没有任何意义
        return 0;
    }
    if (rtk->opt.mode==PMODE_INS_TGNSS&&rtk->opt.insopt.tc<=INSTC_DGPS) {//紧组合且不等于RTK不需要求整周模糊度
        return 0;
    }
    rtk->sol.ratio=0.0;

    /* reset number of dd-ambiguity */
    *namb=0;
    rtk->bias.inherit=0;

    /* tc-mode flag */
    tc=rtk->opt.mode==PMODE_INS_TGNSS;

    tc?x=rtk->ins.x:x=rtk->x;
    tc?P=rtk->ins.P:P=rtk->P;
    tc?Pa=rtk->ins.Pb:Pa=rtk->Pa;
    tc?xb=rtk->ins.xb:xb=rtk->xa;

    /* tc-mode */
    if (opt->mode==PMODE_INS_TGNSS) {
        nx=ins->nx;
        na=ins->nb;
    }
    /* single to double-difference transformation matrix (D') */
    D=zeros(nx,nx);
    if ((nb=ddmat(rtk,D,ddsat,vflg,nv,fixnolock))<=0) {

        trace(2,"no valid double-difference\n");
        free(D);
        return 0;
    }
    //na=9即位置 速度 加速度
    ny=na+nb; y=mat(ny,1); Qy=mat(ny,ny); DP=mat(ny,nx);
    b=mat(nb,2); db=mat(nb,1);
    Qb=mat(nb,nb); Qab=mat(na,nb); QQ=mat(na,nb);
    DD=zeros(ny,ny); index=imat(nb,2);

    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    matmul("TN",ny, 1,nx,1.0,D ,x,0.0,y );
    matmul("TN",ny,nx,nx,1.0,D ,P,0.0,DP);
    matmul("NN",ny,ny,nx,1.0,DP,D,0.0,Qy);

    /* WL ambiguity constraint */
    if (opt->modear==ARMODE_WLNL) {

        /* L1/L2 ambiguity to WL-ambiguity transform matrix */
        if (fixwl&&(*nw=ddmat_WL(na,nb,ddsat,DD,index,wlsat))) {

            /* fix WL ambiguity */
            *fixwl=resamb_WL(rtk,nav,obs,Qy,y,ny,index,DD,nw,ddsat,wlsat,dxwl,Qywl,wlb,nny);
        }
    }
    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb [i+j*nb]=Qy[na+i+(na+j)*ny];
    for (i=0;i<na;i++) for (j=0;j<nb;j++) Qab[i+j*na]=Qy[   i+(na+j)*ny];

    trace(4,"N(0)=");
    tracemat(4,y+na,1,nb,10,3);

    /* lambda/mlambda integer least-square estimation */
    if (!(info=lambda(nb,2,y+na,Qb,b,s))) {

        trace(4,"N(1)="); tracemat(4,b   ,1,nb,10,3);
        trace(4,"N(2)="); tracemat(4,b+nb,1,nb,10,3);

        rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;

        /* validation by popular ratio-test */
        if ((s[0]>0.0&&s[1]/s[0]>=opt->thresar[0])||(inherit=inheritamb(rtk,ddsat,y+na,Qy,na,nb,b,s))) {

            if (inherit) {
                rtk->sol.ratio=(float)(s[1]/s[0]);
                rtk->bias.inherit=1;
            }
            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i=0;i<na;i++) {
                xb[i]=x[i];
                for (j=0;j<na;j++) Pa[i+j*na]=P[i+j*nx];
            }
            /* close loop estimated states */
            if (tc) {
                for (i=0;i<xnCl(insopt);i++) xb[i]=1E-20;
            }
            for (i=0;i<nb;i++) {
                /* fix ambiguity */
                bias[i]=b[i];

                if (fabs(y[na+i]-b[i])>THRES_HOLDAMB) {
                    y[na+i]=1E-5;
                }
                else {
                    /* ambiguity residuals */
                    y[na+i]-=b[i];
                }
            }
            if (!matinv(Qb,nb)) {

                matmul("NN",nb,1,nb, 1.0,Qb ,y+na,0.0,db);
                matmul("NN",na,1,nb,-1.0,Qab,db  ,1.0,xb);

                /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                matmul("NN",na,nb,nb, 1.0,Qab,Qb ,0.0,QQ);
                matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,Pa);

                trace(3,"validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                      nb,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);

                /* restore single-difference ambiguity */
                restamb(rtk,bias,ddsat,nb,xa);
            }
            else nb=0;
        }
        else {
            /* validation failed */
            trace(2,"ambiguity validation fail (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                  nb,s[1]/s[0],s[0],s[1]);
            nb=0;
        }
    }
    else {
        if (inheritamb(rtk,ddsat,y+na,Qy,na,nb,b,s)) {
            rtk->sol.ratio=(float)(s[1]/s[0]);
            rtk->bias.inherit=1;

            rtk->sol.stat=SOLQ_INHERIT;

            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i=0;i<na;i++) {
                xb[i]=x[i];
                for (j=0;j<na;j++) Pa[i+j*na]=P[i+j*nx];
            }
            /* close loop estimated states */
            if (tc) {
                for (i=0;i<xnCl(insopt);i++) xb[i]=0.0;
            }
            for (i=0;i<nb;i++) {
                /* fix ambiguity */
                bias[i]=b[i];

                /* ambiguity residuals */
                y[na+i]-=b[i];
            }
            if (!matinv(Qb,nb)) {

                matmul("NN",nb,1,nb, 1.0,Qb ,y+na,0.0,db);
                matmul("NN",na,1,nb,-1.0,Qab,db  ,1.0,xb);

                /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                matmul("NN",na,nb,nb, 1.0,Qab,Qb ,0.0,QQ);
                matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,Pa);

                trace(3,"validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                      nb,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);

                /* restore single-difference ambiguity */
                restamb(rtk,bias,ddsat,nb,xa);
            }
            else nb=0;
        }
        else {
            trace(3,"lambda error (info=%d)\n",info);
            nb=0;
        }
    }
    free(D ); free(y ); if (flag) free(P);
    free(Qy); free(DP); free(b);
    free(db); free(Qb); free(Qab);
    free(QQ); free(DD); free(index);

    /* number of ambiguities */
    *namb=nb; return nb;
}
static int resamb_LAMBDA2(rtk_t *rtk,const nav_t *nav,const obsd_t *obs,double *bias, double *xa,
                          ddsat_t *ddsat,int *namb,
                          ddsat_t *wlsat,int *nw, double *dxwl,double *Qywl, double *wlb,
                          const int *vflg,
                          int nv,int *nny,int *fixwl)
{
    int nb;

    trace(3,"resamb_LAMBDA2:\n");
#if 1
    if ((nb=resamb_LAMBDA(rtk,nav,obs,bias,xa,ddsat,namb,wlsat,nw,dxwl,Qywl,wlb,vflg,nv,1,nny,fixwl))<MINFIXC) {//锁定宽巷模糊度
        nb=resamb_LAMBDA(rtk,nav,obs,bias,xa,ddsat,namb,wlsat,nw,dxwl,Qywl,wlb,vflg,nv,0,NULL,NULL);//不锁定宽巷模糊度
    }
    return nb;
#else
    return resamb_LAMBDA(rtk,nav,bias,xa,ddsat,namb,wlsat,nw,dxwl,Qywl,vflg,nv,1,nny,1);
#endif
}
/* collect non-ref sats with clear large DD residuals (|v|/sig > thres) ------*/
static int collect_outlier_sats(const double *v, const double *R, const int *vflg,
                                int nv, double thres, int *osat, double *oscore,
                                int maxn)
{
    double score[MAXSAT+1]={0},sig,s;
    int i,j,sat2,n=0,idx[MAXSAT],tmp;

    for (i=0;i<nv;i++) {
        sat2=(vflg[i]>>8)&0xFF;
        if (sat2<=0||sat2>MAXSAT) continue;
        sig=SQRT(R[i+i*nv]);
        if (sig<=0.0) continue;
        s=fabs(v[i])/sig;
        if (s<=thres) continue;          /* only clear outliers */
        if (s>score[sat2]) score[sat2]=s;
    }
    for (i=1;i<=MAXSAT;i++) {
        if (score[i]<=thres) continue;
        idx[n++]=i;
    }
    /* sort outlier sats by score descending (simple insert) */
    for (i=1;i<n;i++) {
        tmp=idx[i];
        for (j=i;j>0&&score[idx[j-1]]<score[tmp];j--) idx[j]=idx[j-1];
        idx[j]=tmp;
    }
    if (n>maxn) n=maxn;
    for (i=0;i<n;i++) {
        osat[i]=idx[i];
        oscore[i]=score[idx[i]];
    }
    return n;
}
/* validation of solution ----------------------------------------------------*/
static int valpos(rtk_t *rtk, const double *v, const double *R, const int *vflg,
                  int nv, int *m,double thres)
{
    prcopt_t *opt=&rtk->opt;
    double vv=0.0;
    double fact=thres*thres;
    int i,stat=1,sat1,sat2,type,freq,np,n,nchi=0,flag=0;
    const char *stype;

    trace(3,"valpos  : nv=%d thres=%.1f\n",nv,thres);

    np=rtk->opt.mode==PMODE_INS_TGNSS?xnP(&opt->insopt):NP(opt);
    flag=(opt->mode==PMODE_INS_TGNSS&&opt->insopt.tc==INSTC_RTK)||opt->mode==PMODE_KINEMA;

    /* post-fit residual test */
    for (*m=0,n=0,i=0;i<nv;i++) {
        type=(vflg[i]>>4)&0xF;
        freq=vflg[i]&0xF;
        if (type!=0&&freq==0) *m=*m+1;
        if (v[i]*v[i]<=fact*R[i+i*nv]) continue;
        sat1=(vflg[i]>>16)&0xFF;
        sat2=(vflg[i]>> 8)&0xFF;
        n++;
        stype=type==0?"L":"C";
        trace(2,"large residual (sat=%2d-%2d %s%d v=%6.3f sig=%.3f)\n",
              sat1,sat2,stype,freq+1,v[i],SQRT(R[i+i*nv]));
    }
    /* rtk position valid solution */
    if (flag&&stat&&nv>np) {

        /* chi-square validation (exclude cycle-slip sats from vv) */
        for (i=0;i<nv;i++) {
            sat1=(vflg[i]>>16)&0xFF;
            sat2=(vflg[i]>> 8)&0xFF;
            freq=vflg[i]&0xF;
            if (!rtk->ssat[sat1-1].vsat[0]||!rtk->ssat[sat2-1].vsat[0]) continue;
            if ((rtk->ssat[sat1-1].slip[freq]&0x03)||(rtk->ssat[sat2-1].slip[freq]&0x03)) continue;
            vv+=v[i]*v[i];
            nchi++;
        }
        if (nchi<=np) {
            trace(3,"validation skip chi-square (nchi=%d np=%d, slip-excluded)\n",
                  nchi,np);
            rtk->failc=0;
        }
        else if (vv>chisqr[nchi-np]*thres) {
            trace(2,"residuals validation fail (nv=%d nchi=%d np=%d vv=%.2f cs=%.2f)\n",
                  nv,nchi,np,vv,chisqr[nchi-np]);
            stat=0;
            rtk->failc++;
        }
        else {
            trace(3,"validation ok (%s nv=%d nchi=%d np=%d vv=%.2f cs=%.2f)\n",
                  time_str(rtk->sol.time,2),nv,nchi,np,vv,chisqr[nchi-np]);
            rtk->failc=0;
        }
    }else{
        stat=0;
    }
    return stat;
}
/* validation of ins solutions------------------------------------------------*/
static int valins(const prcopt_t *opt,const double *x,const insstate_t *ins)
{
    const insopt_t *insopt=&opt->insopt;
    int nba=0,iba=0,nbg=0,ibg=0;

    trace(3,"valins:\n");

    nba=xnBa(insopt); iba=xiBa(insopt);
    nbg=xnBg(insopt); ibg=xiBg(insopt);

    /* check estimated states */
    if (norm(x,3)>5.0*D2R||(nba?norm(x+iba,3)>1E4*Mg2M:false)
        ||(nbg?norm(x+ibg,3)>5.0*D2R:false)) {
        trace(2,"too large estimated state error\n");
        return 0;
    }
    return 1;
}
/* IGG3 variance inflation (Yang 2006): lambda on R, not on predicted P -----*/
static double igg3_lambda(double vt, double k0, double k1)
{
    double a,d;
    vt=fabs(vt);
    if (vt<=k0) return 1.0;
    if (k1<=k0||vt>=k1) return 1e6;
    d=k1-vt;
    if (d<1e-6) return 1e6;
    a=(k1-k0)/d;
    return (vt/k0)*a*a;
}
/* diagonal of Qvv=H'*P*H+R (IGG3 only uses Qii). H columns are sparse. ------*/
static void innov_qvv(const double *H, const double *P, const double *R,
                      int n, int m, double *Q)
{
    int i,j,a,b,nnz,*ix;
    const double *hi;
    double s,ha;

    if (n<=0||m<=0) return;
    ix=imat(n,1);
    for (i=0;i<m;i++) {
        hi=H+i*n;
        nnz=0;
        for (j=0;j<n;j++) if (hi[j]!=0.0) ix[nnz++]=j;
        s=R[i+i*m];
        for (a=0;a<nnz;a++) {
            ha=hi[ix[a]];
            for (b=0;b<nnz;b++) s+=ha*P[ix[a]+ix[b]*n]*hi[ix[b]];
        }
        Q[i+i*m]=s;
    }
    free(ix);
}
static int rkf_satix(const int *sat, int ns, int s)
{
    int i;
    for (i=0;i<ns;i++) if (sat[i]==s) return i;
    return -1;
}
static int rkf_reset_amb(rtk_t *rtk, const obsd_t *obs, const nav_t *nav,
                         const int *sat, const int *iu, const int *ir, int ns,
                         int satno, int f, double *xp)
{
    prcopt_t *opt=&rtk->opt;
    int tc,ib,k;
    double cp,pr,lam,bias,cp1,cp2,pr1,pr2,lam1,lam2,C1,C2;

    k=rkf_satix(sat,ns,satno);
    if (k<0) return 0;
    tc=opt->mode==PMODE_INS_TGNSS?1:0;
    ib=tc?xiBs(&opt->insopt,satno,f):IB(satno,f,opt);

    if (opt->ionoopt==IONOOPT_IFLC) {
        cp1=sdobs(obs,iu[k],ir[k],0);
        cp2=sdobs(obs,iu[k],ir[k],1);
        pr1=sdobs(obs,iu[k],ir[k],NFREQ);
        pr2=sdobs(obs,iu[k],ir[k],NFREQ+1);
        lam1=nav->lam[satno-1][0];
        lam2=nav->lam[satno-1][1];
        if (cp1==0.0||cp2==0.0||pr1==0.0||pr2==0.0||lam1<=0.0||lam2<=0.0) return 0;
        C1= SQR(lam2)/(SQR(lam2)-SQR(lam1));
        C2=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
        bias=(C1*lam1*cp1+C2*lam2*cp2)-(C1*pr1+C2*pr2);
    }
    else {
        cp=sdobs(obs,iu[k],ir[k],f);
        pr=sdobs(obs,iu[k],ir[k],f+NFREQ);
        lam=nav->lam[satno-1][f];
        if (cp==0.0||pr==0.0||lam<=0.0) return 0;
        bias=cp-pr/lam;
    }
    if (tc) insinitx(&rtk->ins,bias,SQR(opt->std[0]),ib);
    else    initx(rtk,bias,SQR(opt->std[0]),ib);
    xp[ib]=bias;
    rtk->ssat[satno-1].lock[f]=-opt->minlock;
    trace(2,"rkf reset amb sat=%3d L%d bias=%.3f\n",satno,f+1,bias);
    return 1;
}
/* inflate Rj by IGG3; phase |vt|>k1 re-inits ambiguity. returns nreset ------*/
static int rkf_reweight(rtk_t *rtk, const obsd_t *obs, const nav_t *nav,
                        const int *sat, const int *iu, const int *ir, int ns,
                        double *xp, const double *Ppred, const double *H,
                        const double *v, double *R, double *Ri, double *Rj,
                        const int *nb, int b, const int *vflg, int nx, int nv,
                        int do_reset)
{
    prcopt_t *opt=&rtk->opt;
    double k0=opt->rkf_k0,k1=opt->rkf_k1,vt,lam,vtmax,*Q;
    int i,nreset=0,sat2,type,freq,nph,nph_k1,imax,mass;

    if (nv<=0||k0<=0.0) return 0;
    if (k1<=k0) k1=k0+4.0;
    Q=mat(nv,nv);
    innov_qvv(H,Ppred,R,nx,nv,Q);

    /* count phase outliers: many large vt => INS/model error, not cycle slips */
    nph=nph_k1=0; imax=-1; vtmax=0.0;
    for (i=0;i<nv;i++) {
        type=(vflg[i]>>4)&0xF;
        if (type!=0||Q[i+i*nv]<=0.0) continue;
        nph++;
        vt=fabs(v[i])/SQRT(Q[i+i*nv]);
        if (vt>k1) nph_k1++;
        if (vt>vtmax) { vtmax=vt; imax=i; }
    }
    mass=(nph_k1>=3)||(nph>0&&nph_k1*2>=nph);

    for (i=0;i<nv;i++) {
        type=(vflg[i]>>4)&0xF;
        freq=vflg[i]&0xF;
        sat2=(vflg[i]>>8)&0xFF;
        if (type>1||Q[i+i*nv]<=0.0) continue;
        vt=fabs(v[i])/SQRT(Q[i+i*nv]);
        lam=igg3_lambda(vt,k0,k1);
        if (type==0) {
            /* phase: reset at most the worst sat, and only if not a mass event */
            if (do_reset&&!mass&&i==imax&&vt>k1) {
                nreset+=rkf_reset_amb(rtk,obs,nav,sat,iu,ir,ns,sat2,freq,xp);
            }
            else if (!mass&&vt>k0&&lam<1e6) {
                Rj[i]*=lam;
                trace(3,"rkf phase inflate sat=%3d L%d vt=%.2f lam=%.3f\n",
                      sat2,freq+1,vt,lam);
            }
        }
        else if (lam!=1.0) { /* code: inflate / reject */
            Rj[i]*=lam;
            trace(3,"rkf code  inflate sat=%3d P%d vt=%.2f lam=%.3g\n",
                  sat2,freq+1,vt,lam);
        }
    }
    ddcov(nb,b,Ri,Rj,nv,R);
    free(Q);
    return nreset;
}
/* iterate filter from the same predicted state; drop worst posterior residual */
static int rkf_iter_filter(double *xp, double *Pp, const double *Ppred,
                           const double *H, const double *v, double *R,
                           double *Ri, double *Rj, const int *nb, int b,
                           const int *vflg, int nx, int nv, int niter, double k0)
{
    double *xp0,*Q,*Hv,*dx,vt,vmax;
    int i,it,info=0,imax,type,allok,nexc=0;

    if (niter<1) niter=1;
    xp0=mat(nx,1); Q=mat(nv,nv); Hv=mat(nv,1); dx=mat(nx,1);
    matcpy(xp0,xp,nx,1);

    for (it=0;it<niter;it++) {
        matcpy(xp,xp0,nx,1);
        matcpy(Pp,Ppred,nx,nx);
        if ((info=filter(xp,Pp,H,v,R,nx,nv))) break;

        for (i=0;i<nx;i++) dx[i]=xp[i]-xp0[i];
        matmul("TN",nv,1,nx,1.0,H,dx,0.0,Hv);
        innov_qvv(H,Pp,R,nx,nv,Q);

        imax=-1; vmax=0.0; allok=1;
        for (i=0;i<nv;i++) {
            type=(vflg[i]>>4)&0xF;
            if (type!=1||Q[i+i*nv]<=0.0) continue; /* iterate-exclude code only */
            vt=fabs(v[i]-Hv[i])/SQRT(Q[i+i*nv]);
            if (vt>k0) {
                allok=0;
                if (vt>vmax) { vmax=vt; imax=i; }
            }
        }
        if (allok||imax<0) break;
        if (nexc>=nv-4) break; /* keep a minimum set */
        Rj[imax]*=1e6;
        ddcov(nb,b,Ri,Rj,nv,R);
        nexc++;
        trace(2,"rkf iter exclude i=%d type=%d vt=%.2f nexc=%d\n",
              imax,(vflg[imax]>>4)&0xF,vmax,nexc);
    }
    free(xp0); free(Q); free(Hv); free(dx);
    return info;
}
/* measurement update: IGG3 R + phase-amb reset + residual iteration ---------*/
static int filter_rkf(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, double dt,
                      double *xp, double *Pp, const double *Ppred,
                      const int *sat, double *y, double *e, double *azel,
                      const int *iu, const int *ir, int ns,
                      double *v, double *H, double *R, double *Ri, double *Rj,
                      int *nb, int *b, int *vflg, const double *rr,
                      int refsat[NUMSYS][2*NFREQ], int usedualfrq, int nx, int *pnv)
{
    int nv=*pnv,nreset,niter;
    double k0;

    if (!rtk->opt.rkf) {
        matcpy(Pp,Ppred,nx,nx);
        return filter(xp,Pp,H,v,R,nx,nv);
    }
    nreset=rkf_reweight(rtk,obs,nav,sat,iu,ir,ns,xp,Ppred,H,v,R,Ri,Rj,nb,*b,
                        vflg,nx,nv,1);
    if (nreset) {
        nv=ddres(rtk,nav,obs,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,vflg,rr,
                 Ri,Rj,nb,b,refsat,usedualfrq);
        if (nv<1) { *pnv=nv; return -1; }
        rkf_reweight(rtk,obs,nav,sat,iu,ir,ns,xp,Ppred,H,v,R,Ri,Rj,nb,*b,
                     vflg,nx,nv,0);
        *pnv=nv;
    }
    niter=rtk->opt.rkf_niter;
    k0=rtk->opt.rkf_k0;
    if (niter<=0) niter=1;
    if (k0<=0.0) k0=2.5;
    return rkf_iter_filter(xp,Pp,Ppred,H,v,R,Ri,Rj,nb,*b,vflg,nx,*pnv,niter,k0);
}
/* relative positioning ------------------------------------------------------*/
static int relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,const nav_t *nav,
                  int degrade,int usedualfrq)
{
    prcopt_t *opt=&rtk->opt;
    insopt_t *insopt=&rtk->opt.insopt;
    insstate_t *ins=&rtk->ins;
    gtime_t time=obs[0].time;
    ddsat_t ddsat[MAXSAT]={{0}};
    ddsat_t wlsat[MAXSAT]={{0}};
    static insstate_t insp={0},instmp={0};
    static int refsat[NUMSYS][2*NFREQ]={0};
    double *Ri,*Rj,dr[3]={0};
    double *rs,*dts,*var,*y,*e,*azel;
    double *v,*H,*R,*xp,*Pp,*xa,*bias,dt,*x,*P,rr[3],*Pa,*dx;
    double *dxwl,*Qywl,*wlb;
    int i,j,k,f,n=nu+nr,ns,ny,nv=0,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT],niter,nx,na;
    int info,vflg[MAXOBS*NFREQ*2+1],svh[MAXOBS*2];
    int stat=rtk->opt.mode<=PMODE_DGPS?SOLQ_DGPS:SOLQ_FLOAT;
    int nf=opt->ionoopt==IONOOPT_IFLC?1:opt->nf,tc;
    int ix,iy,iz,ivx,ivy,ivz,namb=0,nny=0,nw=0,fixwl=0;
    int nb[NFREQ*4*2+2]={0},b=0,m;
    int iba,ibg,nba,nbg;

    iba=xiBa(insopt); ibg=xiBg(insopt);
    nba=xnBa(insopt); nbg=xnBg(insopt);

    /* tc=0: common rtk position mode
     * tc=1: tightly-coupled mode
     * */
    tc=opt->mode==PMODE_INS_TGNSS?1:0;
    nx=tc?ins->nx:rtk->nx;
    na=tc?ins->nb:rtk->na;

    x=tc?ins->x:rtk->x;
    P=tc?ins->P:rtk->P;

    trace(3,"relpos: nx=%d nu=%d nr=%d\n",nx,nu,nr);

    dt=timediff(time,obs[nu].time);

    rs=mat(6,n); dts=mat(2,n);
    var=mat(1,n); y=mat(nf*2,n); e=mat(3,n);
    azel=zeros(2,n);
    Ri=mat(n*nf*2+2,1); Rj=mat(n*nf*2+2,1);

    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=(unsigned char)satsys(i+1,NULL);
        for (j=0;j<NFREQ;j++) rtk->ssat[i].vsat [j]=0;
        for (j=1;j<NFREQ;j++) rtk->ssat[i].snr  [j]=0;
        for (j=0;j<NFREQ;j++) rtk->ssat[i].vsatc[j]=0;
    }
    /* satellite positions/clocks */
    satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);

    /* undifferenced residuals for base station */
    if (!zdres(1,obs+nu,nr,rs+nu*6,dts+nu*2,svh+nu,nav,rtk->rb,opt,1,
               y+nu*nf*2,e+nu*3,azel+nu*2)) {
        errmsg(rtk,"initial base station position error\n");

        free(rs); free(dts); free(var);
        free(y); free(e); free(azel);
        free(Ri); free(Rj);
        return 0;
    }
    /* time-interpolation插值 of residuals (for post-processing) */
    if (opt->intpref) {
        dt=intpres(time,obs+nu,nr,nav,rtk,y+nu*nf*2);
    }
    /* select common satellites between rover and base-station */
    if ((ns=selsat(obs,azel,nu,nr,opt,sat,iu,ir))<=0) {
        errmsg(rtk,"no common satellite\n");

        free(rs); free(dts); free(var);
        free(y); free(e); free(azel);
        free(Ri); free(Rj);
        return 0;
    }
    /* temporal update of states */
    udstate(rtk,obs,sat,rs,iu,ir,ns,nav);

    xp=zeros(nx,1); xa=zeros(nx,1); Pp=zeros(nx,nx);
    dx=zeros(nx,1);

    dxwl=mat(nx,1); Qywl=mat(nx,nx);
    wlb=mat(nx,1);

    /* backup estimated estates */
    matcpy(xp,x,nx,1);
    if (tc) {
        for (i=0;i<xnCl(insopt);i++) xp[i]=1E-20;
    }
    ny=ns*nf*2+2;
    v=mat(ny,1); H=zeros(nx,ny);
    R=mat(ny,ny); bias=mat(nx,1);

    /* backup ins states */
    memcpy(&insp,ins,sizeof(insstate_t));

    /* add 2 iterations for baseline-constraint moving-base */
    niter=opt->niter+(opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0?2:0);

    for (i=0;i<niter;i++) {
        tc?insp2antp(&insp,rr):matcpy(rr,xp,1,3);

        /* undifference residuals for rover */
        if (!zdres(0,obs,nu,rs,dts,svh,nav,rr,opt,0,y,e,azel)) {
            errmsg(rtk,"rover initial position error\n");
            stat=SOLQ_NONE;
            break;
        }
        /* double-difference residuals and partial derivatives */
        if ((nv=ddres(rtk,nav,obs,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,vflg,rr,
                      Ri,Rj,nb,&b,
                      refsat,usedualfrq))<1) {
            errmsg(rtk,"no double-difference residual\n");
            stat=SOLQ_NONE;
            break;
        }
        /* kalman filter measurement update (optional IGG3 robust R) */
        if ((info=filter_rkf(rtk,obs,nav,dt,xp,Pp,P,sat,y,e,azel,iu,ir,ns,
                             v,H,R,Ri,Rj,nb,&b,vflg,rr,refsat,usedualfrq,nx,&nv))) {
            errmsg(rtk,"filter error (info=%d)\n",info);
            stat=SOLQ_NONE;
            break;
        }
        /* close loop for ins state */
        if (tc) {
            if (i==0) matcpy(dx,xp,1,xnCl(insopt));
#if 1
            if (insopt->tc==INSTC_DGPS&&degrade) {
                for (j=iba;j<iba+nba;j++) xp[j]*=0.666;
                for (j=ibg;j<ibg+nbg;j++) xp[j]*=0.666;
            }
#endif
            clp(&insp,insopt,xp);
            for (j=0;j<xnCl(insopt);j++) {
                xp[j]=1E-20;
            }
        }
    }
    tc?insp2antp(&insp,rr):matcpy(rr,xp,1,3);

    if (stat!=SOLQ_NONE&&zdres(0,obs,nu,rs,dts,svh,nav,rr,opt,0,y,e,azel)) {

        /* post-fit residuals for float solution */
        nv=ddres(rtk,nav,obs,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,NULL,
                 R,vflg,rr,Ri,Rj,nb,&b,
                 refsat,usedualfrq);

        /* validation of float solution */
        if (nv&&valpos(rtk,v,R,vflg,nv,&m,4.0)&&(tc?valins(opt,dx,ins):true)) {
#if CHKNUMERIC
            if (!chk_xp(xp,Pp,nx)) {
                fprintf(stderr,"check numeric error: nan or inf\n");
                errmsg(rtk,"check numeric error: nan or inf\n");
                stat=SOLQ_NONE;
            }
            else
#endif
            {
            /* update state and covariance matrix */
            matcpy(P,Pp,nx,nx);
            matcpy(x,xp,nx, 1);

            matcpy(dxwl,xp,1,nx);

            /* update ambiguity control struct */
            rtk->sol.ns=0;
            for (i=0;i<ns;i++) for (f=0;f<nf;f++) {
                    if (!rtk->ssat[sat[i]-1].vsat[f]) continue;
                    rtk->ssat[sat[i]-1].lock[f]++;
                    rtk->ssat[sat[i]-1].outc[f]=0;
                    if (f==0) rtk->sol.ns++; /* valid satellite count by L1 */
                }
            /* lack of valid satellites */
            if (rtk->sol.ns<3) {

                /* valid satellite count by L1-pseudorange */
                for (rtk->sol.ns=0,i=0;i<ns;i++) {
                    for (f=0;f<nf;f++) if (rtk->ssat[sat[i]-1].vsatc[f]) break;
                    rtk->sol.ns++;
                }
                /* tightly-coupled is available though lack satellite */
                if (tc) stat=SOLQ_FLOAT; else stat=SOLQ_NONE;

                /* degrade to pseudorange-dgps */
                if (m>=4&&stat==SOLQ_NONE) {
                    rtk->sol.ns=(unsigned char)m;
                    stat=SOLQ_DGPS;
                }
                rtk->sol.stat=(unsigned char)stat;
            }
            else {
                rtk->sol.stat=SOLQ_FLOAT;
            }
            }
        }
        else {
            /* strict RAIM: only if posopt[4] and clear >=6-sigma outliers exist */
            int raim_ok=0,nos=0,ntry,iex,ns0,sat0[MAXSAT],iu0[MAXSAT],ir0[MAXSAT];
            int osat[4];
            double oscore[4];
            char id[8];

            if (opt->posopt[4]&&nv>0&&
                (nos=collect_outlier_sats(v,R,vflg,nv,6.0,osat,oscore,4))>0&&
                ns>4) {

                ns0=ns;
                for (i=0;i<ns0;i++) {
                    sat0[i]=sat[i]; iu0[i]=iu[i]; ir0[i]=ir[i];
                }
                /* exclude 1..min(2,nos) largest-outlier non-ref sats only */
                for (ntry=1;ntry<=nos&&ntry<=2;ntry++) {
                    memcpy(&insp,ins,sizeof(insstate_t));
                    matcpy(xp,x,nx,1);
                    if (tc) {
                        for (i=0;i<xnCl(insopt);i++) xp[i]=1E-20;
                    }
                    for (i=0;i<MAXSAT;i++) {
                        for (j=0;j<NFREQ;j++) rtk->ssat[i].vsat [j]=0;
                        for (j=0;j<NFREQ;j++) rtk->ssat[i].vsatc[j]=0;
                    }
                    stat=rtk->opt.mode<=PMODE_DGPS?SOLQ_DGPS:SOLQ_FLOAT;

                    for (i=ns=0;i<ns0;i++) {
                        for (iex=0;iex<ntry;iex++) if (sat0[i]==osat[iex]) break;
                        if (iex<ntry) {
                            satno2id(sat0[i],id);
                            trace(2,"relpos RAIM exclude sat=%s (%.1f-sig) ntry=%d\n",
                                  id,oscore[iex],ntry);
                            continue;
                        }
                        sat[ns]=sat0[i]; iu[ns]=iu0[i]; ir[ns]=ir0[i];
                        ns++;
                    }
                    if (ns<4) break;

                    for (i=0;i<niter;i++) {
                        tc?insp2antp(&insp,rr):matcpy(rr,xp,1,3);
                        if (!zdres(0,obs,nu,rs,dts,svh,nav,rr,opt,0,y,e,azel)) {
                            stat=SOLQ_NONE; break;
                        }
                        if ((nv=ddres(rtk,nav,obs,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,
                                      vflg,rr,Ri,Rj,nb,&b,refsat,usedualfrq))<1) {
                            stat=SOLQ_NONE; break;
                        }
                        if ((info=filter_rkf(rtk,obs,nav,dt,xp,Pp,P,sat,y,e,azel,iu,ir,ns,
                                             v,H,R,Ri,Rj,nb,&b,vflg,rr,refsat,usedualfrq,
                                             nx,&nv))) {
                            stat=SOLQ_NONE; break;
                        }
                        if (tc) {
                            if (i==0) matcpy(dx,xp,1,xnCl(insopt));
                            if (insopt->tc==INSTC_DGPS&&degrade) {
                                for (j=iba;j<iba+nba;j++) xp[j]*=0.666;
                                for (j=ibg;j<ibg+nbg;j++) xp[j]*=0.666;
                            }
                            clp(&insp,insopt,xp);
                            for (j=0;j<xnCl(insopt);j++) xp[j]=1E-20;
                        }
                    }
                    if (stat==SOLQ_NONE) continue;

                    tc?insp2antp(&insp,rr):matcpy(rr,xp,1,3);
                    if (!zdres(0,obs,nu,rs,dts,svh,nav,rr,opt,0,y,e,azel)) {
                        stat=SOLQ_NONE; continue;
                    }
                    nv=ddres(rtk,nav,obs,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,NULL,
                             R,vflg,rr,Ri,Rj,nb,&b,refsat,usedualfrq);
                    if (!(nv&&valpos(rtk,v,R,vflg,nv,&m,4.0)&&
                          (tc?valins(opt,dx,ins):true))) {
                        stat=SOLQ_NONE;
                        continue;
                    }

                    /* accept only after re-validation */
#if CHKNUMERIC
                    if (!chk_xp(xp,Pp,nx)) {
                        fprintf(stderr,"check numeric error: nan or inf\n");
                        errmsg(rtk,"check numeric error: nan or inf\n");
                        stat=SOLQ_NONE;
                        continue;
                    }
#endif
                    matcpy(P,Pp,nx,nx);
                    matcpy(x,xp,nx,1);
                    matcpy(dxwl,xp,1,nx);
                    rtk->sol.ns=0;
                    for (i=0;i<ns;i++) for (f=0;f<nf;f++) {
                            if (!rtk->ssat[sat[i]-1].vsat[f]) continue;
                            rtk->ssat[sat[i]-1].lock[f]++;
                            rtk->ssat[sat[i]-1].outc[f]=0;
                            if (f==0) rtk->sol.ns++;
                        }
                    if (rtk->sol.ns<3) {
                        for (rtk->sol.ns=0,i=0;i<ns;i++) {
                            for (f=0;f<nf;f++) if (rtk->ssat[sat[i]-1].vsatc[f]) break;
                            rtk->sol.ns++;
                        }
                        if (tc) stat=SOLQ_FLOAT; else stat=SOLQ_NONE;
                        if (m>=4&&stat==SOLQ_NONE) {
                            rtk->sol.ns=(unsigned char)m;
                            stat=SOLQ_DGPS;
                        }
                        rtk->sol.stat=(unsigned char)stat;
                    }
                    else {
                        rtk->sol.stat=SOLQ_FLOAT;
                    }
                    trace(2,"relpos RAIM ok exclude %d sat(s) ns=%d\n",ntry,ns);
                    raim_ok=1;
                    break;
                }
            }
            if (!raim_ok) stat=SOLQ_NONE;
        }
    }
    /* resolve integer ambiguity by WL-NL */
    if (stat!=SOLQ_NONE&&rtk->opt.modear==ARMODE_WLNLC) {

        if (resamb_WLNL(rtk,obs,sat,iu,ir,ns,nav,azel)) {
            stat=SOLQ_FIX;
        }
    }
        /* resolve integer ambiguity by TCAR */
    else if (stat!=SOLQ_NONE&&rtk->opt.modear==ARMODE_TCARC) {

        if (resamb_TCAR(rtk,obs,sat,iu,ir,ns,nav,azel)) {
            stat=SOLQ_FIX;
        }
    }
    //         /* resolve integer ambiguity by LAMBDA */
    //    else if (stat!=SOLQ_NONE&&resamb_LAMBDA2(rtk,nav,obs,bias,xa,ddsat,&namb,wlsat,&nw,dxwl,Qywl,wlb,vflg,nv,&nny,
    //                                             &fixwl)>2) {
    // else if (stat!=SOLQ_NONE&&PAR_resamb_LAMBDA(rtk,nav,obs,bias,xa,ddsat,&namb,wlsat,&nw,dxwl,Qywl,wlb,vflg,nv,&nny,
    //                                          &fixwl)>2) {
    else if(stat!=SOLQ_NONE&&resamb_parlambda(rtk,nav,obs,bias,xa,ddsat,&namb,wlsat,&nw,dxwl,Qywl,wlb,vflg,nv,&nny,
                                             &fixwl)>2) {

        /* tigtly coupled */
        if (tc) {
            /* adjust corrections */
            if (namb<3) {
                for (j=iba;j<iba+nba;j++) xa[j]*=0.666;//加速度计bias不断缩小，防止错误固定解导致估算错误的bias
                for (j=ibg;j<ibg+nbg;j++) xa[j]*=0.666;//陀螺仪bias不断缩小，防止错误固定解导致估算错误的bias
            }
            /* close loop for ins states */
            clp(&insp,insopt,xa);

            /* convert to gps antenna position */
            //imu位置转gps天线位置
            insp2antp(&insp,rr);
        }
        else {
            matcpy(rr,xa,1,3);
        }
        /* check solutions */
        if (zdres(0,obs,nu,rs,dts,svh,nav,rr,opt,0,y,e,azel)) {

            /* post-fit residuals for fixed solution */
            nv=ddres(rtk,nav,obs,dt,xa,NULL,sat,y,e,azel,iu,ir,ns,v,NULL,R,vflg,rr,NULL,NULL,
                     NULL,NULL,refsat,usedualfrq);

            /* validation of fixed solution */
            if (nv&&valpos(rtk,v,R,vflg,nv,&m,4.0)&&(tc?valins(opt,xa,ins):true)) {
                int accept_fix=1;

#if 0
                /* Height-jump gate. Qt accepts a fix once valpos passes. */
                if (!tc) {
                    static int last_fix_ok=0;
                    static double last_fix_rr[3]={0};
                    if (last_fix_ok) {
                        double pos[3],denu[3],d[3];
                        int k;
                        for (k=0;k<3;k++) d[k]=xa[k]-last_fix_rr[k];
                        ecef2pos(last_fix_rr,pos);
                        ecef2enu(pos,d,denu);
                        if (fabs(denu[2])>0.70) {
                            trace(2,"reject fix height jump dU=%.3f ratio=%.2f\n",
                                  denu[2],rtk->sol.ratio);
                            rtk->nfix=0;
                            accept_fix=0;
                        }
                    }
                    if (accept_fix) {
                        for (i=0;i<3;i++) last_fix_rr[i]=xa[i];
                        last_fix_ok=1;
                    }
                }
#endif
                if (accept_fix) {

                /* hold integer ambiguity */
                if (++rtk->nfix>=rtk->opt.minfix&&
                    rtk->opt.modear>=ARMODE_FIXHOLD) {
                    holdamb(rtk,&insp,xa,ddsat,namb);
                }
                /* store integer ambiguity */
                storeddamb(rtk,ddsat,namb,bias);

                int *flag=imat(namb,1);
                for (i=0;i<namb;i++) flag[i]=1;

                /* search WL interger ambiguity */
                for (k=0,i=0;i<namb;i++) {
                    if (ddsat[i].flag==1) continue;

                    for (j=0;j<namb;j++) {
                        if (ddsat[j].flag==1) continue;
                        if (ddsat[i].sat1==ddsat[j].sat1&&ddsat[i].sat2==ddsat[j].sat2
                            &&flag[j]&&i!=j&&ddsat[j].f==1) {

                            wlsat[k].sat1=ddsat[j].sat1;
                            wlsat[k].sat2=ddsat[j].sat2;
                            wlsat[k].f=0;
                            wlb[k]=bias[i]-bias[j];

                            flag[i]=flag[j]=0;
                            k++;
                        }
                    }
                }
                free(flag);
                if (k) {
                    trace(3,"store WL ambiguity=\n");
                    for (i=0;i<k;i++) {
                        trace(3,"sat=%3d-%3d: WL=%8.2lf flag=%d\n",
                              wlsat[i].sat1,wlsat[i].sat2,
                              wlb[i],wlsat[i].flag);
                    }
                    /* store WL interger ambiguity */
                    storeambwl(rtk,wlsat,k,wlb);
                }
                stat=SOLQ_FIX;
                }
            }
        }
    }
    /* WL interger ambiguity fix update solution */
    else if (stat!=SOLQ_FIX&&fixwl&nny&&nw>2) {

        memcpy(&insp,ins,sizeof(insstate_t));

        /* WL fix to update rover position */
        if (tc) {
            /* close loop for ins states */
            clp(&instmp,insopt,dxwl);

            /* convert to gps antenna position */
            //imu位置转gps天线位置
            insp2antp(&instmp,rr);
        }
        else {
            matcpy(rr,dxwl,1,3);
        }
        /* check solutions */
        if (zdres(0,obs,nu,rs,dts,svh,nav,rr,opt,0,y,e,azel)) {

            /* post-fit residuals for fixed solution */
            nv=ddres(rtk,nav,obs,dt,dxwl,NULL,sat,y,e,azel,iu,ir,ns,v,NULL,R,vflg,rr,NULL,NULL,NULL,NULL,
                     refsat,usedualfrq);

            trace(3,"Qywl=\n");
            tracemat(3,Qywl,nny,nny,12,5);

            /* validation of fixed solution */
            if (nv&&valpos(rtk,v,R,vflg,nv,&m,4.0)&&(tc?valins(opt,xa,ins):true)) {

                if (!tc) {
                    matcpy(rtk->xa,dxwl,1,na);
                    for (i=0;i<na;i++) {
                        //Qywl是通过宽巷模糊度作为观测值修正得到的卡尔曼滤波后验协方差
                        for (j=0;j<na;j++) rtk->Pa[i+j*na]=Qywl[i+j*nny];
                    }
                }
                else {
                    memcpy(&insp,&instmp,sizeof(insstate_t));
                }
                /* set fix flag */
                rtk->sol.ratio=rtk->sol.wlratio;
                stat=SOLQ_FIX;

                /* store WL ambiguity */
                storeambwl(rtk,wlsat,nw,wlb);
            }
        }
    }
#if UPDNEWSAT
    /* update state by new satellites */
    for (k=0,i=0;i<nv;i++) {
        Rj[i]=SQR(100.0);

        /* check observation type */
        if (((vflg[i]>>4)&0xF)==1) continue;
        j=( vflg[i]>>8)&0xFF;
        f=((vflg[i]   )&0xF);
        if (rtk->ssat[j-1].news[f]) {

            /* new satellite */
            Rj[i]=varerr(j,rtk->ssat[j-1].sys,rtk->ssat[j-1].azel[1],baseline(rr,rtk->rb,dr),dt,f,
                         (f>=0&&f<NFREQ)?rtk->ssat[j-1].snr[f]*0.25:0.0,opt);
            k++;
            continue;
        }
    }
    if (k) {
        /* measurement error covariance */
        ddcov(nb,b,Ri,Rj,nv,R);
        if (tc) {
            for (j=0;j<xnCl(insopt);j++) {
                x[j]=1E-20;
            }
        }
        /* kalman filter */
        if ((info=filter(x,P,H,v,R,nx,nv))) {
            errmsg(rtk,"filter error (info=%d)\n",info);
            stat=SOLQ_NONE;
        }
        if (tc) {
            /* close loop for ins states */
            clp(&insp,insopt,x);
        }
    }
#endif
    if (tc&&stat!=SOLQ_NONE) {
        /* update ins states if solution is ok */
        matcpy(ins->re,insp.re,1,3);
        matcpy(ins->ve,insp.ve,1,3);
        matcpy(ins->ae,insp.ae,1,3);
        matcpy(ins->ba,insp.ba,1,3);
        matcpy(ins->bg,insp.bg,1,3);
        matcpy(ins->Ma,insp.Ma,3,3);
        matcpy(ins->Mg,insp.Mg,3,3);

        matcpy(ins->Cbe,insp.Cbe,3,3);
        matcpy(ins->lever,insp.lever,1,3);
    }
    /* fixed states covariance */
    Pa=tc?ins->Pb:rtk->Pa;

    /* index of rover station position states */
    ix=tc?xiP(insopt)+0:0; iy=tc?xiP(insopt)+1:1; iz=tc?xiP(insopt)+2:2;

    /* index of rover station velocity states */
    ivx=tc?xiV(insopt)+0:3;
    ivy=tc?xiV(insopt)+1:4;
    ivz=tc?xiV(insopt)+2:5;

    trace(2,"rtk processing : %s Q=%d\n",time_str(time,0),stat);
    /* save solution status */
    if (stat==SOLQ_FIX) {

        /* rover station position */
        tc?matcpy(rr,ins->re,1,3):matcpy(rr,rtk->xa,1,3);

        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rr[i];
            rtk->sol.qr[i]=(float)Pa[(ix+i)+(ix+i)*na];
        }
        rtk->sol.qr[3]=(float)Pa[ix+iy*na];
        rtk->sol.qr[4]=(float)Pa[iy+iz*na];
        rtk->sol.qr[5]=(float)Pa[iz+ix*na];

        /* velocity and covariance */
        if (rtk->opt.dynamics) {

            /* rover station velocity */
            tc?matcpy(rr,ins->ve,1,3):vectcpy(rr,rtk->xa,3,6);

            for (i=0;i<3;i++) rtk->sol.rr[3+i]=rr[i];

            rtk->sol.qv[0]=(float)Pa[ivx+ivx*na];
            rtk->sol.qv[1]=(float)Pa[ivy+ivy*na];
            rtk->sol.qv[2]=(float)Pa[ivz+ivz*na];

            rtk->sol.qv[3]=(float)Pa[ivx+ivy*na];
            rtk->sol.qv[4]=(float)Pa[ivy+ivz*na];
            rtk->sol.qv[5]=(float)Pa[ivz+ivx*na];
        }
    }
    else if (stat==SOLQ_FLOAT||stat==SOLQ_DGPS) {

        /* rover station position */
        tc?matcpy(rr,ins->re,1,3):matcpy(rr,rtk->x,1,3);

        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rr[i];
            rtk->sol.qr[i]=(float)P[(ix+i)+(ix+i)*nx];
        }
        rtk->sol.qr[3]=(float)P[ix+iy*nx];
        rtk->sol.qr[4]=(float)P[iy+iz*nx];
        rtk->sol.qr[5]=(float)P[iz+ix*nx];

        /* velocity and covariance */
        if (rtk->opt.dynamics) {

            /* rover station velocity */
            tc?matcpy(rr,ins->ve,1,3):vectcpy(rr,rtk->xa,3,6);

            for (i=0;i<3;i++) rtk->sol.rr[3+i]=rr[i];

            rtk->sol.qv[0]=(float)P[ivx+ivx*nx];
            rtk->sol.qv[1]=(float)P[ivy+ivy*nx];
            rtk->sol.qv[2]=(float)P[ivz+ivz*nx];

            rtk->sol.qv[3]=(float)P[ivx+ivy*nx];
            rtk->sol.qv[4]=(float)P[ivy+ivz*nx];
            rtk->sol.qv[5]=(float)P[ivz+ivx*nx];
        }
        rtk->nfix=0;
    }
    /* update reference satellites */
    for (i=0;i<NUMSYS;i++) {
        for (j=0;j<2*NFREQ;j++) rtk->refsat[i][j]=refsat[i][j];
    }
    /* backup last epoch double difference satellite */
    for (rtk->ns=0,i=0;i<nv;i++) {
        rtk->sat[rtk->ns].f   =(vflg[i]    )&0xF;
        rtk->sat[rtk->ns].sat1=(vflg[i]>>16)&0xFF;
        rtk->sat[rtk->ns].sat2=(vflg[i]>> 8)&0xFF;
        rtk->ns++;
    }
    for (i=0;i<n;i++) {
        int rcv0=obs[i].rcv-1;
        double freq;
        if (rcv0<0||rcv0>1) continue;
        if (obs[i].L[0]!=0.0&&rs) {
            for (j=0;j<3;j++) rtk->ssat[obs[i].sat-1].prs[rcv0][j]=rs[i*6+j];
        }
        for (j=0;j<nf;j++) {
            if (obs[i].L[j]==0.0) continue;
            freq=obs_freq(obs+i,j,nav);
            rtk->ssat[obs[i].sat-1].pt[rcv0][j]=obs[i].time;
            /* store Js-corrected phase/code for next-epoch clkjmp / SD detect */
            rtk->ssat[obs[i].sat-1].ph[rcv0][j]=obs[i].L[j]
                +1.0e-3*(double)rtk->opt.Js_h[rcv0]*freq;
            rtk->ssat[obs[i].sat-1].pD[rcv0][j]=obs[i].D[j];
            if (obs[i].P[j]!=0.0) {
                rtk->ssat[obs[i].sat-1].pP[rcv0][j]=obs[i].P[j]
                    +1.0e-3*(double)rtk->opt.Js_P[rcv0]*CLIGHT;
            }
        }
    }
    for (i=0;i<ns;i++) for (j=0;j<nf;j++) {

            /* output snr of rover receiver */
            rtk->ssat[sat[i]-1].snr[j]=obs[iu[i]].SNR[j];
        }
    /* update rover az/el for next epoch (relpos previously never wrote ssat.azel) */
    for (i=0;i<ns;i++) {
        rtk->ssat[sat[i]-1].azel[0]=azel[  iu[i]*2];
        rtk->ssat[sat[i]-1].azel[1]=azel[1+iu[i]*2];
    }
    log_refsat(rtk,refsat,time);
    log_cn0var(rtk,obs,sat,iu,azel,ns,time);
    for (i=0;i<MAXSAT;i++) for (j=0;j<nf;j++) {
            if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
            if (rtk->ssat[i].slip[j]&1) rtk->ssat[i].slipc[j]++;
        }
    if (stat!=SOLQ_NONE) {
        rtk->sol.stat =(unsigned char)stat;
        rtk->ins.ns   =(unsigned char)rtk->sol.ns;
        rtk->ins.gstat=(unsigned char)stat;
    }
    free(rs); free(dts); free(var); free(dx);
    free(xp); free(Pp ); free(xa);
    free(y ); free(e  ); free(azel);
    free(v ); free(H  ); free(bias);
    free(R ); free(Ri ); free(Rj);
    free(dxwl); free(Qywl); free(wlb);
    return stat!=SOLQ_NONE;
}
/* initialize rtk control ----------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*----------------------------------------------------------------------------*/
extern void rtkinit(rtk_t *rtk, const prcopt_t *opt)
{
    sol_t sol0={{0}};
    ambc_t ambc0={{{0}}};
    ssat_t ssat0={0};
    int i = 0,j=0;

    trace(3,"rtkinit :\n");

    rtk->sol=sol0;
    for (i=0;i<6;i++) rtk->rb[i]=0.0;
    rtk->nx=opt->mode<=PMODE_FIXED?NX(opt):pppnx(opt);
    rtk->na=opt->mode<=PMODE_FIXED?NR(opt):pppnx(opt);
    rtk->tt=0.0;
    rtk->x=zeros(rtk->nx,1);
    rtk->P=zeros(rtk->nx,rtk->nx);
    rtk->xa=zeros(rtk->na,1);
    rtk->Pa=zeros(rtk->na,rtk->na);
    rtk->nfix=rtk->neb=0;
    for (i=0;i<MAXSAT;i++) {
        rtk->ambc[i]=ambc0;
        rtk->ssat[i]=ssat0;
    }
    for (i=0;i<MAXERRMSG;i++) rtk->errbuf[i]=0;
    for (i=0;i<NUMSYS;i++) {
        for (j=0;j<NFREQ*2;j++) rtk->refsat[i][j]=0;
    }
    rtk->opt=*opt;
    rtk->opt.nnavsys=countnnavsys(rtk->opt.navsys);
    if (rtk->opt.rkf_k0<=0.0) rtk->opt.rkf_k0=2.5;
    if (rtk->opt.rkf_k1<=0.0) rtk->opt.rkf_k1=8.0;
    if (rtk->opt.rkf_niter<=0) rtk->opt.rkf_niter=1;
    if (rtk->opt.rkf_niter>8) rtk->opt.rkf_niter=8;
    trace(3,"rtkinit: navsys=%d nnavsys=%d\n",rtk->opt.navsys,rtk->opt.nnavsys);
    rtk->opt.insopt.soltype=rtk->opt.soltype;
    rtk->opt.insopt.gopt=&rtk->opt;

    for (i=0;i<2;i++) for (j=0;j<7;j++) rtk->opt.sind[i][j]=opt->sind[i][j];
    if (opt->mode>=PMODE_INS_UPDATE&&
        opt->mode<=PMODE_INS_TGNSS) {

        initlc(&rtk->opt.insopt,&rtk->ins);
        initodo(&rtk->opt.insopt.odopt,&rtk->ins);
        if (opt->mode==PMODE_INS_LGNSS&&opt->insopt.lcopt==IGCOM_USEOBS) {

            /* ins-gnss loosely coupled options */
            switch (opt->insopt.lc) {
                case INSLC_SINGLE: rtk->opt.mode=PMODE_SINGLE    ; break;
                case INSLC_PPK   : rtk->opt.mode=PMODE_PPP_KINEMA; break;
                case INSLC_DGPS  : rtk->opt.mode=PMODE_DGPS      ; break;
                case INSLC_RTK   : rtk->opt.mode=PMODE_KINEMA    ; break;
            }
        }
        // if (opt->mode==PMODE_INS_TGNSS) {
        //     fprintf(stdout,"best choice for dual frequency\n");
        // }
        rtk->ins.rtkp=rtk;
    }
    else {
        /* for backward solution */
        rtk->ins.rtkp=NULL; rtk->ins.x =NULL;
        rtk->ins.xa  =NULL; rtk->ins.xb=NULL;
        rtk->ins.P   =NULL; rtk->ins.Pa=NULL; rtk->ins.Pb=NULL;
    }
}
/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkfree(rtk_t *rtk)
{
    trace(3,"rtkfree :\n");

    rtk->nx=rtk->na=rtk->ins.nx=rtk->ins.nb=0;
    if (rtk->x ) free(rtk->x ); rtk->x =NULL;
    if (rtk->P ) free(rtk->P ); rtk->P =NULL;
    if (rtk->xa) free(rtk->xa); rtk->xa=NULL;
    if (rtk->Pa) free(rtk->Pa); rtk->Pa=NULL;

    if (rtk->ins.x ) free(rtk->ins.x ); rtk->ins.x =NULL;
    if (rtk->ins.P ) free(rtk->ins.P ); rtk->ins.P =NULL;
    if (rtk->ins.Pa) free(rtk->ins.Pa); rtk->ins.Pa=NULL;
    if (rtk->ins.Pb) free(rtk->ins.Pb); rtk->ins.Pb=NULL;
    if (rtk->ins.xb) free(rtk->ins.xb); rtk->ins.xb=NULL;
    if (rtk->ins.xa) free(rtk->ins.xa); rtk->ins.xa=NULL;

    if (rtk->ins.gmeas.data) {
        free(rtk->ins.gmeas.data);
    }
    rtk->ins.gmeas.data=NULL;
    rtk->ins.nx=rtk->ins.nb=0;
    rtk->ins.gmeas.n=rtk->ins.gmeas.nmax=0;

    if (rtk->bias.amb) {
        free(rtk->bias.amb);
    }
    rtk->bias.amb=NULL;
}
/* precise positioning ---------------------------------------------------------
* input observation data and navigation message, compute rover position by
* precise positioning
* args   : rtk_t *rtk       IO  rtk control/result struct
*            rtk->sol       IO  solution
*                .time      O   solution time
*                .rr[]      IO  rover position/velocity
*                               (I:fixed mode,O:single mode)
*                .dtr[0]    O   receiver clock bias (s)
*                .dtr[1]    O   receiver glonass-gps time offset (s)
*                .Qr[]      O   rover position covarinace
*                .stat      O   solution status (SOLQ_???)
*                .ns        O   number of valid satellites
*                .age       O   age of differential (s)
*                .ratio     O   ratio factor for ambiguity validation
*            rtk->rb[]      IO  base station position/velocity
*                               (I:relative mode,O:moving-base mode)
*            rtk->nx        I   number of all states
*            rtk->na        I   number of integer states
*            rtk->ns        O   number of valid satellite
*            rtk->tt        O   time difference between current and previous (s)
*            rtk->x[]       IO  float states pre-filter and post-filter
*            rtk->P[]       IO  float covariance pre-filter and post-filter
*            rtk->xa[]      O   fixed states after AR
*            rtk->Pa[]      O   fixed covariance after AR
*            rtk->ssat[s]   IO  sat(s+1) status
*                .sys       O   system (SYS_???)
*                .az   [r]  O   azimuth angle   (rad) (r=0:rover,1:base)
*                .el   [r]  O   elevation angle (rad) (r=0:rover,1:base)
*                .vs   [r]  O   data valid single     (r=0:rover,1:base)
*                .resp [f]  O   freq(f+1) pseudorange residual (m)
*                .resc [f]  O   freq(f+1) carrier-phase residual (m)
*                .vsat [f]  O   freq(f+1) data vaild (0:invalid,1:valid)
*                .fix  [f]  O   freq(f+1) ambiguity flag
*                               (0:nodata,1:float,2:fix,3:hold)
*                .slip [f]  O   freq(f+1) slip flag
*                               (bit8-7:rcv1 LLI, bit6-5:rcv2 LLI,
*                                bit2:parity unknown, bit1:slip)
*                .lock [f]  IO  freq(f+1) carrier lock count
*                .outc [f]  IO  freq(f+1) carrier outage count
*                .slipc[f]  IO  freq(f+1) cycle slip count
*                .rejc [f]  IO  freq(f+1) data reject count
*                .gf        IO  geometry-free phase (L1-L2) (m)
*                .gf2       IO  geometry-free phase (L1-L5) (m)
*            rtk->nfix      IO  number of continuous fixes of ambiguity
*            rtk->neb       IO  bytes of error message buffer
*            rtk->errbuf    IO  error message buffer
*            rtk->tstr      O   time string for debug
*            rtk->opt       I   processing options
*            rtk->ins       IO  ins/gnss coupled state
*          obsd_t *obs      I   observation data for an epoch
*                               obs[i].rcv=1:rover,2:reference
*                               sorted by receiver and satellte
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation messages
* return : status (0:no solution,1:valid solution)
* notes  : before calling function, base station position rtk->sol.rb[] should
*          be properly set for relative mode except for moving-baseline
*-----------------------------------------------------------------------------*/
extern int rtkpos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    insopt_t *insopt=&opt->insopt;
    sol_t solb={{0}};
    gtime_t time;
    static obsd_t obsd[MAXOBS];
    int fi=0,fj=1,fk=2;
    int i,j,nu,nr,stat=0,tcs=0,tcp=0;
    char msg[128]="";

    trace(3,"rtkpos  : time=%s n=%d\n",time_str(obs[0].time,3),n);
    trace(3,"obs=\n"); traceobs(3,obs,n);

    double sow;
    sow=time2gpst(obs[0].time,NULL);

    /* check tc-mode */
    tcs=opt->mode==PMODE_INS_TGNSS&&insopt->tc==INSTC_SINGLE;
    tcp=opt->mode==PMODE_INS_TGNSS&&insopt->tc==INSTC_PPK;

    if (n<=0) return 0;

    /* set base staion position */
    if (opt->refpos<=POSOPT_RINEX&&opt->mode!=PMODE_SINGLE&&opt->mode!=PMODE_MOVEB) {
        for (i=0;i<6;i++) {
            rtk->rb[i]=i<3?opt->rb[i]:0.0;
        }
    }
    /* count rover/base station observations */
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ;

    /* for rover and base observation data */
    for (i=0;i<nu+nr&&opt->adjobs;i++) {
        memcpy(&obsd[i],&obs[i],sizeof(obsd_t));

        /* adjust observation */
        if (adjsind(opt,&obs[i],&fi,&fj,&fk)) {
            trace(4,"adjust observation data signal index ok\n");
        }
        /* here just adjust three frequency */
        for (j=0;j<3;j++) {
            int src=j==0?fi:j==1?fj:fk;
            obsd[i].LLI[j]=obs[i].LLI[src];
            obsd[i].SNR[j]=obs[i].SNR[src];
            obsd[i].code[j]=obs[i].code[src];

            obsd[i].P[j]=obs[i].P[src];
            obsd[i].L[j]=obs[i].L[src];
            obsd[i].D[j]=obs[i].D[src];
        }
        /* index of frequency */
        fi=0;
        fj=1;
        fk=2;
    }
    if (opt->adjobs) {
        trace(3,"adjust obs=\n");
        traceobs(3,obsd,n);
    }
    /* Slot wavelength is GPS L2 on index 1. Galileo E5b is FREQ7.
       Take the frequency from the observation code, kinematic and tc. */
    {
        const obsd_t *obsw=opt->adjobs?obsd:obs;
        nav_t *navw=(nav_t *)nav;
        int k,f,sys,nf=NF(opt);
        for (k=0;k<nu+nr;k++) {
            sys=satsys(obsw[k].sat,NULL);
            for (f=0;f<nf;f++) {
                double freq=code2freq(sys,obsw[k].code[f],0);
                if (freq>0.0) navw->lam[obsw[k].sat-1][f]=CLIGHT/freq;
            }
        }
    }
    /* previous epoch */
    time=rtk->sol.time;

    /* rover position by single point positioning */
    if (!pntpos(opt->adjobs?obsd:obs,nu,nav,&rtk->opt,&rtk->sol,tcs?&rtk->ins:NULL,
                NULL,rtk->ssat,msg)) {
        trace(3,"rtk %d\n",stat);
        errmsg(rtk,"point pos error (%s)\n",msg);

        if (!rtk->opt.dynamics&&opt->mode<=PMODE_PPP_FIXED) {
            outsolstat(rtk);
            return 0;
        }
        /* retry single-tc mode */
        if (opt->mode==PMODE_INS_TGNSS&&insopt->tc>INSTC_SINGLE) {
            stat=pntpos(opt->adjobs?obsd:obs,nu,nav,
                        &rtk->opt,
                        &rtk->sol,
                        &rtk->ins,
                        NULL,NULL,msg);
            goto exit;
        }
    }
    if (time.time!=0) rtk->tt=timediff(rtk->sol.time,time);
    if (fabs(rtk->tt)<DTTOL&&opt->mode<=PMODE_FIXED) return stat;

    /* single point positioning */
    if (opt->mode==PMODE_SINGLE||insopt->tc==INSTC_SINGLE) {
        outsolstat(rtk);
        return 1;
    }
    /* suppress output of single solution */
    if (!opt->outsingle) {
        rtk->sol.stat=SOLQ_NONE;
    }
    /* precise point positioning */
    if ((opt->mode>=PMODE_PPP_KINEMA&&opt->mode<PMODE_INS_UPDATE)||tcp) {
        pppos(rtk,opt->adjobs?obsd:obs,nu,nav);
        outsolstat(rtk);
        return 1;
    }
    /* check number of data of base station and age of differential */
    if (nr==0) {
        errmsg(rtk,"no base station observation data for rtk\n");
        outsolstat(rtk);
        return 1;
    }
    /*  moving baseline */
    if (opt->mode==PMODE_MOVEB) {

        /* estimate position/velocity of base station */
        if (!pntpos(opt->adjobs?obsd:obs+nu,nr,nav,&rtk->opt,&solb,NULL,NULL,NULL,msg)) {
            errmsg(rtk,"base station position error (%s)\n",msg);
            return 0;
        }
        rtk->sol.age=(float)timediff(rtk->sol.time,solb.time);

        if (fabs(rtk->sol.age)>TTOL_MOVEB) {
            errmsg(rtk,"time sync error for moving-base (age=%.1f)\n",rtk->sol.age);
            return 0;
        }
        for (i=0;i<6;i++) {
            rtk->rb[i]=solb.rr[i];
        }
        /* time-synchronized position of base station */
        for (i=0;i<3;i++) {
            rtk->rb[i]+=rtk->rb[i+3]*rtk->sol.age;
        }
    }
    else {
        rtk->sol.age=(float)timediff(obs[0].time,obs[nu].time);

        if (fabs(rtk->sol.age)>opt->maxtdiff) {
            errmsg(rtk,"age of differential error (age=%.1f)\n",rtk->sol.age);
            outsolstat(rtk);
            return 1;
        }
    }
    /* relative positioning */
    stat=relpos(rtk,opt->adjobs?obsd:obs,nu,nr,nav,0,1);
    if (!stat) {
        stat=relpos(rtk,opt->adjobs?obsd:obs,nu,nr,nav,0,0);
    }
#if 1
    /* degrade to dgps-tc mode if rtk-tc fail */
    if (stat==0&&opt->mode==PMODE_INS_TGNSS) {
        insopt->tc=INSTC_DGPS;
        stat=relpos(rtk,opt->adjobs?obsd:obs,nu,nr,nav,1,1);
        insopt->tc=INSTC_RTK;
        if (stat) goto exit;
    }
    /* degrade to single-tc mode if dgps-tc fail */
    if (stat==0&&opt->mode==PMODE_INS_TGNSS) {
        stat=pntpos(opt->adjobs?obsd:obs,nu,nav,&rtk->opt,
                    &rtk->sol,
                    &rtk->ins,
                    NULL,NULL,msg);
        goto exit;
    }
exit:
#endif
    outsolstat(rtk);
    return stat;
}
