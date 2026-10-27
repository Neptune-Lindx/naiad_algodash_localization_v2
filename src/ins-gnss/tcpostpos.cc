/*-----------------------------------------------------------------------------
* tcpostpos.cc : ins-gnss tightly coupled post-processing app.
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2019/04/14 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants-----------------------------------------------------------------*/
#define INTKEEPALIVE 1000             /* keep alive interval (ms) */
#define OUTSOLFRQ    50               /* frequency of output ins solutions */
#define SOL_OUTPUT_FILE 1             /* output solution to file */
#define MINVEL        5.0             /* min velocity for initial ins states */
#define MAXGYRO       (30.0*D2R)      /* max rotation speed value for initial */
#define MAXVAR_POSE   (5.0*D2R)       /* max variance of pose measurement */
#define MAXDIFF       10.0            /* max time difference between solution */
#define ADJOBS        0               /* adjust observation data */
#define MAXINFILE   1000         /* max number of input files */

/* constants/global variables -----------------------------------------------*/
static int nepoch=0;                  /* number of observation epochs */
static int nimu=0;                    /* number of imu measurements epochs */
static int iobsu =0;                  /* current rover observation data index */
static int iobsr =0;                  /* current reference observation data index */
static int iimu  =0;                  /* current imu measurement data */
static int keepalive=1;               /* keep alive flag */
static int timeout  =10000;           /* timeout time (ms) */
static int reconnect=10000;           /* reconnect interval (ms) */
static int week=0;                    /* GPS week */
static obs_t obss={0};                /* observation data */
static nav_t navs={0};                /* navigation data */
static imu_t imus={0};                /* imu measurement data */
static stream_t moni={0};             /* monitor stream */
static stream_t frst={0};             /* solution result file */
static sta_t stas[MAXRCV];      /* station infomation */
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    trace(3,"freeobsnav:\n");

    free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}
/* adjust imu measurement data-----------------------------------------------*/
static void adj_imudata(const prcopt_t *opt,imu_t *imu)
{
    int i;
    trace(3,"adjustimu:\n");
    for (i=0;i<imu->n;i++) adjustimu(opt,&imu->data[i]);
}
/* read imu measurement data-------------------------------------------------*/
static int readimu(const char *file,int type,const prcopt_t *prcopt,imu_t *imu)
{
    int nimu=0;
    switch (type) {
        case STRFMT_M39:
            readimub(file,imu,prcopt->insopt.imudecfmt,
                     prcopt->insopt.imuformat,
                     prcopt->insopt.imucoors,
                     prcopt->insopt.imuvalfmt);
            break;
        case STRFMT_STIM300:
            readstim300(file,imu);
            break;
        case STRFMT_CSVIMU:
            readimucsv(file,imu);
            break;
        default:
            trace(2,"no matched type\n");
            break;
    }
    /* sort imu measurement data */
    nimu=sortimudata(imu);

    /* adjust imu measurement data */
    adj_imudata(prcopt,imu);
    trace(3,"imu start:%d end:%d\n",imu->data[0].time.time,imu->data[imu->n - 1].time.time);
    return nimu;
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, char **infile,
                      const int *index, int n, prcopt_t *prcopt,
                      obs_t *obs, nav_t *nav, sta_t *sta)
{
    int i,j,ind=0,nobs=0,rcv=1;

    trace(3,"readobsnav: ts=%s n=%d\n",time_str(ts,0),n);

    obs->data=NULL; obs->n =obs->nmax =0;
    nav->eph =NULL; nav->n =nav->nmax =0;
    nav->geph=NULL; nav->ng=nav->ngmax=0;
    nav->seph=NULL; nav->ns=nav->nsmax=0;
    nepoch=0;

    for (i=0;i<n;i++) {
//        if (checkbrk("")) return 0;

        if (index[i]!=ind) {
            if (obs->n>nobs) rcv++;
            ind=index[i]; nobs=obs->n;
        }
        /* read rinex obs and nav file */
        if (readrnxt(infile[i],rcv,ts,te,ti,prcopt->rnxopt[rcv<=1?0:1],obs,nav,
                     rcv<=2?sta+rcv-1:NULL)<0) {
            trace(3,"error : insufficient memory");
            trace(1,"insufficient memory\n");
            return 0;
        }
    }
    if (obs->n<=0) {
        trace(3,"error : no obs data");
        trace(1,"\n");
        return 0;
    }
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        trace(3,"error : no nav data");
        trace(1,"\n");
        return 0;
    }
    /* sort observation data */
    nepoch=sortobs(obs);

    /* observation signal index for rover and base */
    for (i=0;i<2;i++) {
        for (j=0;j<7;j++) prcopt->sind[i][j]=obs->sind[i][j];
        for (j=0;j<7;j++) nav->sind[i][j]=obs->sind[i][j];
    }
    /* delete duplicated ephemeris */
    uniqnav(nav);

    /* set time span for progress display */
    if (ts.time==0||te.time==0) {
        for (i=0;   i<obs->n;i++) if (obs->data[i].rcv==1) break;
        for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
        if (i<j) {
            if (ts.time==0) ts=obs->data[i].time;
            if (te.time==0) te=obs->data[j].time;
            settspan(ts,te);
        }
    }
    return 1;
}
/* read obs and nav data -----------------------------------------------------*/
//static int readobsnav(const char *rfile, const char *bfile, const char *gpsnavf,
//                      const char *bdsnavf,const char *glonavf,const char *mixnavf,
//                      prcopt_t *prcopt,obs_t *obs,
//                      nav_t *nav, sta_t *sta, int *nepoch)
//{
//    int i,j;
//    gtime_t ts={0},te={0};

//    obs->data=NULL; obs->n =obs->nmax =0;
//    nav->eph =NULL; nav->n =nav->nmax =0;
//    nav->geph=NULL; nav->ng=nav->ngmax=0;
//    nav->seph=NULL; nav->ns=nav->nsmax=0;
//    *nepoch=0;

//    /* read all data*/
//    if (  rfile&&readrnxt(rfile  ,1,ts,te,0.0,prcopt->rnxopt[0],obs,nav,NULL)<0) return 0;
//    if (  bfile&&readrnxt(bfile  ,2,ts,te,0.0,prcopt->rnxopt[1],obs,nav,NULL)<0) return 0;
//    if (gpsnavf&&readrnxt(gpsnavf,0,ts,te,0.0,prcopt->rnxopt[0],obs,nav,NULL)<0) return 0;
////    if (bdsnavf&&readrnxt(bdsnavf,0,ts,te,0.0,prcopt->rnxopt[0],obs,nav,NULL)<0) return 0;
////    if (glonavf&&readrnxt(glonavf,0,ts,te,0.0,prcopt->rnxopt[0],obs,nav,NULL)<0) return 0;
////    if (mixnavf&&readrnxt(mixnavf,0,ts,te,0.0,prcopt->rnxopt[0],obs,nav,NULL)<0) return 0;
//    if (obs->n<=0) {
//        trace(1,"\n");
//        return 0;
//    }
//    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
//        trace(1,"\n");
//        return 0;
//    }
//    /* sort observation data */
//    *nepoch=sortobs(obs);

//    /* observation signal index for rover and base */
//    for (i=0;i<2;i++) {
//        for (j=0;j<7;j++) prcopt->sind[i][j]=obs->sind[i][j];
//        for (j=0;j<7;j++) nav->sind[i][j]=obs->sind[i][j];
//    }
//    /* delete duplicated ephemeris */
//    uniqnav(nav);
//    return *nepoch>0;
//}
/* initialization position mode/ionosphere and troposphere option------------*/
static void initrtkpos(rtk_t *rtk,prcopt_t *prcopt)
{
    prcopt->mode   =PMODE_KINEMA;
    prcopt->ionoopt=IONOOPT_BRDC;
    prcopt->tropopt=TROPOPT_SAAS;
#if ADJOBS
    prcopt->adjobs=1;
#endif
    rtkinit(rtk,prcopt);
}
/* solution convert to velocity----------------------------------------------*/
static void sol2vel(const sol_t *sol1,const sol_t *sol2,double *v)
{
    v[0]=(sol1->rr[0]-sol2->rr[0])/timediff(sol1->time,sol2->time);
    v[1]=(sol1->rr[1]-sol2->rr[1])/timediff(sol1->time,sol2->time);
    v[2]=(sol1->rr[2]-sol2->rr[2])/timediff(sol1->time,sol2->time);
}
/* initial ins states--------------------------------------------------------*/
static int init_ins(const imud_t *imu,const obsd_t *obs, int n, const nav_t *nav,
                    const prcopt_t *prcopt,
                    insstate_t *ins,
                    int ins_mode)
{
    /* global variables for rtk positioning */
    static int first=1,i,minsol=5;
    static prcopt_t popt=*prcopt;
    static rtk_t rtk={0};
    static sol_t sols[5]={0};
    double vr[3]={0};

    if (n<=0) {
        trace(2,"no observation data to initial\n");
        return 0;
    }
    /* initial gps position options */
    if (first) {
        initrtkpos(&rtk,&popt); first=0;
    }
    rtkpos(&rtk,obs,n,nav);

    /* save position solution to buffer */
    for (i=0;i<minsol-1;i++)sols[i]=sols[i+1]; sols[i]=rtk.sol;
    for (i=0;i<minsol;i++) {
        if (sols[i].stat>popt.insopt.iisu||sols[i].stat==SOLQ_NONE) {
            trace(2,"check solution status fail,%d\n",popt.insopt.iisu);
            return 0;
        }
    }
    for (i=0;i<minsol-1;i++) {
        if (timediff(sols[i+1].time,sols[i].time)>3.0) {
            return 0;
        }
    }
    /* compute velocity from solutions */
    matcpy(vr,sols[minsol-1].rr+3,1,3);
    if (norm(vr,3)==0.0) {
        sol2vel(sols+minsol-1,sols+minsol-2,vr);
    }
    if(ins_mode == PMODE_INS_TGNSS){
        /* check velocity ok? */
        trace(3,"norm(vr,3):%lf\n",norm(vr,3));
        if (norm(vr,3)<MINVEL||norm(imu->gyro,3)>MAXGYRO) return 0;

        for (i=0;i<minsol-1;i++) {
            if (timediff(sols[i+1].time,sols[i].time)>MAXDIFF) {
                trace(2,"large time difference of solution\n");
                return 0;
            }
            if (fabs(timediff(sols[i+1].time,sols[i].time))<1E-5) {
                trace(2,"duplicate gps measurement\n");
                return 0;
            }
        }
        /* initialize ins states */
        if (!ant2inins(sols[minsol-1].time,sols[minsol-1].rr,vr,&popt.insopt,
                       NULL,ins,NULL)) {
            return 0;
        }
        /* update ins state in n-frame */
        updinsn(ins);
        ins->time=sols[minsol-1].time;

        /* reset rtk position options */
        rtkfree(&rtk);
        first=1;

        trace(3,"initial ins state ok\n");
        return 1;
    }else{
        trace(3,"initial rtk ok\n");
        return 1;
    }
}
/* open output solution file--------------------------------------------------*/
static int open_solfile(const char *file)
{
    if (!stropen(&frst,STR_FILE,STR_MODE_W,file)) return 0;
    return 1;
}
/* thread to send keep alive for monitor port --------------------------------*/
static void *sendkeepalive(void *arg)
{
    trace(3,"sendkeepalive: start\n");

    while (keepalive) {
        strwrite(&moni,(unsigned char *)"\r",1);
        sleepms(INTKEEPALIVE);
    }
    trace(3,"sendkeepalive: stop\n");
    return NULL;
}
/* open monitor port ---------------------------------------------------------*/
static int openmoni(int port)
{
    pthread_t thread;
    char path[64];

    trace(3,"openmomi: port=%d\n",port);

    sprintf(path,":%d",port);

    if (!stropen(&moni,STR_TCPSVR,STR_MODE_RW,path)) return 0;
    strsettimeout(&moni,timeout,reconnect);
    keepalive=1;
    pthread_create(&thread,NULL,sendkeepalive,NULL);
    return 1;
}
/* find observation index for tightly coupling--------------------------------*/
static int fnobs(gtime_t imut,int *iobs,const imu_t *imu,const obs_t *obs)
{
    double sow1,sow2;
    int i,info=0;

    for (i=*iobs-100<0?0:*iobs-10;i<obs->n;i++) {
        sow1=time2gpst(obs->data[i].time,&week);
        sow2=time2gpst(imut,NULL);
        if (fabs(sow1-sow2)<DTTOL) {
            *iobs=i; info=1;
            break;
        }
        else if (sow1-sow2>2.0*DTTOL) {
            info=0; break;
        }
    }
    return info;
}
/* input imu measurement data-------------------------------------------------*/
static int inputimu(imud_t *imudata,const prcopt_t* opt,imud_t *imuz,int ws)
{
    int i,k;
    if (iimu<0||iimu>=imus.n) {
        return 0; /* no imu measurement data */
    }
    /* prepare imu data for static detect 静止采集imu，用于做静止检测*/
    for (k=0,i=iimu;k<ws&&i>=0&&i<imus.n&&imuz;k++,i++) {
        imuz[k]=imus.data[i];
    }
    /* forward input */
    *imudata=imus.data[iimu++];
//    imudata->time=timeadd(imudata->time,week*604800.0);
//     trace(3,"iimu:%d %d %d %d\n",iimu,imus.n,imudata->time.time,week);
    return 1;
}
/* search next observation data index ----------------------------------------*/
extern int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs,int solq, const prcopt_t *popt)
{
    int i,nu,nr,n=0;

    if (0<iobsu&&iobsu>=obss.n) {
        return 0;
    }
    if ((nu=nextobsf(&obss,&iobsu,1))<=0) return -1;
    if (popt->intpref) {
        for (;(nr=nextobsf(&obss,&iobsr,2))>0;iobsr+=nr)
            if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)>-DTTOL) break;
    }
    else {
        for (i=iobsr;(nr=nextobsf(&obss,&i,2))>0;iobsr=i,i+=nr)
            if (timediff(obss.data[i].time,obss.data[iobsu].time)>DTTOL) break;
    }
    nr=nextobsf(&obss,&iobsr,2);
    if (nr<=0) {
        nr=nextobsf(&obss,&iobsr,2);
    }
    for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu+i];
    for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr+i];
    trace(3,"su:%d\n",iobsu+i);
    iobsu+=nu;
    return n;
}
/* write solution header to output stream ------------------------------------*/
static void writesolhead(stream_t *stream, const solopt_t *solopt)
{
    unsigned char buff[1024];
    int n;

    n=outsolheads(buff,solopt);
    strwrite(stream,buff,n);
}
/* write solution status output stream ---------------------------------------*/
static int wrt_solution(rtk_t *rtk,const solopt_t *solopt)
{
    unsigned char buff[1024];
    static int c=0,out_head=0;
    int n=0;

    if (rtk->sol.time.time==0) return 0;
    if (!out_head) {
        writesolhead(&frst,solopt);
        out_head=1;
    }
#if SOL_OUTPUT_FILE
    /* output solution to file */
    n=outsols(buff,&rtk->sol,rtk->rb,solopt,&rtk->ins,&rtk->opt.insopt,0);
    if (frst.port) {
        strwrite(&frst,buff,n);
    }
#endif
    /* output solution to monitor */
    if (moni.port) {
        n=outsols(buff,&rtk->sol,rtk->rb,solopt,&rtk->ins,&rtk->opt.insopt,1);
        if (c++>OUTSOLFRQ) {
            strwrite(&moni,buff,n);
            c=0;
        }
    }
    return n;
}
/* update solution status----------------------------------------------------*/
static void update_stat(const gmea_t *gm,insstate_t *ins)
{
    /* update ins solution status */
    if (ins->stat==INSS_LCUD&&gm) {
        ins->ns=gm->ns; ins->gstat=gm->stat;
    }
    else if (ins->stat==INSS_MECH) {
        ins->ns=0;
        ins->gstat=SOLQ_NONE;
    }
}
/* tightly coupled------------------------------------------------------------*/
static int tcfilt(const prcopt_t *popt,const solopt_t *solopt)
{
    static gtime_t last_time;
    double pos[3];
    int i,flag=0,nobs,n=0,nc=0,init=0,ws=10,zf;
    obsd_t obs[MAXOBS*2]={{0}};
    imud_t imud={0},*imuz;
    rtk_t rtk={0};
    rtkinit(&rtk,popt);

    /* static detect window size */
    ws=rtk.opt.insopt.zvopt.ws<=0?5:rtk.opt.insopt.zvopt.ws;
    imuz=(imud_t*)malloc(sizeof(imud_t)*ws);

    /* tightly coupled process */
    while (inputimu(&imud,popt,imuz,0)) {

        /* match observation for imu measurement data */
        flag=fnobs(imud.time,&iobsu,&imus,&obss);
        if (flag) {

            nobs=inputobs(obs,rtk.sol.stat,popt);
            //暂且先这么改，目前会间隔不更新数据
            if(obs[0].time.time == last_time.time){
                nobs = 0;
            }
            last_time = obs[0].time;
            if (nobs) {
                /* get all gnss measurement data */
                for (i=n=0;i<nobs;i++) {
                    if ((satsys(obs[i].sat,NULL)&popt->navsys)&&popt->exsats[obs[i].sat-1]!=1) {
                        obs[n++]=obs[i];
                    }
                }
                if (n<=0) continue;

                /* carrier-phase bias correction */
                if (navs.nf>0) {
                    corr_phase_bias_fcb(obs,n,&navs);
                }
                else if (!strstr(popt->pppopt,"-DIS_FCB")) {
                    corr_phase_bias_ssr(obs,n,&navs);
                }
                /* initial ins states */
                if (init==0) {
                    if (init_ins(&imud,obs,n,&navs,popt,&rtk.ins,rtk.opt.mode)) init=1;
                    continue;
                }

                /* tightly coupled */
                tcigpos(&rtk.opt,obs,n,&navs,&imud,&rtk,&rtk.ins,INSUPD_MEAS);

                /* doppler measurement aid */
                if (popt->insopt.dopp) {
                    doppler(obs,nobs,&navs,&rtk.opt,&rtk.ins);
                }
            }
        }
        else {
            if(rtk.opt.mode == PMODE_INS_TGNSS){
                /* ins still initialing */
                if (init==0) continue;
                /* ins mechanization */
                tcigpos(&rtk.opt,NULL,n,&navs,&imud,&rtk,&rtk.ins,INSUPD_TIME);
            }
        }
        if(rtk.opt.mode == PMODE_INS_TGNSS){
            /* non-holonomic constraint 非完整性约束*/
            if (popt->insopt.nhc
                &&(nc++>rtk.opt.insopt.nhz?nc=0,true:false)) {
                nhc(&rtk.ins,&rtk.opt.insopt,&imud);
            }
            /* odometry velocity aid */
            if (popt->insopt.odo) {
                odo(&rtk.opt.insopt,&imud,&imud.odo,&rtk.ins);
            }
            /* static detector */
            if (flag) {
                /* zero velocity/zero angular rate update */
                ecef2pos(rtk.ins.re,pos);

                zf=detstc(imuz,ws,&rtk.opt.insopt,pos);
                if (popt->insopt.odo) {
                    zf|=detstatic_ODO(&rtk.opt.insopt,&imud.odo);
                }
                /* zero velocity update */
                if (zf&&popt->insopt.zvu) {
                    zvu(&rtk.ins,&rtk.opt.insopt,imuz,1);
                }
                /* zero angular rate update */
                if (zf&&popt->insopt.zaru) {
                    zaru(&rtk.ins,&rtk.opt.insopt,imuz,1);
                }
            }
            /* solution. */
            ins2sol(&rtk.ins,&rtk.opt.insopt,&rtk.sol);
        }
        /* write solution */
        if (!wrt_solution(&rtk,solopt)) continue;
    }
    rtkfree(&rtk);
    free(imuz);
    return init;
}
/* close monitor port---------------------------------------------------------*/
static int close_moni(stream_t *moni)
{
    trace(3,"closemoni:\n");
    keepalive=0;

    /* send disconnect message */
    strwrite(moni,(unsigned char *)MSG_DISCONN,strlen(MSG_DISCONN));

    /* wait fin from clients */
    sleepms(1000);
    strclose(moni);
}
/* average of single position ------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  const prcopt_t *opt)
{
    obsd_t data[MAXOBS];
    gtime_t ts={0};
    sol_t sol={{0}};
    int i,j,n=0,m,iobs;
    char msg[128];

    trace(3,"avepos: rcv=%d obs.n=%d\n",rcv,obs->n);

    for (i=0;i<3;i++) ra[i]=0.0;

    for (iobs=0;(m=nextobsf(obs,&iobs,rcv))>0;iobs+=m) {

        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 hz */

        if (!pntpos(data,j,nav,opt,&sol,NULL,NULL,NULL,msg)) continue;

        for (i=0;i<3;i++) ra[i]+=sol.rr[i];
        n++;
    }
    if (n<=0) {
        trace(1,"no average of base station position\n");
        return 0;
    }
    for (i=0;i<3;i++) ra[i]/=n;
    return 1;
}
/* station position from file ------------------------------------------------*/
static int getstapos(const char *file, char *name, double *r)
{
    FILE *fp;
    char buff[256],sname[256],*p,*q;
    double pos[3];

    trace(3,"getstapos: file=%s name=%s\n",file,name);

    if (!(fp=fopen(file,"r"))) {
        trace(1,"station position file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if ((p=strchr(buff,'%'))) *p='\0';

        if (sscanf(buff,"%lf %lf %lf %s",pos,pos+1,pos+2,sname)<4) continue;

        for (p=sname,q=name;*p&&*q;p++,q++) {
            if (toupper((int)*p)!=toupper((int)*q)) break;
        }
        if (!*p) {
            pos[0]*=D2R;
            pos[1]*=D2R;
            pos2ecef(pos,r);
            fclose(fp);
            return 1;
        }
    }
    fclose(fp);
    trace(1,"no station position: %s %s\n",name,file);
    return 0;
}
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile)
{
    double *rr=rcvno==1?opt->ru:opt->rb,del[3],pos[3],dr[3]={0};
    int i,postype=rcvno==1?opt->rovpos:opt->refpos;
    char *name;

    trace(3,"antpos  : rcvno=%d\n",rcvno);

    if (postype==POSOPT_SINGLE) { /* average of single position */
        if (!avepos(rr,rcvno,obs,nav,opt)) {
            showmsg("error : station pos computation");
            return 0;
        }
    }
    else if (postype==POSOPT_FILE) { /* read from position file */
        name=stas[rcvno==1?0:1].name;
        if (!getstapos(posfile,name,rr)) {
            showmsg("error : no position of %s in %s",name,posfile);
            return 0;
        }
    }
    else if (postype==POSOPT_RINEX) { /* get from rinex header */
        if (norm(stas[rcvno==1?0:1].pos,3)<=0.0) {
            showmsg("error : no position in rinex header");
            trace(1,"no position position in rinex header\n");
            return 0;
        }
        /* antenna delta */
        if (stas[rcvno==1?0:1].deltype==0) { /* enu */
            for (i=0;i<3;i++) del[i]=stas[rcvno==1?0:1].del[i];
            del[2]+=stas[rcvno==1?0:1].hgt;
            ecef2pos(stas[rcvno==1?0:1].pos,pos);
            enu2ecef(pos,del,dr);
        }
        else { /* xyz */
            for (i=0;i<3;i++) dr[i]=stas[rcvno==1?0:1].del[i];
        }
        for (i=0;i<3;i++) rr[i]=stas[rcvno==1?0:1].pos[i]+dr[i];
    }
    return 1;
}
/* ins-gnss tightly coupled post-processing ---------------------------------
 * args:
 *        prcopt_t *popt    I  options
 *        solopt_t *solopt  I  solution options
 *        int port          I  port for solution monitor (0: none)
 *        char  *outfile    I  file path of output solution (NULL: no output)
 *        char **infiles    I  [0]: rover observation data file (RINEX)
 *                             [1]: base observation data file (RINEX)
 *                             [2]: gps-navigation file
 *                             [3]: bds-navigation file
 *                             [4]: glo-navigation file
 *                             [5]: mixed-navigation file
 *                             [6]: imu measurement data file
 * return: status (1: ok,0: fail)
 * ---------------------------------------------------------------------------*/
extern int tcpostpos(prcopt_t *popt, const solopt_t *solopt, const filopt_t *fopt,
                     int port, const char *outfile, char** infiles)
{
    int flag=1;
    gtime_t ts,te;
    char tracefile[1024];
    ts.time = 0;te.time = 0;
    ts.sec = 0;te.sec = 0;
    double ti = 0;
    int n  = 3;//3个文件，后期得修改
    int index[MAXINFILE]={0};
    for (int i=0;i<n;i++) index[i]=i;

    /* open debug trace (same rule as postpos) */
    if (solopt->trace>0) {
        if (outfile&&*outfile) {
            strcpy(tracefile,outfile);
            strcat(tracefile,".trace");
        }
        else if (fopt&&*fopt->trace) {
            strcpy(tracefile,fopt->trace);
        }
        else {
            tracefile[0]='\0';
        }
        if (tracefile[0]) {
            traceclose();
            traceopen(tracefile);
            tracelevel(solopt->trace);
        }
    }

    trace(3,"tcpostpos: port=%d  file=%s\n",port,outfile);

    /* open result file */
    if (outfile&&!open_solfile(outfile)) {
        trace(2,"open file fail\n");
        flag=0;
        goto exit;
    }
    /* open monitor */
    if (port&&!openmoni(port)) {
        trace(2,"open monitor fail\n");
        flag=0;
        goto exit;
    }
 
    /* read imu data 如果有INS模式,才读取imu数据*/
    if(popt->mode >=9){
        if (!(nimu=readimu(infiles[6],STRFMT_CSVIMU,popt,&imus))) {
            trace(2,"read imu data fail\n");
            flag=0;
            goto exit;
        }
    }
    /* read nav/obs data */
    if(!readobsnav(ts,te,ti,infiles,index,n,popt,&obss,&navs,stas)){
//    if (!readobsnav(infiles[0],infiles[1],infiles[2],infiles[3],infiles[4],infiles[5],popt,
//                    &obss,&navs,stas,&nepoch)) {
        trace(2,"read obs fail\n");
        flag=0;
        goto exit;
    }
    /* read dcb parameters */
//    if (*fopt->dcb) {
//        reppath(fopt->dcb,path,ts,"","");
//        readdcb(path,&navs,stas);
//    }
    /* set antenna paramters */
//    if (popt_.mode!=PMODE_SINGLE) {
//        setpcv(obss.n>0?obss.data[0].time:timeget(),&popt_,&navs,&pcvss,&pcvsr,
//               stas);
//    }
    /* read ocean tide loading parameters */
//    if (popt_.mode>PMODE_SINGLE&&*fopt->blq) {
//        readotl(&popt_,fopt->blq,stas);
//    }
    /* rover/reference fixed position */
//    if (popt_.mode==PMODE_FIXED) {
//        if (!antpos(&popt_,1,&obss,&navs,stas,fopt->stapos)) {
//            freeobsnav(&obss,&navs);
//            return 0;
//        }
//    }
//    else
//    if (PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_STATIC) {
        if (!antpos(popt,2,&obss,&navs,NULL,NULL)) {
            freeobsnav(&obss,&navs);
            return 0;
        }
//    }
    /* tightly coupled */
    if (!tcfilt(popt,solopt)) {
        trace(2,"tightly coupled fail\n");
        flag=0;
        goto exit;
    }
exit:
    /* close monitor/file */
    if (port) close_moni(&moni);
    if (outfile) strclose(&frst);
    if (imus.n) {
        free(imus.data);
    }
    freeobsnav(&obss,&navs);
    if (solopt->trace>0) traceclose();
    return flag;
}
