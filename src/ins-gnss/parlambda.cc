#include "parlambda.h"

// find minimum in s[] and return its index
static int find_min(double *s, int n)
{
	int i, index = -1;
	double MAX = 1000000000.0;
	if (s == NULL)
	{
		return -1;
	}
	for (i = 0; i < n; i++)
	{
		if (s[i] < MAX)
		{
			MAX = s[i];
			index = i;
		}

	}
	return index;
}

/* matrix diagonal -> val (n x 1) */
static void diag(const double *matrix, int n, double *val)
{
	int i;
	if (matrix==NULL||val==NULL) return;
	for (i=0;i<n;i++) val[i]=matrix[i*n+i];
}

/* find maximum in s[] and return its index */
static int find_max(double *s, int n)
{
	int i, index=-1;
	double vmin=0.0;
	if (s==NULL) return -1;
	for (i=0;i<n;i++) {
		if (s[i]>vmin) {
			vmin=s[i];
			index=i;
		}
	}
	return index;
}
/* max decorrelated variance index in Qahat via Qz=Z'*Q*Z (PAD-style) */
static int find_max_qz_index(const double *Qahat, int n, double *Z, double *ZTQ,
							double *Qzhat, double *Qvar)
{
	int i, ind, Qahat_index;
	double amax, ai;

	if (Qahat==NULL||n<=0||Z==NULL||ZTQ==NULL||Qzhat==NULL||Qvar==NULL) return -1;
	if (lambda_reduction(n,Qahat,Z)) return -1;
	matmul("TN",n,n,n,1.0,Z,Qahat,0.0,ZTQ); /* Z'*Qahat */
	matmul("NN",n,n,n,1.0,ZTQ,Z,0.0,Qzhat);  /* Qzhat=Z'*Qahat*Z */
	diag(Qzhat,n,Qvar);
	ind=find_max(Qvar,n);
	if (ind<0||ind>=n) return -1;
	/* map z-index to ahat: dominant |Z(:,ind)| (approx. PAR record[]) */
	Qahat_index=0;
	amax=fabs(Z[ind*n]);
	for (i=1;i<n;i++) {
		ai=fabs(Z[i+ind*n]);
		if (ai>amax) { amax=ai; Qahat_index=i; }
	}
	return Qahat_index;
}

static int find_satadd(int *pafixvalid, int n) {
	int i,j=0;
	for(i=0;i<n;i++) {
		if (pafixvalid[i]) {
			j++;
		}
	}
	if(j < PAR_MIN_SAT) return -1;

	for (i = 0; i < n; i++) {
		if (!pafixvalid[i]) {
			return i;
		}
	}
	return -1;
}

static int find_min_ele(double *sat_azi,double *sat_ele, int n,int area) {
	double *s;
	int i,k;
	s = mat(n,1);
	memcpy(s,sat_ele, n * sizeof(double));
	for (i = 0; i < n; i++) {
		if ((sat_azi[i] < (area*90.0*D2R)) || (sat_azi[i] >= ((area + 1)*90.0*D2R))) { 
			s[i] = 90.0*D2R; 
		}
	}
	k = find_min(s, n);
	free(s);
	return k;
}

static int find_min_i(int *s, int n)
{
	int i, index = -1;
	int MAX = 1000000000;
	if (s == NULL)
	{
		return -1;
	}
	for (i = 0; i < n; i++)
	{
		if (s[i] < MAX)
		{
			MAX = s[i];
			index = i;
		}

	}
	return index;
}

/* lock PAR only after tracking is mature and lock counts spread enough */
#ifndef PAR_LOCK_READY
#define PAR_LOCK_READY  5   /* min lock epochs to treat a DD as mature */
#endif
#ifndef PAR_LOCK_GAP
#define PAR_LOCK_GAP    3   /* min (max_lock-min_lock) to allow lock-based drop */
#endif
static int lock_par_ready(const rtk_t *rtk, int n)
{
	int i,f,lk,nready=0,lmin=1000000000,lmax=-1000000000;
	int ready_thres;

	if (rtk==NULL||n<=0) return 0;
	ready_thres=rtk->opt.minlock>0?rtk->opt.minlock:1;
	if (ready_thres<PAR_LOCK_READY) ready_thres=PAR_LOCK_READY;

	for (i=0;i<n;i++) {
		f=rtk->ddfreq[i];
		if (f<0||f>=NFREQ) continue;
		lk=rtk->ssat[rtk->ddsat2[i]-1].lock[f];
		if (lk>lmax) lmax=lk;
		if (lk<lmin) lmin=lk;
		if (lk>=ready_thres) nready++;
	}
	if (nready<PAR_MIN_SAT) return 0;
	if (lmax<lmin||(lmax-lmin)<PAR_LOCK_GAP) return 0;
	return 1;
}

// remove one row and one col at "remove_index" in Q 
static void remove_mat(double *Q, int n, int remove_index)
{
    int i, j, row=0;
    int col;
	double *newQ;

	if (Q==NULL || n<=1 || remove_index<0 || remove_index>=n) return;
	newQ = mat(n-1,n-1);
	if (newQ==NULL) return;
    for (i = 0; i < n; i++)
    {
        if (i!=remove_index)
        {
            col = 0;
            for (j = 0; j < n; j++)
            {
                if (j!=remove_index)
                {
                    newQ[row*(n-1)+col]=Q[i*n+j];
                    col++;
                }
            }
            row++;
        }

    }
	memcpy(Q, newQ, (n-1)*(n-1)*sizeof(double));
	free(newQ);

}

static void remove_Qab(double *Qab,int na,int nb,int remove_index)
{
	int i,j,row = 0;
	if (Qab == NULL || remove_index >= na)return;
	for (i = 0; i < na; i++){
		for (j = 0; j < nb; j++) {
			Qab[row*nb + j] = Qab[i*nb + j];
		}
		if (i != remove_index)row++;
	}
	for (j = 0; j < nb; j++)Qab[row*nb + j] = 0;
}
// remove one row at "remove_index" in x and get new x
static void remove_vec(double *x, int n, int remove_index)
{
	int i, row = 0;;
	if (x == NULL|| remove_index < 0 || remove_index >= n)return;
	for (i = 0; i < n; i++)
	{
		if (i != remove_index)x[row++] = x[i];
	}
	x[row] = 0;
}
// remove one row at "remove_index" in x and get new x
static void remove_sat(int *sat, int n, int remove_index)
{
	int i, row = 0;;
	if (sat == NULL || remove_index < 0 || remove_index >= n)return;
	for (i = 0; i < n; i++)
	{
		if (i != remove_index)sat[row++] = sat[i];
	}
	sat[row] = 0;
}
/* first DD-entry index that has ambslip flag (1 flag per DD entry) */
static int find_ambslip_sat(const int *ambslip, int nsat)
{
	int i;
	if (ambslip==NULL) return -1;
	for (i=0;i<nsat;i++) {
		if (ambslip[i]) return i;
	}
	return -1;
}

static void mat_split(double *Qz, double *Qz1, double *Qz2, double *Qz1z2, double *Qz2z1, int n, int a)
{
	int j, i;
	for (j = 0; j < a; j++) for (i = 0; i < a; i++) {
		Qz1[a*j + i] = Qz[n*j + i];
	}
	for (j = 0; j < n - a; j++) for (i = 0; i < n - a; i++) {
		Qz2[(n - a)*j + i] = Qz[n*(j + a) + i + a];
	}
	for (j = 0; j < a; j++) for (i = 0; i < n - a; i++) {
		Qz1z2[(n - a)*j + i] = Qz[n*j + i + a];
	}
	for (j = 0; j < n - a; j++) for (i = 0; i < a; i++) {
		Qz2z1[a*j + i] = Qz[n*(j + a) + i];
	}
}
// find integer val in s[] and return its index
static int findi(const int *s, int n, int val)
{
	int i;
	if (s == NULL)
	{
		return -1;
	}
	for (i = 0; i < n; i++)
	{
		if (s[i] == val)
		{
			return i;
		}

	}
	return -1;
}

/* DD system index: 1=GLO, 2=GAL, 3=BDS, 4=GPS/QZS/SBS (same as ddres) */
static int par_sysind(int sys)
{
	switch (sys) {
		case SYS_GPS: return 4;
		case SYS_QZS: return 4;
		case SYS_SBS: return 4;
		case SYS_GLO: return 1;
		case SYS_GAL: return 2;
		case SYS_CMP: return 3;
	}
	return 4;
}
/* whether DD system-group m is enabled in navsys */
static int par_group_enabled(int navsys, int m)
{
	switch (m) {
		case 1: return (navsys&SYS_GLO)!=0;
		case 2: return (navsys&SYS_GAL)!=0;
		case 3: return (navsys&SYS_CMP)!=0;
		case 4: return (navsys&(SYS_GPS|SYS_QZS|SYS_SBS))!=0;
	}
	return 0;
}
/* recount per-group fixed sats from fixsat[] (unique sat, multi-freq once) */
static void update_nfixsat_s(rtk_t *rtk)
{
	int i,m,sat;
	int seen[MAXSAT];

	for (i=0;i<6;i++) rtk->nfixsat_s[i]=0;
	memset(seen,0,sizeof(seen));

	for (i=0;i<rtk->nfixsat;i++) {
		sat=rtk->fixsat[i];
		if (sat<=0||sat>MAXSAT) continue;
		if (seen[sat-1]) continue;
		seen[sat-1]=1;
		m=par_sysind(rtk->ssat[sat-1].sys);
		if (m>=0&&m<6) rtk->nfixsat_s[m]++;
	}
}
/* true if every enabled DD group has fewer than PAR_MIN_SAT fixed sats */
static int all_groups_lack_fix(rtk_t *rtk)
{
	int i,lack=1;

	for (i=0;i<6;i++) {
		if (!par_group_enabled(rtk->opt.navsys,i)) continue;
		lack=lack&&(rtk->nfixsat_s[i]<PAR_MIN_SAT);
	}
	return lack;
}

//判断当前历元的卫星是否全部包含锁定卫星
static int contain_fixsat(rtk_t *rtk,int n)
{
	int i,j,f,contain_num=0;
	bool lackfix_flag = false;
	if(n < PAR_MIN_SAT) return 0;

	lackfix_flag = rtk->nfixsat < (PAR_MIN_SAT+rtk->opt.nnavsys);
	lackfix_flag = lackfix_flag && all_groups_lack_fix(rtk);
	if(lackfix_flag) return 0;

	for (f = 0; f < NFREQ; f++) {
		for (i = 0; i < rtk->nfixsat; i++) for (j = 0; j < n; j++) {
			//在不发生周跳和半周跳的前提下，锁定卫星包含在全集卫星中，则contain_num加1
			if (rtk->fixsat[i] == rtk->ddsat2[j] && (!(rtk->ssat[rtk->fixsat[i] - 1].slip[f] & 0x03)))
				contain_num++;
		}
	}
	if (contain_num == rtk->nfixsat) {
		trace(3, "contain fixsat ok\n");
		return 1;
	}else {
		trace(3, "contain fixsat failed\n");
		//trace(3, "fix stat:"); for (i = 0; i < rtk->nfixsat; i++)trace(3, "%d,", rtk->fixsat[i]); trace(3, "\n");
		//trace(3, "full stat:"); for (i = 0; i < n; i++)trace(3, "%d,", rtk->ddsat2[i]); trace(3,"\n");
 		return 0;
	}
}

//部分模糊度映射到全局模糊度
static void par2fullamb(const double *fa,const double *Qa_,double *afixhat,int fixed,int n,int m,int *delete_record,double *afix)
{
	int unfixed, inc = 0, i, j;
	double *record_Z_row, *Qa, *ZTQ, *Qzhat_record, *Za_record;
	double *ZQaZT1, *ZQaZT2, *ZQaZT12, *ZQaZT21;
	double *b_record, *QQ, *Za_unfix, *F_join;

	if (n<=0||m<=0||fixed<=0||fixed>n) return;

	unfixed = n - fixed;
	/* all fixed: no conditional update; keep lambda layout (column-major n x m) */
	if (unfixed==0) {
		memcpy(afix, afixhat, sizeof(double)*n*m);
		return;
	}
	record_Z_row = mat(n, n);
	Qa = mat(n, n);
	ZTQ = mat(n, n);
	Qzhat_record = mat(n, n);
	Za_record = mat(n, 1);
	F_join = mat(n, m);
	ZQaZT1 = mat(unfixed, unfixed);
	ZQaZT2 = mat(fixed, fixed);
	ZQaZT12 = mat(unfixed, fixed);
	ZQaZT21 = mat(fixed, unfixed);
	b_record = mat(fixed, m);
	QQ = mat(unfixed, fixed);
	Za_unfix = mat(unfixed, m);
	memset(record_Z_row, 0, n*n * sizeof(double));
	//先排无法固定的模糊度，再排已经固定的模糊度
	for (i = 0; i < unfixed; i++)
	{
		if (delete_record[i]<0||delete_record[i]>=n) {
			free(record_Z_row); free(Qa); free(ZTQ); free(Qzhat_record); free(Za_record);
			free(ZQaZT1); free(ZQaZT2); free(ZQaZT12); free(ZQaZT21);
			free(b_record); free(QQ); free(Za_unfix); free(F_join);
			return;
		}
		// unfixed ambiguities
		record_Z_row[delete_record[i] * n + inc] = 1.0;
		inc++;
	}
	for (i = 0; i < n; i++)
	{
		// fixed ambiguities
		if (findi(delete_record, unfixed, i) < 0)
		{
			if (inc>=n) break;
			record_Z_row[n*i + inc] = 1.0;
			inc++;
		}
	}
	//ZQbZT=Zp*Qb*Zp'; 协方差整数变换后，重新排列未固定模糊度和固定模糊度的等到新的协方差ZQZT
	memcpy(Qa, Qa_, n*n*sizeof(double));
	matmul("NN", n, n, n, 1.0, record_Z_row, Qa, 0.0, ZTQ);
	matmul("NT", n, n, n, 1.0, ZTQ, record_Z_row, 0.0, Qzhat_record);
	//将浮点解分为可固定的和不可固定两部分
	matmul("NN", n, 1, n, 1, record_Z_row, fa, 0, Za_record);
	//《A partial Carrier phase  integer ambiguity fixing algorithm for combinatiorial optimization between network rtk
	// reference  stations》公式20 
	mat_split(Qzhat_record, ZQaZT1, ZQaZT2, ZQaZT12, ZQaZT21, n, unfixed);
	if (!matinv(ZQaZT2, fixed)) {
		for (i = 0; i < m; i++) {
			for (j = 0; j < unfixed; j++) {
				Za_unfix[i*unfixed + j] = Za_record[j];
			}
		}
		for (i = 0; i < m; i++) {
			for (j = 0; j < fixed; j++) {
				b_record[i*(fixed)+j] = Za_record[j + unfixed] - afixhat[i*fixed + j];
			}
			matmul("NN", unfixed, fixed, fixed, 1.0, ZQaZT21, ZQaZT2, 0.0, QQ);
			matmul("NN", unfixed, 1, fixed, -1.0, QQ, b_record + i * fixed, 1.0, Za_unfix + i * unfixed);
		}
		for (i = 0; i < m; i++) {
			for (j = 0; j < unfixed; j++) {
				F_join[i*n + j] = Za_unfix[i*unfixed + j];
			}
		}
		for (i = 0; i < m; i++) {
			for (j = 0; j < fixed; j++) {
				F_join[i*n + j + unfixed] = afixhat[i*fixed + j];
			}
		}
		matmul("TN", n, m, n, 1, record_Z_row, F_join, 0, afix);
	}
	free(record_Z_row); free(Qa); free(ZTQ); free(Qzhat_record); free(Za_record);
	free(ZQaZT1); free(ZQaZT2); free(ZQaZT12); free(ZQaZT21);
	free(b_record); free(QQ); free(Za_unfix); free(F_join);
}

/* split float ambiguities into integer + fractional parts */
static void mat_int_frac(const double *full, int *int_part, double *frac_part, int n)
{
	int i;
	for (i=0; i<n; i++) {
		int_part[i] = ROUND(full[i]);
		frac_part[i] = full[i]-int_part[i];
	}
}
/* restore full ambiguities: full = int_part + frac_part */
static void mat_add(const int *int_part, const double *frac_part, double *full, int n)
{
	int i;
	for (i=0; i<n; i++) {
		full[i] = int_part[i] + frac_part[i];
	}
}

extern int parlambda(rtk_t *rtk,int n, int m, const double *a, const double *Q, double *F,
	double *s,int *sat_record)
{
	int del_num = 0;
	double ratio = 0;
	double ratio_mask = 0;
	//注：锁定卫星个数表示持续多个历元共用一些固定的卫星解算的卫星个数
	//    锁定卫星跟踪数表示每颗卫星有效的历元个数，rtk->ssat[i].lock[f]
	int i, j,f, k, valid_fixsat = 0; int Qahat_index = -1;
	int info = -1, area = 0, *sat_lock, *delete_record,*record_sat,*ambslip, *pafixvalid,*orig_idx;
	int *int_part;
	double *sat_azi,*sat_ele,*Qahat,*ahat,*afixhat,*frac_part;
	sat_lock=imat(n,1);delete_record=imat(n, 1);
	record_sat=imat(n, 1);ambslip=imat(n,1);pafixvalid=imat(n, 1);
	orig_idx=imat(n,1);int_part=imat(n,1);
	Qahat = mat(n, n);ahat = mat(n, 1);sat_azi = mat(n, 1);sat_ele = mat(n, 1);
	afixhat = mat(n,m);frac_part = mat(n, 1);

	if(!checkeig(n,Q))return -1; //协方差矩阵不可逆

	/* search on fractional parts only; restore integers after success */
	mat_int_frac(a, int_part, frac_part, n);
	//如果已经锁定卫星个数,且锁定个数超过一半以上的可视卫星，就直接进行部分模糊度解算
	memset(sat_azi, 0, n*1*sizeof(double));
	memset(sat_ele, 0, n*1*sizeof(double));
	memset(ambslip, 0, n*sizeof(int));
	memset(pafixvalid, 0, n*1*sizeof(int));
	
	if (rtk->nfixsat < rtk->opt.minfix)trace(3, "lack of fix sat:%d\n", rtk->nfixsat);
	
	ratio_mask = FindMinRatio(0.95, n);

 	if (contain_fixsat(rtk,n)&&(rtk->nfixsat >= rtk->opt.minfix)) {
		//部分模糊度
		del_num = 0;
		for (i = 0; i < n;i++) {
			if (!del_num) {
				memcpy(ahat, a, n*sizeof(double));
				memcpy(Qahat, Q, n*n*sizeof(double));
				/* record_sat is n x 1; do not copy full MAXSAT ddsat2 */
				memcpy(record_sat, rtk->ddsat2, n*sizeof(int));
			}
			if (findi(rtk->fixsat,rtk->nfixsat,rtk->ddsat2[i])<0) {
				/* find remaining copy (multi-freq keeps one entry per DD) */
				Qahat_index = findi(record_sat,n - del_num,rtk->ddsat2[i]);
				if (Qahat_index<0) continue; /* avoid remove_*(-1) heap corruption */
				delete_record[del_num++] = i;
				remove_sat(record_sat, n - del_num + 1, Qahat_index);
				remove_vec(ahat, n - del_num + 1, Qahat_index);
				remove_mat(Qahat, n - del_num + 1, Qahat_index);
			}
		}
		info = lambda(n - del_num, 2, ahat, Qahat, afixhat, s);
		par2fullamb(a, Q, afixhat, n - del_num, n, 2, delete_record, F);
		//end部分模糊度解算
		if (s[0] <= 0.0 || s[1] / s[0] >= ratio_mask) valid_fixsat = 1;
		else valid_fixsat = 0;
		trace(3, "contain fixsat,full sat=%d,fix sat=%d,s1/s0=%f\n", n, rtk->nfixsat, s[1] / s[0]);
	}
	else valid_fixsat = 0;
	//如果锁定卫星无法解算模糊度时，直接重新锁定卫星个数,并循环到模糊度搜索成功固定
	if (!valid_fixsat) {
		del_num = 0;
		for (i=0;i<n;i++) orig_idx[i]=i;
		for (k = 0; k <= n - PAR_MIN_SAT; k++) {//(3+rtk->opt.nvsysum)
			if (k == 0) {//全集模糊度
				info=lambda(n, 2, frac_part, Q, F, s);
				rtk->nfixsat = n;
				for (i = 0; i < rtk->nfixsat; i++)rtk->fixsat[i] = rtk->ddsat2[i];
				update_nfixsat_s(rtk);
				/* clear ambslip; compare each DD vs same-freq previous fix only */
				for (f = 0; f < NFREQ; f++)for (i = 0; i < MAXSAT; i++) rtk->ssat[i].ambslip[f] = 0;
				for (i = 0; i < n; i++) {
					f=rtk->ddfreq[i];
					if (f<0||f>=NFREQ) continue;
					if (rtk->ssat[rtk->ddsat2[i] - 1].pafixvalid[f]&&
						fabs(F[i]+int_part[i] - rtk->ssat[rtk->ddsat2[i] - 1].pafix[f]) > 0.99)
						rtk->ssat[rtk->ddsat2[i] - 1].ambslip[f] = 1;
				}
				trace(2,"PAR ambslip set (n=%d):",n);
				for (i=0;i<n;i++) {
					f=rtk->ddfreq[i];
					trace(2," sat%d/L%d[%d]",rtk->ddsat2[i],f+1,
					      (f>=0&&f<NFREQ)?rtk->ssat[rtk->ddsat2[i]-1].ambslip[f]:0);
				}
				trace(2,"\n");
			} else {//部分模糊度
				if (k == 1) {//第一次需要拷贝全局变量
					memcpy(ahat, frac_part, n*sizeof(double));
					memcpy(Qahat, Q, n*n*sizeof(double));
					for (i = 0; i < n; i++) {
						orig_idx[i]=i;
						f=rtk->ddfreq[i];
						sat_azi[i] = rtk->ssat[rtk->ddsat2[i] - 1].azel[0];
						sat_ele[i] = rtk->ssat[rtk->ddsat2[i] - 1].azel[1];
						/* one flag per DD entry (this entry's frequency only) */
						ambslip[i]=(f>=0&&f<NFREQ)?rtk->ssat[rtk->ddsat2[i]-1].ambslip[f]:0;
						pafixvalid[i]=(f>=0&&f<NFREQ)?rtk->ssat[rtk->ddsat2[i]-1].pafixvalid[f]:0;
					}
					trace(2,"PAR copy arrays (n=%d):\n",n);
					for (i=0;i<n;i++) {
						trace(2,"  i=%2d sat=%3d L%d orig=%2d pafix=%d ambslip=%d azi=%.1f ele=%.1f ahat=%.3f\n",
						      i,rtk->ddsat2[i],rtk->ddfreq[i]+1,orig_idx[i],pafixvalid[i],
						      ambslip[i],sat_azi[i]*R2D,sat_ele[i]*R2D,ahat[i]);
					}
				}
				//如果当前模糊度跟上一模糊度不一样，优先剔除该卫星进行部分模糊度解算
				// //否则剔除最低高度角法进行部分模糊度解算
				trace(2,"PAR remain k=%d ns=%d area=%d:",k,n-k+1,area);
				for (i=0;i<n-k+1;i++) {
					trace(2," [%d]sat=%d pafix=%d slip=%d azi=%.0f ele=%.0f",
					      i,rtk->fixsat[i],pafixvalid[i],ambslip[i],
					      sat_azi[i]*R2D,sat_ele[i]*R2D);
				}
				trace(2,"\n");
				Qahat_index = -1;
				if(Qahat_index<0)Qahat_index = find_ambslip_sat(ambslip, n - k + 1);
				trace(2,"PAR after ambslip: idx=%d\n",Qahat_index);
				if (Qahat_index < 0)Qahat_index = find_satadd(pafixvalid, n - k + 1);
				trace(2,"PAR after satadd: idx=%d\n",Qahat_index);
				if(Qahat_index<0)Qahat_index = find_min_ele(sat_azi,sat_ele, n - k + 1,area);
				trace(2,"PAR after min_ele: idx=%d\n",Qahat_index);
				area=(area+1) % 4;//将卫星按方位角分为4个区域
				if (Qahat_index<0||Qahat_index>=n-k+1) {
					trace(2,"PAR select fail idx=%d ns=%d, break\n",Qahat_index,n-k+1);
					break;
				}
				/* keep original amb index (not satno); multi-freq satno would collide */
				trace(2,"PAR drop idx=%d sat=%d orig_idx=%d del_num=%d\n",
				      Qahat_index,rtk->fixsat[Qahat_index],orig_idx[Qahat_index],del_num);
				delete_record[del_num++] = orig_idx[Qahat_index];
				rtk->nfixsat = n-k;
				remove_sat(orig_idx, n - k + 1, Qahat_index);
				remove_vec(sat_azi, n - k + 1, Qahat_index);
				remove_vec(sat_ele, n - k + 1, Qahat_index);
				remove_vec(ahat, n - k + 1, Qahat_index);
				remove_sat(ambslip, n - k + 1, Qahat_index);
				remove_sat(pafixvalid, n - k + 1, Qahat_index);
				remove_sat(rtk->fixsat, rtk->nfixsat + 1, Qahat_index);
				remove_mat(Qahat, n - k + 1, Qahat_index);
				update_nfixsat_s(rtk);
				info = lambda(n-k, 2, ahat, Qahat, afixhat, s);
				par2fullamb(frac_part,Q, afixhat, n - k, n,2,delete_record,F);
				//end部分模糊度解算
				ratio = s[1]/s[0];
			}
			if (s[0] <= 0.0 || s[1] / s[0] >= ratio_mask)break;
		}
		trace(3, "resize fixsat,full sat=%d,fix sat=%d,s1/s0=%f\n", n, rtk->nfixsat, s[1] / s[0]);

		/* 方差剔除法 (from PAR_LAMBDA): drop max decorrelated var, then lambda */		
		bool lackfix_flag = rtk->nfixsat < (PAR_MIN_SAT+rtk->opt.nnavsys);
		lackfix_flag = lackfix_flag && all_groups_lack_fix(rtk);
		if (lackfix_flag && (s[0] > 0.0 && s[1] / s[0] < ratio_mask))
		{
			double *Z,*ZTQ,*Qzhat,*Qvar;
			del_num = 0;
			Z=mat(n,n); ZTQ=mat(n,n); Qzhat=mat(n,n); Qvar=mat(n,1);
			for (k = 1; k <= n - PAR_MIN_SAT; k++) {
				if (k == 1) {//第一次需要拷贝全局变量
					memcpy(ahat, frac_part, n*sizeof(double));
					memcpy(Qahat, Q, n*n*sizeof(double));
					for (i=0;i<n;i++) orig_idx[i]=i;
					rtk->nfixsat = n;
					for (i = 0; i < rtk->nfixsat; i++)rtk->fixsat[i] = rtk->ddsat2[i];
					update_nfixsat_s(rtk);
				}
				/* Qzhat=Z'*Qahat*Z, remove max decorrelated variance */
				Qahat_index = find_max_qz_index(Qahat, n - k + 1, Z, ZTQ, Qzhat, Qvar);
				if (Qahat_index<0||Qahat_index>=n-k+1) break;
				delete_record[del_num++] = orig_idx[Qahat_index];
				rtk->nfixsat = n - k;
				remove_sat(orig_idx, n - k + 1, Qahat_index);
				remove_sat(rtk->fixsat, rtk->nfixsat + 1, Qahat_index);
				remove_vec(ahat, n - k + 1, Qahat_index);
				remove_mat(Qahat, n - k + 1, Qahat_index);
				update_nfixsat_s(rtk);
				info = lambda(n - k, 2, ahat, Qahat, afixhat, s);
				par2fullamb(frac_part,Q, afixhat, n - k, n,2,delete_record, F);
				/* end部分模糊度解算 */
				if (s[0] <= 0.0 || s[1] / s[0] >= ratio_mask)break;
			}
			free(Z); free(ZTQ); free(Qzhat); free(Qvar);
			trace(3, "var PAR,full sat=%d,fix sat=%d,s1/s0=%f\n", n, rtk->nfixsat, s[1] / s[0]);
		}
		
		//重新用lock锁定卫星跟踪数来求取部分模糊度
		lackfix_flag = rtk->nfixsat < (PAR_MIN_SAT+rtk->opt.nnavsys);
		lackfix_flag = lackfix_flag && all_groups_lack_fix(rtk);
		if (lackfix_flag && (s[0] > 0.0 && s[1] / s[0] < ratio_mask) &&
			lock_par_ready(rtk, n))
		{
			del_num = 0;
			for (k = 1; k <= n - PAR_MIN_SAT; k++) {//(3+rtk->opt.nvsysum)
				if (k == 1) {//第一次需要拷贝全局变量
					memcpy(ahat, frac_part, n*sizeof(double));
					memcpy(Qahat, Q, n*n*sizeof(double));
					for (i=0;i<n;i++) {
						orig_idx[i]=i;
						/* one lock count per DD entry (this entry's frequency only) */
						f=rtk->ddfreq[i];
						sat_lock[i]=(f>=0&&f<NFREQ)?rtk->ssat[rtk->ddsat2[i]-1].lock[f]:0;
					}
					rtk->nfixsat = n;
					for (i = 0; i < rtk->nfixsat; i++)rtk->fixsat[i] = rtk->ddsat2[i];
					update_nfixsat_s(rtk);
					trace(2,"PAR lock copy (n=%d):\n",n);
					for (i=0;i<n;i++) {
						f=rtk->ddfreq[i];
						trace(2,"  i=%2d sat=%3d L%d orig=%2d lock=%d (ssat.lock=%d/%d/%d)\n",
						      i,rtk->ddsat2[i],f+1,orig_idx[i],sat_lock[i],
						      rtk->ssat[rtk->ddsat2[i]-1].lock[0],
						      rtk->ssat[rtk->ddsat2[i]-1].lock[1],
						      rtk->ssat[rtk->ddsat2[i]-1].lock[2]);
					}
				}
				//部分模糊度：剔除本条目频点 lock 最小的 DD
				trace(2,"PAR lock remain k=%d ns=%d\n",k,n-k+1);
				for (i=0;i<n-k+1;i++) {
					trace(2,"  [%d]sat=%d lock=%d\n",i,rtk->fixsat[i],sat_lock[i]);
				}
				Qahat_index = find_min_i(sat_lock, n - k + 1);
				trace(2,"PAR after lock_min: idx=%d sat=%d lock=%d\n",
				      Qahat_index,
				      (Qahat_index>=0&&Qahat_index<n-k+1)?rtk->fixsat[Qahat_index]:-1,
				      (Qahat_index>=0&&Qahat_index<n-k+1)?sat_lock[Qahat_index]:-1);
				if (Qahat_index<0||Qahat_index>=n-k+1) break;
				delete_record[del_num++] = orig_idx[Qahat_index];			// store all of the deleted ambiguitiy indexs
				rtk->nfixsat = n - k;
				remove_sat(orig_idx, n - k + 1, Qahat_index);
				remove_sat(rtk->fixsat, rtk->nfixsat + 1, Qahat_index);
				remove_sat(sat_lock, n - k + 1, Qahat_index);
				remove_vec(ahat, n - k + 1, Qahat_index);
				remove_mat(Qahat, n - k + 1, Qahat_index);
				update_nfixsat_s(rtk);
				info = lambda(n - k, 2, ahat, Qahat, afixhat, s);
				par2fullamb(frac_part,Q, afixhat, n - k, n,2,delete_record, F);
				//end部分模糊度解算
				if (s[0] <= 0.0 || s[1] / s[0] >= ratio_mask)break;
			}
			trace(3, "lock PAR,full sat=%d,fix sat=%d,s1/s0=%f\n", n, rtk->nfixsat, s[1] / s[0]);
		}
		else if (lackfix_flag && (s[0] > 0.0 && s[1] / s[0] < ratio_mask)) {
			trace(2,"PAR lock skip: lock not ready (need mature+spread)\n");
		}
		
	}

	/* restore full ambiguities: F = int_part + fractional integers */
	for (i = 0; i < m; i++) {
		mat_add(int_part, F+i*n, F+i*n, n);
	}
	
	/* record previous fix per DD frequency (do not broadcast F[i] to all freqs) */
	update_nfixsat_s(rtk);
	for (i = 0; i < MAXSAT; i++) for (f = 0; f < NFREQ; f++) rtk->ssat[i].pafixvalid[f] = 0;
	if (s[0] <= 0.0 || s[1] / s[0] >= ratio_mask) {
		for (i = 0; i < n; i++) {
			f=rtk->ddfreq[i];
			if (f<0||f>=NFREQ) continue;
			rtk->ssat[rtk->ddsat2[i] - 1].pafix[f] = F[i];
			rtk->ssat[rtk->ddsat2[i] - 1].pafixvalid[f] = 1;
		}
	}

	//for (i = 0; i < 6; i++)if (rtk->ref_sat[i] != 0)printf("%d ", rtk->ref_sat[i]);
	//printf("%d %d %f\n", n, rtk->nfixsat,qs[1]/qs[0]);
	//if(rtk->nfixsat<n)for (i = 0; i < rtk->nfixsat; i++)printf("%d %.1f,",rtk->fixsat[i],afixhat[i]);
	//else for (i = 0; i < rtk->nfixsat; i++)printf("%d %.1f,", rtk->fixsat[i], afix[i]);
	//for (f = 0; f < NFREQ; f++)for (i = 0; i < n; i++) printf("%d %.1f,", rtk->ddsat2[i], rtk->ssat[rtk->ddsat2[i] - 1].pafix[f]);
	//printf("%f\n",qs[1]/qs[0]);
    //for (i = 0; i < n; i++)printf("%d %.1f,", rtk->ddsat2[i], afix[i]);
    //printf("\n");
	trace(3,"parlambda,full sat=%d,par sat=%d\n", n, rtk->nfixsat);

	free(Qahat);free(ahat);free(sat_lock);free(delete_record);free(record_sat);free(ambslip);
	free(pafixvalid);free(sat_azi);free(sat_ele);free(afixhat);free(orig_idx);
	free(int_part);free(frac_part);
	return info;
}

/* ===== resamb interface (same as PAR_resamb_LAMBDA in PAR_LAMBDA.cc) ===== */

#ifndef THRES_HOLDAMB
#define THRES_HOLDAMB   3.0
#endif
#define NF_(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP_(opt)     ((opt)->dynamics==0?3:9)
#define NI_(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT_(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL_(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB_(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF_(opt))
#define NR_(opt)     (NP_(opt)+NI_(opt)+NT_(opt)+NL_(opt))
#define IB_(s,f,opt) (NR_(opt)+MAXSAT*(f)+(s)-1)

/* single to double-difference transformation matrix (D') --------------------*/
static int par_ddmat(rtk_t *rtk, double *D, ddsat_t *ddsat, const int *vflg,
					 int nv, int flag)
{
	int i,j,k,m,f,tc,sat1,sat2;
	int nb=0,nx=rtk->nx,na=rtk->na;
	static unsigned char fixflag[NFREQ]={0};
	double *x=rtk->x;

	trace(3,"par_ddmat :\n");

	tc=rtk->opt.mode==PMODE_INS_TGNSS;

	if (tc) {nx=rtk->ins.nx; na=rtk->ins.nb;}
	if (tc) x=rtk->ins.x;

	for (i=0;i<MAXSAT;i++) {
		memcpy(rtk->ssat[i].fix,fixflag,sizeof(unsigned char)*NFREQ);
	}
	for (i=0;i<na;i++) D[i+i*nx]=1.0;

	for (i=0;i<nv;i++) {
		if (((vflg[i]>>4)&0xF)==1) continue;

		sat1=(vflg[i]>>16)&0xFF;
		sat2=(vflg[i]>> 8)&0xFF;
		f=vflg[i]&0xF;

		m=par_sysind(rtk->ssat[sat2-1].sys);

		if ((m==1&&rtk->opt.glomodear==0)||
			(m==3&&rtk->opt.bdsmodear==0)) {
			continue;
		}
		tc?k=xiBs(&rtk->opt.insopt,sat1,f):k=IB_(sat1,f,&rtk->opt);
		tc?j=xiBs(&rtk->opt.insopt,sat2,f):j=IB_(sat2,f,&rtk->opt);

		rtk->ssat[sat1-1].fix[f]=2;

		if (x[j]==0.0||x[k]==0.0) continue;
		if (rtk->ssat[sat2-1].slip[f]) continue;

		if ((rtk->ssat[sat2-1].lock[f]>=0||flag)&&
			rtk->ssat[sat2-1].vsat[f]&&
			rtk->ssat[sat2-1].azel[1]>=rtk->opt.elmaskar) {

			D[k+(na+nb)*nx]= 1.0;
			D[j+(na+nb)*nx]=-1.0;

			ddsat[nb].f=f;
			ddsat[nb].sat1=sat1;
			ddsat[nb].sat2=sat2;

			rtk->ssat[sat2-1].fix[f]=2;
			nb++;
		}
		else rtk->ssat[sat2-1].fix[f]=1;
	}
	return nb;
}
/* restore single-differenced ambiguity --------------------------------------*/
static void par_restamb(rtk_t *rtk, const double *bias, ddsat_t *ddsat,
						int nb, double *xa)
{
	int i,f,sat1,sat2,ib1,ib2,tc,nx,na;

	trace(3,"par_restamb :nb=%d\n",nb);

	tc=rtk->opt.mode==PMODE_INS_TGNSS;

	tc?nx=rtk->ins.nx:nx=rtk->nx;
	tc?na=rtk->ins.nb:na=rtk->na;

	for (i=0;i<nx;i++) tc?xa[i]=rtk->ins.x[i]:xa[i]=rtk->x[i];
	for (i=0;i<na;i++) tc?xa[i]=rtk->ins.xb[i]:xa[i]=rtk->xa[i];

	for (i=0;i<nb;i++) {
		sat1=ddsat[i].sat1;
		sat2=ddsat[i].sat2;
		f=ddsat[i].f;

		ib1=tc?xiBs(&rtk->opt.insopt,sat1,f):IB_(sat1,f,&rtk->opt);
		ib2=tc?xiBs(&rtk->opt.insopt,sat2,f):IB_(sat2,f,&rtk->opt);

		xa[ib2]=xa[ib1]-bias[i];
	}
}
/* resolve integer ambiguity by parlambda (internal) -------------------------*/
static int resamb_parlambda_i(rtk_t *rtk, const nav_t *nav, const obsd_t *obs,
							  double *bias, double *xa,
							  ddsat_t *ddsat, int *namb,
							  ddsat_t *wlsat, int *nw,
							  double *dxwl, double *Qywl, double *wlb,
							  const int *vflg, int nv, int fixnolock,
							  int *nny, int *fixwl)
{
	prcopt_t *opt=&rtk->opt;
	insstate_t *ins=&rtk->ins;
	insopt_t *insopt=&rtk->opt.insopt;
	int i,j,ny,nb,info,nx=rtk->nx,na=rtk->na,tc,*sat_record;
	double *D,*DP,*y,*Qy,*b,*db,*Qb,*Qab,*QQ,s[2]={0};
	double *x,*P,*Pa,*xb;

	(void)nav; (void)obs; (void)wlsat; (void)nw;
	(void)dxwl; (void)Qywl; (void)wlb; (void)nny; (void)fixwl;

	trace(3,"resamb_parlambda_i : nx=%d\n",nx);

	if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||
		rtk->opt.thresar[0]<1.0) {
		trace(2,"resamb_parlambda_i early-return1: mode=%d(<=DGPS?%d) modear=%d(OFF?%d) thresar0=%.3f(<1?%d) fixnolock=%d nv=%d\n",
			  rtk->opt.mode,rtk->opt.mode<=PMODE_DGPS,
			  rtk->opt.modear,rtk->opt.modear==ARMODE_OFF,
			  rtk->opt.thresar[0],rtk->opt.thresar[0]<1.0,
			  fixnolock,nv);
		if (namb) *namb=0;
		return 0;
	}
	if (rtk->opt.mode==PMODE_INS_TGNSS&&rtk->opt.insopt.tc<=INSTC_DGPS) {
		trace(2,"resamb_parlambda_i early-return2: mode=%d tc=%d(<=DGPS?%d INSTC_DGPS=%d INSTC_RTK=%d) fixnolock=%d nv=%d\n",
			  rtk->opt.mode,rtk->opt.insopt.tc,
			  rtk->opt.insopt.tc<=INSTC_DGPS,INSTC_DGPS,INSTC_RTK,
			  fixnolock,nv);
		if (namb) *namb=0;
		return 0;
	}
	rtk->sol.ratio=0.0;
	if (namb) *namb=0;
	rtk->bias.inherit=0;

	tc=rtk->opt.mode==PMODE_INS_TGNSS;

	tc?x=rtk->ins.x:x=rtk->x;
	tc?P=rtk->ins.P:P=rtk->P;
	tc?Pa=rtk->ins.Pb:Pa=rtk->Pa;
	tc?xb=rtk->ins.xb:xb=rtk->xa;

	if (opt->mode==PMODE_INS_TGNSS) {
		nx=ins->nx;
		na=ins->nb;
	}
	D=zeros(nx,nx);
	if ((nb=par_ddmat(rtk,D,ddsat,vflg,nv,fixnolock))<=0) {
		trace(2,"no valid double-difference\n");
		free(D);
		return 0;
	}
	ny=na+nb;
	y=mat(ny,1); Qy=mat(ny,ny); DP=mat(ny,nx);
	b=mat(nb,2); db=mat(nb,1);
	Qb=mat(nb,nb); Qab=mat(na,nb); QQ=mat(na,nb);
	sat_record=imat(nb,1);
	memset(rtk->ddsat2,0,sizeof(rtk->ddsat2));
	memset(rtk->ddfreq,0,sizeof(rtk->ddfreq));
	for (i=0;i<nb;i++) {
		/* parlambda uses 1-based sat no in ddsat2; keep sat_record for API */
		sat_record[i]=ddsat[i].sat2;
		rtk->ddsat2[i]=ddsat[i].sat2;
		rtk->ddfreq[i]=ddsat[i].f;
	}

	matmul("TN",ny, 1,nx,1.0,D ,x,0.0,y );
	matmul("TN",ny,nx,nx,1.0,D ,P,0.0,DP);
	matmul("NN",ny,ny,nx,1.0,DP,D,0.0,Qy);

	for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb[i+j*nb]=Qy[na+i+(na+j)*ny];
	for (i=0;i<na;i++) for (j=0;j<nb;j++) Qab[i+j*na]=Qy[i+(na+j)*ny];

	trace(4,"N(0)=");
	tracemat(4,y+na,1,nb,10,3);

	if (!(info=parlambda(rtk,nb,2,y+na,Qb,b,s,sat_record))) {

		trace(4,"N(1)="); tracemat(4,b   ,1,nb,10,3);
		trace(4,"N(2)="); tracemat(4,b+nb,1,nb,10,3);

		rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
		if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;

		double ratio_mask = FindMinRatio(0.95, nb);
		if (s[0]>0.0&&s[1]/s[0]>=ratio_mask) {

			for (i=0;i<na;i++) {
				xb[i]=x[i];
				for (j=0;j<na;j++) Pa[i+j*na]=P[i+j*nx];
			}
			if (tc) {
				for (i=0;i<xnCl(insopt);i++) xb[i]=1E-20;
			}
			for (i=0;i<nb;i++) {
				bias[i]=b[i];
				y[na+i]-=b[i];
			}
			if (!matinv(Qb,nb)) {
				matmul("NN",nb,1,nb, 1.0,Qb ,y+na,0.0,db);
				matmul("NN",na,1,nb,-1.0,Qab,db  ,1.0,xb);

				matmul("NN",na,nb,nb, 1.0,Qab,Qb ,0.0,QQ);
				matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,Pa);

				trace(3,"validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
					  nb,s[1]/s[0],s[0],s[1]);

				par_restamb(rtk,bias,ddsat,nb,xa);
			}
			else nb=0;
		}
		else {
			trace(2,"ambiguity validation fail (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
				  nb,s[0]>0.0?s[1]/s[0]:0.0,s[0],s[1]);
			nb=0;
		}
	}
	else {
		trace(2,"parlambda error (info=%d)\n",info);
		nb=0;
	}
	free(D); free(y); free(Qy); free(DP);
	free(b); free(db); free(Qb); free(Qab); free(QQ);
	free(sat_record);
	if (namb) *namb=nb;
	return nb;
}
/* resolve integer ambiguity by parlambda ------------------------------------
* same interface as PAR_resamb_LAMBDA / resamb_LAMBDA2
*----------------------------------------------------------------------------*/
extern int resamb_parlambda(rtk_t *rtk, const nav_t *nav, const obsd_t *obs,
							double *bias, double *xa,
							ddsat_t *ddsat, int *namb,
							ddsat_t *wlsat, int *nw,
							double *dxwl, double *Qywl, double *wlb,
							const int *vflg, int nv,
							int *nny, int *fixwl)
{
	int nb;

	trace(3,"resamb_parlambda:\n");
	/* fixnolock=0: par_ddmat requires lock>=0, so new sats wait minlock epochs */
	nb=resamb_parlambda_i(rtk,nav,obs,bias,xa,ddsat,namb,wlsat,nw,dxwl,Qywl,wlb,
						 vflg,nv,0,nny,fixwl);
	return nb;
}


static double FS[7][60] = {
	{ 405284.07,999.00,141.11,53.44,29.75,20.03,15.02,12.05,10.11,8.75,
		7.76,7.00,6.41,5.93,5.54,5.20,4.92,4.68,4.47,4.29,
		4.13,3.98,3.85,3.74,3.63,3.53,3.44,3.36,3.29,3.22,
		3.15,3.09,3.04,2.98,2.93,2.89,2.84,2.80,2.76,2.73,
		2.69,2.66,2.63,2.60,2.57,2.54,2.51,2.49,2.46,2.44,
		2.42,2.40,2.38,2.36,2.34,2.32,2.30,2.28,2.27,2.25 },

	{ 16210.72, 199.00, 47.47, 23.15, 14.94, 11.07, 8.89, 7.50, 6.54, 5.85,
		5.32, 4.91, 4.57, 4.30, 4.07, 3.87, 3.71, 3.56, 3.43, 3.32,
		3.22, 3.12, 3.04, 2.97, 2.90, 2.84, 2.78, 2.72, 2.67, 2.63,
		2.58, 2.54, 2.51, 2.47, 2.44, 2.41, 2.38, 2.35, 2.32, 2.30,
		2.27, 2.25, 2.23, 2.21, 2.19, 2.17, 2.15, 2.13, 2.11, 2.10,
		2.08, 2.07, 2.05, 2.04, 2.02, 2.01, 2.00, 1.99, 1.97, 1.96},

	{ 4052.18, 99.00, 29.46, 15.98, 10.97, 8.47, 6.99, 6.03, 5.35, 4.85,
		4.46, 4.16, 3.91, 3.70, 3.52, 3.37, 3.24, 3.13, 3.03, 2.94,
		2.86, 2.78, 2.72, 2.66, 2.60, 2.55, 2.51, 2.46, 2.42, 2.39,
		2.35, 2.32, 2.29, 2.26, 2.23, 2.21, 2.18, 2.16, 2.14, 2.11,
		2.09, 2.08, 2.06, 2.04, 2.02, 2.01, 1.99, 1.98, 1.96, 1.95,
		1.94, 1.92, 1.91, 1.90, 1.89, 1.88, 1.87, 1.86, 1.85, 1.84},

	{ 647.79, 39.00, 15.44, 9.60, 7.15, 5.82, 4.99, 4.43, 4.03, 3.72,
		3.48, 3.29, 3.12, 2.98, 2.87, 2.76, 2.68, 2.60, 2.53, 2.47,
		2.41, 2.36, 2.31, 2.27, 2.23, 2.20, 2.16, 2.13, 2.10, 2.07,
		2.05, 2.03, 2.00, 1.98, 1.96, 1.94, 1.92, 1.91, 1.89, 1.88,
		1.86, 1.85, 1.83, 1.82, 1.81, 1.80, 1.78, 1.77, 1.76, 1.75,
		1.74, 1.73, 1.72, 1.72, 1.71, 1.70, 1.69, 1.68, 1.67, 1.67},

	{ 161.45, 19.00, 9.28, 6.39, 5.05, 4.28, 3.79, 3.44, 3.18, 2.98,
		2.82, 2.69, 2.58, 2.48, 2.40, 2.33, 2.27, 2.22, 2.17, 2.12,
		2.08, 2.05, 2.01, 1.98, 1.96, 1.93, 1.90, 1.88, 1.86, 1.84,
		1.82, 1.80, 1.79, 1.77, 1.76, 1.74, 1.73, 1.72, 1.70, 1.69,
		1.68, 1.67, 1.66, 1.65, 1.64, 1.63, 1.62, 1.62, 1.61, 1.60,
		1.59, 1.58, 1.58, 1.57, 1.56, 1.56, 1.55, 1.55, 1.54, 1.53},

	{ 39.86, 9.00, 5.39, 4.11, 3.45, 3.05, 2.78, 2.59, 2.44, 2.32,
		2.23, 2.15, 2.08, 2.02, 1.97, 1.93, 1.89, 1.85, 1.82, 1.79,
		1.77, 1.74, 1.72, 1.70, 1.68, 1.67, 1.65, 1.63, 1.62, 1.61,
		1.59, 1.58, 1.57, 1.56, 1.55, 1.54, 1.53, 1.52, 1.51, 1.51,
		1.50, 1.49, 1.48, 1.48, 1.47, 1.46, 1.46, 1.45, 1.45, 1.44,
		1.44, 1.43, 1.43, 1.42, 1.42, 1.41, 1.41, 1.40, 1.40, 1.40},

	{ 9.47, 4.00, 2.94, 2.48, 2.23, 2.06, 1.94, 1.86, 1.79, 1.73,
		1.69, 1.65, 1.61, 1.58, 1.56, 1.54, 1.52, 1.50, 1.48, 1.47,
		1.45, 1.44, 1.43, 1.42, 1.41, 1.40, 1.39, 1.38, 1.37, 1.36,
		1.36, 1.35, 1.34, 1.34, 1.33, 1.33, 1.32, 1.32, 1.31, 1.31,
		1.30, 1.30, 1.30, 1.29, 1.29, 1.28, 1.28, 1.28, 1.27, 1.27,
		1.27, 1.26, 1.26, 1.26, 1.26, 1.25, 1.25, 1.25, 1.25, 1.24}

};


double FindMinRatio(double cl, int npar)
{
	int n = npar;

	if (n == 0)
	{
		n = 1;
	}
	else if (n < 60)
	{
		n = npar;
	}
	else
	{
		n = 60;
	}

	if (cl == 1.0)
	{
		return 3.0;
	}
	else if (cl >= 0.999)
	{
		return FS[0][n - 1];
	}
	else if (cl >= 0.995)
	{
		return FS[1][n - 1];
	}
	else if (cl >= 0.99)
	{
		return FS[2][n - 1];
	}
	else if (cl >= 0.975)
	{
		return FS[3][n - 1];
	}
	else if (cl >= 0.95)
	{
		return FS[4][n - 1];
	}
	else if (cl >= 0.9)
	{
		return FS[5][n - 1];
	}
	else
	{
		return FS[6][n - 1];
	}
}

/* Householder reduction of symmetric Q to tridiagonal (column-major n x n)
 * q: n x n work (Householder / eigenvectors), b: diagonal, c: subdiagonal */
static void eastrq(int n, const double *Q, double *q, double *b, double *c)
{
	int i,j,k;
	double h,f,g,h2,c0;

	for (i=0;i<n;i++) for (j=0;j<n;j++) q[i+j*n]=Q[i+j*n];

	for (i=n-1;i>=1;i--) {
		h=0.0;
		if (i>1) for (k=0;k<=i-1;k++) h+=q[i+k*n]*q[i+k*n];

		if (h+1.0==1.0) {
			c[i-1]=0.0;
			if (i==1) c[i-1]=q[i+(i-1)*n];
			b[i]=0.0;
		}
		else {
			c[i-1]=sqrt(h);
			if (q[i+(i-1)*n]>0.0) c[i-1]=-c[i-1];
			h=h-q[i+(i-1)*n]*c[i-1];
			q[i+(i-1)*n]-=c[i-1];
			f=0.0;
			for (j=0;j<=i-1;j++) {
				q[j+i*n]=q[i+j*n]/h;
				g=0.0;
				for (k=0;k<=j;k++) g+=q[j+k*n]*q[i+k*n];
				if (j+1<=i-1)
					for (k=j+1;k<=i-1;k++) g+=q[k+j*n]*q[i+k*n];
				if (j==0) c0=g/h;
				else      c[j-1]=g/h;
				f+=g*q[j+i*n];
			}
			h2=f/(h+h);
			for (j=0;j<=i-1;j++) {
				f=q[i+j*n];
				if (j==0) { g=c0-h2*f; c0=g; }
				else      { g=c[j-1]-h2*f; c[j-1]=g; }
				for (k=0;k<=j;k++) {
					if (k==0) q[j+k*n]-=f*c0+g*q[i+k*n];
					else      q[j+k*n]-=f*c[k-1]+g*q[i+k*n];
				}
			}
			b[i]=h;
		}
	}
	b[0]=0.0;
	for (i=0;i<n;i++) {
		if (b[i]!=0.0&&i-1>=0) {
			for (j=0;j<=i-1;j++) {
				g=0.0;
				for (k=0;k<=i-1;k++) g+=q[i+k*n]*q[k+j*n];
				for (k=0;k<=i-1;k++) q[k+j*n]-=g*q[k+i*n];
			}
		}
		b[i]=q[i+i*n];
		q[i+i*n]=1.0;
		if (i-1>=0) for (j=0;j<=i-1;j++) {
			q[i+j*n]=0.0;
			q[j+i*n]=0.0;
		}
	}
}

/* QR eigenvalues of tridiagonal; eigenvalues in b[]; return 1=ok, 0=fail */
static int ebstq(int n, double *q, double *b, double *c, double eps, int lmax)
{
	int i,j,k,m,it;
	double d,f,h,g,p,r,e,s;

	c[n-1]=0.0;
	d=0.0; f=0.0;
	for (j=0;j<n;j++) {
		it=0;
		h=eps*(fabs(b[j])+fabs(c[j]));
		if (h>d) d=h;
		m=j;
		while (m<=n-1&&fabs(c[m])>d) m++;
		if (m!=j) {
			do {
				if (it==lmax) return 0;
				it++;
				g=b[j];
				p=(b[j+1]-g)/(2.0*c[j]);
				r=sqrt(p*p+1.0);
				if (p>=0.0) b[j]=c[j]/(p+r);
				else        b[j]=c[j]/(p-r);
				h=g-b[j];
				for (i=j+1;i<n;i++) b[i]-=h;
				f+=h;
				p=b[m]; e=1.0; s=0.0;
				for (i=m-1;i>=j;i--) {
					g=e*c[i]; h=e*p;
					if (fabs(p)>=fabs(c[i])) {
						e=c[i]/p; r=sqrt(e*e+1.0);
						c[i+1]=s*p*r; s=e/r; e=1.0/r;
					}
					else {
						e=p/c[i]; r=sqrt(e*e+1.0);
						c[i+1]=s*c[i]*r; s=1.0/r; e=e/r;
					}
					p=e*b[i]-s*g;
					b[i+1]=h+s*(e*g+s*b[i]);
					for (k=0;k<n;k++) {
						h=q[k+(i+1)*n];
						q[k+(i+1)*n]=s*q[k+i*n]+e*h;
						q[k+i*n]=e*q[k+i*n]-s*h;
					}
				}
				c[j]=s*p; b[j]=e*p;
			} while (fabs(c[j])>d);
		}
		b[j]+=f;
	}
	for (i=0;i<n;i++) {
		k=i; p=b[i];
		for (j=i+1;j<n;j++) if (b[j]<=p) { k=j; p=b[j]; }
		if (k!=i) {
			b[k]=b[i]; b[i]=p;
			for (j=0;j<n;j++) {
				p=q[j+i*n]; q[j+i*n]=q[j+k*n]; q[j+k*n]=p;
			}
		}
		}
		return 1;
}

/* check positive-definite of symmetric covariance Q (column-major n x n)
 * same layout as Qb/Qy/lambda: Q[i+j*n] = row i, col j
 * return: 1=PD, 0=not PD or eig failed */
extern int checkeig(int n, const double *Q)
{
	int i;
	double *q,*b,*c;

	if (n<=0||Q==NULL) return 0;

	q=mat(n,n); b=mat(n,1); c=mat(n,1);
	if (q==NULL||b==NULL||c==NULL) {
		free(q); free(b); free(c);
		return 0;
	}
	eastrq(n,Q,q,b,c);
	if (!ebstq(n,q,b,c,1e-12,100)) {
		free(q); free(b); free(c);
		return 0;
	}
	for (i=0;i<n;i++) {
		if (b[i]<0.0) {
			free(q); free(b); free(c);
			return 0;
		}
	}
	free(q); free(b); free(c);
	return 1;
}

//
//static int ratioinv(/*const double Pf_FIX,*/const float PfILS, const int N, double *mu)
//{
//	//double kPf=round(Pf_FIX*1000);
//	int n=N,i;
//	const int nrow=31;
//	const int ncol=64;
//	//matrix_t table;
//	float x[31],y[31];
//
//	*mu=0.0;
//	//table.nrow=31;
//	//table.ncol=64;
//	//table.data=tabledata;
//
//	//if(kPf==1) {
//	//	table.nrow=31;
//	//	table.ncol=64;
//	//	table.data=tabledata1;
//	//} else if(kPf==10) {
//	//	table.nrow=31;
//	//	table.ncol=41;
//	//	table.data=tabledata10;
//	//} else return 0;
//
//	if(n<1) return 0;
//	if(n>ncol-1) n=ncol-1;
//	for(i=0;i<nrow;i++) {
//		x[i]=tabledata[i*ncol];
//		y[i]=tabledata[i*ncol+n];
//	}
//	if(!interp1q(x,y,nrow,PfILS,mu)) return 0;
//	return 1;
//}
