/*------------------------------------------------------------------------------
* cn0_varerr.cc : C/N0 departure from the elevation template
*
* Yuan Li et al., Meas. Sci. Technol. 34 (2023) 015011, eqs (7)(9).
*   C/N0*(θ) = a1 + a2 θ + a3 θ^2 + a4 θ^3
*   STD*(θ)  = exp(b1 + b2 θ + b3 θ^2 + b4 θ^3)
*   dcn0     = |C/N0 - C/N0*|
*   ncn0     = |C/N0 - C/N0*| / STD*
* θ and the template are in dB-Hz, θ in degrees. b* are log-cubic
* coefficients, so the exponential keeps STD* in dB-Hz. The elevation
* fed to the variance is eqs (9)(12)(13):
*   |cn0-cn0*| <= ω STD*            keep θ
*   ω STD* < |cn0-cn0*| <= φ STD*  walk θ̄ by 1° until (9) holds
*   |cn0-cn0*| > φ STD*             θ̄ = 0 (delete)
* θ̄ is clamped to [1°, 90°]. ω, φ are cn0_w, cn0_phi per frequency.
* A group with a* or b* all zero has no template.
* m=0 GPS/QZS/SBS/GAL, m=1 BDS. GLONASS has no coefficients.
*-----------------------------------------------------------------------------*/
#include <math.h>
#include "carvig.h"

static double poly3(const double *p, double th)
{
    return p[0]+th*(p[1]+th*(p[2]+th*p[3]));
}
static int coeff_set(const double *p)
{
    return p[0]!=0.0||p[1]!=0.0||p[2]!=0.0||p[3]!=0.0;
}
/* C/N0* and STD* at elevation th (deg). 1 if this system/freq has both. --*/
static int template_at(const prcopt_t *opt, int sys, int freq, double th,
                      double *cn0s, double *stds)
{
    int m;
    const double *a,*b;

    if (sys&(SYS_GPS|SYS_QZS|SYS_SBS|SYS_GAL)) m=0;
    else if (sys==SYS_CMP) m=1;
    else return 0;
    if (!opt||freq<0||freq>=NFREQ) return 0;
    a=opt->cn0_a[m][freq];
    b=opt->cn0_b[m][freq];
    if (!coeff_set(a)||!coeff_set(b)) return 0;
    *cn0s=poly3(a,th);
    *stds=exp(poly3(b,th));
    return *stds>0.0;
}
/* dcn0 = |cn0-C/N0*|, ncn0 = dcn0/STD* at elevation el (rad).
 * cn0 is the measured C/N0 (dB-Hz). 1 if this system/freq has a template. -*/
extern int cn0_varerr(const prcopt_t *opt, int sys, int freq, double el,
                      double cn0, double *dcn0, double *ncn0)
{
    double th,cn0s,stds;

    if (dcn0) *dcn0=0.0;
    if (ncn0) *ncn0=0.0;
    if (!opt||freq<0||freq>=NFREQ||el<=0.0||cn0<=0.0) return 0;
    th=el*R2D;
    if (!template_at(opt,sys,freq,th,&cn0s,&stds)) return 0;
    if (dcn0) *dcn0=fabs(cn0-cn0s);
    if (ncn0) *ncn0=fabs(cn0-cn0s)/stds;
    return 1;
}
/* eq (12): step θ̄ by 1° toward the template until eq (9) holds at θ̄.
 * Measured C/N0 is fixed; C/N0* and STD* are recomputed at θ̄. ----------*/
static double virtual_el(const prcopt_t *opt, int sys, int freq, double cn0,
                        double th, double w)
{
    double thb=th,cn0s,stds;
    int k,dir=0,step;

    for (k=0;k<90;k++) {
        if (!template_at(opt,sys,freq,thb,&cn0s,&stds)) break;
        if (fabs(cn0-cn0s)<=w*stds) break;
        step=(cn0<=cn0s)?-1:1;
        if (dir&&step!=dir) break;          /* 1° step crossed the ω band */
        dir=step;
        thb+=step;
        if (thb<1.0) { thb=1.0; break; }    /* paper clamp γ=1°, δ=90° */
        if (thb>90.0) { thb=90.0; break; }
    }
    return thb*D2R;
}
/* elevation for the elevation variance (rad). eqs (9)(12)(13).
 * ω<=0 or no template leaves el unchanged. φ<=ω has no middle band. ----*/
extern double cn0_el(const prcopt_t *opt, int sys, int freq, double el,
                     double cn0)
{
    double dcn0,ncn0,w,phi;

    if (!opt||freq<0||freq>=NFREQ) return el;
    w=opt->cn0_w[freq];
    if (!(w>0.0)) return el;
    if (!cn0_varerr(opt,sys,freq,el,cn0,&dcn0,&ncn0)) return el;
    if (ncn0<=w) return el;
    phi=opt->cn0_phi[freq];
    if (!(phi>w)||ncn0>phi) return 0.0;
    return virtual_el(opt,sys,freq,cn0,el*R2D,w);
}
