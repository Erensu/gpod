/*------------------------------------------------------------------------------
 * podlib.h: PODLIB constants, types and function prototypes
 *
 * version : $Revision:$ $Date:$
 * Copyright(c) 2023-2025 by sujinglan, all rights reserved
 * history : 2024/10/17 1.0  new
 *-----------------------------------------------------------------------------*/
#ifndef PODLIB_H
#define PODLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#include "rtklib.h"
#include "sofa.h"
#include "log.h"

/* constants -----------------------------------------------------------------*/
#define GRAVITY_MODEL_GGM03C    0      /* earth's gravity field model: GGM03C */
#define GRAVITY_MODEL_EGM2008   1      /* earth's gravity field model: EGM2008 */

#define GRAVITY_MODEL_TIDE_ZERO 0      /* tide system of earth's gravity field model: tide-zero */
#define GRAVITY_MODEL_TIDE_FREE 1      /* tide system of earth's gravity field model: tide-free */

#define SHADOW_CYLINDRICAL      0      /* shadow model: cylindrical */
#define SHADOW_GEOMETRICAL      1      /* shadow model: geometrical */
#define SHADOW_CONICAL          2      /* shadow model: conical */

#define ATMODENSITY_HPM         0      /* atmospheric density model: Harris-Priester */
#define ATMODENSITY_NRL         1      /* atmospheric density model: nrlmsise-00 */

#define OCEANTIDECORR_DEFAULT   1      /* ocean tide corrections for normalized earth potential coefficients: default */
#define OCEANTIDECORR_OTMODEL   2      /* ocean tide corrections for normalized earth potential coefficients: ocean tide model */

#define ORBITINT_GSLODE_RK8PD   1      /* orbit numerical integration method: GSL-RK8PD */
#define ORBITINT_GSLODE_RKF45   2      /* orbit numerical integration method: GSL-RKF45 */
#define ORBITINT_GSLODE_MSADAMS 3      /* orbit numerical integration method: GSL-MSADAMS */
#define ORBITINT_GSLODE_MSBDF   4      /* orbit numerical integration method: GSL-MSBDF */
#define ORBITINT_RKF78M         5      /* orbit numerical integration method: RKF78 */
#define ORBITINT_ODEINTV2_ADAMSADAMS 6 /* orbit numerical integration method: odeint-v2-adaptive_adams_bashforth_moulton */
#define ORBITINT_ODEINTV2_CONMSADAMS 7 /* orbit numerical integration method: odeint-v2-controlled_adams_bashforth_moulton */
#define ORBITINT_ODEINTV2_RKF78      8  /* orbit numerical integration method: odeint-v2-runge_kutta_fehlberg78 */
#define ORBITINT_ODEINTV2_RK_CASH_KARP54 9 /* orbit numerical integration method: odeint-v2-runge_kutta_cash_karp54 */

#define PMUT1_OCEAN_CORR        0      /* diurnal/subdiurnal tidal effects on polar motion("),UT1(s) and LOD(s) in time domain, IERS conventions 2010 */

#define SOLRAD_MODEL_DEFAULT    0      /* solar radiation pressure model: default */
#define SOLRAD_MODEL_ECOM1_5    1      /* solar radiation pressure model: ECOM1(5) */
#define SOLRAD_MODEL_ECOM1_9    2      /* solar radiation pressure model: ECOM1(9) */
#define SOLRAD_MODEL_ECOM2      3      /* solar radiation pressure model: ECOM2 */
#define SOLRAD_MODEL_ECOMC      4      /* solar radiation pressure model: ECOMC */
#define SOLRAD_MODEL_CODE_5     5      /* solar radiation pressure model: CODE(5) */
#define SOLRAD_MODEL_CODE_9     6      /* solar radiation pressure model: CODE(9) */

#define UDTRANSMATRIX_INTGRTN   1      /* update transition matrix using numerical integration */
#define UDTRANSMATRIX_NUMDIF_2P 2      /* update transition matrix using numerical differential: two-points */
#define UDTRANSMATRIX_NUMDIF_3P 3      /* update transition matrix using numerical differential: three-points */
#define ACCLGRAD_ANALYSIS       0      /* partial differential of acceleration wrt to parameters using analysis */  
#define ACCLGRAD_NUMDIFF_ALL    1      /* partial differential of acceleration wrt to parameters using numerical differential */
#define ACCLGRAD_NUMDIFF_SEP    2

#define GI_MERCURY  (0)                /* planets */
#define GI_VENUS    (1)
#define GI_MARS     (2)
#define GI_JUPITER  (3)
#define GI_SATURN   (4)
#define GI_URANUS   (5)
#define GI_NEPTUNE  (6)
#define GI_PLUTO    (7)
#define GI_MOON     (8)
#define GI_SUN      (9)
#define GI_SUNSSB   (10)
#define GI_EARTH    (11)

#define GNP        (3)                /* satellite position */
#define GNV        (3)                /* satellite velocity */
#define GNS        (15)               /* SRP parameters */
#define GNC        (2)                /* satellite clock */
#define GNDCB      (1)                /* satellite DCB */
#define GIP        (0)
#define GIV        (3)
#define GIS        (6)
#define GIS_D0     (6)                 /* srp parameters: D0,Dc,Ds,Dc2,Ds2,Dc4,Ds4,Y0,Yc,Ys,B0,Bc,Bs */
#define GIS_DC     (7)
#define GIS_DS     (8)
#define GIS_DC2    (9)
#define GIS_DS2    (10)
#define GIS_DC4    (11)
#define GIS_DS4    (12)
#define GIS_Y0     (13)
#define GIS_YC     (14)
#define GIS_YS     (15)
#define GIS_B0     (16)
#define GIS_BC     (17)
#define GIS_BS     (18)
#define GIC        (21)
#define GICR       (22)  

#define GNRCVS     (135)               /* max numbers of receivers */

#define GNRP       (3)                 /* receiver position */
#define GNRC       (6)                 /* receiver clock bias */
#define GNRCR      (1)                 /* receiver clock drift */
#define GNRTR      (3)                 /* receiver ZTD parameter */

#define GIRP       (0)
#define GIRC       (GNRP)
#define GIRD       (GNRP+GNRC)
#define GITR       (GNRP+GNRC+GNRCR)

#define GNX        (GNP+GNV+GNS+GNC)
#define GNRX       (GNRP+GNRC+GNRCR+GNRTR)

#define GPOD_OBSS_TYPE_SATPOS    0
#define GPOD_OBSS_TYPE_SATOBS    1

#define GPOD_FIXRCVPOS_ESTSATPOS 0
#define GPOD_FIXSATPOS_ESTRCVPOS 1

/* type definitions ----------------------------------------------------------*/
typedef struct {
    int year,mon,day;
    int bsrn,nd,avgap,c9,isn,q,aparr[8],kparr[8],sumkp;
    double adjf10,adjctrf81,adjlstf81,obsf10,obsctrf81,obslstf81,cp;
} spwd_t;

typedef struct {
    int n,nmax;
    spwd_t *data;
} spw_t;

typedef struct {
    double area;
    double mass;
    double CR,P0,Au;
    double D0,Y0,B0,Bc,Bs,Yc,Ys,Dc,Ds;
    double Dc2,Ds2,Dc4,Ds4;
    int shmopt;
    int solradmdl;
    char blocktype[32];
    struct force_model *fmdl;
} solrad_pressure_t;

typedef struct {
    int atmdensopt;
    double area;
    double mass;
    double CD;
    struct force_model *fmdl;
} atmospheric_drag_t;

typedef struct {
    double rr[32][4];
    struct force_model *fmdl;
} nbody_perturbation_t;

typedef struct {
    char name[32];
    int N,M;
    int solidtide;
    int oceantide;
    int poletide;
    int anelasticearth;
    int rotcorr;
    int tidesystem;
    double *C,*S;
    double RE,GM;
    struct force_model *fmdl;
} gravity_model_t;

typedef struct {
    struct force_model *fmdl;
} relatvty_mode_t;

typedef struct {
    char darw[16];
    int l,m;
    double dn[6];
    double dCp;
    double dSp;
    double dCm;
    double dSm;
} oceantide_model_coeffs_t;

typedef struct {
    int maxl,maxm;
    int n,nmax;
    oceantide_model_coeffs_t *mdldata;
} oceantide_model_t;

typedef struct {
    int relativity;
    int gravity;
    int nbody;
    int atmdrag;
    int solrad;
    int shmopt;
    int atmdensopt;
    int solidtide;
    int oceantide;
    int poletide;
    int anelasticearth;
    int solradmdl;
    int gN,gM;
    int oL,oM;
    int eopfiletype;
    int udtransmatrix;
    int udacclgrad;
    int usesamegravitymdl;
    int usesameeop;
    int usesamespw;
    int usesameotm;
    char jplde440file[MAXSTRPATH];
    char spwfile[MAXSTRPATH];
    char eopfile[MAXSTRPATH];
    char gfcfile[MAXSTRPATH];
    char octfile[MAXSTRPATH];
    double solradparas[8];
    double atmdrgparas[8];
    double CR,P0,Au;
    double D0,Y0,B0,Bc,Bs,Yc,Ys,Dc,Ds;
    double Dc2,Ds2,Dc4,Ds4;
} force_model_opt_t;

typedef struct force_model {
    erp_t erp;
    spw_t spw;
    gravity_model_t gm;
    atmospheric_drag_t ad;
    nbody_perturbation_t np;
    solrad_pressure_t sp;
    force_model_opt_t opt;
    oceantide_model_t ot;
    relatvty_mode_t rt;

    int noneedfreeerp;
    int noneedfreespw;
    int noneedfreegry;
    int noneedfreeotm;
} force_model_t;

typedef struct pod_rcv {
    int id;
    int clkref;
    char name[8];
    double pos[3];
    sta_t sta;
} pod_rcv_t;

typedef struct {
    char id[8];
    char blocktype[32];
    char yaw_bias[8];
    double mass;
    double yaw_rate;
} pod_satinfo_t;

typedef struct podobs {
    int type;
    int nrcv;
    int nobs[GNRCVS];
    int ircv[GNRCVS];
    char name[GNRCVS][8];
    obsd_t obs[GNRCVS][MAXOBS];
    double pos[MAXSAT][6];
    double vel[MAXSAT][6];
    gtime_t tutc;
} podobs_t;

typedef struct podobss {
    int n,nmax;
    podobs_t *data;
    obs_t obss[GNRCVS];
    sta_t stas[GNRCVS];
} podobss_t;

typedef struct satorbit {
    gtime_t tutc;
    int sat;
    int state;
    double x[GNX];
    double F[GNX*GNX];
    force_model_t fmdl;
} satorbit_t;

typedef struct pod_flt {
    gtime_t time;
    ssat_t ssat[GNRCVS][MAXSAT];
    sol_t sols[GNRCVS];
    int nx,state;
    int ismap[MAXSAT];
    int ns;
    int ibmap[GNRCVS][MAXSAT];
    int nbias[GNRCVS];
    double dt;
    double *x,*P;
    double *xp,*Pp;
    double *x0,*P0;
    double clkref[NSYS];
} pod_flt_t;

typedef struct {
    force_model_opt_t fmdlopts[MAXSAT];
    prcopt_t estopt;
    sta_t podrcvs[GNRCVS];
    double satvar[GNX];
    double satprn[GNX];
    double rcvvar[GNRX];
    char navfile[MAXSTRPATH];
    char clkref[MAXSTRPATH];
    char satantp[MAXSTRPATH];
    char rcvantp[MAXSTRPATH];
    char dcb[MAXSTRPATH];
    char blq[MAXSTRPATH];
    char staposfile[MAXSTRPATH];
    char snxfile[MAXSTRPATH];
    int udorbitint;
    int estmode;
    int onlyusecode;
    int outldetexcs;
    int satclkmode;
} pod_opt_t;

typedef struct pod_sol {
    gtime_t tutc;
    int ns[MAXSAT];
    int stat;
    double xs[MAXSAT][GNX];
    double xr[GNRCVS][GNRX];
    double Qs[MAXSAT][GNX];
    double Qr[GNRCVS][GNRX];
    double dops[MAXSAT];
} pod_sol_t;

typedef struct pod {
    pod_flt_t flt;
    pod_sol_t sol;
    pod_opt_t opt;
    pod_rcv_t rcv[GNRCVS];
    satorbit_t orbits[MAXSAT];
    nav_t nav;
    pcvs_t pcvss;
    pcvs_t pcvsr;
} pod_t;

/* satellite informations---------------------------------------------------*/
extern pod_satinfo_t gsatinfos[256];
extern pod_satinfo_t *podgetsatinfo(int sat);

/* estimate state index define----------------------------------------------*/
extern int GIXSAT(pod_t *pod, int sat);
extern int GIXRCV(pod_t *pod, int rcv);
extern int GIXRCV_POS(pod_t *pod, int rcv);
extern int GIXRCV_CLK(pod_t *pod, int rcv, int sys);
extern int GIXRCV_CKR(pod_t *pod, int rcv);
extern int GIXRCV_TRP(pod_t *pod, int rcv);
extern int GIXSAT_POS(pod_t *pod, int sat);
extern int GIXSAT_VEL(pod_t *pod, int sat);
extern int GIXSAT_SRP(pod_t *pod, int sat, int i);
extern int GIXSAT_CLK(pod_t *pod, int sat);
extern int GIXSAT_CKR(pod_t *pod, int sat);
extern int GIXSAT_BIAS(pod_t *pod, int rcv, int sat);
extern int GIXDCB(pod_t *pod);

/* function definitions-------------------------------------------------------*/
typedef int (*deqfunc)(double t, const double *y, double *dydt, void *params);

extern double cal2mjd(const double *dt);
extern void mjd2cal(const double mjd, double *dt);
extern double mjdtdb(const double mjd, const double *erpv);
extern double mjdtt(const double mjd, const double *erpv);
extern double mjdut1(const double mjd, const double *erpv);
extern void reftimediff(const double ut1_utc, const double tai_utc, double *rdt);
extern double doy(const double *ep);
extern double time2mjd(gtime_t tutc);
extern double time2mjdtt(gtime_t tutc, const double *erpv);
extern double time2mjdut1(gtime_t tutc, const double *erpv);

extern int rkf78m(deqfunc func, void *param, const double *y0, double x, double h, double xmax, double tol, double *y, int n);
extern void ecef2satf(const double *rs, double *E);
extern void ecsf2ecef(gtime_t tutc, const double *erpv, double *U, double *dU, double *gmst);
extern void eci2ecef2(gtime_t tutc, const erp_t *erp, double *U, double *dU, double *gmst);
extern void satfixed(const double *rsat, const double *rsun, double *ex, double *ey, double *ez);
extern void satf2ecsf(const double *rsat, const double *rsun, double *E);
extern void ecsf2satf(const double *rs, double *E);
extern double relcorr(const double *rsat, const double *rrcv, int opt);
extern void pmut1_ocean(gtime_t tutc, double *cor);
extern void utlibr(gtime_t tutc, double *cor);
extern void pmsdnut2(gtime_t tutc, double *cor);
extern int readjplde440coeffs(const char *file);
extern void rcvposeci(gtime_t tutc, const double *rr_ecef, const erp_t *erp, double dt, double *rr_eci, double *T);

extern void doodsonargs(gtime_t tutc, const erpd_t *erp, double *beta, double *funt);
extern void legendre(int N, int M, double fi, double *P, double *P1, double *P2);

extern int jpleph_de440(gtime_t tutc, const double *erpv, double *mercury, double *venus,
                        double *earth, double *mars, double *jupiter, double *saturn,
                        double *uranus, double *neptune, double *pluto, double *moon,
                        double *sun, double *sunssb);

extern int readgfc(const char *file, double **C, double **S, double *R, double *GM, char *name, int *tidesyetem);
extern int readeop(const char *file, erp_t *erp, int type);
extern int readspw(const char *file, spw_t *spw);
extern int readotfile(const char *file, oceantide_model_t *otmdl);
extern int geteop(const erp_t *erp, gtime_t time, erpd_t *erpd);

extern void forcemdlinit(force_model_t *fmdl, const force_model_opt_t *opt);
extern void forcemdlfree(force_model_t *fmdl);
extern void gravityforce(gtime_t tutc, const erpd_t *erp, const gravity_model_t *gm, const double *rs, double *fa, double *dadx, int uddadxmethod);
extern void solradpforce(gtime_t tutc, const erpd_t *erp, const solrad_pressure_t *sp, const double *rs, double *fa, double *dadx, int uddadxmethod);
extern void atmdragforce(gtime_t tutc, const erpd_t *erp, const atmospheric_drag_t *ad, const double *rs, double *fa, double *dadx, int uddadxmethod);
extern void nbpertnforce(gtime_t tutc, const erpd_t *erp, nbody_perturbation_t *np, const double *rs, double *fa, double *dadx, int uddadxmethod);
extern void relattyforce(gtime_t tutc, const erpd_t *erp, const relatvty_mode_t *rt, const double *rs, double *fa, double *dadx, int uddadxmethod);
extern void satforce(force_model_t *fmdl, gtime_t tutc, const erpd_t *erp, const double *rs, double *fa, double *dadx);

extern int satposeci(gtime_t tutc, const erp_t *erp, const double *rs_ecef, double *rs_eci);
extern int satposecef(gtime_t tutc, const erp_t *erp, const double *rs_eci, double *rs_ecef);

extern int readpodobs_sp3(const char **files, int n, podobss_t *obss);
extern int readpodobs_rnx(const char **files, int n, podobss_t *obss);
extern void readpodobs_rnxe(const char **files, int n, podobss_t *obss);
extern void readpodobs_rnxet(const char **files, int n, podobss_t *obss, gtime_t ts, gtime_t te, double tint);
extern int sortpodobs(podobss_t *obss);
extern int podobssget(const podobss_t *obss, gtime_t tc, podobs_t *podobs);

extern int satorbit(int method, satorbit_t *orbit, double dt);
extern void satorbitinit(satorbit_t *orbit, const force_model_opt_t *opt, const double *x0, gtime_t tutc0);
extern void satorbitfree(satorbit_t *orbit);

extern void orbintj2_U(gtime_t tutc, const double *U, double dt, const double *x0, double *x, double *Phi);

extern void orbitgetx(const satorbit_t *orbit, double *x);
extern void orbitsetx(const double *x, satorbit_t *orbit);

extern void podinit(pod_t *pod, const pod_opt_t *opt);
extern void podfree(pod_t *pod);
extern void podsatinit(pod_t *pod, int sat, const double *x0, gtime_t tutc0);
extern void podudrcvsnx(pod_t *pod, const sta_t *stas, int nsta, const char *snxfile);
extern void podsetclkref(pod_t *pod, const char *staname);

extern int podflt(pod_t *pod, const podobs_t *obs);

#ifdef __cplusplus
}
#endif
#endif