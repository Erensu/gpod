/*------------------------------------------------------------------------------
 * gpod.c : GNSS orbit estimate functions

 * author  : sujinglan
 * version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
 * Copyright(c) 2023-2025 by sujinglan, all rights reserved
 * history : 2024/10/17 1.0  new
 *----------------------------------------------------------------------------*/
#include "podlib.h"

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

#define THRES_REJECT 4.0            /* reject threshold of posfit-res (sigma) */
#define THRES_MW_JUMP 10.0

#define CONSTEPH_POSVEL 0
#define CONSTEPH_CLK    0
#define ORBINTJ2     1
#define MAXDOPS      10.0

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_VEL     SQR(10.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */

#define NF(opt)    ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)

#define ZD_SATI(val) ((val>>16)&0xFF)
#define ZD_RCVI(val) ((val>> 8)&0xFF)
#define ZD_TYPE(val) ((val>> 4)&0x0F)
#define ZD_IRES(val) ((val    )&0x0F)

/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
                   const pcvs_t *pcvr, const sta_t *sta)
{
    pcv_t *pcv,pcv0={0};
    double pos[3],del[3];
    int i,j,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    char id[64];

    /* set satellite antenna parameters */
    for (i=0;i<MAXSAT;i++) {
        nav->pcvs[i]=pcv0;
        if (!(satsys(i+1,NULL)&popt->navsys)) continue;
        if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
            satno2id(i+1,id);
            log_trace(3,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->pcvs[i]=*pcv;
    }
    popt->pcvr[0]=pcv0;
    strcpy(popt->anttype[0],sta->antdes);

    if (sta->deltype==1) { /* xyz */
        if (norm(sta->pos,3)>0.0) {
            ecef2pos(sta->pos,pos);
            ecef2enu(pos,sta->del,del);
            for (j=0;j<3;j++) popt->antdel[0][j]=del[j];
        }
    }
    else { /* enu */
        for (j=0;j<3;j++) popt->antdel[0][j]=sta->del[j];
    }
    if (!(pcv=searchpcv(0,popt->anttype[0],time,pcvr))) {
        log_trace(2,"no receiver antenna pcv: %s\n",popt->anttype[0]);
        *popt->anttype[0]='\0';
        return;
    }
    strcpy(popt->anttype[0],pcv->type);
    popt->pcvr[0]=*pcv;
}
/* set receiver for POD------------------------------------------------------*/
static void setrcvs(pod_t *pod, const pod_opt_t *opt)
{
    int i;

    for (i=0;i<GNRCVS;i++) {
        strncpy(pod->rcv[i].name,opt->podrcvs[i].name,4);
        pod->rcv[i].id=i+1;
        pod->rcv[i].clkref=0;
    }
}
/* read ocean tide loading parameters ---------------------------------------*/
static void readotl(prcopt_t *popt, const char *file, const sta_t *sta)
{
    static double odisp[MAXRCV][6*11]={0};

    if (odisp[sta->staid][0]) {
        matcpy(popt->odisp[0],odisp[sta->staid],1,6*11);
        return;
    }
    readblq(file,sta->name,popt->odisp[0]);
    matcpy(odisp[sta->staid],popt->odisp[0],1,6*11);
}
/* read station position from file-------------------------------------------*/
static void readstapos_snx(pod_t *pod, const pod_opt_t *opt)
{
    podudrcvsnx(pod,opt->podrcvs,GNRCVS,opt->snxfile);
}
static void readstapos_xyz(pod_t *pod, const pod_opt_t *opt)
{
    FILE *fp;
    char buff[1024];
    int i;

    for (i=0;i<GNRCVS;i++) {
        if (strcmp(pod->rcv[i].name,"")==0) continue;
        if (norm(pod->rcv[i].pos,3)>0.0) continue;
        if (!(fp=fopen(opt->staposfile,"r"))) return;

        double pos[3];
        char name[16];
        while (fgets(buff,sizeof(buff),fp)) {
            if (sscanf(buff,"%s %lf %lf %lf\n",name,pos,pos+1,pos+2)<4) continue;

            if (strcmp(name,pod->rcv[i].name)==0) {
                if (norm(pos,3)) matcpy(pod->rcv[i].pos,pos,1,3);
                break;
            }
        }
        fclose(fp);
    }
}
static void readstapos(pod_t *pod, const pod_opt_t *opt)
{
    readstapos_snx(pod,opt);
    readstapos_xyz(pod,opt);
}
/* initial orbit estimate-----------------------------------------------------
 * args:      pod_t *pod      IO  orbit estimate struct
 *            pod_opt_t *opt  I   orbit estimate options
 * return : none
 *---------------------------------------------------------------------------*/
extern void podinit(pod_t *pod, const pod_opt_t *opt)
{
    int i,j,nx=MAXSAT*GNX+GNRCVS*GNRX+GNRCVS*MAXSAT+GNDCB;
    double x0[6]={0};
    gtime_t t0={0};

    /* set receiver options */
    setrcvs(pod,opt);

    /* receiver positions */
    readstapos(pod,opt);

    /* set reference clock */
    podsetclkref(pod,opt->clkref);

    /* read navigation date */
    readrnx(opt->navfile,0,"",NULL,&pod->nav,NULL);
    readrnx(opt->navfile2,0,"",NULL,&pod->nav,NULL);
    uniqnav(&pod->nav);

    /* read satellite antenna parameters */
    readpcv(opt->satantp,&pod->pcvss);

    /* read receiver antenna parameters */
    readpcv(opt->rcvantp,&pod->pcvsr);

    /* read dcb parameters */
    readdcb(opt->dcb,&pod->nav,NULL);

    for (i=0;i<MAXSAT;i++) {
        satorbitinit(&pod->orbits[i],&opt->fmdlopts[i],x0,t0);

        pod_satinfo_t *satinfo=podgetsatinfo(i+1);
        if (satinfo==NULL) continue;
        pod->orbits[i].fmdl.sp.mass=satinfo->mass;
        strcpy(pod->orbits[i].fmdl.sp.blocktype,satinfo->blocktype);
    }
    pod->flt.nx=nx;
    pod->flt.x =zeros( 1,nx);
    pod->flt.P =zeros(nx,nx);
    pod->flt.xp=zeros( 1,nx);
    pod->flt.Pp=zeros(nx,nx);
    pod->opt=*opt;

    if (pod->flt.x ==NULL||pod->flt.P ==NULL||
        pod->flt.xp==NULL||pod->flt.Pp==NULL) {
        fprintf(stderr,"no enough memory\n");
        exit(0);
    }
}
/* free orbit estimate-------------------------------------------------------*/
extern void podfree(pod_t *pod)
{
    if (pod->flt.x ) free(pod->flt.x );
    if (pod->flt.P ) free(pod->flt.P );
    if (pod->flt.xp) free(pod->flt.xp);
    if (pod->flt.Pp) free(pod->flt.Pp);
    if (pod->flt.x0) free(pod->flt.x0);
    if (pod->flt.P0) free(pod->flt.P0);

    freenav(&pod->nav,0xFF);

    /* free antenna parameters */
    free(pod->pcvss.pcv); pod->pcvss.pcv=NULL; pod->pcvss.n=pod->pcvss.nmax=0;
    free(pod->pcvsr.pcv); pod->pcvsr.pcv=NULL; pod->pcvsr.n=pod->pcvsr.nmax=0;

    int i;
    for (i=0;i<MAXSAT;i++) {
        satorbitfree(&pod->orbits[i]);
    }
}
/* set reference clock station------------------------------------------------*/
extern void podsetclkref(pod_t *pod, const char *staname)
{
    int i;
    for (i=0;i<GNRCVS;i++) {
        if (strcmp(staname,pod->rcv[i].name)==0) {
            pod->rcv[i].clkref=1;
        }
    }
}
/* satellite position(t+dt) in ECI--------------------------------------------*/
static void satposdt(const double *state, double dt, double *rs, double *F)
{
    memset(F,0,sizeof(double)*9);
    F[0]=F[4]=F[8]=1.0;
    if (norm(state,3)<=0.0) return;

    const double GM_E=398600.4415E9;
    double r=norm(state,3),k=1.0-GM_E/(2.0*r*r*r)*dt*dt,r2=r*r,g;
    rs[0]=state[0]*k+state[3]*dt;
    rs[1]=state[1]*k+state[4]*dt;
    rs[2]=state[2]*k+state[5]*dt;

    /* state transition matrix(t to t+dt) */
    g=GM_E/(2.0*r2*r2*r)*dt*dt;
    F[0]=1.0+g*(3.0*state[0]*state[0]-r2);
    F[4]=1.0+g*(3.0*state[1]*state[1]-r2);
    F[8]=1.0+g*(3.0*state[2]*state[2]-r2);
    F[1]=F[3]=g*3.0*state[0]*state[1];
    F[2]=F[6]=g*3.0*state[0]*state[2];
    F[5]=F[7]=g*3.0*state[1]*state[2];
}
/* station position(t+dt) in ECI----------------------------------------------*/
static void staposdt(gtime_t tutc, const double *state, const erp_t *erp, double dt, double *rr, double *T)
{

#define Rz(t,X) do {                                     \
    (X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0;         \
    (X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)

    double erpv[6]={0},U[9],R[9]={0};

    /* earth rotation parameter values */
    geterp(erp,utc2gpst(tutc),erpv);

    /* ECI to ECEF transformation matrix */
    ecsf2ecef(tutc,erpv,U,NULL,NULL);

    /* correct earth rotation effect */
    Rz(-OMGE*dt,R);

    /* station position in ECI */
    matmul33("TNN",U,R,state,3,3,3,1,rr);

    if (T) {
        matmul("TN",3,3,3,1.0,U,R,0.0,T);
    }
}
extern void rcvposeci(gtime_t tutc, const double *rr_ecef, const erp_t *erp, double dt, double *rr_eci, double *T)
{
    staposdt(tutc,rr_ecef,erp,dt,rr_eci,T);
}
/* satellite-station range----------------------------------------------------
 * args   : gtime_t tutc  I  time in utc
 *          erp_t *erp    I  earth rotation parameters
 *          double dtr    I  station receiver clock bias (s)
 *          double dts    I  satellite clock bias (s)
 *          double *rr    I  station position (ECEF/m)
 *          double *rs    I  satellite position (ECI/m)
 *          double *rstx  O  satellite position on transmit signal (ECI/m)
 *          duoble *rrrx  O  station position on receive signal (ECEF/m)
 *          double *drds  O  partial derivatives of range by satellite state
 * return: satellite-station range (m)
 *----------------------------------------------------------------------------*/
extern double ssrange(gtime_t tutc, const erp_t *erp, double dtr, double dts, const double *rr, const double *rs,
                      double *rstx, double *rrrx, double *drds)
{
    double rrt[3],rst[3]={0},dt,rng=0.0,F[9]={0},rsr[3],rngt;
    int i,j;

    staposdt(tutc,rr,erp,-dtr,rrt,NULL);

    for (i=0;i<16;i++) {
        dt=dtr+rng/CLIGHT;
        satposdt(rs,-dt,rst,F);
        for (j=0;j<3;j++) rsr[j]=rst[j]-rrt[j];
        rngt=rng;
        rng=norm(rsr,3);
        if (fabs(rng-rngt)<1E-4) break;
    }
    if (drds) {
        for (i=0;i<3;i++) {
            rsr[i]/=rng;
        }
        matmul("TN",3,1,3,1.0,F,rsr,0.0,drds);
        for (i=0;i<3;i++) {
            drds[3+i]=-rsr[i]*dt;
        }
    }
    if (rstx) matcpy(rstx,rst,1,3);
    if (rrrx) matcpy(rrrx,rrt,1,3);
    return rng;
}
/* estimate orbit using satellite position/velocity measurement---------------
 * args:      pod_t *pod    IO  orbit estimate struct
 *            int sat       I   satellite no
 *            double *rs    I   satellite position/velocity measurement (ECEF,m,m/s)
 *            double *var   I   satellite position/velocity measurement variance (ECEF,m^2,m^2/s^2)
 *            gtime_t tutc  I   time of satellite position
 * return : status (0:fail,1:ok)
 * --------------------------------------------------------------------------*/
extern int podsatposflt(pod_t *pod, int sat, const double *rs, const double *var, gtime_t tutc)
{
    satorbit_t *orbit=&pod->orbits[sat-1];
    pod_opt_t *opt=&pod->opt;
    int i,j,state=1,nx=GNX,nxa=pod->flt.nx;
    double erpv[6]={0},U[9],dU[9],rs_eci[3]={0},vs_eci[3]={0},dt=timediff(tutc,orbit->tutc);
    double *xp,*Pp,*P0,*R,*v,*H;
    double v1[3],v2[3];

    /* invalid update time */
    if (!tutc.time) return 0;

    /* initial estimated states and covariance */
    if (norm(pod->flt.x+GIXSAT(pod,sat),nx)<=0.0) {

        for (i=0;i<nx;i++) {
            pod->flt.x[GIXSAT(pod,sat)+i]=orbit->x[i];
            for (i=0;i<nx;i++) {
                pod->flt.P[i+GIXSAT(pod,sat)+(i+GIXSAT(pod,sat))*pod->flt.nx]=opt->satvar[i];
            }
        }
    }
    if (fabs(dt)<1E-10) return 0;

    matcpy(pod->flt.x0,pod->flt.x,  1,nxa);
    matcpy(pod->flt.P0,pod->flt.P,nxa,nxa);

    /* predict satellite orbit */
    satorbit(opt->udorbitint,orbit,dt);
    pod->flt.dt=dt;

    xp=mat(1,nx); Pp=mat(nx,nx); P0=mat(nx,nx); v=zeros(1,6); H=zeros(6,nx); R=zeros(6,6);

    /* predict estimated states and covariance */
    for (i=0;i<nx;i++) {
        xp[i]=orbit->x[i];
        for (j=0;j<nx;j++) {
            P0[i+j*nx]=pod->flt.P[i+GIXSAT(pod,sat)+(j+GIXSAT(pod,sat))*pod->flt.nx];
        }
    }
    matmul33("NNT",orbit->F,P0,orbit->F,nx,nx,nx,nx,Pp);
    for (i=0;i<nx;i++) Pp[i+i*nx]+=SQR(opt->satprn[i])*fabs(dt);

    for (i=0;i<nx;i++) {
        pod->flt.x [i+GIXSAT(pod,sat)]=xp[i];
        pod->flt.xp[i+GIXSAT(pod,sat)]=xp[i];
        for (j=0;j<nx;j++) {
            pod->flt.P [i+GIXSAT(pod,sat)+(j+GIXSAT(pod,sat))*pod->flt.nx]=Pp[i+j*nx];
            pod->flt.Pp[i+GIXSAT(pod,sat)+(j+GIXSAT(pod,sat))*pod->flt.nx]=Pp[i+j*nx];
        }
    }
    /* earth rotation parameter values */
    geterp(&orbit->fmdl.erp,utc2gpst(tutc),erpv);
    if (norm(erpv,4)<=0.0) {

        free(xp); free(Pp); free(P0); free(v); free(H); free(R);
        return 0;
    }
    /* invalid satellite position measurement */
    if (norm(rs,6)<=0.0||norm(var,6)<=0.0) {

        free(xp); free(Pp); free(P0); free(v); free(H); free(R);
        return 0;
    }
    /* compute satellite position in ECI(m) */
    ecsf2ecef(tutc,erpv,U,dU,NULL);
    matmul("TN",3,1,3,1.0,U,rs,0.0,rs_eci);

    /* satellite velocity in ECI(m/s) */
    if (norm(rs+3,3)>0.0) {
        matmul("TN",3,1,3,1.0,U,rs+3,0.0,v1);
        matmul("TN",3,1,3,1.0,dU,rs,0.0,v2);
        for (i=0;i<3;i++) vs_eci[i]=v1[i]+v2[i];
    }
    /* update estimated states and covariance */
    for (i=0;i<6;i++) R[i+i*6]=SQR(100.0);

    /* satellite position measurement */
    for (i=0;i<3;i++) {
        if (!rs_eci[i]||var[i]<=0.0) continue;
        v[i]=rs_eci[i]-xp[GIP+i];
        R[i+i*6]=var[i];
        H[GIP+i+i*GNX]=1.0;
    }
    /* satellite velocity measurement */
    for (i=0;i<3;i++) {
        if (!vs_eci[i]||var[i+3]<=0.0) continue;
        v[i+3]=vs_eci[i]-xp[GIV+i];
        R[i+3+(i+3)*6]=var[i+3];
        H[GIV+i+(i+3)*GNX]=1.0;
    }
    /* EKF satellite position measurement update */
    if (filter(xp,Pp,H,v,R,nx,6,NULL,0)) state=0;

    /* update satellite states and covariance */
    if (state) {
        pod->flt.time=tutc;
        orbitsetx(xp,orbit);
        for (i=0;i<nx;i++) {
            pod->flt.x[i+GIXSAT(pod,sat)]=xp[i];
            for (j=0;j<nx;j++) {
                pod->flt.P[i+GIXSAT(pod,sat)+(j+GIXSAT(pod,sat))*pod->flt.nx]=Pp[i+j*nx];
            }
        }
    }
    free(xp); free(Pp); free(P0); free(v); free(H); free(R);
    return state;
}
/* satellite orbit solution data--------------------------------------------*/
typedef struct orbfitsol {
    gtime_t tutc;
    gtime_t tutcp;
    double x [GNX];
    double x0[GNX];
    double xp[GNX];
    double F [GNX*GNX];
    double P [GNX*GNX];
    double P0[GNX*GNX];
    double Pp[GNX*GNX];
} orbfitsol_t;

/* fit satellite orbit using precise ephemeris by FBS-----------------------*/
static int podsatorbfitfbsprc(pod_t *pod, const podobss_t *obss, int sat, orbfitsol_t *fitsols)
{
    double rss[6],var[6],x0[GNX];
    int i,j,k,state;

    if (obss->n<=0) return 0;

    orbfitsol_t* solsf=calloc(obss->n,sizeof(orbfitsol_t));
    orbfitsol_t* solsb=calloc(obss->n,sizeof(orbfitsol_t));
    int nsolf=0,nsolb=0;
    int nsols=0;
    pod_flt_t *flt=&pod->flt;

    for (j=0;j<3;j++) {
        rss[j+0]=obss->data[0].pos[sat-1][j];
        rss[j+3]=obss->data[0].vel[sat-1][j];
    }
    /* satellite position/velocity in ECI */
    satposeci(obss->data[0].tutc,&pod->orbits[sat-1].fmdl.erp,rss,x0);
    podsatinit(pod,sat,x0,obss->data[0].tutc);

    for (i=0;i<obss->n;i++) {
        for (j=0;j<3;j++) {
            rss[j+0]=obss->data[i].pos[sat-1][j+0];
            rss[j+3]=obss->data[i].vel[sat-1][j+0];
            var[j+0]=obss->data[i].pos[sat-1][j+3];
            var[j+3]=obss->data[i].vel[sat-1][j+3];
        }
        /* satellite orbit estimate */
        if (!(state=podsatposflt(pod,sat,rss,var,obss->data[i].tutc))) continue;

        solsf[nsolf].tutc=obss->data[i].tutc;
        matcpy(solsf[nsolf].x,&flt->x[GIXSAT(pod,sat)],1,GNX);
        for (j=0;j<GNX;j++) {
            for (k=0;k<GNX;k++) {
                solsf[nsolf].P[j+k*GNX]=flt->P[GIXSAT(pod,sat)+j+(GIXSAT(pod,sat)+k)*flt->nx];
            }
        }
        nsolf++;
    }
    for (j=0;j<3;j++) {
        rss[j+0]=obss->data[obss->n-1].pos[sat-1][j];
        rss[j+3]=obss->data[obss->n-1].vel[sat-1][j];
    }
    /* satellite position/velocity in ECI */
    satposeci(obss->data[obss->n-1].tutc,&pod->orbits[sat-1].fmdl.erp,rss,x0);
    podsatinit(pod,sat,x0,obss->data[obss->n-1].tutc);

    for (i=obss->n-1;i>=0;i--) {
        for (j=0;j<3;j++) {
            rss[j+0]=obss->data[i].pos[sat-1][j+0];
            rss[j+3]=obss->data[i].vel[sat-1][j+0];
            var[j+0]=obss->data[i].pos[sat-1][j+3];
            var[j+3]=obss->data[i].vel[sat-1][j+3];
        }
        /* satellite orbit estimate */
        if (!(state=podsatposflt(pod,sat,rss,var,obss->data[i].tutc))) continue;

        solsb[nsolb].tutc=obss->data[i].tutc;
        matcpy(solsb[nsolb].x,&flt->x[GIXSAT(pod,sat)],1,GNX);
        for (j=0;j<GNX;j++) {
            for (k=0;k<GNX;k++) {
                solsb[nsolb].P[j+k*GNX]=flt->P[GIXSAT(pod,sat)+j+(GIXSAT(pod,sat)+k)*flt->nx];
            }
        }
        nsolb++;
    }
    for (i=0;i<nsolf;i++) {
        orbfitsol_t sol=solsf[i];

        for (k=-1,j=0;j<nsolb;j++) {
            if (fabs(timediff(solsf[i].tutc,solsb[j].tutc))<1E-5) {
                k=j; break;
            }
        }
        if (k<0) continue;

        double xs[GNX],Ps[GNX*GNX];
        if (smoother(solsf[i].x,solsf[i].P,solsb[k].x,solsb[k].P,GNX,xs,Ps)) continue;

        fitsols[nsols].tutc=solsf[i].tutc;
        matcpy(fitsols[nsols].x,xs,  1,GNX);
        matcpy(fitsols[nsols].P,Ps,GNX,GNX);
        nsols++;
    }
    free(solsf);
    free(solsb);
    return nsols;
}
static int podsatorbfitfbs(pod_t *pod, const podobss_t *obss, int sat, orbfitsol_t **fitsols, int *nsols)
{
    *fitsols=calloc(obss->n,sizeof(orbfitsol_t));
    *nsols=podsatorbfitfbsprc(pod,obss,sat,*fitsols);
    return *nsols;
}
/* RTS smoother----------------------------------------------------------------*/
static int rts(const double *xf_k_1, const double *Pf_k_1, const double *Ff_k_1,
               const double *xpf_k, const double *Ppf_k, const double *xs_k, const double *Ps_k,
               double *xs_k_1, double *Ps_k_1, int nx)
{
    double *dx=mat(1,nx),*g=mat(1,nx),*dP=mat(nx,nx),*Pg=mat(nx,nx);
    double *G=mat(nx,nx),*Pi=mat(nx,nx);
    int i,stat=0;

    matcpy(Pi,Ppf_k,nx,nx);
    if (!matinv(Pi,nx)) {
        /* G=Pf_k_1*Ff_k_1'*Pi */
        matmul33("NTN",Pf_k_1,Ff_k_1,Pi,nx,nx,nx,nx,G);

        /* xs_k_1=xf_k_1+G*(xs_k-xpf_k) */
        for (i=0;i<nx;i++) dx[i]=xs_k[i]-xpf_k[i];
        matmul("NN",nx,1,nx,1.0,G,dx,0.0,g);
        for (i=0;i<nx;i++) xs_k_1[i]=xf_k_1[i]+g[i];

        /* Ps_k_1=Pf_k_1+G*(Ps_k-Ppf_k)*G' */
        for (i=0;i<nx*nx;i++) dP[i]=Ps_k[i]-Ppf_k[i];
        matmul33("NNT",G,dP,G,nx,nx,nx,nx,Pg);
        for (i=0;i<nx*nx;i++) Ps_k_1[i]=Pf_k_1[i]+Pg[i];
        stat=1;
    }
    free(G); free(Pi); free(dx);
    free(g); free(dP); free(Pg);
    return stat;
}
/* fit satellite orbit using precise ephemeris by RTS-----------------------*/
static int podsatorbfitrtsprc(pod_t *pod, const podobss_t *obss, int sat, orbfitsol_t *fitsols)
{
    double rss[6],var[6],x0[GNX];
    int i,j,k,nx=GNX;

    if (obss->n<=0) return 0;

    orbfitsol_t* solsf=calloc(obss->n,sizeof(orbfitsol_t));
    pod_flt_t *flt=&pod->flt;
    pod_opt_t *opt=&pod->opt;
    satorbit_t *orbit=&pod->orbits[sat-1];
    int nsolf=0,nsols=0;

    for (j=0;j<3;j++) {
        rss[j+0]=obss->data[0].pos[sat-1][j];
        rss[j+3]=obss->data[0].vel[sat-1][j];
    }
    satposeci(obss->data[0].tutc,&orbit->fmdl.erp,rss,x0);
    podsatinit(pod,sat,x0,obss->data[0].tutc);

    for (i=0;i<obss->n;i++) {
        for (j=0;j<3;j++) {
            rss[j+0]=obss->data[i].pos[sat-1][j+0];
            rss[j+3]=obss->data[i].vel[sat-1][j+0];
            var[j+0]=obss->data[i].pos[sat-1][j+3];
            var[j+3]=obss->data[i].vel[sat-1][j+3];
        }
        solsf[nsolf].tutcp=flt->time;
        if (!podsatposflt(pod,sat,rss,var,obss->data[i].tutc)) continue;

        solsf[nsolf].tutc=flt->time;
        for (j=0;j<nx;j++) {
            solsf[nsolf].x [j]=flt->x [GIXSAT(pod,sat)+j];
            solsf[nsolf].x0[j]=flt->x0[GIXSAT(pod,sat)+j];
            solsf[nsolf].xp[j]=flt->xp[GIXSAT(pod,sat)+j];

            for (k=0;k<nx;k++) {
                solsf[nsolf].P0[j+k*nx]=flt->P0[GIXSAT(pod,sat)+j+(GIXSAT(pod,sat)+k)*flt->nx];
                solsf[nsolf].P [j+k*nx]=flt->P [GIXSAT(pod,sat)+j+(GIXSAT(pod,sat)+k)*flt->nx];
                solsf[nsolf].Pp[j+k*nx]=flt->Pp[GIXSAT(pod,sat)+j+(GIXSAT(pod,sat)+k)*flt->nx];
            }
        }
        for (j=0;j<nx;j++) {
            for (k=0;k<nx;k++) {
                solsf[nsolf].F[j+k*nx]=orbit->F[j+k*nx];
            }
        }
        nsolf++;
    }
    double *xsp=zeros( 1,nx);
    double *xss=zeros( 1,nx);
    double *Psp=zeros(nx,nx);
    double *Pss=zeros(nx,nx);

    matcpy(xsp,solsf[nsolf-1].x, 1,nx);
    matcpy(Psp,solsf[nsolf-1].P,nx,nx);

    int nsol=0;
    for (i=nsolf-1;i>=1;i--) {
        if (!rts(solsf[i].x0,solsf[i].P0,solsf[i].F,solsf[i].xp,solsf[i].Pp,xsp,Psp,xss,Pss,nx)) continue;

        matcpy(xsp,xss, 1,nx);
        matcpy(Psp,Pss,nx,nx);

        matcpy(fitsols[nsol].x,xss, 1,nx);
        matcpy(fitsols[nsol].P,Pss,nx,nx);
        fitsols[nsol].tutc=solsf[i].tutcp;
        nsol++;
    }
    free(solsf);
    free(xsp); free(xss);
    free(Psp); free(Pss);
    return nsol;
}
static int podsatorbfitrts(pod_t *pod, const podobss_t *obss, int sat, orbfitsol_t **fitsols, int *nsols)
{
    int i;

    for (i=0;i<MAXSAT;i++) {
        if (sat>0&&i+1!=sat) continue;
        fitsols[i]=calloc(obss->n,sizeof(orbfitsol_t));
        nsols[i]=podsatorbfitrtsprc(pod,obss,i+1,fitsols[i]);
    }
    return 1;
}
/* fit satellite orbit using precise ephemeris by RTS/FBS--------------------
 * args:      int type        I   fit method type (0:FBS, 1:RTS)
 *            pod_opt_t *opt  I   orbit estimate options
 *            int sat         I   satellite no. (<0: all satellites)
 *            podobss_t *obss I   precise ephemeris observation data
 *            char *solfile   I   output solution data file path
 * return : status (0:fail,1:ok)
 * --------------------------------------------------------------------------*/
extern int podsatorbfit(int type, const pod_opt_t *opt, const podobss_t *obss, int sat, const char *solfile)
{
    FILE *fp=fopen(solfile,"w");
    if (fp==NULL) return 0;

    orbfitsol_t *fitsols=NULL;
    int i,j,state,nsols=0;
    pod_t pod={0};

    podinit(&pod,opt);

    if (type==1) {
        int nx=pod.flt.nx;
        pod.flt.x0=zeros( 1,nx);
        pod.flt.P0=zeros(nx,nx);
    }
    switch (type) {
        case 0: state=podsatorbfitfbs(&pod,obss,sat,&fitsols,&nsols); break;
        case 1: state=podsatorbfitrts(&pod,obss,sat,&fitsols,&nsols); break;
        default:
            state=podsatorbfitfbs(&pod,obss,sat,&fitsols,&nsols); break;
    }
    if (state) {
        char prn[8];
        satno2id(sat,prn);

        if (nsols) {
            fprintf(fp,"%4s\n",prn);
        }
        for (j=0;j<nsols;j++) {
            fprintf(fp,"%s POS: %15.5lf %15.5lf %15.5lf %15.5lf %15.5lf %15.5lf "
                       "SRP: D0=%8.3lf Ds=%8.3lf Dc=%8.3lf B0=%8.3lf Bs=%8.3lf Bc=%8.3lf Y0=%8.3lf Ys=%8.3lf Yc=%8.3lf\n",
                    time_str(fitsols[j].tutc,3),
                    fitsols[j].x[0],fitsols[j].x[1],fitsols[j].x[2],
                    fitsols[j].x[3],fitsols[j].x[4],fitsols[j].x[5],
                    fitsols[j].x[GIS_D0],fitsols[j].x[GIS_DS],fitsols[j].x[GIS_DC],
                    fitsols[j].x[GIS_B0],fitsols[j].x[GIS_BS],fitsols[j].x[GIS_BC],
                    fitsols[j].x[GIS_Y0],fitsols[j].x[GIS_YS],fitsols[j].x[GIS_YC]);
            fflush(fp);
        }
    }
    podfree(&pod);
    fclose(fp);
    if (fitsols) free(fitsols);
    return state;
}
/* initial satellite orbit estimator------------------------------------------*/
extern void podsatinit(pod_t *pod, int sat, const double *x0, gtime_t tutc0)
{
    int i,j,nx=GNX;

    satorbitinit(&pod->orbits[sat-1],&pod->opt.fmdlopts[sat-1],x0,tutc0);
    pod->orbits[sat-1].sat=sat;

    for (i=0;i<nx;i++) pod->flt.x[GIXSAT(pod,sat)+i]=0.0;
    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) {
            pod->flt.P[GIXSAT(pod,sat)+i+(GIXSAT(pod,sat)+j)*pod->flt.nx]=0.0;
        }
    }
}
/* index of satellite system (m=0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN) -----*/
static int indsys(int sys)
{
    switch (sys) {
        case SYS_GPS: return 0;
        case SYS_SBS: return 0;
        case SYS_GLO: return 1;
        case SYS_GAL: return 2;
        case SYS_CMP: return 3;
        case SYS_QZS: return 4;
        case SYS_IRN: return 5;
        default:
            return 0;
    }
    return 0;
}
/* satellite orbit estimator using EKF and satellite position-----------------*/
static int podflt_satpos(pod_t *pod, const podobs_t *obs)
{
    double rss[6];
    double var[6];
    int i,j,state;

    for (state=i=0;i<MAXSAT;i++) {
        for (j=0;j<3;j++) {
            rss[j+0]=obs->pos[i][j+0];
            rss[j+3]=obs->vel[i][j+0];
            var[j+0]=obs->pos[i][j+3];
            var[j+3]=obs->vel[i][j+3];
        }
        state+=podsatposflt(pod,i+1,rss,var,obs->tutc);
    }
    return state;
}
/* initialize state and covariance -------------------------------------------*/
static void initx(pod_t *pod, double xi, double var, int i)
{
    int j;
    pod->flt.x[i]=xi;
    for (j=0;j<pod->flt.nx;j++) {
        pod->flt.P[i+j*pod->flt.nx]=pod->flt.P[j+i*pod->flt.nx]=i==j?var:0.0;
    }
}
/* initial satellite orbit position/velocity using navigation data------------*/
static void initsatorbclksrp_nav(pod_t *pod, int sat, gtime_t tutc)
{
    double rs[6],dts[2],x0[6],var,relclk,dant[3];
    int i,svh;

    if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
        if (norm(&pod->flt.x[GIXSAT_POS(pod,sat)],3)&&
            norm(&pod->flt.x[GIXSAT_VEL(pod,sat)],3)&&
            norm(&pod->flt.x[GIXSAT_CLK(pod,sat)],2)) {

            double vars=0.0,maxvars=SQR(10.0);
            for (i=0;i<3;i++) vars+=pod->flt.P[i+GIXSAT_POS(pod,sat)+(i+GIXSAT_POS(pod,sat))*pod->flt.nx];
            if (vars<maxvars) return;
        }
    }
    gtime_t time=utc2gpst(tutc);
    if (!satpos(time,time,sat,pod->opt.estopt.sateph,&pod->nav,rs,dts,&var,&svh)) {
        return;
    }
    if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
        satantoff(time,rs,sat,&pod->nav,dant);
        for (i=0;i<3;i++) rs[i]-=dant[i];
    }
    satposeci(tutc,&pod->orbits[sat-1].fmdl.erp,rs,x0);
    relclk=relcorr(x0,NULL,0);
    podsatinit(pod,sat,x0,tutc);

    if (pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS) relclk=0.0;

    initx(pod,dts[0]*CLIGHT+relclk,pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS?0.0:pod->opt.satvar[GIC+0],GIXSAT_CLK(pod,sat));
    initx(pod,dts[1]*CLIGHT       ,pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS?0.0:pod->opt.satvar[GIC+1],GIXSAT_CKR(pod,sat));

    for (i=0;i<3;i++) {
        initx(pod,x0[i+0],pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS?0.0:pod->opt.satvar[GIP+i],GIXSAT_POS(pod,sat)+i);
        initx(pod,x0[i+3],pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS?0.0:pod->opt.satvar[GIV+i],GIXSAT_VEL(pod,sat)+i);
    }
    if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
        for (i=0;i<GNS;i++) {
            initx(pod,pod->orbits[sat-1].x[i+GIS],pod->opt.satvar[GIS+i],GIXSAT_SRP(pod,sat,i));
        }
    }
}
/* initial receiver position--------------------------------------------------*/
static void initrcvpos(pod_t *pod, int rcv, gtime_t tutc)
{
    if (rcv<0||rcv>MAXRCV) return;

    double dt=pod->flt.dt;
    int i;

    if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {

        if (pod->opt.estopt.mode==PMODE_FIXED) {
            for (i=0;i<3;i++) {
                initx(pod,pod->rcv[rcv-1].pos[i],0.0,GIXRCV_POS(pod,rcv)+i);
            }
            return;
        }
        /* initialize position for first epoch */
        if (norm(pod->flt.x+GIXRCV_POS(pod,rcv),3)<=0.0) {
            for (i=0;i<3;i++) {
                initx(pod,pod->rcv[rcv-1].pos[i],SQR(1E-2),GIXRCV_POS(pod,rcv)+i);
            }
        }
        /* static ppp mode */
        for (i=0;i<3;i++) {
            pod->flt.P[(i+GIXRCV_POS(pod,rcv))*(1+pod->flt.nx)]+=SQR(pod->opt.estopt.prn[5])*fabs(dt);
        }
    }
    else {
        /* initialize position for first epoch */
        if (norm(pod->flt.x+GIXRCV_POS(pod,rcv),3)<=0.0) {
            for (i=0;i<3;i++) {
                initx(pod,pod->flt.sols[rcv].rr[i],VAR_POS,GIXRCV_POS(pod,rcv)+i);
            }
        }
        /* static ppp mode */
        for (i=0;i<3;i++) {
            pod->flt.P[(i+GIXRCV_POS(pod,rcv))*(1+pod->flt.nx)]+=SQR(pod->opt.estopt.prn[5])*fabs(dt);
        }
    }
}
/* initial receiver clock-----------------------------------------------------*/
static void initrcvclk(pod_t *pod, int rcv, const obsd_t *obs, int nobs)
{
    prcopt_t opt=pod->opt.estopt;
    sol_t *sol=&pod->flt.sols[rcv];
    char msg[32];

    if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
        opt.mode=PMODE_FIXED;
        matcpy(opt.ru,pod->rcv[rcv-1].pos,1,3);
    }
    /* receiver single point positioning */
    if (!pntpos(obs,nobs,&pod->nav,&opt,sol,NULL,pod->flt.ssat[rcv],msg)) {
        return;
    }
    if (sol->stat!=SOLQ_SINGLE) return;

    double dtr;
    int i;

    /* initialize receiver clock (white noise) */
    for (i=0;i<NSYS;i++) {
        if (opt.sateph==EPHOPT_PREC) {
            /* prec ephemeris is based gpst neglect receiver inter-system bias  */
            dtr=sol->dtr[0];
        }
        else {
            /* update receiver clock */
            dtr=i==0?sol->dtr[0]:sol->dtr[0]+sol->dtr[i];
        }
        initx(pod,CLIGHT*dtr,pod->opt.rcvvar[GIRC],GIXRCV_CLK(pod,rcv,i));
    }
    /* initialize receiver clock drift */
    if (!pod->flt.x[GIXRCV_CKR(pod,rcv)]) {
        initx(pod,sol->dtrr,pod->opt.rcvvar[GIRD],GIXRCV_CKR(pod,rcv));
    }
    /* set reference clock */
    if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
        if (pod->rcv[rcv-1].clkref) {
            for (i=0;i<NSYS;i++) {
                dtr=i==0?sol->dtr[0]:sol->dtr[0]+sol->dtr[i];
                if (pod->flt.clkref[i]==0.0) pod->flt.clkref[i]=dtr*CLIGHT;
            }
        }
    }
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, pod_t *pod)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa,rs[3];
    int i,j;
    const char *type;

    log_trace(3,"testeclipse:\n");

    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
    normv3(rsun,esun);

    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;

        matcpy(rs,pod->flt.x+GIXSAT_POS(pod,obs[i].sat),1,3);
        if ((r=norm(rs,3))<=0.0) continue;

        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;

        /* sun-earth-satellite angle */
        cosa=dot(rs,esun,3)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        ang=acos(cosa);

        /* test eclipse */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;

        log_trace(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),obs[i].sat);

        for (j=0;j<3;j++) {
            pod->flt.x[GIXSAT_POS(pod,obs[i].sat)+i]=0.0;
        }
    }
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
                     double *yaw)
{
    *yaw=yaw_nominal(beta,mu);
    return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
                   const double *rs, double *exs, double *eys)
{
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0};
    int i;

    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL);

    /* beta and orbit angle */
    matcpy(ri,rs,6,1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    cross3(ri,ri+3,n);
    cross3(rsun,n,p);
    if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||
        !normv3(p,ep)) return 0;
    beta=PI/2.0-acos(dot(esun,en,3));
    E=acos(dot(es,ep,3));
    mu=PI/2.0+(dot(es,esun,3)<=0?-E:E);
    if      (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;

    /* yaw-angle of satellite */
    if (!yaw_angle(sat,type,opt,beta,mu,&yaw)) return 0;

    /* satellite fixed x,y-vector */
    cross3(en,es,ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0;i<3;i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* phase windup model --------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char *type, int opt,
                     const double *rs, const double *rr, double *phw)
{
    double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
    double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph;
    int i;

    if (opt<=0) return 1; /* no phase windup */

    /* satellite yaw attitude model */
    if (!sat_yaw(time,sat,type,opt,rs,exs,eys)) return 0;

    /* unit vector satellite to receiver */
    for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
    if (!normv3(r,ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr,pos);
    xyz2enu(pos,E);
    exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x = north */
    eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y = west  */

    /* phase windup effect */
    cross3(ek,eys,eks);
    cross3(ek,eyr,ekr);
    for (i=0;i<3;i++) {
        ds[i]=exs[i]-ek[i]*dot(ek,exs,3)-eks[i];
        dr[i]=exr[i]-ek[i]*dot(ek,exr,3)+ekr[i];
    }
    cosp=dot(ds,dr,3)/norm(ds,3)/norm(dr,3);
    if      (cosp<-1.0) cosp=-1.0;
    else if (cosp> 1.0) cosp= 1.0;
    ph=acos(cosp)/2.0/PI;
    cross3(ds,dr,drs);
    if (dot(ek,drs,3)<0.0) ph=-ph;

    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
    return 1;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, int idx, int type,
                     const prcopt_t *opt)
{
    double fact=1.0,sinel=sin(el);

    if (type==1) fact*=opt->eratio[idx==0?0:1];
    fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);

    if (sys==SYS_GPS||sys==SYS_QZS) {
        if (idx==2) fact*=EFACT_GPS_L5; /* GPS/QZS L5 error factor */
    }
    if (opt->ionoopt==IONOOPT_IFLC) fact*=3.0;
    return SQR(fact*opt->err[1])+SQR(fact*opt->err[2]/sinel);
}
/* geometry-free phase measurement -------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav)
{
    double freq1,freq2;

    freq1=sat2freq(obs->sat,obs->code[0],nav);
    freq2=sat2freq(obs->sat,obs->code[1],nav);
    if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[1]==0.0) return 0.0;
    return (obs->L[0]/freq1-obs->L[1]/freq2)*CLIGHT;
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav)
{
    double freq1,freq2;

    freq1=sat2freq(obs->sat,obs->code[0],nav);
    freq2=sat2freq(obs->sat,obs->code[1],nav);

    if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[1]==0.0||
        obs->P[0]==0.0||obs->P[1]==0.0) return 0.0;
    return (obs->L[0]-obs->L[1])*CLIGHT/(freq1-freq2)-
           (freq1*obs->P[0]+freq2*obs->P[1])/(freq1+freq2);
}
/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                      const prcopt_t *opt, const double *dantr,
                      const double *dants, double phw, double *L, double *P,
                      double *Lc, double *Pc)
{
    double freq[NFREQ]={0},C1,C2;
    int i,sys=satsys(obs->sat,NULL);

    for (i=0;i<NFREQ;i++) {
        L[i]=P[i]=0.0;
        freq[i]=sat2freq(obs->sat,obs->code[i],nav);
        if (freq[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) continue;

        /* antenna phase center and phase windup correction */
        L[i]=obs->L[i]*CLIGHT/freq[i]-dants[i]-dantr[i]-phw*CLIGHT/freq[i];
        P[i]=obs->P[i]-dants[i]-dantr[i];

        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
        if (sys==SYS_GPS||sys==SYS_GLO) {
            if (obs->code[i]==CODE_L1C) P[i]+=nav->cbias[obs->sat-1][1];
            if (obs->code[i]==CODE_L2C) P[i]+=nav->cbias[obs->sat-1][2];
        }
    }
    /* iono-free LC */
    *Lc=*Pc=0.0;
    if (freq[0]==0.0||freq[1]==0.0) return;
    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
    C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));

    if (L[0]!=0.0&&L[1]!=0.0) *Lc=C1*L[0]+C2*L[1];
    if (P[0]!=0.0&&P[1]!=0.0) *Pc=C1*P[0]+C2*P[1];
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(pod_t *pod, int rcv, const obsd_t *obs, int n)
{
    int i,j;

    for (i=0;i<n&&i<MAXOBS;i++) {
        for (j=0;j<pod->opt.estopt.nf;j++) {
            if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&3)) continue;

            log_trace(3,"detslp_ll: slip detected sat=%2d f=%d\n",obs[i].sat,j+1);
            pod->flt.ssat[rcv][obs[i].sat-1].slip[j]=1;
        }
    }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(pod_t *pod, int rcv, const obsd_t *obs, int n, const nav_t *nav)
{
    double g0,g1;
    int i,j;

    for (i=0;i<n&&i<MAXOBS;i++) {

        if ((g1=gfmeas(obs+i,nav))==0.0) continue;

        g0=pod->flt.ssat[rcv][obs[i].sat-1].gf[0];
        pod->flt.ssat[rcv][obs[i].sat-1].gf[0]=g1;

        if (g0!=0.0&&fabs(g1-g0)>pod->opt.estopt.thresslip) {
            log_trace(3,"detslip_gf: slip detected sat=%2d gf=%8.3f->%8.3f\n",obs[i].sat,g0,g1);

            for (j=0;j<pod->opt.estopt.nf;j++) {
                pod->flt.ssat[rcv][obs[i].sat-1].slip[j]|=1;
            }
        }
    }
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
static void detslp_mw(pod_t *pod, int rcv, const obsd_t *obs, int n, const nav_t *nav)
{
    double w0,w1;
    int i,j;

    for (i=0;i<n&&i<MAXOBS;i++) {
        if ((w1=mwmeas(obs+i,nav))==0.0) continue;

        w0=pod->flt.ssat[rcv][obs[i].sat-1].mw[0];
        pod->flt.ssat[rcv][obs[i].sat-1].mw[0]=w1;

        if (w0!=0.0&&fabs(w1-w0)>THRES_MW_JUMP) {
            log_trace(3,"detslip_mw: slip detected sat=%2d mw=%8.3f->%8.3f\n",obs[i].sat,w0,w1);

            for (j=0;j<pod->opt.estopt.nf;j++) {
                pod->flt.ssat[rcv][obs[i].sat-1].slip[j]|=1;
            }
        }
    }
}
/* temporal update of position -----------------------------------------------*/
static void udpos_pod(pod_t *pod, int rcv, gtime_t tutc)
{
    initrcvpos(pod,rcv,tutc);
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_pod(pod_t *pod, int rcv, const obsd_t *obs, int nobs)
{
    initrcvclk(pod,rcv,obs,nobs);
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_pod(pod_t *pod, gtime_t tutc, int rcv)
{
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    double dt=pod->flt.dt;
    int i=GIXRCV_TRP(pod,rcv),j;

    if (pod->flt.x[i]==0.0) {
        ecef2pos(pod->flt.x+GIXRCV_POS(pod,rcv),pos);
        ztd=sbstropcorr(pod->flt.sols[rcv].time,pos,azel,&var);
        initx(pod,ztd,var,i);

        if (pod->opt.estopt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) initx(pod,1E-6,VAR_GRA,j);
        }
    }
    else {
        pod->flt.P[i+i*pod->flt.nx]+=SQR(pod->opt.estopt.prn[2])*fabs(dt);

        if (pod->opt.estopt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) {
                pod->flt.P[j+j*pod->flt.nx]+=SQR(pod->opt.estopt.prn[2]*0.1)*fabs(dt);
            }
        }
    }
}
/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_pod(pod_t *pod)
{
    int i=GIXDCB(pod);

    if (pod->flt.x[i]==0.0) {
        initx(pod,1E-6,VAR_DCB,i);
    }
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_pod(pod_t *pod, gtime_t tutc, int rcv, const obsd_t *obs, int n, const nav_t *nav)
{
    double L[NFREQ],P[NFREQ],Lc,Pc,bias[MAXOBS],offset=0.0,pos[3]={0};
    double freq1,freq2,ion,dantr[NFREQ]={0},dants[NFREQ]={0},dt=pod->flt.dt;
    int i,j,k,f,sat,slip[MAXOBS]={0},clk_jump=0;

    /* handle day-boundary clock jump */
    if (pod->opt.estopt.posopt[5]) {
        clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;
    }
    for (i=0;i<MAXSAT;i++) {
        for (j=0;j<pod->opt.estopt.nf;j++) pod->flt.ssat[rcv][i].slip[j]=0;
    }
    /* detect cycle slip by LLI */
    detslp_ll(pod,rcv,obs,n);

    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(pod,rcv,obs,n,nav);

    /* detect slip by Melbourne-Wubbena linear combination jump */
    detslp_mw(pod,rcv,obs,n,nav);

    for (f=0;f<NF(&pod->opt.estopt);f++) {

        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) {
            if (clk_jump||++pod->flt.ssat[rcv][i].outc[f]>(uint32_t)pod->opt.estopt.maxout) {
                if (pod->flt.x[GIXSAT_BIAS(pod,rcv,i+1)]) {
                    initx(pod,0.0,0.0,GIXSAT_BIAS(pod,rcv,i+1));
                }
            }
        }
        for (i=k=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=GIXSAT_BIAS(pod,rcv,sat);
            corr_meas(obs+i,nav,pod->flt.ssat[rcv][sat-1].azel,&pod->opt.estopt,dantr,dants,0.0,L,P,&Lc,&Pc);

            bias[i]=Lc-Pc;
            slip[i]=pod->flt.ssat[rcv][sat-1].slip[0]||pod->flt.ssat[rcv][sat-1].slip[1];

            if (pod->flt.x[j]==0.0||slip[i]||bias[i]==0.0) continue;

            offset+=bias[i]-pod->flt.x[j];
            k++;
        }
        /* correct phase-code jump to ensure phase-code coherency */
        if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) {
            for (i=0;i<MAXSAT;i++) {
                j=GIXSAT_BIAS(pod,rcv,i+1);
                if (pod->flt.x[j]!=0.0) pod->flt.x[j]+=offset/k;
            }
            log_trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",time_str(obs[0].time,0),k,offset/k/CLIGHT);
        }
        for (i=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=GIXSAT_BIAS(pod,rcv,sat);

            pod->flt.P[j+j*pod->flt.nx]+=SQR(pod->opt.estopt.prn[0])*fabs(dt);

            if (bias[i]==0.0||(pod->flt.x[j]!=0.0&&!slip[i])) continue;

            /* reinitialize phase-bias if detecting cycle slip */
            initx(pod,bias[i],VAR_BIAS,GIXSAT_BIAS(pod,rcv,sat));
        }
    }
}
/* temporal update of satellite position/velocity/clock/SRP---------------------*/
static void udsatorbclksrp_pod(pod_t *pod, gtime_t tutc)
{
    satorbit_t *orbit;
    int i,j;

    for (i=0;i<MAXSAT;i++) {

        /* unselected sat sys */
        if (!(satsys(i+1,NULL)&pod->opt.estopt.navsys)) continue;

        initsatorbclksrp_nav(pod,i+1,tutc);

        if (pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS) continue;
        orbit=&pod->orbits[i];

        if (orbit->sat<=0) continue;

        double dt=timediff(tutc,orbit->tutc);
        if (fabs(dt)<1E-8) continue;

        /* update satellite position/velocity */
        satorbit(pod->opt.udorbitint,&pod->orbits[i],dt);

        /* satellite clock transition matrix */
        if (pod->opt.satclkmode==0) {
            orbit->F[GICR+GICR*GNX]=1.0;
            orbit->F[GIC +GICR*GNX]=0.0;
        }
        /* update covariance matrix of satellite position/velocity */
        int ix,jx,sat=i+1;
        double *Ps=zeros(GNX,pod->flt.nx);
        double *Pk=zeros(GNX,pod->flt.nx);

        for (ix=0;ix<GNX;ix++) {
            for (jx=0;jx<pod->flt.nx;jx++) Ps[ix+jx*GNX]=pod->flt.P[ix+GIXSAT(pod,sat)+jx*pod->flt.nx];
        }
        matmul("NN",GNX,pod->flt.nx,GNX,1.0,orbit->F,Ps,0.0,Pk);

        for (ix=0;ix<GNX;ix++) {
            for (jx=0;jx<pod->flt.nx;jx++) pod->flt.P[ix+GIXSAT(pod,sat)+jx*pod->flt.nx]=Pk[ix+jx*GNX];
        }
        for (ix=0;ix<pod->flt.nx;ix++) {
            for (jx=0;jx<GNX;jx++) Ps[ix+jx*pod->flt.nx]=pod->flt.P[ix+(jx+GIXSAT(pod,sat))*pod->flt.nx];
        }
        matmul("NT",pod->flt.nx,GNX,GNX,1.0,Ps,orbit->F,0.0,Pk);

        for (ix=0;ix<pod->flt.nx;ix++) {
            for (jx=0;jx<GNX;jx++) pod->flt.P[ix+(jx+GIXSAT(pod,sat))*pod->flt.nx]=Pk[ix+jx*pod->flt.nx];
        }
        for (ix=0;ix<GNX;ix++) {
            if (pod->flt.x[ix+GIXSAT(pod,sat)]) {

                if (ix==GIC) {
                    double ddclk=0.0,dts[2];

                    /* current satellite clock */
                    ephclk(utc2gpst(orbit->tutc),utc2gpst(orbit->tutc),sat,&pod->nav,dts);

                    /* satellite clock variance */
                    ddclk=3.0*(dts[0]*CLIGHT-pod->flt.x[GIXSAT_CLK(pod,sat)]);

                    pod->flt.P[GIXSAT_CLK(pod,sat)+GIXSAT_CLK(pod,sat)*pod->flt.nx]+=SQR(ddclk?ddclk:pod->opt.satprn[ix]);
                }
                else {
                    pod->flt.P[ix+GIXSAT(pod,sat)+(ix+GIXSAT(pod,sat))*pod->flt.nx]+=SQR(pod->opt.satprn[ix])*fabs(dt);
                }
            }
        }
        free(Ps); free(Pk);

        /* update POD states */
        for (ix=0;ix<GNX-GNC;ix++) {
            pod->flt.x[GIXSAT(pod,sat)+ix]=orbit->x[ix];
        }
        /* update satellite clock */
        if (pod->opt.satclkmode==0) {
            double dts[2],relclk;

            /* white noise */
            ephclk(utc2gpst(orbit->tutc),utc2gpst(orbit->tutc),sat,&pod->nav,dts);
            relclk=relcorr(orbit->x,NULL,0);
            initx(pod,dts[0]*CLIGHT+relclk,VAR_CLK,GIXSAT_CLK(pod,sat));
        }
        else {
            double corrtime=3600.0*6.0;

            /* 1st order gauss-marcov */
            int jc=GIXSAT_CLK(pod,sat);
            int jr=GIXSAT_CKR(pod,sat);
            pod->flt.x[jc]+=(corrtime-corrtime*exp(-dt/corrtime))*pod->flt.x[jr];
            pod->flt.x[jr]*=exp(-dt/corrtime);
        }
    }
}
/* temporal update of states --------------------------------------------------*/
static void udstate_pod(pod_t *pod, gtime_t tutc, int rcv, const obsd_t *obs, int n)
{
    log_trace(3,"udstate_pod: n=%d\n",n);

    /* temporal update of clock */
    udclk_pod(pod,rcv,obs,n);

    /* temporal update of position */
    udpos_pod(pod,rcv,gpst2utc(obs[0].time));

    /* temporal update of tropospheric parameters */
    if (pod->opt.estopt.tropopt==TROPOPT_EST||pod->opt.estopt.tropopt==TROPOPT_ESTG) {
        udtrop_pod(pod,tutc,rcv);
    }
    /* temporal update of L5-receiver-dcb parameters */
    if (pod->opt.estopt.nf>=3) {
        uddcb_pod(pod);
    }
    /* temporal update of phase-bias */
    udbias_pod(pod,tutc,rcv,obs,n,&pod->nav);
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
                      double *dant)
{
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;

    for (i=0;i<3;i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return;

    cosa=dot(eu,ez,3);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);

    antmodel_s(pcv,nadir,dant);
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double *pos,
                              const double *azel, const double *x, double *dtdx,
                              double *var)
{
    const double zazel[]={0.0,PI/2.0};
    double zhd,m_h,m_w,cotz,grad_n,grad_e;

    /* zenith hydrostatic delay */
    zhd=tropmodel(time,pos,zazel,0.0);

    /* mapping function */
    m_h=tropmapf(time,pos,azel,&m_w);

    if (azel[1]>0.0) {

        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[1]+grad_e*x[2];
        dtdx[1]=grad_n*(x[0]-zhd);
        dtdx[2]=grad_e*(x[0]-zhd);
    }
    dtdx[0]=m_w;
    *var=SQR(0.01);
    return m_h*zhd+m_w*(x[0]-zhd);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(pod_t *pod, int rcv, gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, const double *x, double *dtdx,
                      const nav_t *nav, double *dtrp, double *var)
{
    double trp[3]={0};

    if (opt->tropopt==TROPOPT_SAAS) {
        *dtrp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS);
        return 1;
    }
    if (opt->tropopt==TROPOPT_SBAS) {
        *dtrp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
        matcpy(trp,x+GIXRCV_TRP(pod,rcv),opt->tropopt==TROPOPT_EST?1:3,1);
        *dtrp=trop_model_prec(time,pos,azel,trp,dtdx,var);
        return 1;
    }
    return 0;
}

/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, int sat, const double *x,
                      const nav_t *nav, double *dion, double *var)
{
    if (opt->ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,dion,var);
    }
    if (opt->ionoopt==IONOOPT_BRDC) {
        *dion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*dion*ERR_BRDCI);
        return 1;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        *dion=*var=0.0;
        return 1;
    }
    return 0;
}
/* phase and code residuals --------------------------------------------------*/
static int ppp_res(int post, int rcv, const obsd_t *obs, int n,
                   const double *dr, const nav_t *nav,
                   const double *x, const double *Pp, pod_t *pod, double *v, double *H, double *var, double *azel, int *exc, int *vflag)
{
    prcopt_t *opt=&pod->opt.estopt;
    double y,r,cdtr,cdts,bias,relcor=0.0,C=0.0,rr[3],rr_eci[3],pos[3],e[3],dtdx[3],L[NFREQ],P[NFREQ],Lc,Pc;
    double dtrp=0.0,vart=0.0,dcb=0.0,freq;
    double dantr[NFREQ]={0},dants[NFREQ]={0},satdant[3],satdant_eci[3],azelt[2*MAXOBS]={0};
    double ve[MAXOBS*2*NFREQ]={0},vmax=0;
    char str[32];
    int ne=0,obsi[MAXOBS*2*NFREQ]={0},frqi[MAXOBS*2*NFREQ],maxobs,maxfrq,rej;
    int i,j,k,sat,sys,nv=0,nx=pod->flt.nx,stat=1;

    time2str(obs[0].time,str,2);

    for (i=0;i<MAXSAT;i++) {
        for (j=0;j<opt->nf;j++) pod->flt.ssat[rcv][i].vsat[j]=0;
        for (j=0;j<opt->nf;j++) pod->flt.ssat[rcv][i].resp[j]=0;
        for (j=0;j<opt->nf;j++) pod->flt.ssat[rcv][i].resc[j]=0;
    }
    for (i=0;i<3;i++) rr[i]=x[GIXRCV_POS(pod,rcv)+i]+dr[i];
    ecef2pos(rr,pos);
    if (norm(rr,3)<=0.0) return 0;

    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;

        /* check DOPS of satellite */
        if (pod->sol.dops[sat-1]>MAXDOPS) continue;

        gtime_t ts;
        satorbit_t *orb=&pod->orbits[sat-1];
        double pr,U[9],erpv[6],tt,rs_ecef[6],dion=0.0,vari=0.0,dt,rs[6],dts[2],vars,rs_eci[3],rng_eci,e_eci[3],dr_eci[3],E[9];
        int svh;

        for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;
        if (pr==0.0) continue;

        /* transmission time by satellite clock */
        ts=timeadd(obs[i].time,-pr/CLIGHT);

        /* satellite clock bias by broadcast ephemeris */
        if (!ephclk(ts,obs[i].time,sat,nav,&dt)) continue;
        ts=timeadd(ts,-dt);

        /* satellite position and clock at transmission time */
        if (!satpos(ts,obs[i].time,obs[i].sat,opt->sateph,nav,rs,dts,&vars,&svh)) continue;
        ts=timeadd(obs[i].time,-pr/CLIGHT-x[GIXSAT_CLK(pod,sat)]/CLIGHT);
        tt=timediff(gpst2utc(ts),orb->tutc);

        /* transform matrix of ECI to ECEF */
        eci2ecef2(gpst2utc(ts),&orb->fmdl.erp,U,NULL,NULL);

        /* satellite orbit propagator process */
        double xs[GNX]={0},Fs[GNX*GNX]={0};
        if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
#if ORBINTJ2
            double xs_[6],Fs_[6*6];
            orbintj2_U(gpst2utc(ts),U,tt,orb->x,xs_,Fs_);

            for (j=0;j<GNX;j++) Fs[j+j*GNX]=1.0;
            for (j=0;j<6;j++) {
                for (k=0;k<6;k++) xs[k]=xs_[k];
                for (k=0;k<6;k++) Fs[j+k*GNX]=Fs_[j+k*6];
            }
#else
            satorbite(ORBITINT_ODEINTV2_RKF78,orb,tt,xs,Fs);
#endif
        }
        /* satellite position in ECEF */
        satposecef(gpst2utc(ts),&orb->fmdl.erp,xs+GIP,rs_ecef);

        /* satellite antenna phase center offset in ECI */
        satantoff(obs[i].time,rs,sat,nav,satdant);
        matmul("TN",3,1,3,1.0,U,satdant,0.0,satdant_eci);
        for (j=0;j<3;j++) rs_eci[j]=xs[GIP+j]+satdant_eci[j];

        /* receiver position in ECI */
        switch (sys) {
            case SYS_GLO: k=1; break;
            case SYS_GAL: k=2; break;
            case SYS_CMP: k=3; break;
            case SYS_IRN: k=4; break;
            default:      k=0; break;
        }
        double Tr[9];
        rcvposeci(gpst2utc(obs[i].time),rr,&orb->fmdl.erp,-x[GIXRCV_CLK(pod,rcv,k)]/CLIGHT,rr_eci,Tr);

        /* unit vector of satellite to receiver */
        for (j=0;j<3;j++) dr_eci[j]=rs_eci[j]-rr_eci[j];
        rng_eci=norm(dr_eci,3);
        for (j=0;j<3;j++) e_eci[j]=dr_eci[j]/rng_eci;

        /* eci->radial/alongtrk/crosstrk */
        ecef2satf(xs,E);
        matmul("NN",3,1,3,1.0,E,e_eci,0.0,e);
        azel[0+2*i]=atan2(-e[2],e[1]);
        azel[1+2*i]=asin(e[0]/norm(e,3));

        /* geometric distance in ECEF */
        if ((r=geodist(rs,rr,e))<=0.0||satazel(pos,e,azelt+i*2)<opt->elmin) {
            exc[i]=1;
            continue;
        }
        if (!(sys=satsys(sat,NULL))||!pod->flt.ssat[rcv][sat-1].vs||exc[i]||
            satexclude(obs[i].sat,pod->flt.ssat[rcv][sat-1].var,pod->flt.ssat[rcv][sat-1].svh,opt)) {
            exc[i]=1;
            continue;
        }
        /* tropospheric and ionospheric model */
        if (!model_trop(pod,rcv,obs[i].time,pos,azelt+i*2,opt,x,dtdx,nav,&dtrp,&vart)||
            !model_iono(obs[i].time,pos,azelt+i*2,opt,sat,x,nav,&dion,&vari)) continue;

        /* satellite and receiver antenna model */
        if (opt->posopt[0]) satantpcv(rs,rr,nav->pcvs+sat-1,dants);
        antmodel(opt->pcvr,opt->antdel[0],azelt+i*2,opt->posopt[1],dantr);

        /* phase windup model */
        if (!model_phw(obs[i].time,sat,nav->pcvs[sat-1].type,opt->posopt[2]?2:0,rs,rr,&pod->flt.ssat[rcv][sat-1].phw)) {
            continue;
        }
        /* corrected phase and code measurements */
        corr_meas(obs+i,nav,azelt+i*2,&pod->opt.estopt,dantr,dants,pod->flt.ssat[rcv][sat-1].phw,L,P,&Lc,&Pc);

        /* partial derivatives of range by satellite state */
        double ee[GNX]={0},es[GNX]={0},er[3];
        matcpy(ee,e_eci,1,3);
        matmul("NN",1,GNX,GNX, 1.0,ee,Fs,0.0,es);
        matmul("NT",1,GNP,GNP,-1.0,ee,Tr,0.0,er);

        /* stack phase and code residuals {L1,P1,L2,P2,...} */
        for (j=0;j<2*NF(opt);j++) {

            if (pod->opt.onlyusecode==1) if (j%2==0) continue;

            if (pod->flt.ssat[rcv][sat-1].rejc[j%NFREQ]) continue;

            dcb=bias=cdtr=cdts=relcor=0.0;

            if (opt->ionoopt==IONOOPT_IFLC) {
                if ((y=j%2==0?Lc:Pc)==0.0) continue;
            }
            else {
                if ((y=j%2==0?L[j/2]:P[j/2])==0.0) continue;

                if ((freq=sat2freq(sat,obs[i].code[j/2],nav))==0.0) continue;
                C=SQR(FREQ1/freq)*ionmapf(pos,azelt+i*2)*(j%2==0?-1.0:1.0);
            }
            memset(&H[nx*nv],0,sizeof(double)*nx);

            /* satellite/receiver position */
            if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
                for (k=0;k<GNX;k++) H[GIXSAT(pod,sat)+k+nx*nv]=es[k];
                for (k=0;k<GNP;k++) H[GIXRCV(pod,rcv)+k+nx*nv]=er[k];
            }
            else if (pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS) {
                for (k=0;k<GNP;k++) H[GIXRCV_POS(pod,rcv)+k+nx*nv]=-e[k];
            }
            else continue;

            /* troposphere delay */
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                    H[GIXRCV_TRP(pod,rcv)+k+nx*nv]=dtdx[k];
                }
            }
            /* L5-receiver-dcb */
            if (j/2==2&&j%2==1) {
                dcb+=pod->flt.x[GIXDCB(pod)];
                H[GIXDCB(pod)+nx*nv]=1.0;
            }
            /* satellite-receiver phase bias */
            if (j%2==0) {
                if ((bias=x[GIXSAT_BIAS(pod,rcv,sat)])==0.0) continue;
                H[GIXSAT_BIAS(pod,rcv,sat)+nx*nv]=1.0;
            }
            /* relativity correction */
            if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) relcor=relcorr(&x[GIXSAT_POS(pod,sat)],rr_eci,0);

            /* receiver clock */
            switch (sys) {
                case SYS_GLO: k=1; break;
                case SYS_GAL: k=2; break;
                case SYS_CMP: k=3; break;
                case SYS_IRN: k=4; break;
                default:      k=0; break;
            }
            cdtr=x[GIXRCV_CLK(pod,rcv,k)]-pod->flt.clkref[k];
            H[GIXRCV_CLK(pod,rcv,k)+nx*nv]=1.0;

            /* set reference clock */
            if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
                if (pod->rcv[rcv-1].clkref==1) {
                    H[GIXRCV_CLK(pod,rcv,k)+nx*nv]=0.0;
                }
            }
            /* satellite clock bias */
            cdts=x[GIXSAT_CLK(pod,sat)]-pod->flt.clkref[k];
            if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) H[GIXSAT_CLK(pod,sat)+nx*nv]=-1.0;

            /* phase/code measurement residual */
            if (pod->opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) v[nv]=y-(rng_eci+cdtr-cdts+dtrp+dcb+bias+relcor+C*dion);
            else v[nv]=y-(r+cdtr-cdts+dtrp+dcb+bias+relcor+C*dion);

            if (j%2==0) pod->flt.ssat[rcv][sat-1].resc[j/2]=v[nv];
            else        pod->flt.ssat[rcv][sat-1].resp[j/2]=v[nv];

            /* phase/code measurement variance */
            var[nv]=varerr(obs[i].sat,sys,azelt[1+i*2],j/2,j%2,opt)+vart;
            if (pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS) var[nv]+=vars;

            if (sys==SYS_GLO&&j%2==1) var[nv]+=VAR_GLO_IFB;

            char prn[8];
            satno2id(sat,prn);
            log_trace(3,"%4s sat=%2d %s%d res=%9.4f sig=%9.4f el=%4.1f cdtr=%10.3lf cdts=%10.3lf\n",prn,sat,j%2?"P":"L",j/2+1,v[nv],sqrt(var[nv]),azelt[1+i*2]*R2D,cdtr,cdts);

            /* reject satellite by pre-fit residuals */
            if (fabs(v[nv])>SQRT(VAR_CLK)*THRES_REJECT) {
                log_trace(2,"outlier (%d) rejected sat=%2d %s%d res=%9.4f el=%4.1f\n",post,sat,j%2?"P":"L",j/2+1,v[nv],azelt[1+i*2]*R2D);
                exc[i]=1;
                pod->flt.ssat[rcv][sat-1].rejc[j%NFREQ]++;
                continue;
            }
            if (j%2==0) pod->flt.ssat[rcv][sat-1].vsat[j/2]=1;
            vflag[nv]=(sat<<16)|(rcv<<8)|((j%NF(opt)?1:0)<<4)|nv;
            nv++;
        }
    }
    return nv;
}
/* detect outliers through mad------------------------------------------------*/
static int maddetoutl(const double *v, int nv, int *outl_ind, double thres)
{
    int i,outl_n=0;
    double m,med;

    if (nv<3) return 0;

    med=median(v,nv); m=mad(v,nv);
    for (i=0;i<nv;i++) {
        if (fabs((v[i]-med)/(1.4826*m))<thres) continue;
        if (outl_ind) outl_ind[outl_n++]=i;
    }
    return outl_n;
}
/* detect measurement outliers ------------------------------------------------*/
static int measdetoutl(const double *v, int nv, int *outl_ind, int *flag, double thresmad)
{
    int i,outl_n=0;
    for (i=0;flag&&i<nv;i++) flag[i]=1;
    if (thresmad) {
        outl_n+=maddetoutl(v,nv,outl_ind+outl_n,thresmad);
    }
    for (i=0;flag&&i<outl_n;i++) flag[outl_ind[i]]=0;
    return outl_n;
}
/* POD measurement for given receiver ----------------------------------------*/
static int podmeasresrcv(pod_t *pod, int post, gtime_t tutc, int rcv, int rind, const obsd_t *obs, int n, const double *xp, const double *Pp, double *v, double *H, double *var, int *exc)
{
    nav_t *nav=&pod->nav;
    const prcopt_t *opt=&pod->opt.estopt;
    double azel[MAXOBS*2]={0},dr[3]={0},std[3];
    char str[32];
    int i,j,nv,info,nx=pod->flt.nx,vflag[MAXOBS*2];

    time2str(obs[0].time,str,2);
    log_trace(3,"ppposprc   : time=%s nx=%d n=%d rcv=%2d\n",str,pod->flt.nx,n,rcv);

    /* exclude measurements of eclipsing satellite (block IIA) */
    if (pod->opt.estopt.posopt[3]) {
        testeclipse(obs,n,nav,pod);
    }
    /* earth tides correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),xp+GIXRCV_POS(pod,rcv),opt->tidecorr==1?1:7,&nav->erp,opt->odisp[0],dr);
    }
    /* prefit residuals */
    if (!(nv=ppp_res(post,rcv,obs,n,dr,nav,xp,Pp,pod,v,H,var,azel,exc,vflag))) {
        log_trace(2,"%s no valid obs data\n",str);
    }
    return nv;
}
/* update POD solution status-------------------------------------------------*/
static void udpossolstat(pod_t *pod, const podobs_t *obss, int stat)
{
    pod->sol.stat=stat;
    pod->sol.tutc=obss->tutc;
    pod->flt.time=obss->tutc;

    int i,j,k;

    for (i=0;i<obss->nrcv;i++) {
        const prcopt_t *opt=&pod->opt.estopt;

        int rcv=obss->ircv[i];
        int n=obss->nobs[i];
        const obsd_t *obs=obss->obs[i];

        for (j=0;j<GNRX;j++) {
            pod->sol.xr[rcv][j]=pod->flt.x[GIXRCV(pod,rcv)+j];
            pod->sol.Qr[rcv][j]=pod->flt.P[GIXRCV(pod,rcv)+j+(GIXSAT(pod,rcv)+j)*pod->flt.nx];
        }
        for (j=0;j<n;j++) {
            for (k=0;k<GNX;k++) pod->sol.xs[obs[j].sat-1][k]=pod->flt.x[GIXSAT(pod,obs[j].sat)+k];
            for (k=0;k<GNX;k++) pod->sol.Qs[obs[j].sat-1][k]=pod->flt.P[GIXSAT(pod,obs[j].sat)+k+(GIXSAT(pod,obs[j].sat)+k)*pod->flt.nx];
        }
        /* test # of valid satellites */
        pod->flt.sols[rcv].ns=0;
        for (j=0;j<n&&j<MAXOBS;j++) {
            for (k=0;k<opt->nf;k++) {
                if (!pod->flt.ssat[rcv][obs[j].sat-1].vsat[k]) continue;
                pod->flt.ssat[rcv][obs[j].sat-1].lock[k]++;
                pod->flt.ssat[rcv][obs[j].sat-1].outc[k]=0;
                if (k==0) pod->flt.sols[rcv].ns++;
            }
        }
        for (j=0;j<n&&j<MAXOBS;j++) {
            for (k=0;k<opt->nf;k++) {
                pod->flt.ssat[rcv][obs[j].sat-1].snr[k]=obs->SNR[k];
            }
        }
        for (j=0;j<MAXSAT;j++) {
            for (k=0;k<opt->nf;k++) {
                if (pod->flt.ssat[rcv][j].slip[k]&3) pod->flt.ssat[rcv][j].slipc[k]++;
            }
        }
        pod->flt.sols[rcv].stat=pod->flt.sols[rcv].ns<5?SOLQ_NONE:stat;

        for (k=0;k<3;k++) {
            pod->flt.sols[rcv].rr[k]=pod->flt.x[GIXRCV_POS(pod,rcv)+k];
            pod->flt.sols[rcv].qr[k]=(float)pod->flt.P[GIXRCV_POS(pod,rcv)+k+(GIXRCV_POS(pod,rcv)+k)*pod->flt.nx];
        }
        pod->flt.sols[rcv].qr[3]=(float)pod->flt.P[GIXRCV_POS(pod,rcv)+1];
        pod->flt.sols[rcv].qr[4]=(float)pod->flt.P[GIXRCV_POS(pod,rcv)+2+pod->flt.nx];
        pod->flt.sols[rcv].qr[5]=(float)pod->flt.P[GIXRCV_POS(pod,rcv)+2];

        for (j=0;j<n&&j<MAXOBS;j++) {
            satorbit_t *orb=&pod->orbits[obs[j].sat-1];
            for (k=0;k<GNP;k++) orb->x[GIP+k]=pod->flt.x[GIXSAT_POS(pod,obs[j].sat)+k];
            for (k=0;k<GNV;k++) orb->x[GIV+k]=pod->flt.x[GIXSAT_VEL(pod,obs[j].sat)+k];
            for (k=0;k<GNC;k++) orb->x[GIC+k]=pod->flt.x[GIXSAT_CLK(pod,obs[j].sat)+k];
            for (k=0;k<GNS;k++) orb->x[GIS+k]=pod->flt.x[GIXSAT_SRP(pod,obs[j].sat,k)];

            /* srp parameters: D0,Dc,Ds,Dc2,Ds2,Dc4,Ds4,Y0,Yc,Ys,B0,Bc,Bs */
            orb->fmdl.sp.D0 =orb->x[GIS_D0 ];
            orb->fmdl.sp.Dc =orb->x[GIS_DC ];
            orb->fmdl.sp.Ds =orb->x[GIS_DS ];
            orb->fmdl.sp.Dc2=orb->x[GIS_DC2];
            orb->fmdl.sp.Ds2=orb->x[GIS_DS2];
            orb->fmdl.sp.Dc4=orb->x[GIS_DC4];
            orb->fmdl.sp.Ds4=orb->x[GIS_DS4];
            orb->fmdl.sp.Y0 =orb->x[GIS_Y0 ];
            orb->fmdl.sp.Yc =orb->x[GIS_YC ];
            orb->fmdl.sp.Ys =orb->x[GIS_YS ];
            orb->fmdl.sp.B0 =orb->x[GIS_B0 ];
            orb->fmdl.sp.Bc =orb->x[GIS_BC ];
            orb->fmdl.sp.Bs =orb->x[GIS_BS ];
        }
    }
}
/* constraint of satellite position/velocity/clock using eph------------------*/
static int constsatposveleph(pod_t *pod, int post, const podobs_t *obss, const double *xp, double *v, double *H, double *var, int *outlflag)
{
    int i,j,k,nv,udflg[MAXSAT]={0},nx=pod->flt.nx;
    double varpos_const=SQR(5.0);
    double varvel_const=SQR(1E-3);

    if (pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS) return 0;

    for (nv=i=0;i<obss->nrcv;i++) {
        for (j=0;j<obss->nobs[i];j++) {

            double rs[6],dts[2],vars,rs_eci[6],dants[3],dants_eci[3],U[9];
            int sat=obss->obs[i][j].sat,svh;

            if (udflg[sat-1]) continue;
            satorbit_t *orb=&pod->orbits[sat-1];
            gtime_t gpst=utc2gpst(orb->tutc);

            /* satellite position/velocity */
            if (!satpos(gpst,gpst,sat,pod->opt.estopt.sateph,&pod->nav,rs,dts,&vars,&svh)) continue;

            /* transform matrix of ECI to ECEF */
            satposeci(orb->tutc,&orb->fmdl.erp,rs,rs_eci);

            /* transform matrix of ECI to ECEF */
            eci2ecef2(gpst2utc(gpst),&orb->fmdl.erp,U,NULL,NULL);

            /* satellite antenna phase center offset in ECI */
            satantoff(gpst,rs,sat,&pod->nav,dants);
            matmul("TN",3,1,3,1.0,U,dants,0.0,dants_eci);
            for (k=0;k<3;k++) rs_eci[k]=rs_eci[k]-dants_eci[k];

            /* satellite position/velocity measurement residuals */
            for (k=0;k<3;k++) {
                memset(&H[nx*nv],0,sizeof(double)*nx);
                v[nv]=rs_eci[k]-xp[GIXSAT_POS(pod,sat)+k];
                H[GIXSAT_POS(pod,sat)+k+nv*nx]=1.0;
                var[nv]=varpos_const;
                nv++;
            }
            for (k=0;k<3;k++) {
                memset(&H[nx*nv],0,sizeof(double)*nx);
                v[nv]=rs_eci[k+3]-xp[GIXSAT_VEL(pod,sat)+k];
                H[GIXSAT_VEL(pod,sat)+k+nv*nx]=1.0;
                var[nv]=varvel_const;
                nv++;
            }
            udflg[sat-1]=1;
        }
    }
    return nv;
}
static int constsatclk(pod_t *pod, int post, const podobs_t *obss, const double *xp, double *v, double *H, double *var, int *outlflag)
{
    int i,j,k,nv,udflg[MAXSAT]={0},nx=pod->flt.nx;
    double varclk_const=SQR(5.0);
    double varclkr_const=SQR(1E-5);

    if (pod->opt.estmode==GPOD_FIXSATPOS_ESTRCVPOS) return 0;

    for (nv=i=0;i<obss->nrcv;i++) {
        for (j=0;j<obss->nobs[i];j++) {

            double rs[6],dts[2],vars,rs_eci[6],dants[3],dants_eci[3],U[9];
            int sat=obss->obs[i][j].sat,svh;

            if (udflg[sat-1]) continue;
            satorbit_t *orb=&pod->orbits[sat-1];
            gtime_t gpst=utc2gpst(orb->tutc);

            /* satellite clock */
            if (!satpos(gpst,gpst,sat,pod->opt.estopt.sateph,&pod->nav,rs,dts,&vars,&svh)) continue;

            /* satellite clock measurement residuals */
            memset(&H[nx*nv],0,sizeof(double)*nx);
            v[nv]=dts[0]*CLIGHT+relcorr(&xp[GIXSAT_POS(pod,sat)],NULL,0)-xp[GIXSAT_CLK(pod,sat)];
            H[GIXSAT_CLK(pod,sat)+nv*nx]=1.0;
            var[nv]=varclk_const;
            nv++;

            memset(&H[nx*nv],0,sizeof(double)*nx);
            v[nv]=dts[1]*CLIGHT-xp[GIXSAT_CKR(pod,sat)];
            H[GIXSAT_CKR(pod,sat)+nv*nx]=1.0;
            var[nv]=varclkr_const;
            nv++;
            udflg[sat-1]=1;
        }
    }
    return nv;
}
/* POD measurement residuals--------------------------------------------------*/
static int podmeasres(pod_t *pod, int post, const podobs_t *obss, const double *xp, const double *Pp, double *v, double *H, double *var, int *exc)
{
    pod_opt_t *opt=&pod->opt;
    int i,j,stat=0,nv,nx=pod->flt.nx;

    for (nv=i=0;i<obss->nrcv;i++) {
        sta_t *sta=&pod->rcv[obss->ircv[i]-1].sta;

        /* set antenna parameters */
        setpcv(utc2gpst(obss->tutc),&opt->estopt,&pod->nav,&pod->pcvss,&pod->pcvsr,sta);

        /* read ocean tide loading parameters */
        readotl(&opt->estopt,opt->blq,sta);

        /* satellite orbit estimator of PPP */
        nv+=podmeasresrcv(pod,post,obss->tutc,obss->ircv[i],i,obss->obs[i],obss->nobs[i],xp,Pp,v+nv,H+nv*nx,var+nv,exc+i*MAXOBS);
    }
#if CONSTEPH_POSVEL
    /* constraint of satellite position/velocity using eph */
    nv+=constsatposveleph(pod,post,obss,xp,v+nv,H+nv*nx,var+nv,outlflag);
#endif

#if CONSTEPH_CLK
    /* constraint of satellite clock using eph */
    nv+=constsatclk(pod,post,obss,xp,v+nv,H+nv*nx,var+nv,outlflag);
#endif
    return nv;
}
/* adjust observation data----------------------------------------------------*/
static void adjobss(pod_t *pod, podobs_t *obss)
{
    int i,j,k,flag;

    for (i=0;i<obss->nrcv;i++) {
        for (flag=j=0;j<GNRCVS;j++) {
            if (strcmp(pod->rcv[j].name,obss->name[i])!=0) continue;
            obss->ircv[i]=pod->rcv[j].id;
            flag=1;
            break;
        }
        if (flag==0) {
            obss->nobs[i]=0;
            obss->ircv[i]=-1;
        }
    }
    for (k=i=0;i<obss->nrcv;i++) {
        if (obss->nobs[i]==0||obss->ircv[i]==-1) continue;
        obss->ircv[k]=obss->ircv[i];
        obss->nobs[k]=obss->nobs[i];
        memcpy(obss->obs[k],obss->obs[i],sizeof(obss->obs[i]));
        strcpy(obss->name[k],obss->name[i]);
        k++;
    }
    obss->nrcv=k;
}
/* detect measurement outliers------------------------------------------------*/
static int poddetmeasoutl(pod_t *pod, const podobs_t *obss)
{
    int i,j,nv,*outli=imat(obss->nrcv,2),*vind=imat(obss->nrcv,2),m=0;
    double *v=mat(obss->nrcv,2);

    for (i=0;i<MAXSAT;i++) {
        for (nv=j=0;j<obss->nrcv;j++) {
            if (obss->ircv[j]<0) continue;
            if (pod->flt.ssat[obss->ircv[j]][i].vsat[0]==0) continue;
            v[nv]=pod->flt.ssat[obss->ircv[j]][i].resc[0]; vind[nv++]=obss->ircv[j];
            v[nv]=pod->flt.ssat[obss->ircv[j]][i].resp[0]; vind[nv++]=obss->ircv[j];
        }
        if (nv<=0) continue;

        /* measurement outliers detect */
        int outln=measdetoutl(v,nv,outli,NULL,4.0);
        m+=outln;

        for (j=0;j<outln;j++) {
            pod->flt.ssat[vind[outli[j]]][i].rejc[0]=1;
            pod->flt.ssat[vind[outli[j]]][i].rejc[1]=1;
        }
    }
    free(outli); free(vind); free(v);
    return m;
}
/* satellite orbit estimator using EKF and satellite observation--------------*/
static int podflt_satobs(pod_t *pod, const podobs_t *obss)
{
    pod_opt_t *opt=&pod->opt;
    int i,j,k,stat=0,nv,nx=pod->flt.nx,iter,*exc,maxobsn=-1;
    double *v,*H,*R,*xp,*Pp,*var;
    podobs_t obss_=*obss;

    pod->sol.stat=SOLQ_NONE;

    /* adjust observation data */
    adjobss(pod,&obss_);

    for (i=0;i<obss->nrcv;i++) {
        for (k=0,j=0;j<obss->nobs[i];j++) {

            /* unselected sat sys */
            if (!(satsys(obss->obs[i][j].sat,NULL)&pod->opt.estopt.navsys)) continue;
            obss_.obs[i][k++]=obss->obs[i][j];
        }
        obss_.nobs[i]=k;
        if (maxobsn<0||k>maxobsn) maxobsn=k;
    }
    if (pod->flt.time.time) {
        pod->flt.dt=timediff(obss_.tutc,pod->flt.time);
    }
    nv=obss_.nrcv*maxobsn*2+MAXSAT*6;
    v=mat(nv,1); H=zeros(nv,nx); var=mat(nv,1); R=zeros(nv,nv);
    exc=imat(1,MAXOBS*obss->nrcv);

    for (k=0;k<obss_.nrcv;k++) {
        for (i=0;i<MAXSAT;i++) {
            for (j=0;j<opt->estopt.nf;j++) pod->flt.ssat[obss_.ircv[k]-1][i].rejc[j]=0;
        }
    }
    /* temporal update of ekf states */
    for (k=0;k<obss_.nrcv;k++) {
        sta_t *sta=&pod->rcv[obss_.ircv[k]-1].sta;

        /* set antenna parameters */
        setpcv(utc2gpst(obss_.tutc),&opt->estopt,&pod->nav,&pod->pcvss,&pod->pcvsr,sta);

        /* update of POD states */
        udstate_pod(pod,obss_.tutc,obss_.ircv[k],obss_.obs[k],obss_.nobs[k]);
    }
    /* update of satellite position/velocity/clock/SRP */
    udsatorbclksrp_pod(pod,obss_.tutc);

    /* prefit measurement residuals */
    nv=podmeasres(pod,0,&obss_,pod->flt.x,pod->flt.P,v,H,var,exc);

    /* detect measurement outliers */
    if (pod->opt.outldetexcs) {
        if (poddetmeasoutl(pod,&obss_)>0) {
            nv=podmeasres(pod,0,&obss_,pod->flt.x,pod->flt.P,v,H,var,exc);
        }
    }
    /* measurement variance */
    for (memset(R,0,nv*nv*sizeof(double)),i=0;i<nv;i++) R[i+i*nv]=var[i];

    /* measurement update of ekf states */
    if (filter(pod->flt.x,pod->flt.P,H,v,R,nx,nv,NULL,0)) {
        log_trace(2,"EKF filter error\n");
        stat=SOLQ_NONE;
    }
    else stat=SOLQ_PPP;

    /* update POD solution status */
    udpossolstat(pod,&obss_,stat);

    free(v); free(H); free(R); free(var); free(exc);
    return stat;
}
/* satellite orbit estimator using EKF----------------------------------------*/
extern int podflt(pod_t *pod, const podobs_t *obs)
{
    switch (obs->type) {
        case GPOD_OBSS_TYPE_SATPOS: return podflt_satpos(pod,obs);
        case GPOD_OBSS_TYPE_SATOBS: return podflt_satobs(pod,obs);
    }
    return 0;
}
/* forward/backward combined filter-------------------------------------------*/
extern int podfbcmb(const pod_opt_t *opt, pod_sol_t **sols)
{
    FILE *fpobs=fopen(opt->obsbinfile,"rb");
    pod_t podf={0};
    pod_t podb={0};
    pod_sol_t *solf;
    pod_sol_t *solb;
    int nsolf=0,nmaxsolf=1024;
    int nsolb=0,nmaxsolb=1024;
    int i,j,k,x;

    if (fpobs==NULL) return -1;

    podinit(&podf,opt);
    solf=calloc(nmaxsolf,sizeof(pod_sol_t));

    /* POD forward filter */
    for (;;) {
        podobs_t podobs={0};
        if (!fread(&podobs,sizeof(podobs),1,fpobs)) break;

        podflt(&podf,&podobs);

        fprintf(stderr,"POD forward processing: %s nrcv=%3d\n",time_str(utc2gpst(podobs.tutc),2),podobs.nrcv);

        if (podf.sol.stat>0) {
            pod_sol_t *solft;
            if (nsolf>nmaxsolf) nmaxsolf*=2,solft=realloc(solf,nmaxsolf*sizeof(pod_sol_t)),solf=solft;
            solf[nsolf++]=podf.sol;
        }
    }
    podfree(&podf);

    podinit(&podb,opt);
    solb=calloc(nmaxsolb,sizeof(pod_sol_t));

    fseeko64(fpobs,0,SEEK_END);

    /* POD backward filter */
    for (;;) {
        podobs_t podobs={0};

        if (fseeko64(fpobs,(_off64_t)-sizeof(podobs),SEEK_CUR)<0) break;
        if (!fread(&podobs,sizeof(podobs),1,fpobs)) break;
        if (fseeko64(fpobs,(_off64_t)-sizeof(podobs),SEEK_CUR)<0) break;

        fprintf(stderr,"POD backward processing: %s nrcv=%3d\n",time_str(utc2gpst(podobs.tutc),2),podobs.nrcv);

        podflt(&podb,&podobs);

        if (podb.sol.stat>0) {
            pod_sol_t *solbt;
            if (nsolf>nmaxsolf) nmaxsolb*=2,solbt=realloc(solf,nmaxsolf*sizeof(pod_sol_t)),solb=solbt;
            solb[nsolb++]=podb.sol;
        }
    }
    podfree(&podb);

    *sols=calloc(max(nsolf,nsolb),sizeof(pod_sol_t));
    int nsols=0;

    /* POD forward/backward combined filter */
    for (i=0;i<nsolf;i++) {
        for (j=0;j<nsolb;j++) {
            if (fabs(timediff(solf[i].tutc,solb[j].tutc))<1E-2) {
                (*sols)[nsols]=solf[i];

                for (k=0;k<MAXSAT;k++) {
                    for (x=0;x<GNX;x++) {
                        if (fabs(solf[i].xs[k][x])<1E-30||fabs(solb[j].xs[k][x])<1E-30) continue;
                        (*sols)[nsols].Qs[k][x]=1.0/(1.0/solf[i].Qs[k][x]+1.0/solb[j].Qs[k][x]);
                        (*sols)[nsols].xs[k][x]=(*sols)[nsols].Qs[k][x]*(1.0/solf[i].Qs[k][x]*solf[i].xs[k][x]+1.0/solb[j].Qs[k][x]*solb[j].xs[k][x]);
                    }
                }
                nsols++; break;
            }
        }
    }
    fclose(fpobs);
    free(solf);
    free(solb);
    return nsols;
}
extern int podfbcmbe(const pod_opt_t *opt, pod_sol_t **sols)
{
    FILE *fpobs=fopen(opt->obsbinfile,"rb");
    pod_t pod={0};
    pod_sol_t *solf;
    int nsolf=0,nmaxsolf=1024;
    int i,j,k,x;

    if (fpobs==NULL) return -1;

    podinit(&pod,opt);
    solf=calloc(nmaxsolf,sizeof(pod_sol_t));

    /* POD forward filter */
    for (;;) {
        podobs_t podobs={0};
        if (!fread(&podobs,sizeof(podobs),1,fpobs)) break;

        podflt(&pod,&podobs);

        fprintf(stderr,"POD forward processing: %s nrcv=%3d\n",time_str(utc2gpst(podobs.tutc),2),podobs.nrcv);

        if (pod.sol.stat>0) {
            pod_sol_t *solft;
            if (nsolf>nmaxsolf) nmaxsolf*=2,solft=realloc(solf,nmaxsolf*sizeof(pod_sol_t)),solf=solft;
            solf[nsolf++]=pod.sol;
        }
    }
    *sols=calloc(nsolf*2,sizeof(pod_sol_t));
    int nsols=0;

    /* POD backward filter */
    for (;;) {
        podobs_t podobs={0};

        if (fseeko64(fpobs,(_off64_t)-sizeof(podobs),SEEK_CUR)<0) break;
        if (!fread(&podobs,sizeof(podobs),1,fpobs)) break;
        if (fseeko64(fpobs,(_off64_t)-sizeof(podobs),SEEK_CUR)<0) break;

        fprintf(stderr,"POD backward processing: %s nrcv=%3d\n",time_str(utc2gpst(podobs.tutc),2),podobs.nrcv);

        podflt(&pod,&podobs);

        if (pod.sol.stat>0) {
            (*sols)[nsols++]=pod.sol;
        }
    }
    fclose(fpobs);
    free(solf);
    return nsols;
}