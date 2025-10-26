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

#define THRES_MW_JUMP 10.0

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

/* initial satellite orbit position/velocity using navigation data------------*/
extern void initsatorbclksrp_nav(pod_t *pod, int sat, gtime_t tutc);

/* initialize state and covariance -------------------------------------------*/
extern void initx(pod_t *pod, double xi, double var, int i);

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
    return 0;
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
extern int model_phw(gtime_t time, int sat, const char *type, int opt,
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
extern double podvarerr(int sat, int sys, double el, int idx, int type,
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
extern void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
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
extern void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
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
extern int model_trop(pod_t *pod, int rcv, gtime_t time, const double *pos, const double *azel,
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
extern int model_iono(gtime_t time, const double *pos, const double *azel,
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
extern int ppp_res(int post, int rcv, const obsd_t *obs, int n,
                   const double *dr, const nav_t *nav,
                   const double *x, const double *Pp, pod_t *pod, double *v, double *H, double *var, double *azel, int *exc, int *vflag);

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