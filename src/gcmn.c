/*------------------------------------------------------------------------------
 * orbit.c : GNSS orbit functions

 * author  : sujinglan
 * version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
 * Copyright(c) 2023-2025 by sujinglan, all rights reserved
 * history : 2024/10/17 1.0  new
 *-----------------------------------------------------------------------------*/
#include "podlib.h"
#include "eopspw.h"

/* satellite informations------------------------------------------------------*/
extern pod_satinfo_t gsatinfos[256]={
        {"G01","BLOCK IIIA" ,"N",2161.0,0.11},
        {"G02","BLOCK IIR-B","U",1080.0,0.20},
        {"G03","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G04","BLOCK IIIA" ,"N",2161.0,0.11},
        {"G05","BLOCK IIR-M","U",1080.0,0.20},
        {"G06","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G07","BLOCK IIR-M","U",1088.0,0.20},
        {"G08","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G09","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G10","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G11","BLOCK IIIA" ,"N",2161.0,0.11},
        {"G12","BLOCK IIR-M","U",1088.0,0.20},
        {"G13","BLOCK IIR-A","U",1088.0,0.20},
        {"G14","BLOCK IIIA" ,"N",2161.0,0.11},
        {"G15","BLOCK IIR-M","U",1088.0,0.20},
        {"G16","BLOCK IIR-A","U",1088.0,0.20},
        {"G17","BLOCK IIR-M","U",1088.0,0.20},
        {"G18","BLOCK IIIA" ,"N",2161.0,0.11},
        {"G19","BLOCK IIR-B","U",1088.0,0.20},
        {"G20","BLOCK IIR-A","U",1088.0,0.20},
        {"G21","BLOCK IIIA" ,"N",2161.0,0.11},
        {"G22","BLOCK IIR-A","U",1088.0,0.20},
        {"G23","BLOCK IIIA" ,"N",2161.0,0.11},
        {"G24","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G25","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G26","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G27","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G28","BLOCK IIIA" ,"N",2161.0,0.11},
        {"G29","BLOCK IIR-M","N",1088.0,0.20},
        {"G30","BLOCK IIF"  ,"N",1633.0,0.11},
        {"G31","BLOCK IIR-M","U",1088.0,0.20},
        {"G32","BLOCK IIF"  ,"N",1633.0,0.11},

        {"R01","GLONASS-M"  ,"U",1415.0,0.25},
        {"R02","GLONASS-M"  ,"U",1415.0,0.25},
        {"R03","GLONASS-M"  ,"U",1415.0,0.25},
        {"R04","GLONASS-M+" ,"U",1415.0,0.25},
        {"R05","GLONASS-M+" ,"U",1415.0,0.25},
        {"R06","GLONASS-M"  ,"U",1415.0,0.25},
        {"R07","GLONASS-M"  ,"U",1415.0,0.25},
        {"R08","GLONASS-M"  ,"U",1415.0,0.25},
        {"R09","GLONASS-K1B","U", 995.0,0.25},
        {"R10","GLONASS-M"  ,"U",1415.0,0.25},
        {"R11","GLONASS-K1B","U", 995.0,0.25},
        {"R12","GLONASS-M+" ,"U",1415.0,0.25},
        {"R13","GLONASS-M"  ,"U",1415.0,0.25},
        {"R14","GLONASS-M"  ,"U",1415.0,0.25},
        {"R15","GLONASS-M"  ,"U",1415.0,0.25},
        {"R16","GLONASS-M+" ,"U",1415.0,0.25},
        {"R17","GLONASS-M"  ,"U",1415.0,0.25},
        {"R18","GLONASS-M"  ,"U",1415.0,0.25},
        {"R19","GLONASS-K1+","U", 995.0,0.25},
        {"R20","GLONASS-M"  ,"U",1415.0,0.25},
        {"R21","GLONASS-M+" ,"U",1415.0,0.25},
        {"R22","GLONASS-K1B","U", 995.0,0.25},
        {"R23","GLONASS-M"  ,"U",1415.0,0.25},
        {"R24","GLONASS-M+" ,"U",1415.0,0.25},
        {"R25","GLONASS-M"  ,"U",1415.0,0.25},
        {"R26","GLONASS-K2" ,"U",1645.0,0.25},
        {"R27","GLONASS-K2" ,"U",1645.0,0.25},

        {"E01","GALILEO-2","U",706673/1E3,0.40},
        {"E02","GALILEO-2","U",709858/1E3,0.40},
        {"E03","GALILEO-2","U",707634/1E3,0.40},
        {"E04","GALILEO-2","U",708765/1E3,0.40},
        {"E05","GALILEO-2","U",709755/1E3,0.40},
        {"E06","GALILEO-2","U",708301/1E3,0.40},
        {"E07","GALILEO-2","U",706643/1E3,0.40},
        {"E08","GALILEO-2","U",709135/1E3,0.40},
        {"E09","GALILEO-2","U",707879/1E3,0.40},
        {"E10","GALILEO-2","U",710670/1E3,0.40},
        {"E11","GALILEO-1","U",696802/1E3,0.40},
        {"E12","GALILEO-1","U",695318/1E3,0.40},
        {"E13","GALILEO-2","U",710795/1E3,0.40},
        {"E14","GALILEO-2","U",662141/1E3,0.40},
        {"E15","GALILEO-2","U",710703/1E3,0.40},
        {"E16","GALILEO-2","U",706384/1E3,0.40},
        {"E18","GALILEO-2","U",660977/1E3,0.40},
        {"E19","GALILEO-1","U",697632/1E3,0.40},
        {"E20","GALILEO-1","U",693842/1E3,0.40},
        {"E21","GALILEO-2","U",708500/1E3,0.40},
        {"E22","GALILEO-2","U",697701/1E3,0.40},
        {"E23","GALILEO-2","U",711103/1E3,0.40},
        {"E24","GALILEO-2","U",708783/1E3,0.40},
        {"E25","GALILEO-2","U",710475/1E3,0.40},
        {"E26","GALILEO-2","U",705685/1E3,0.40},
        {"E27","GALILEO-2","U",710244/1E3,0.40},
        {"E29","GALILEO-2","U",710649/1E3,0.40},
        {"E30","GALILEO-2","U",707734/1E3,0.40},
        {"E31","GALILEO-2","U",711514/1E3,0.40},
        {"E33","GALILEO-2","U",709500/1E3,0.40},
        {"E34","GALILEO-2","U",709500/1E3,0.40},
        {"E36","GALILEO-2","U",711410/1E3,0.40},

        {"C01","BEIDOU-2G"       ,"U",1500.0,0.30},
        {"C02","BEIDOU-2G"       ,"U",1551.0,0.30},
        {"C03","BEIDOU-2G"       ,"U",1486.0,0.30},
        {"C04","BEIDOU-2G"       ,"U",1536.0,0.30},
        {"C05","BEIDOU-2G"       ,"U",1502.0,0.30},
        {"C06","BEIDOU-2I"       ,"U",1284.0,0.30},
        {"C07","BEIDOU-2I"       ,"U",1284.0,0.30},
        {"C08","BEIDOU-2I"       ,"U",1284.0,0.30},
        {"C09","BEIDOU-2I"       ,"U",1284.0,0.30},
        {"C10","BEIDOU-2I"       ,"U",1278.0,0.30},
        {"C11","BEIDOU-2M"       ,"U",1193.0,0.30},
        {"C12","BEIDOU-2M"       ,"U",1176.0,0.30},
        {"C13","BEIDOU-2I"       ,"U",1283.0,0.30},
        {"C14","BEIDOU-2M"       ,"U",1184.0,0.30},
        {"C19","BEIDOU-3M-CAST"  ,"U", 943.0,0.30},
        {"C20","BEIDOU-3M-CAST"  ,"U", 942.0,0.30},
        {"C21","BEIDOU-3M-CAST"  ,"U", 942.0,0.30},
        {"C22","BEIDOU-3M-CAST"  ,"U", 941.0,0.30},
        {"C23","BEIDOU-3M-CAST"  ,"U", 945.0,0.30},
        {"C24","BEIDOU-3M-CAST"  ,"U", 946.0,0.30},
        {"C25","BEIDOU-3M-SECM-A","U",1043.0,0.30},
        {"C26","BEIDOU-3M-SECM-A","U",1041.0,0.30},
        {"C27","BEIDOU-3M-SECM-A","U",1018.0,0.30},
        {"C28","BEIDOU-3M-SECM-A","U",1014.0,0.30},
        {"C29","BEIDOU-3M-SECM-A","U",1010.0,0.30},
        {"C30","BEIDOU-3M-SECM-A","U",1008.0,0.30},
        {"C31","BEIDOU-3SI-SECM" ,"U", 848.0,0.30},
        {"C32","BEIDOU-3SI-CAST" ,"U",1007.0,0.30},
        {"C33","BEIDOU-3SI-CAST" ,"U",1007.0,0.30},
        {"C34","BEIDOU-3M-SECM-A","U",1046.0,0.30},
        {"C35","BEIDOU-3M-SECM-A","U",1045.0,0.30},
        {"C36","BEIDOU-3M-CAST"  ,"U",1061.0,0.30},
        {"C37","BEIDOU-3M-CAST"  ,"U",1061.0,0.30},
        {"C38","BEIDOU-3I"       ,"U",2952.0,0.30},
        {"C39","BEIDOU-3I"       ,"U",2949.0,0.30},
        {"C40","BEIDOU-3I"       ,"U",2870.0,0.30},
        {"C41","BEIDOU-3M-CAST"  ,"U",1059.0,0.30},
        {"C42","BEIDOU-3M-CAST"  ,"U",1059.0,0.30},
        {"C43","BEIDOU-3M-SECM-B","U",1078.0,0.30},
        {"C44","BEIDOU-3M-SECM-B","U",1075.0,0.30},
        {"C45","BEIDOU-3M-CAST"  ,"U",1058.0,0.30},
        {"C46","BEIDOU-3M-CAST"  ,"U",1059.0,0.30},
        {"C47","BEIDOU-3M-SECM-B","U",1075.0,0.30},
        {"C48","BEIDOU-3M-CAST"  ,"U",1060.0,0.30},
        {"C49","BEIDOU-3M-SECM-B","U",1075.0,0.30},
        {"C50","BEIDOU-3M-CAST"  ,"U",1060.0,0.30},
        {"C56","BEIDOU-3SI-CAST" ,"U",2800.0,0.30},
        {"C57","BEIDOU-3SM-CAST" ,"U",1014.0,0.30},
        {"C58","BEIDOU-3SM-CAST" ,"U",1014.0,0.30},
        {"C59","BEIDOU-3G"       ,"U",2968.0,0.30},
        {"C60","BEIDOU-3G"       ,"U",2968.0,0.30},
        {"C61","BEIDOU-3G"       ,"U",2968.0,0.30},
        {"C62","BEIDOU-3G"       ,"U",2968.0,0.30},

        {"J01","QZSS-1" ,"U",2207.0,0.00},
        {"J02","QZSS-2I","U",2212.0,0.00},
        {"J03","QZSS-2I","U",2215.0,0.00},
        {"J04","QZSS-2A","U",2297.0,0.00},
        {"J07","QZSS-2G","U",2376.0,0.00},
        {"J08","QZSS-3G","U",2700.0,0.00},

        {"I01","IRNSS-1I","U", 700.0,0.00},
        {"I02","IRNSS-1I","U", 700.0,0.00},
        {"I03","IRNSS-1G","U", 700.0,0.00},
        {"I11","IRNSS-2I","U",2250.0,0.00},
};
/* get satellite information-------------------------------------------------*/
extern pod_satinfo_t *podgetsatinfo(int sat)
{
    int i;

    for (i=0;i<sizeof(gsatinfos)/sizeof(pod_satinfo_t);i++) {
        if (satid2no(gsatinfos[i].id)!=sat) continue;
        return &gsatinfos[i];
    }
    return NULL;
}
/* eci/ecsf to ecef transformation matrix -------------------------------------
 * compute eci to ecef transformation matrix
 * args   : gtime_t tutc     I   time in utc
 *          double *erpv     I   erp values {xp,yp,ut1_utc,lod,tai-utc} (rad,rad,s,s/d,s)
 *          double *U        O   eci to ecef transformation matrix (3 x 3)
 *          double *dU       O   dot of eci to ecef transformation matrix (3 x 3)
 *          double *gmst     IO  greenwich mean sidereal time (rad)
 *                               (NULL: no output)
 * return : none
 *-----------------------------------------------------------------------------*/
extern void ecsf2ecef(gtime_t tutc, const double *erpv, double *U, double *dU, double *gmst)
{
    const static double mjd0=2400000.5;
    static double U_[9],dU_[9],gmst_;
    static gtime_t tutc_;
    double mjd_utc,mjd_ut1,mjd_tt;
    double ep[6],dt[6],T[3][3],I[3][3]={0.0},W[3][3],gast,omg,S[9]={0},Rz[9],O[9],NP[9],PT[9];
    double corr[6]={0};
    int i,j;

    if (tutc.time&&fabs(timediff(tutc,tutc_))<0.01) { /* read cache */
        if (U) matcpy(U,U_,3,3);
        if (dU) matcpy(dU,dU_,3,3);
        if (gmst) *gmst=gmst_;
        return;
    }
#if PMUT1_OCEAN_CORR
    pmut1_ocean(tutc,corr);
#endif
    tutc_=tutc;
    time2epoch(tutc_,ep);
    mjd_utc=cal2mjd(ep);

    reftimediff(erpv[2]+corr[2],erpv[4],dt);
    mjd_ut1=mjd_utc+(erpv[2]+corr[2])/86400.0;
    mjd_tt=mjd_utc+dt[3]/86400.0;

    /* iau precession and nutation */
    iauPnm06a(mjd0,mjd_tt,T);
    gast=iauGst06(mjd0,mjd_ut1,mjd0,mjd_tt,T);
    iauPom00(erpv[0]+corr[0]*AS2R,erpv[1]+corr[1]*AS2R,iauSp00(mjd0,mjd_tt),W);

    /* earth rotation matrix */
    I[0][0]=I[1][1]=I[2][2]=1.0;
    iauRz(gast,I);
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) Rz[i+j*3]=I[i][j];
    }
    /* derivative of earth rotation matrix */
    S[1]=-1.0; S[3]=1.0;
    omg=OMGE-0.843994809*1E-9*(erpv[3]+corr[3]);
    matmul("NN",3,3,3,omg,S,Rz,0.0,O);

    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) NP[i+j*3]=T[i][j];
        for (j=0;j<3;j++) PT[i+j*3]=W[i][j];
    }
    matmul("NN",3,3,3,1.0,PT,Rz,0.0,S); /* PT=Ry(-xp)*Rx(-yp) */
    matmul("NN",3,3,3,1.0,S ,NP,0.0,U_); /* U=PT*Rz(gast)*NP */

    matmul("NN",3,3,3,1.0,PT,O,0.0,S); /* PT=Ry(-xp)*Rx(-yp) */
    matmul("NN",3,3,3,1.0,S,NP,0.0,dU_); /* dU=PT*O(gast)*NP */
    gmst_=iauGmst06(mjd0,mjd_ut1,mjd0,mjd_tt);

    if (U) matcpy(U,U_,3,3);
    if (dU) matcpy(dU,dU_,3,3);
    if (gmst) *gmst=gmst_;

    log_trace(4,"gmst=%.12f gast=%.12f\n",gmst_,gast);
}
extern void ecsf2ecef2(gtime_t tutc, const double *erpv, double *U, double *dU, double *gmst)
{
    static double U_[9],dU_[9],gmst_;
    static gtime_t tutc_;
    const double mjd0=2400000.5;
    const double omge=2.0*PI*1.00273781191135448/86400.0;
    double mjd_utc,mjd_ut1,mjd_tt,ep[6],dt[6],x,y,s,sp,era,dera;
    double Q[3][3]={0},R[3][3]={0},W[3][3]={0};
    double Q_[9],R_[9],W_[9],T[9],dR[9]={0},corr[6]={0};
    int i,j;

    if (tutc.time&&fabs(timediff(tutc,tutc_))<1E-6) { /* read cache */
        if (U) matcpy(U,U_,3,3);
        if (dU) matcpy(dU,dU_,3,3);
        if (gmst) *gmst=gmst_;
        return;
    }
#if PMUT1_OCEAN_CORR
    pmut1_ocean(tutc,corr);
#endif
    time2epoch(tutc,ep);
    mjd_utc=cal2mjd(ep);

    reftimediff(erpv[2]+corr[2],erpv[4],dt);
    mjd_ut1=mjd_utc+(erpv[2]+corr[2])/86400.0;
    mjd_tt=mjd_utc+dt[3]/86400.0;

    iauXy06(mjd0,mjd_tt,&x,&y);
    s=iauS06(mjd0,mjd_tt,x,y);
    iauC2ixys(x,y,s,Q);
    era=iauEra00(mjd0,mjd_ut1);

    R[0][0]=R[1][1]=R[2][2]=1.0;
    iauRz(era,R);
    sp=iauSp00(mjd0,mjd_tt);
    iauPom00(erpv[0]+corr[0]*AS2R,erpv[1]+corr[1]*AS2R,sp,W);

    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) W_[i+j*3]=W[i][j];
        for (j=0;j<3;j++) R_[i+j*3]=R[i][j];
        for (j=0;j<3;j++) Q_[i+j*3]=Q[i][j];
    }
    matmul("NN",3,3,3,1.0,W_,R_,0.0,T);
    matmul("NN",3,3,3,1.0,T,Q_,0.0,U_);

    dera=omge-0.843994809*1E-9*(erpv[3]+corr[3]);
    dR[0+0*3]=-sin(era)*dera; dR[0+1*3]= cos(era)*dera;
    dR[1+0*3]=-cos(era)*dera; dR[1+1*3]=-sin(era)*dera;

    matmul("NN",3,3,3,1.0,W_,dR,0.0,T);
    matmul("NN",3,3,3,1.0,T,Q_,0.0,dU_);

    gmst_=iauGmst06(mjd0,mjd_ut1,mjd0,mjd_tt);

    if (U) matcpy(U,U_,3,3);
    if (dU) matcpy(dU,dU_,3,3);
    if (gmst) *gmst=gmst_;

    log_trace(3,"gmst=%.12f\n",gmst_);
}
extern void eci2ecef2(gtime_t tutc, const erp_t *erp, double *U, double *dU, double *gmst)
{
    double erpv[8];

    /* earth rotation parameter values */
    geterp(erp,utc2gpst(tutc),erpv);

    /* ECI to ECEF transformation matrix */
    ecsf2ecef(tutc,erpv,U,NULL,NULL);
}
/*  the relation between different time system follows as:
 *
 *           -14s
 *     -----------------> BDT(Compass Time)
 *     |
 *     |         +19s             +32.184s           +rel.effects
 *    GPS -------------> TAI ----------------> TT -----------------> TDB
 *                       T |
 *            -(UT1-TAI) | |    -leap seconds
 *    UT1 ---------------| |--------------------> UTC
 *     |
 *     |   earth rotation
 *     ---------------------> GAST
 *
 *  Most of the methods are modified from the code provided by
 *  Montenbruck, Oliver, 2001, as the attachement of the book
 *  "Satellite Orbits: Models, Methods and applications".
 * */
/* reference system time difference-------------------------------------------
 * args   : double ut1_utc   I   UT1-UTC (s)
 *          double tai_utc   I   TAI-UTC (s)
 *          double *rdt      O   reference system time difference (s)
 *                               (0:UT1-TAI,1:UTC-GPS,2:UT1-GPS,3:TT-UTC,4:GPS-UTC)
 * return : none
 *-----------------------------------------------------------------------------*/
extern void reftimediff(const double ut1_utc, const double tai_utc, double *rdt)
{
    const static double tt_tai=32.184;
    const static double gps_tai=-19.0;
    rdt[0]= ut1_utc-tai_utc; /* UT1-TAI */
    rdt[1]=-tai_utc-gps_tai; /* UTC-GPS */
    rdt[2]=rdt[0]-gps_tai;   /* UT1-GPS */
    rdt[3]=tt_tai+tai_utc;   /* TT-UTC */
    rdt[4]=gps_tai+tai_utc;  /* GPS-UTC */
}
/* chebyshev approximation of 3D vectors--------------------------------------*/
static void cheb3d(int n, double t, double ta, double tb, const double *Cx,
                   const double *Cy, const double *Cz, double *rv)
{
    int i;
    double tau,f1[3]={0},f2[3]={0},tmp[3];

    tau=(2.0*t-ta-tb)/(tb-ta);
    for (i=n-1;i>=1;i--) {
        matcpy(tmp,f1,1,3);
        f1[0]=2.0*tau*f1[0]-f2[0]+Cx[i];
        f1[1]=2.0*tau*f1[1]-f2[1]+Cy[i];
        f1[2]=2.0*tau*f1[2]-f2[2]+Cz[i];
        matcpy(f2,tmp,1,3);
    }
    rv[0]=tau*f1[0]-f2[0]+Cx[0];
    rv[1]=tau*f1[1]-f2[1]+Cy[0];
    rv[2]=tau*f1[2]-f2[2]+Cz[0];
}
/* sun,moon and nine major planets' equatorial position using JPL ephemerides
* referred to the international celestial reference frame (ICRF)
*----------------------------------------------------------------------------*/
static double jpleph_de440coeff[12566][1018];
static int jpleph_udflag=0;
extern int jpleph_de440(gtime_t tutc, const double *erpv, double *mercury, double *venus,
                        double *earth, double *mars, double *jupiter, double *saturn,
                        double *uranus, double *neptune, double *pluto, double *moon,
                        double *sun, double *sunssb)
{
    const static double mjd0=2400000.5;
    double mjd_utc,mjd_tdb,ep[6],jd,t1,dt,Cmx[13*32],Cmy[13*32],Cmz[13*32],mjd0_;
    int i,j,k=-1,tmp[13];

    time2epoch(tutc,ep);
    mjd_utc=cal2mjd(ep);
    mjd_tdb=mjdtdb(mjd_utc,erpv);
    jd=mjd_tdb+mjd0;

    for (i=0;i<12566;i++) {
        if (jpleph_de440coeff[i][0]<=jd&&jpleph_de440coeff[i][1]>=jd) {k=i; break;}
    }
    if (k<0) return 0;

    t1=jpleph_de440coeff[k][0]-mjd0;
    dt=mjd_tdb-t1;
    for (j=0,i=230;i<=269;i+=13) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,13);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,13);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,13);
    matcpy(Cmx+13,&jpleph_de440coeff[k][tmp[0]+39],1,13);
    matcpy(Cmy+13,&jpleph_de440coeff[k][tmp[1]+39],1,13);
    matcpy(Cmz+13,&jpleph_de440coeff[k][tmp[2]+39],1,13);

    if (dt>=0.0&&dt<=16.0) {
        j=0;
        mjd0_=t1;
    }
    else if (16.0<dt&&dt<=32.0) {
        j=1;
        mjd0_=t1+16*j;
    }
    double r_EarthMoon[3];
    cheb3d(13,mjd_tdb,mjd0_,mjd0_+16.0,&Cmx[13*j],&Cmy[13*j],&Cmz[13*j],r_EarthMoon);
    r_EarthMoon[0]*=1E3;
    r_EarthMoon[1]*=1E3;
    r_EarthMoon[2]*=1E3;

    for (j=0,i=440;i<=479;i+=13) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,13);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,13);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,13);

    for (i=1;i<=7;i++) {
        matcpy(Cmx+13*i,&jpleph_de440coeff[k][tmp[0]+39*i],1,13);
        matcpy(Cmy+13*i,&jpleph_de440coeff[k][tmp[1]+39*i],1,13);
        matcpy(Cmz+13*i,&jpleph_de440coeff[k][tmp[2]+39*i],1,13);
    }
    if (dt>=0.0&&dt<=4.0) {
        j=0;
        mjd0_=t1;
    }
    else if (4.0<dt&&dt<=8.0) {
        j=1;
        mjd0_=t1+4.0*j;
    }
    else if (8.0<dt&&dt<=12.0) {
        j=2;
        mjd0_=t1+4.0*j;
    }
    else if (12.0<dt&&dt<=16.0) {
        j=3;
        mjd0_=t1+4.0*j;
    }
    else if (16.0<dt&&dt<=20.0) {
        j=4;
        mjd0_=t1+4.0*j;
    }
    else if (20.0<dt&&dt<=24.0) {
        j=5;
        mjd0_=t1+4.0*j;
    }
    else if (24.0<dt&&dt<=28.0) {
        j=6;
        mjd0_=t1+4.0*j;
    }
    else if (28.0<dt&&dt<=32.0) {
        j=7;
        mjd0_=t1+4.0*j;
    }
    double r_Moon[3];
    cheb3d(13,mjd_tdb,mjd0_,mjd0_+4.0,&Cmx[13*j],&Cmy[13*j],&Cmz[13*j],r_Moon);
    r_Moon[0]*=1E3;
    r_Moon[1]*=1E3;
    r_Moon[2]*=1E3;

    for (j=0,i=752;i<=785;i+=11) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,11);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,11);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,11);
    matcpy(Cmx+11,&jpleph_de440coeff[k][tmp[0]+33],1,11);
    matcpy(Cmy+11,&jpleph_de440coeff[k][tmp[1]+33],1,11);
    matcpy(Cmz+11,&jpleph_de440coeff[k][tmp[2]+33],1,11);

    if (dt>=0.0&&dt<=16.0) {
        j=0;
        mjd0_=t1;
    }
    else if (16.0<dt&&dt<=32.0) {
        j=1;
        mjd0_=t1+16*j;
    }
    double r_Sun[3];
    cheb3d(11,mjd_tdb,mjd0_,mjd0_+16.0,&Cmx[11*j],&Cmy[11*j],&Cmz[11*j],r_Sun);
    r_Sun[0]*=1E3;
    r_Sun[1]*=1E3;
    r_Sun[2]*=1E3;
    for (j=0,i=2;i<=44;i+=14) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,14);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,14);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,14);
    for (i=1;i<=3;i++) {
        matcpy(Cmx+14*i,&jpleph_de440coeff[k][tmp[0]+42*i],1,14);
        matcpy(Cmy+14*i,&jpleph_de440coeff[k][tmp[1]+42*i],1,14);
        matcpy(Cmz+14*i,&jpleph_de440coeff[k][tmp[2]+42*i],1,14);
    }
    if (dt>=0.0&&dt<=8.0) {
        j=0;
        mjd0_=t1;
    }
    else if (8.0<dt&&dt<=16.0) {
        j=1;
        mjd0_=t1+8.0*j;
    }
    else if (16.0<dt&&dt<=24.0) {
        j=2;
        mjd0_=t1+8.0*j;
    }
    else if (24.0<dt&&dt<=32.0) {
        j=2;
        mjd0_=t1+8.0*j;
    }
    double r_Mercury[3];
    cheb3d(14,mjd_tdb,mjd0_,mjd0_+8.0,&Cmx[14*j],&Cmy[14*j],&Cmz[14*j],r_Mercury);
    r_Mercury[0]*=1E3;
    r_Mercury[1]*=1E3;
    r_Mercury[2]*=1E3;

    for (j=0,i=170;i<=200;i+=10) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,10);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,10);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,10);
    matcpy(Cmx+10,&jpleph_de440coeff[k][tmp[0]+30],1,10);
    matcpy(Cmy+10,&jpleph_de440coeff[k][tmp[1]+30],1,10);
    matcpy(Cmz+10,&jpleph_de440coeff[k][tmp[2]+30],1,10);

    if (dt>=0.0&&dt<=16.0) {
        j=0;
        mjd0_=t1;
    }
    else if (16.0<dt&&dt<=32.0) {
        j=1;
        mjd0_=t1+16*j;
    }
    double r_Venus[3];
    cheb3d(10,mjd_tdb,mjd0_,mjd0_+16.0,&Cmx[10*j],&Cmy[10*j],&Cmz[10*j],r_Venus);
    r_Venus[0]*=1E3;
    r_Venus[1]*=1E3;
    r_Venus[2]*=1E3;

    for (j=0,i=308;i<=341;i+=11) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,11);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,11);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,11);
    j=0; mjd0_=t1;

    double r_Mars[3];
    cheb3d(11,mjd_tdb,mjd0_,mjd0_+32.0,&Cmx[11*j],&Cmy[11*j],&Cmz[11*j],r_Mars);
    r_Mars[0]*=1E3;
    r_Mars[1]*=1E3;
    r_Mars[2]*=1E3;

    for (j=0,i=341;i<=365;i+=8) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,8);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,8);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,8);
    j=0; mjd0_=t1;

    double r_Jupiter[3];
    cheb3d(8,mjd_tdb,mjd0_,mjd0_+32.0,&Cmx[8*j],&Cmy[8*j],&Cmz[8*j],r_Jupiter);
    r_Jupiter[0]*=1E3;
    r_Jupiter[1]*=1E3;
    r_Jupiter[2]*=1E3;

    for (j=0,i=365;i<=386;i+=7) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,7);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,7);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,7);
    j=0; mjd0_=t1;

    double r_Satrun[3];
    cheb3d(7,mjd_tdb,mjd0_,mjd0_+32.0,&Cmx[7*j],&Cmy[7*j],&Cmz[7*j],r_Satrun);
    r_Satrun[0]*=1E3;
    r_Satrun[1]*=1E3;
    r_Satrun[2]*=1E3;

    for (j=0,i=386;i<=404;i+=6) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,6);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,6);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,6);
    j=0; mjd0_=t1;

    double r_Uranus[3];
    cheb3d(6,mjd_tdb,mjd0_,mjd0_+32.0,&Cmx[6*j],&Cmy[6*j],&Cmz[6*j],r_Uranus);
    r_Uranus[0]*=1E3;
    r_Uranus[1]*=1E3;
    r_Uranus[2]*=1E3;

    for (j=0,i=404;i<=422;i+=6) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,6);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,6);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,6);
    j=0; mjd0_=t1;

    double r_Neptune[3];
    cheb3d(6,mjd_tdb,mjd0_,mjd0_+32.0,&Cmx[6*j],&Cmy[6*j],&Cmz[6*j],r_Neptune);
    r_Neptune[0]*=1E3;
    r_Neptune[1]*=1E3;
    r_Neptune[2]*=1E3;

    for (j=0,i=422;i<=440;i+=6) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,6);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,6);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,6);
    j=0; mjd0_=t1;

    double r_Pluto[3];
    cheb3d(6,mjd_tdb,mjd0_,mjd0_+32.0,&Cmx[6*j],&Cmy[6*j],&Cmz[6*j],r_Pluto);
    r_Pluto[0]*=1E3;
    r_Pluto[1]*=1E3;
    r_Pluto[2]*=1E3;
    for (j=0,i=818;i<=838;i+=10) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,10);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,10);

    for (i=1;i<=3;i++) {
        matcpy(Cmx+10*i,&jpleph_de440coeff[k][tmp[0]+20*i],1,10);
        matcpy(Cmy+10*i,&jpleph_de440coeff[k][tmp[1]+20*i],1,10);
    }
    memset(Cmz,0,sizeof(Cmz));

    if (dt>=0.0&&dt<=8.0) {
        j=0;
        mjd0_=t1;
    }
    else if (8.0<dt&&dt<=16.0) {
        j=1;
        mjd0_=t1+8.0*j;
    }
    else if (16.0<dt&&dt<=24.0) {
        j=2;
        mjd0_=t1+8.0*j;
    }
    else if (24.0<dt&&dt<=32.0) {
        j=2;
        mjd0_=t1+8.0*j;
    }
    double r_Nutations[3];
    cheb3d(10,mjd_tdb,mjd0_,mjd0_+8.0,&Cmx[10*j],&Cmy[10*j],Cmz,r_Nutations);

    for (j=0,i=898;i<=928;i+=10) {
        tmp[j++]=i;
    }
    matcpy(Cmx,&jpleph_de440coeff[k][tmp[0]],1,10);
    matcpy(Cmy,&jpleph_de440coeff[k][tmp[1]],1,10);
    matcpy(Cmz,&jpleph_de440coeff[k][tmp[2]],1,10);
    for (i=1;i<=3;i++) {
        matcpy(Cmx+10*i,&jpleph_de440coeff[k][tmp[0]+30*i],1,10);
        matcpy(Cmy+10*i,&jpleph_de440coeff[k][tmp[1]+30*i],1,10);
        matcpy(Cmz+10*i,&jpleph_de440coeff[k][tmp[2]+30*i],1,10);
    }
    if (dt>=0.0&&dt<=8.0) {
        j=0;
        mjd0_=t1;
    }
    else if (8.0<dt&&dt<=16.0) {
        j=1;
        mjd0_=t1+8.0*j;
    }
    else if (16.0<dt&&dt<=24.0) {
        j=2;
        mjd0_=t1+8.0*j;
    }
    else if (24.0<dt&&dt<=32.0) {
        j=2;
        mjd0_=t1+8.0*j;
    }
    double Librations[3];
    cheb3d(10,mjd_tdb,mjd0_,mjd0_+8.0,&Cmx[10*j],&Cmy[10*j],&Cmz[10*j],Librations);
    double EMRAT=81.3005682214972154;
    double EMRAT1=1.0/(1.0+EMRAT);
    double r_Earth[3];
    for (i=0;i<3;i++) r_Earth[i]=r_EarthMoon[i]-EMRAT1*r_Moon[i];

    if (mercury) for (i=0;i<3;i++) mercury[i]=-r_Earth[i]+r_Mercury[i];
    if (venus  ) for (i=0;i<3;i++) venus  [i]=-r_Earth[i]+r_Venus [i];
    if (sunssb ) for (i=0;i<3;i++) sunssb [i]=r_Sun[i];
    if (earth  ) for (i=0;i<3;i++) earth  [i]=r_Earth[i];
    if (moon   ) for (i=0;i<3;i++) moon   [i]=r_Moon[i];
    if (mars   ) for (i=0;i<3;i++) mars   [i]=-r_Earth[i]+r_Mars[i];
    if (jupiter) for (i=0;i<3;i++) jupiter[i]=-r_Earth[i]+r_Jupiter[i];
    if (saturn ) for (i=0;i<3;i++) saturn [i]=-r_Earth[i]+r_Satrun[i];
    if (uranus ) for (i=0;i<3;i++) uranus [i]=-r_Earth[i]+r_Uranus[i];
    if (neptune) for (i=0;i<3;i++) neptune[i]=-r_Earth[i]+r_Neptune[i];
    if (pluto  ) for (i=0;i<3;i++) pluto  [i]=-r_Earth[i]+r_Pluto[i];
    if (sun    ) for (i=0;i<3;i++) sun    [i]=-r_Earth[i]+r_Sun[i];
}
/* read jpleph-de440 coeffs-----------------------------------------------------*/
extern int readjplde440coeffs(const char *file)
{
    char buff[64*1018],*p,*q,*val[4096];
    int i,j=0,n;
    FILE *fp;

    if (jpleph_udflag) return j;

    if (!(fp=fopen(file,"r"))) {
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        for (n=0,p=buff;*p&&n<1018;p=q+1) {
            if ((q=strchr(p,','))||(q=strchr(p,'\n'))) {
                val[n++]=p; *q='\0';
            }
            else break;
        }
        if (n<1018) continue;
        for (i=0;i<n;i++) {
            jpleph_de440coeff[j][i]=atof(val[i]);
        }
        j++;
    }
    if (j) jpleph_udflag=1;
    fclose(fp); return j;
}
/* calender date/time convert to MJD--------------------------------------------
 * args   : double    *ep    I    calender date/time[year,month,day(,hour,min,sec)]
 * return : modefied julan date
 * -----------------------------------------------------------------------------*/
extern double cal2mjd(const double *ep)
{
    double djm0,djm;
    int year=ep[0],month=ep[1],day=ep[2];

    iauCal2jd(year,month,day,&djm0,&djm);
    djm+=ep[3]/24.0+ep[4]/1440.0+ep[5]/86400.0;
    return djm;
}
/* MJD convert to calender date/time--------------------------------------------*/
extern void mjd2cal(const double mjd, double *ep)
{
    int year,month,day;
    const static double mjd0=2400000.5;
    double fd;

    iauJd2cal(mjd0,mjd,&year,&month,&day,&fd);
    ep[0]=year;
    ep[1]=month;
    ep[2]=day;
    ep[3]=(int)(fd*86400.0/3600.0);
    ep[4]=(int)((fd*86400.0-ep[3]*3600.0)/60.0);
    ep[5]=fd*86400.0-ep[3]*3600.0-ep[4]*60.0;
}
/* transform ecef to satellite orbit coordinate---------------------------------
 * args:    double   *rs   satellite position/velocity(ecef/m)
 *          double   *E    transformation matrix ecef to satellite orbit coordinate
 *                         [r_radial,r_along-track,r_cross-track]=E*r
 * return none;
 *-----------------------------------------------------------------------------*/
extern void ecef2satf(const double *rs, double *E)
{
    double pos[3],vel[3],w[3]={0,0,OMGE};
    double crt[3],alt[3],a,b,c;
    int i;

    matcpy(pos,rs,1,3);
    cross3(w,pos,vel);
    for (i=0;i<3;i++) vel[i]+=rs[i+3];
    cross3(pos,vel,crt);
    cross3(crt,pos,alt);

    a=norm(pos,3);
    b=norm(alt,3);
    c=norm(crt,3);
    E[0]=pos[0]/a; E[3]=pos[1]/a; E[6]=pos[2]/a;
    E[1]=alt[0]/b; E[4]=alt[1]/b; E[7]=alt[2]/b;
    E[2]=crt[0]/c; E[5]=crt[1]/c; E[8]=crt[2]/c;
}
/* unit vector of satellite fixed coordinate------------------------------------*/
extern void satfixed(const double *rsat, const double *rsun, double *ex, double *ey, double *ez)
{
    double esun[3],ss[3],a;
    int i;

    for (i=0;i<3;i++) ss[i]=rsun[i]-rsat[i];
    for (i=0;i<3;i++) esun[i]=ss[i]/norm(ss,3);
    ez[0]=-rsat[0]/norm(rsat,3);
    ez[1]=-rsat[1]/norm(rsat,3);
    ez[2]=-rsat[2]/norm(rsat,3);
    cross3(ez,esun,ey);
    a=norm(ey,3);
    for (i=0;i<3;i++) ey[i]/=a;
    cross3(ey,ez,ex);
}
/* relativity correction--------------------------------------------------------*/
extern double relcorr(const double *rsat, const double *rrcv, int opt)
{
    const double GME_=3.986004415E+14; /* earth gravitational constant */
    double rels,rs,rr,rrs,rsr[3];

    /* satellite clock correction */
    rels=2.0/CLIGHT*dot(rsat,rsat+3,3);

    /* shapiro time delay correction */
    if (opt==1) {
        rs=norm(rsat,3);
        rr=norm(rrcv,3);
        rsr[0]=rsat[0]-rrcv[0];
        rsr[1]=rsat[1]-rrcv[1];
        rsr[2]=rsat[2]-rrcv[2];
        rrs=norm(rsr,3);
        rels+=2.0*GME_/CLIGHT/CLIGHT*log((rs+rr+rrs)/(rs+rr-rrs));
    }
    return rels;
}
/* MJD for barycentric dynamical time--------------------------------------------*/
extern double mjdtdb(const double mjd, const double *erpv)
{
    const double mjd_j2000=51544.5;
    double t,dt[6],mjd_tt;

    /* compute MJD(TT) */
    reftimediff(erpv[2],erpv[4],dt);
    mjd_tt=mjd+dt[3]/86400.0;
    t=(mjd_tt-mjd_j2000)/36525.0;

    /* compute MJD(TDB) */
    return mjd_tt+(0.001657*sin(628.3076*t+6.2401)+
                   0.000022*sin( 575.3385*t+4.2970)+
                   0.000014*sin(1256.6152*t+6.1969)+
                   0.000005*sin( 606.9777*t+4.0212)+
                   0.000005*sin(  52.9691*t+0.4444)+
                   0.000002*sin(  21.3299*t+5.5431)+
                   0.000010*sin( 628.3076*t+4.2490))/86400.0;
}
/* MJD for TT time---------------------------------------------------------------*/
extern double mjdtt(const double mjd, const double *erpv)
{
    double dt[6],mjd_tt;

    /* compute MJD(TT) */
    reftimediff(erpv[2],erpv[4],dt);
    mjd_tt=mjd+dt[3]/86400.0;
    return mjd_tt;
}
/* MJD for UT1 time--------------------------------------------------------------*/
extern double mjdut1(const double mjd, const double *erpv)
{
    double mjd_ut1;

    /* compute MJD(UT1) */
    mjd_ut1=mjd+erpv[2]/86400.0;
    return mjd_ut1;
}
/* eci to satellite-fixed coordinate transformation matrix -------------------
 * args:    double *rs  I  satellite position/velocity[x;y;z;vx;vy;vz](m) (eci)
 *          double *E   O  eci to satellite-fixed coordinate transformation matrix(3x3)
 *                         [r_radial;r_along-track;r_cross-track]=E*r  
 *----------------------------------------------------------------------------*/
extern void ecsf2satf(const double *rs, double *E)
{
	double crt[3],alt[3],nrad,ncrt,nalt;
	
	cross3(rs,rs+3,crt);
	cross3(crt,rs,alt);
	nrad=norm(rs,3);
	nalt=norm(alt,3);
	ncrt=norm(crt,3);
	E[0]=rs[0]/nrad;
	E[3]=rs[1]/nrad;
	E[6]=rs[2]/nrad;
	E[1]=alt[0]/nalt;
	E[4]=alt[1]/nalt;
	E[7]=alt[2]/nalt;
	E[2]=crt[0]/ncrt;
	E[5]=crt[1]/ncrt;
	E[8]=crt[2]/ncrt;
}
/* days of year-----------------------------------------------------------------*/
extern double doy(const double *ep)
{
    double month[12],days;
    int i;

    for (i=0;i<12;i++) {
        month[i]=31.0;
        if (i==2) month[i]=28.0;
        if (i==4||i==6||i==9||i==11) month[i]=30.0;
    }
    if (((int)ep[0])%4==0) {
        month[1]=29.0;
        if (((int)ep[0])%100==0&&((int)ep[0])%400!=0) month[1]=28.0;
    }
    i=1;
    days=0.0;

    while (i<((int)ep[1])&&i<12) {
        days+=month[i-1]; 
        i++;
    }
    return days+ep[2]+ep[3]/24.0+ep[4]/1440.0+ep[5]/86400.0;
}
/* time convert to MJD--------------------------------------------------------*/
extern double time2mjd(gtime_t tutc)
{
    double ep[6];
    time2epoch(tutc,ep);
    return cal2mjd(ep);
}
/* time convert to MJD(TT)----------------------------------------------------*/
extern double time2mjdtt(gtime_t tutc, const double *erpv)
{
    return mjdtt(time2mjd(tutc),erpv);
}
/* time convert to MJD(UT1)---------------------------------------------------*/
extern double time2mjdut1(gtime_t tutc, const double *erpv)
{
    return mjdut1(time2mjd(tutc),erpv);
}
/* trim space-----------------------------------------------------------------*/
static char *rtrim(char *str)
{
    if (str==NULL||*str=='\0') return str;
    int len=strlen(str); char *p=str+len-1;
    while (p>=str&&isspace(*p)) {*p='\0'; --p;}
    return str;
}
static char *ltrim(char *str)
{
    if (str==NULL||*str=='\0') return str;
    int len=0; char *p=str;
    while (*p!='\0'&&isspace(*p)) {++p; ++len;}
    memmove(str,p,strlen(str)-len+1);
    return str;
}
/* read gfc format earth gravity filed model data-----------------------------*/
extern int readgfc(const char *file, double **C, double **S, double *RE, double *GM, char *name, int *tidesyetem)
{
    FILE *fp;
    if (!(fp=fopen(file,"r"))) return -1;
    *C=*S=NULL;

    char buff[1024],*p,*q,*val[64],*str;
    int n,N=0;
    while (fgets(buff,sizeof(buff),fp)) {
        str=rtrim(buff);
        str=ltrim(buff);
        for (n=0,p=str;*p&&n<1018;p=q+1) {
            if ((q=strchr(p,' '))||(q=strchr(p,'\n'))||(q=strchr(p,'\0'))) {
                val[n]=p; *q='\0';
                if (strcmp(val[n],"")!=0) n++;
            }
            else break;
        }
        if (n<2) continue;

        if      (strstr(val[0],"modelname")) strcpy(name,val[1]);
        else if (strstr(val[0],"earth_gravity_constant")) *GM=atof(val[1]);
        else if (strstr(val[0],"radius")) *RE=atof(val[1]);
        else if (strstr(val[0],"max_degree")) N=atoi(val[1]);
        else if (strstr(val[0],"tide_system")) {
            if      (strstr(val[1],"zero")) *tidesyetem=GRAVITY_MODEL_TIDE_ZERO;
            else if (strstr(val[1],"free")) *tidesyetem=GRAVITY_MODEL_TIDE_FREE;
        }
        if (N>0) {
            if (*C==NULL) *C=zeros(N+1,N+1);
            if (*S==NULL) *S=zeros(N+1,N+1);
        }
        if (!*C||!*S) continue;
        if (strstr(val[0],"gfc")) {
            int i=atoi(val[1]);
            int j=atoi(val[2]);
            (*C)[j+i*(N+1)]=atof(val[3]);
            (*S)[j+i*(N+1)]=atof(val[4]);
        }
    }
    fclose(fp); return N;
}
/* read EOP data from http://celestrak.org/SpaceData/EOP-All.txt--------------*/
static int readeop_celst(const char *file, erp_t *erp)
{
    erpd_t *erp_data;
    FILE *fp;
    double v[14]={0};
    char buff[1024];

    if (!(fp=fopen(file,"r"))) return 0;

    erp->n=0;
    while (fgets(buff,sizeof(buff),fp)) {
        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10,v+11,v+12)<13) {
            continue;
        }
        if (erp->n>=erp->nmax) {
            erp->nmax=erp->nmax<=0?128:erp->nmax*2;
            erp_data=(erpd_t *)realloc(erp->data,sizeof(erpd_t)*erp->nmax);
            if (!erp_data) {
                free(erp->data); erp->data=NULL; erp->n=erp->nmax=0;
                fclose(fp);
                return 0;
            }
            erp->data=erp_data;
        }
        erp->data[erp->n].mjd=v[3];
        erp->data[erp->n].xp=v[4]*AS2R;
        erp->data[erp->n].yp=v[5]*AS2R;
        erp->data[erp->n].ut1_utc=v[6];
        erp->data[erp->n].lod=v[7];
        erp->data[erp->n].tai_utc=v[12];
        erp->n++;
    }
    fclose(fp);
    initeop(file,geoparr,&gjdeopstart);
    return erp->n;
}
/* read EOP data via a IGS file-----------------------------------------------*/
static int readeop_igsv2(const char *file, erp_t *erp)
{
    return readerp(file,erp);
}
/* read EOP data form file----------------------------------------------------*/
extern int readeop(const char *file, erp_t *erp, int type)
{
    switch (type) {
        case 0: return readeop_celst(file,erp);
        case 1: return readeop_igsv2(file,erp);
    }
}
/* get earth rotation parameters----------------------------------------------*/
extern int geteop(const erp_t *erp, gtime_t time, erpd_t *erpd)
{
    double mjd_utc=floor(time2mjd(time));
    double fixf=time2mjd(time)-mjd_utc;
    int i,j,k;

    for (j=-1,k=erp->n-1;j<k-1;) {
        i=(j+k)/2;
        if (mjd_utc<erp->data[i].mjd) k=i;
        else j=i;
    }
    if (j<0) return 0;
    erpd->xp=erp->data[j].xp+(erp->data[MIN(j+1,erp->n-1)].xp-erp->data[j].xp)*fixf;
    erpd->yp=erp->data[j].yp+(erp->data[MIN(j+1,erp->n-1)].yp-erp->data[j].yp)*fixf;
    erpd->ut1_utc=erp->data[j].ut1_utc+(erp->data[MIN(j+1,erp->n-1)].ut1_utc-erp->data[j].ut1_utc)*fixf;
    erpd->lod=erp->data[j].lod+(erp->data[MIN(j+1,erp->n-1)].lod-erp->data[j].lod)*fixf;
    erpd->tai_utc=erp->data[j].tai_utc;
    erpd->mjd=time2mjd(time);
    return 1;
}
/* read space weather data----------------------------------------------------*/
extern int readspw(const char *file, spw_t *spw)
{
    FILE *infile;
    char longstr[140],str[9],blk[30],buff[1024];
    int numrecsobs=0,numrecspred,i,year,mon,day,tmp;
    spwd_t *spw_data;

    if (!(infile=fopen(file,"r"))) return 0;

    while (fgets(buff,sizeof(buff),infile)) {
        if (strstr(buff,"NUM_OBSERVED_POINTS")) {
            sscanf(buff,"%s %i",blk,&numrecsobs);
            break;
        }
    }
    if (numrecsobs<=0) return 0;
    spw->data=calloc(numrecsobs,sizeof(spwd_t));

    /* find epoch date */
    fgets(longstr,140,infile);

    /* process observed records */
    for (i=0;i<numrecsobs;i++) {
        /* use d format for integers with leading 0's */
        fscanf(infile,"%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4lf %2i %4i %6lf %2i %6lf %6lf %6lf %6lf %6lf \n",
               &spw->data[i].year,&spw->data[i].mon,&spw->data[i].day,&spw->data[i].bsrn,&spw->data[i].nd,
               &spw->data[i].kparr[0],&spw->data[i].kparr[1],&spw->data[i].kparr[2],&spw->data[i].kparr[3],
               &spw->data[i].kparr[4],&spw->data[i].kparr[5],&spw->data[i].kparr[6],&spw->data[i].kparr[7],&spw->data[i].sumkp,
               &spw->data[i].aparr[0],&spw->data[i].aparr[1],&spw->data[i].aparr[2],&spw->data[i].aparr[3],
               &spw->data[i].aparr[4],&spw->data[i].aparr[5],&spw->data[i].aparr[6],&spw->data[i].aparr[7],&spw->data[i].avgap,
               &spw->data[i].cp,&spw->data[i].c9,&spw->data[i].isn,&spw->data[i].adjf10,&spw->data[i].q,
               &spw->data[i].adjctrf81,&spw->data[i].adjlstf81,&spw->data[i].obsf10,&spw->data[i].obsctrf81,&spw->data[i].obslstf81);
    }
    fgets(longstr,140,infile);
    fgets(longstr,140,infile);
    fscanf(infile,"%s %i ",blk,&numrecspred);
    fgets(longstr,140,infile);

    spw_data=realloc(spw->data,(numrecsobs+numrecspred)*sizeof(spwd_t));

    if (spw_data==NULL) {
        free(spw->data); spw->n=spw->nmax=0;
        return 0;
    }
    spw->data=spw_data;

    /* process predicted records */
    for (i=numrecsobs;i<numrecsobs+numrecspred;i++) {
        /* use d format for integers with leading 0's */
        fscanf(infile,"%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4lf %2i %4i %6lf %6lf %6lf %6lf %6lf %6lf \n",
               &spw->data[i].year,&spw->data[i].mon,&spw->data[i].day,&spw->data[i].bsrn,&spw->data[i].nd,
               &spw->data[i].kparr[0],&spw->data[i].kparr[1],&spw->data[i].kparr[2],&spw->data[i].kparr[3],
               &spw->data[i].kparr[4],&spw->data[i].kparr[5],&spw->data[i].kparr[6],&spw->data[i].kparr[7],&spw->data[i].sumkp,
               &spw->data[i].aparr[0],&spw->data[i].aparr[1],&spw->data[i].aparr[2],&spw->data[i].aparr[3],
               &spw->data[i].aparr[4],&spw->data[i].aparr[5],&spw->data[i].aparr[6],&spw->data[i].aparr[7],&spw->data[i].avgap,
               &spw->data[i].cp,&spw->data[i].c9,&spw->data[i].isn,&spw->data[i].adjf10,
               &spw->data[i].adjctrf81,&spw->data[i].adjlstf81,&spw->data[i].obsf10,&spw->data[i].obsctrf81,&spw->data[i].obslstf81);
        spw->data[i].q=0;
    }
    spw->n=spw->nmax=numrecsobs+numrecspred;
    fclose(infile);

    initspw(file,gspwarr,&gjdspwstart);
    return spw->n;
}
/* string is number-----------------------------------------------------------*/
static int isdigitstr(char *str)
{
    return (strspn(str,"0123456789")==strlen(str));
}
/* add ocean tide model data -------------------------------------------------*/
static int addotmodeldata(oceantide_model_t *otmdl, const oceantide_model_coeffs_t *data)
{
    oceantide_model_coeffs_t *mdldata;

    if (otmdl->nmax<=otmdl->n) {
        if (otmdl->nmax<=0) otmdl->nmax=1024;
        else otmdl->nmax*=2;
        if (!(mdldata=(oceantide_model_coeffs_t *)realloc(otmdl->mdldata,sizeof(oceantide_model_coeffs_t)*otmdl->nmax))) {
            free(otmdl->mdldata);
            otmdl->mdldata=NULL;
            otmdl->n=otmdl->nmax=0;
            return -1;
        }
        otmdl->mdldata=mdldata;
    }
    otmdl->mdldata[otmdl->n++]=*data;
    return 1;
}
/* read ocean tide model coeffs-----------------------------------------------*/
extern int readotfile(const char *file, oceantide_model_t *otmdl)
{
    otmdl->n=0;
    FILE *fp;
    if (!(fp=fopen(file,"r"))) return -1;

    char buff[1024],*p,*q,*val[64],t[8];
    int len,n;
    while (fgets(buff,sizeof(buff),fp)) {
        for (n=0,p=buff;*p&&n<1018;p=q+1) {
            if ((q=strchr(p,' '))||(q=strchr(p,'\n'))) {
                val[n]=p; *q='\0';
                if (strcmp(val[n],"")!=0) n++;
            }
            else break;
        }
        if (n<8) continue;

        /* doodson number, coded into doodson variable multipliers */
        oceantide_model_coeffs_t coeffs={0};

        if (!isdigitstr(val[2])||!isdigitstr(val[3])) continue;

        /* darwin's symbol */
        strcpy(coeffs.darw,val[1]);

        len=strlen(val[0]);
        if (len<=0) continue;
        sprintf(t,"%c",val[0][len-1]); coeffs.dn[5]=atoi(t);
        sprintf(t,"%c",val[0][len-2]); coeffs.dn[4]=atoi(t);
        sprintf(t,"%c",val[0][len-3]); coeffs.dn[3]=atoi(t);
        sprintf(t,"%c",val[0][len-5]); coeffs.dn[2]=atoi(t);
        sprintf(t,"%c",val[0][len-6]); coeffs.dn[1]=atoi(t);
        if (len>=7) {
            sprintf(t,"%c",val[0][len-7]); coeffs.dn[0]=atoi(t);
        }
        /* degree/order */
        coeffs.l=atoi(val[2]);
        coeffs.m=atoi(val[3]);

        /* coefficients */
        coeffs.dCp=atof(val[4])*1E-12;
        coeffs.dSp=atof(val[5])*1E-12;
        coeffs.dCm=atof(val[6])*1E-12;
        coeffs.dSm=atof(val[7])*1E-12;

        if (coeffs.l<=otmdl->maxl) {
            if (otmdl->maxm>0) if (coeffs.m>otmdl->maxm) continue;
            if (addotmodeldata(otmdl,&coeffs)==-1) continue;
        }
    }
    fclose(fp); return otmdl->n;
}
/* satellite position/velocity in ECI from ECEF-----------------------------*/
extern int satposeci(gtime_t tutc, const erp_t *erp, const double *rs_ecef, double *rs_eci)
{
    double erpv[6],U[9],dU[9];

    /* earth rotation parameter values */
    geterp(erp,utc2gpst(tutc),erpv);
    if (norm(erpv,6)<=0.0) return 0;

    /* satellite position/velocity in ECI(m/s) */
    ecsf2ecef(tutc,erpv,U,dU,NULL);
    matmul("TN",3,1,3,1.0,U,rs_ecef,0.0,rs_eci);

    matmul("TN",3,1,3,1.0,U,rs_ecef+3,0.0,rs_eci+3);
    matmul("TN",3,1,3,1.0,dU,rs_ecef,1.0,rs_eci+3);
    return 1;
}
/* satellite position/velocity in ECEF from ECI-----------------------------*/
extern int satposecef(gtime_t tutc, const erp_t *erp, const double *rs_eci, double *rs_ecef)
{
    double erpv[6],U[9],dU[9];

    /* earth rotation parameter values */
    geterp(erp,utc2gpst(tutc),erpv);
    if (norm(erpv,6)<=0.0) return 0;

    /* satellite position/velocity in ECEF(m/s) */
    ecsf2ecef(tutc,erpv,U,dU,NULL);
    matmul("NN",3,1,3,1.0,U,rs_eci,0.0,rs_ecef);

    matmul("NN",3,1,3,1.0,U,rs_eci+3,0.0,rs_ecef+3);
    matmul("NN",3,1,3,1.0,dU,rs_eci,1.0,rs_ecef+3);
    return 1;
}
/* update POD receiver station position from SNX solutions------------------*/
extern void podudrcvsnx(pod_t *pod, const sta_t *stas, int nsta, const char *snxfile)
{
    if (strcmp(snxfile,"")==0) return;

    int i,j;
    for (i=0;i<nsta;i++) {
        char staname[8]={0};
        sol_t snxsol={0};

        if (strcmp(stas[i].name,"")==0) break;

        for (j=0;j<GNRCVS;j++) {
            strncpy(staname,stas[i].name,4);

            if (strcmp(pod->rcv[j].name,staname)==0) {
                pod->rcv[j].sta=stas[i];
                readsnx(snxfile,staname,&snxsol);
                matcpy(pod->rcv[pod->rcv[j].id-1].pos,snxsol.rr,1,3);
                break;
            }
        }
    }
}
/* satellite fixed coordinate to eci transformation matrix -------------------*/
extern void satf2ecsf(const double *rsat, const double *rsun, double *E)
{
    double esun[3],*ex=E,*ey=E+3,*ez=E+6,en;

    esun[0]=rsun[0]-rsat[0];
    esun[1]=rsun[1]-rsat[1];
    esun[2]=rsun[2]-rsat[2];

    en=norm(esun,3);
    esun[0]/=en;
    esun[1]/=en;
    esun[2]/=en;

    ez[0]=-rsat[0];
    ez[1]=-rsat[1];
    ez[2]=-rsat[2];
    en=norm(ez,3);
    ez[0]/=en;
    ez[1]/=en;
    ez[2]/=en;
    cross3(ez,esun,ey);

    en=norm(ey,3);
    ey[0]/=en;
    ey[1]/=en;
    ey[2]/=en;
    cross3(ey,ez,ex);
}
/* read snx solutions----------------------------------------------------------*/
extern int readsnx(const char *file, const char *staname, sol_t *sol)
{
    FILE *fp;
    sol_t sol0={0};

    *sol=sol0;
    if (!(fp=fopen(file,"r"))) return 0;

    char buff[1024];
    int flag=0;

    while (fgets(buff,sizeof(buff),fp)) {
        if (strstr(buff,"+SOLUTION/ESTIMATE")) flag=1;
        if (strstr(buff,"-SOLUTION/ESTIMATE")) break;

        if (flag) {
            char name[8],type[8],pt[8],refep[32],unit[8];
            double val,std;
            int n,soln,s;

            if (sscanf(buff,"%d %s %s %s %d %s %s %d %lf %lf\n",&n,type,name,pt,&soln,refep,unit,&s,&val,&std)<8) continue;
            if (!strstr(name,staname)) continue;

            if      (strstr(type,"STAX")) sol->rr[0]=val,sol->qr[0]=(float)std;
            else if (strstr(type,"STAY")) sol->rr[1]=val,sol->qr[1]=(float)std;
            else if (strstr(type,"STAZ")) sol->rr[2]=val,sol->qr[2]=(float)std;
        }
    }
    fclose(fp);
    return norm(sol->rr,3)>0.0;
}
/* add POD observation data --------------------------------------------------*/
static int addpodobs(podobss_t *obs, const podobs_t *data)
{
    podobs_t *obs_data;

    if (obs->nmax<=obs->n) {
        if (obs->nmax<=0) obs->nmax=1024; else obs->nmax*=2;
        if (!(obs_data=(podobs_t *)realloc(obs->data,sizeof(podobs_t)*obs->nmax))) {
            free(obs->data);
            obs->data=NULL;
            obs->n=obs->nmax=0;
            return -1;
        }
        obs->data=obs_data;
    }
    obs->data[obs->n++]=*data;
    return 1;
}
/* read POD satellite position observation from SP3 file----------------------*/
extern int readpodobs_sp3(const char **files, int n, podobss_t *obss)
{
    double rs[6],dts[2],var,tint=300.0;
    int i,j,k;

    for (i=0;i<n;i++) {
        nav_t nav={0};

        readsp3(files[i],&nav,0);
        uniqnav(&nav);

        for (j=0;j<nav.ne;j++) {
            for (k=0;k<MAXSAT;k++) {
                if (!peph2pos(nav.peph[j].time,k+1,&nav,0,rs,dts,&var)) continue;
                matcpy(nav.peph[j].vel[k],rs+3,1,3);
            }
            gtime_t t0={0};
            if (!screent(nav.peph[j].time,t0,t0,tint)) continue;

            podobs_t obs={0};
            obs.type=GPOD_OBSS_TYPE_SATPOS;
            for (k=0;k<MAXSAT;k++) {
                obs.tutc=gpst2utc(nav.peph[j].time);
                matcpy(obs.pos[k],nav.peph[j].pos[k],1,3);
                matcpy(obs.vel[k],nav.peph[j].vel[k],1,3);

                obs.pos[k][3]=SQR(nav.peph[j].std[k][0]?nav.peph[j].std[k][0]:0.05);
                obs.pos[k][4]=SQR(nav.peph[j].std[k][1]?nav.peph[j].std[k][1]:0.05);
                obs.pos[k][5]=SQR(nav.peph[j].std[k][2]?nav.peph[j].std[k][2]:0.05);
                obs.vel[k][3]=SQR(nav.peph[j].vst[k][0]?nav.peph[j].vst[k][0]:0.05);
                obs.vel[k][4]=SQR(nav.peph[j].vst[k][1]?nav.peph[j].vst[k][1]:0.05);
                obs.vel[k][5]=SQR(nav.peph[j].vst[k][2]?nav.peph[j].vst[k][2]:0.05);
            }
            addpodobs(obss,&obs);
        }
        freenav(&nav,0xFF);
    }
    return sortpodobs(obss);
}
/* compare POD observation data ----------------------------------------------*/
static int cmppodobs(const void *p1, const void *p2)
{
    podobs_t *q1=(podobs_t*)p1,*q2=(podobs_t*)p2;
    double tt=timediff(q1->tutc,q2->tutc);
    if (fabs(tt)>DTTOL) return tt<0?-1:1;
    else return 0;
}
/* sort and unique POD observation data---------------------------------------*/
extern int sortpodobs(podobss_t *obss)
{
    int i,j,n;

    if (obss->n<=0) return 0;

    qsort(obss->data,obss->n,sizeof(podobs_t),cmppodobs);

    /* delete duplicated data */
    for (i=j=0;i<obss->n;i++) {
        if (timediff(obss->data[i].tutc,obss->data[j].tutc)!=0.0) {
            obss->data[++j]=obss->data[i];
        }
    }
    obss->n=j+1;

    for (i=n=0;i<obss->n;i=j,n++) {
        for (j=i+1;j<obss->n;j++) {
            if (timediff(obss->data[j].tutc,obss->data[i].tutc)>DTTOL) break;
        }
    }
    return n;
}
/* read POD satellite code/phase observation from RINEX file------------------*/
extern int readpodobs_rnx(const char **files, int n, podobss_t *obss)
{
    obs_t obsrnxs[MAXRCV]={0};
    int nepoch[MAXRCV];
    int i,j,k,ti=30;

    for (i=0;i<n&&i<MAXRCV;i++) {
        readrnx(files[i],i+1,"",&obsrnxs[i],NULL,&obss->stas[i]);
        obss->stas[i].staid=i+1;
        nepoch[i]=sortobs(&obsrnxs[i]);
    }
    gtime_t tmin={0};
    gtime_t tmax={0};
    for (i=0;i<n;i++) {
        if (obsrnxs[i].n<=0) continue;
        if (!tmin.time||timediff(tmin,obsrnxs[i].data[0].time)<0.0) {
            tmin=obsrnxs[i].data[0].time;
        }
        if (!tmax.time||timediff(tmax,obsrnxs[i].data[obsrnxs[i].n-1].time)<0.0) {
            tmax=obsrnxs[i].data[obsrnxs[i].n-1].time;
        }
    }
    if (!tmin.time) return 0;
    if (!tmax.time) return 0;

    gtime_t tcur=tmin;
    for (;;) {
        if (timediff(tcur,tmax)>0.0) break;

        podobs_t podobs={0};

        for (j=0;j<n;j++) {
            obsd_t obsd[MAXOBS]={0};
            int nobs=0;

            for (k=0;k<obsrnxs[j].n;k++) {
                if (fabs(timediff(obsrnxs[j].data[k].time,tcur))<1E-2) {
                    obsd[nobs++]=obsrnxs[j].data[k];
                }
            }
            podobs.nobs[j]=nobs;
            for (k=0;k<nobs;k++) {
                podobs.obs[j][k]=obsd[k];
            }
        }
        podobs.tutc=gpst2utc(tcur);
        podobs.type=GPOD_OBSS_TYPE_SATOBS;
        podobs.nrcv=n;
        addpodobs(obss,&podobs);
        tcur=timeadd(tcur,ti);
    }
    for (i=0;i<n;i++) freeobs(&obsrnxs[i]);
    return sortpodobs(obss);
}
extern void readpodobs_rnxe(const char **files, int n, podobss_t *obss)
{
    int i;

    for (i=0;i<n&&i<GNRCVS;i++) {
        readrnx(files[i],i+1,"",&obss->obss[i],NULL,&obss->stas[i]);
        obss->stas[i].staid=i+1;
        sortobs(&obss->obss[i]);
    }
    obss->n=n;
}
extern void readpodobs_rnxet(const char **files, int n, podobss_t *obss, gtime_t ts, gtime_t te, double tint)
{
    int i;

    for (i=0;i<n&&i<GNRCVS;i++) {
        readrnxt(files[i],i+1,ts,te,tint,"",&obss->obss[i],NULL,&obss->stas[i]);
        obss->stas[i].staid=i+1;
        sortobs(&obss->obss[i]);
    }
    obss->n=n;
}
/* get POD satellite code/phase observation---------------------------------*/
extern int podobssget(const podobss_t *obss, gtime_t tc, podobs_t *podobs)
{
    int j,k,obsc=0,nrcv=0;

    memset(podobs,0,sizeof(*podobs));

    for (j=0;j<obss->n;j++) {
        obsd_t obsd[MAXOBS]={0};
        int nobs=0;

        for (k=0;k<obss->obss[j].n;k++) {
            if (fabs(timediff(obss->obss[j].data[k].time,tc))<1E-2) {
                obsd[nobs++]=obss->obss[j].data[k];
            }
        }
        if (nobs<=0) continue;
        obsc+=nobs;

        podobs->nobs[nrcv]=nobs;
        podobs->ircv[nrcv]=obss->stas[j].staid;
        strncpy(podobs->name[nrcv],obss->stas[j].name,4);
        for (k=0;k<nobs;k++) {
            podobs->obs[nrcv][k]=obsd[k];
        }
        nrcv++;
    }
    if (obsc<=0) return 0;
    podobs->tutc=gpst2utc(tc);
    podobs->type=GPOD_OBSS_TYPE_SATOBS;
    podobs->nrcv=nrcv;
    return 1;
}