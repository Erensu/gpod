#include "src/podlib.h"

extern void openmatlib();

int main(int argc, const char *argv[])
{
    openmatlib();

    const char *eopfile="data\\model\\EOP-All.txt";
    const char *spwfile="data\\model\\SW-All.txt";
    const char *gfcfile="data\\model\\GGM03C.gfc";
    const char *jplde440file="data\\model\\de440coeff";
    const char *octfile="data\\model\\fes2004_Cnm-Snm.dat";
    const char *snxfile="data\\model\\igs20p21130.snx";
    const char *navfile="data\\nav\\BRDC00IGS_R_20201870000_01D_MN.rnx";
    const char *navfile1="data\\nav\\BRDC00IGS_R_20201880000_01D_MN.rnx";
    const char *sp3file[3]={
            "D:\\data\\nav\\igs21130.sp3",
            "D:\\data\\nav\\igs21131.sp3",
            "D:\\data\\nav\\igs21132.sp3"
    };
    char *solfile="pode.sol";
    const char *obsbinfile="data\\obs\\podobs_2020187.bin";
    const char *stabinfile="data\\obs\\podsta_2020187.bin";
    FILE *fp=fopen(solfile,"w");
    FILE *fpobs=fopen(obsbinfile,"rb");
    FILE *fpsta=fopen(stabinfile,"rb");
    int nstation=124;
    sta_t stas[256]={0};

    for (int i=0;i<nstation;i++) {
        fread(stas+i,sizeof(sta_t),1,fpsta);
    }
    fclose(fpsta);

    force_model_opt_t fmdlopt={0};
    fmdlopt.gravity=1;
    fmdlopt.nbody=1;
    fmdlopt.solidtide=1;
    fmdlopt.solrad=1;
    fmdlopt.poletide=0;
    fmdlopt.relativity=1;
    fmdlopt.atmdrag=1;
    fmdlopt.shmopt=SHADOW_CYLINDRICAL;
    fmdlopt.atmdensopt=ATMODENSITY_NRL;
    fmdlopt.gM=10;
    fmdlopt.gN=10;
    fmdlopt.oM=10;
    fmdlopt.oL=10;
    fmdlopt.oceantide=OCEANTIDECORR_DEFAULT;
    fmdlopt.anelasticearth=0;
    fmdlopt.udtransmatrix=UDTRANSMATRIX_INTGRTN;
    fmdlopt.usesamegravitymdl=1;
    fmdlopt.usesameeop=1;
    fmdlopt.usesamespw=1;
    fmdlopt.usesameotm=1;

    fmdlopt.atmdrgparas[0]=55.64;
    fmdlopt.atmdrgparas[1]=8000;
    fmdlopt.atmdrgparas[2]=2.7;

    fmdlopt.solradparas[0]=55.64;
    fmdlopt.solradparas[1]=8000;
    fmdlopt.solradparas[2]=1.0;
    fmdlopt.solradparas[3]=4.55982118135874e-06;
    fmdlopt.solradparas[4]=149597870700.000;
    fmdlopt.solradmdl=SOLRAD_MODEL_ECOMC;
    fmdlopt.udacclgrad=ACCLGRAD_ANALYSIS;
    fmdlopt.D0=-1.0;
    fmdlopt.Dc=1E-30;
    fmdlopt.Ds=1E-30;
    fmdlopt.B0=1E-30;
    fmdlopt.Bc=1E-30;
    fmdlopt.Bs=1E-30;
    fmdlopt.Y0=1E-30;
    fmdlopt.Yc=1E-30;
    fmdlopt.Ys=1E-30;

    fmdlopt.Dc2=1E-30;
    fmdlopt.Dc4=1E-30;
    fmdlopt.Ds2=1E-30;
    fmdlopt.Ds4=1E-30;

    strcpy(fmdlopt.eopfile,eopfile);
    strcpy(fmdlopt.spwfile,spwfile);
    strcpy(fmdlopt.gfcfile,gfcfile);
    strcpy(fmdlopt.octfile,octfile);
    strcpy(fmdlopt.jplde440file,jplde440file);

    pod_t pod={0};
    pod_opt_t opt={0};

    int i,j;
    for (i=GIP;i<GIP+GNP;i++) opt.satvar[i]=SQR(5.0);
    for (i=GIV;i<GIV+GNV;i++) opt.satvar[i]=SQR(0.0001);
    for (i=GIS;i<GIS+GNS;i++) opt.satvar[i]=SQR(0.1);
    opt.satvar[GIC+0]=SQR(3.0);
    opt.satvar[GIC+1]=SQR(0.0001);

    for (i=GIP;i<GIP+GNP;i++) opt.satprn[i]=1E-30;
    for (i=GIV;i<GIV+GNV;i++) opt.satprn[i]=1E-15;
    for (i=GIS;i<GIS+GNS;i++) opt.satprn[i]=1E-15;
    opt.satprn[GIC+0]=1.0;
    opt.satprn[GIC+1]=1E-8;

    for (i=GIRP;i<GIRP+GNRP;i++) opt.rcvvar[i]=SQR(1E-3);
    for (i=GIRC;i<GIRC+GNRC;i++) opt.rcvvar[i]=SQR(60.0);
    for (i=GIRD;i<GIRD+GNRCR;i++) opt.rcvvar[i]=SQR(1.0);
    for (i=GITR;i<GITR+GNRTR;i++) opt.rcvvar[i]=SQR(0.5);

    strcpy(opt.satantp,"data\\model\\igs_absolute_14.atx");
    strcpy(opt.rcvantp,"data\\model\\igs_absolute_14.atx");
    strcpy(opt.dcb,"data\\model\\P1C12004.DCB");
    strcpy(opt.blq,"data\\model\\FES2014b.BLQ");
    strcpy(opt.snxfile,snxfile);
    strcpy(opt.staposfile,"");

    log_trace_open("");
    log_set_level(1);

    opt.udorbitint=ORBITINT_ODEINTV2_RKF78;
    strcpy(opt.navfile,navfile);

    for (i=0;i<MAXSAT;i++) {
        opt.fmdlopts[i]=fmdlopt;
    }
    opt.estopt=prcopt_default_ppp;
    opt.estopt.sateph=EPHOPT_BRDC;
    opt.estopt.mode=PMODE_FIXED;
    opt.estopt.ionoopt=IONOOPT_IFLC;
    opt.estopt.tropopt=TROPOPT_EST;
    opt.estopt.nf=2;
    opt.estopt.elmin=10.0*D2R;
    opt.estopt.navsys=SYS_GPS;
    opt.estopt.posopt[0]=1;
    opt.estopt.posopt[2]=1;
    opt.estopt.posopt[3]=1;
    opt.estopt.posopt[4]=1;
    opt.estopt.posopt[5]=1;
    opt.estopt.tidecorr=7;
    opt.estmode=GPOD_FIXRCVPOS_ESTSATPOS;
    opt.outldetexcs=0;
    opt.satclkmode=1;
    strncpy(opt.clkref,stas[3].name,4);

    int backfilter=1;

    for (i=0;i<nstation;i++) {
        opt.podrcvs[i]=stas[i];
        strncpy(opt.podrcvs[i].name,stas[i].name,4);
    }
    podinit(&pod,&opt);
    readrnx(navfile1,0,"",NULL,&pod.nav,NULL);
    uniqnav(&pod.nav);

    for (i=0;i<3;i++) {
        readsp3(sp3file[i],&pod.nav,0);
    }
    log_trace_open(".\\pod_trace.out");
    log_set_level(1);

    static double drtns[3][MAXSAT][3][1024]={{0.0}};
    static int nrtns[3][MAXSAT]={0};

    for (;;) {

        podobs_t podobs={0};

        if (!fread(&podobs,sizeof(podobs),1,fpobs)) break;

        podflt(&pod,&podobs);

        fprintf(stderr,"POD processing: %s nrcv=%3d\n",time_str(utc2gpst(podobs.tutc),2),podobs.nrcv);

        if (pod.opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
            fprintf(fp,"\ntime=%s\n",time_str(utc2gpst(pod.sol.tutc),3));
            for (i=0;i<MAXSAT;i++) {
                char prn[8];
                satno2id(i+1,prn);

                double rs[6],dts[2],rs_eci[6],var,relclk;

                peph2pos(utc2gpst(pod.sol.tutc),i+1,&pod.nav,0,rs,dts,&var);
                satposeci(pod.sol.tutc,&pod.orbits[0].fmdl.erp,rs,rs_eci);

                relclk=relcorr(rs_eci,NULL,0);

                double E[9],dxyz[3],drtn[3];

                dxyz[0]=pod.sol.xs[i][0]?pod.sol.xs[i][0]-rs_eci[0]:0.0;
                dxyz[1]=pod.sol.xs[i][1]?pod.sol.xs[i][1]-rs_eci[1]:0.0;
                dxyz[2]=pod.sol.xs[i][2]?pod.sol.xs[i][2]-rs_eci[2]:0.0;

                ecsf2satf(rs,E);
                matmul("NN",3,1,3,1.0,E,dxyz,0.0,drtn);

                drtns[0][i][0][nrtns[0][i]]=drtn[0];
                drtns[0][i][1][nrtns[0][i]]=drtn[1];
                drtns[0][i][2][nrtns[0][i]]=drtn[2];
                nrtns[0][i]++;

                fprintf(fp,"%s  [%6.2lf,%3d]",prn,pod.sol.dops[i],pod.sol.ns[i]);
                for (j=0;j<6;j++) fprintf(fp,"%15.4lf  [%6.3lf] ",pod.sol.xs[i][j]?pod.sol.xs[i][j]-rs_eci[j]:0.0,sqrt(pod.sol.Qs[i][j]));
                fprintf(fp,"%15.4lf [%6.3lf]  %15.4lf ",
                        pod.sol.xs[i][GIC ]?pod.sol.xs[i][GIC ]-relclk-dts[0]*CLIGHT:0.0,sqrt(pod.sol.Qs[i][GIC]),
                        pod.sol.xs[i][GICR]?pod.sol.xs[i][GICR]-dts[1]*CLIGHT:0.0);

                for (j=6;j<GNX;j++) fprintf(fp,"%15.4lf ",pod.sol.xs[i][j]);
                fprintf(fp,"\n");
                fflush(fp);
            }
        }
        else {
            fprintf(fp,"\ntime=%s\n",time_str(utc2gpst(pod.sol.tutc),3));
            for (i=0;i<31;i++) {

                fprintf(fp,"%s  ",pod.rcv[i+1].name);
                for (j=0;j<3;j++) fprintf(fp,"%15.4lf ",pod.sol.xr[i+1][j]-pod.rcv[i+1].pos[j]);
                fprintf(fp,"\n");
                fflush(fp);
            }
        }
    }
    if (backfilter) {
        fseeko64(fpobs,0,SEEK_END);

        for (;;) {
            podobs_t podobs={0};

            if (fseeko64(fpobs,(_off64_t)-sizeof(podobs),SEEK_CUR)<0) break;

            if (!fread(&podobs,sizeof(podobs),1,fpobs)) break;

            podflt(&pod,&podobs);

            if (fseeko64(fpobs,(_off64_t)-sizeof(podobs),SEEK_CUR)<0) break;

            fprintf(stderr,"POD processing: %s nrcv=%3d\n",time_str(utc2gpst(podobs.tutc),2),podobs.nrcv);

            if (pod.opt.estmode==GPOD_FIXRCVPOS_ESTSATPOS) {
                fprintf(fp,"\ntime=%s\n",time_str(utc2gpst(pod.sol.tutc),3));
                for (i=0;i<MAXSAT;i++) {
                    char prn[8];
                    satno2id(i+1,prn);

                    double rs[6],dts[2],rs_eci[6],var,relclk;

                    peph2pos(utc2gpst(pod.sol.tutc),i+1,&pod.nav,0,rs,dts,&var);
                    satposeci(pod.sol.tutc,&pod.orbits[0].fmdl.erp,rs,rs_eci);

                    relclk=relcorr(rs_eci,NULL,0);

                    double E[9],dxyz[3],drtn[3];

                    dxyz[0]=pod.sol.xs[i][0]?pod.sol.xs[i][0]-rs_eci[0]:0.0;
                    dxyz[1]=pod.sol.xs[i][1]?pod.sol.xs[i][1]-rs_eci[1]:0.0;
                    dxyz[2]=pod.sol.xs[i][2]?pod.sol.xs[i][2]-rs_eci[2]:0.0;

                    ecsf2satf(rs,E);
                    matmul("NN",3,1,3,1.0,E,dxyz,0.0,drtn);

                    drtns[1][i][0][nrtns[1][i]]=drtn[0];
                    drtns[1][i][1][nrtns[1][i]]=drtn[1];
                    drtns[1][i][2][nrtns[1][i]]=drtn[2];
                    nrtns[1][i]++;

                    fprintf(fp,"%s  [%6.2lf,%3d]",prn,pod.sol.dops[i],pod.sol.ns[i]);
                    for (j=0;j<6;j++) fprintf(fp,"%15.4lf  [%6.3lf] ",pod.sol.xs[i][j]?pod.sol.xs[i][j]-rs_eci[j]:0.0,sqrt(pod.sol.Qs[i][j]));
                    fprintf(fp,"%15.4lf  %15.4lf ",
                            pod.sol.xs[i][GIC ]?pod.sol.xs[i][GIC ]-relclk-dts[0]*CLIGHT:0.0,
                            pod.sol.xs[i][GICR]?pod.sol.xs[i][GICR]-dts[1]*CLIGHT:0.0);

                    for (j=6;j<GNX;j++) fprintf(fp,"%15.4lf ",pod.sol.xs[i][j]);
                    fprintf(fp,"\n");
                    fflush(fp);
                }
            }
            else {
                fprintf(fp,"\ntime=%s\n",time_str(utc2gpst(pod.sol.tutc),3));
                for (i=0;i<31;i++) {

                    fprintf(fp,"%s  ",pod.rcv[i+1].name);
                    for (j=0;j<3;j++) fprintf(fp,"%15.4lf ",pod.sol.xr[i+1][j]-pod.rcv[i+1].pos[j]);
                    fprintf(fp,"\n");
                    fflush(fp);
                }
            }
        }
    }
    FILE *fprms_forward=fopen("satrms_forward.out","w");

    for (int c=0;c<nrtns[0][0];c++) {
        for (i=0;i<MAXSAT;i++) {
            for (j=0;j<3;j++) fprintf(fprms_forward,"%6.2lf %6.2lf %6.2lf ",drtns[0][i][0][c],drtns[0][i][1][c],drtns[0][i][2][c]);
        }
        fprintf(fprms_forward,"\n");
        fflush(fprms_forward);
    }
    for (i=0;i<MAXSAT;i++) {
        double rmss[3];
        char prn[32];

        satno2id(i+1,prn);

        rmss[0]=rms(drtns[0][i][0],nrtns[0][i]);
        rmss[1]=rms(drtns[0][i][1],nrtns[0][i]);
        rmss[2]=rms(drtns[0][i][2],nrtns[0][i]);
        fprintf(stderr,"%s RMS[forward]=%6.2lf %6.2lf %6.2lf %6.2lf\n",prn, rmss[0],rmss[1],rmss[2],SQRT(SQR(rmss[0])+SQR(rmss[1])+SQR(rmss[2])));
    }
    fprintf(fprms_forward,"\n\n");
    fflush(fprms_forward);

    FILE *fprms_backword=fopen("satrms_backword.out","w");

    for (int c=0;c<nrtns[1][0];c++) {
        for (i=0;i<MAXSAT;i++) {
            for (j=0;j<3;j++) fprintf(fprms_backword,"%6.2lf %6.2lf %6.2lf ",drtns[1][i][0][c],drtns[1][i][1][c],drtns[1][i][2][c]);
        }
        fprintf(fprms_backword,"\n");
        fflush(fprms_backword);
    }
    for (i=0;i<MAXSAT;i++) {
        double rmss[3];
        char prn[32];

        satno2id(i+1,prn);

        rmss[0]=rms(drtns[1][i][0],nrtns[1][i]);
        rmss[1]=rms(drtns[1][i][1],nrtns[1][i]);
        rmss[2]=rms(drtns[1][i][2],nrtns[1][i]);
        fprintf(stderr,"%s RMS[backword]=%6.2lf %6.2lf %6.2lf %6.2lf\n",prn, rmss[0],rmss[1],rmss[2],SQRT(SQR(rmss[0])+SQR(rmss[1])+SQR(rmss[2])));
    }
    fprintf(fprms_backword,"\n\n");
    fflush(fprms_backword);

    fclose(fprms_forward);
    fclose(fprms_backword);

    podfree(&pod);
    fclose(fp);
    log_trace_close();
    return 1;
}