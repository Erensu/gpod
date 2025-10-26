/* ---------------------------------------------------------------------
*
*                              eopspw.cpp
*
*  this file contains routines to read the eop and space weather data
*  from the cssi files on celestrak. it also includes routines to
*  interpolate and convert between the values.
*
*                          companion code for
*             fundamentals of astrodynamics and applications
*                                  2007
*                            by david vallado
*
*     (w) 719-573-2600, email dvallado@agi.com
*     *****************************************************************
*  current :
*            21 mar 08  david vallado
*                           misc fixes
*  changes :
*            14 dec 05  david vallado
*                           misc fixes
*            21 oct 05  david vallado
*                           original version
*       ----------------------------------------------------------------      */
#include "eopspw.h"

spwdata gspwarr[spwsize];
eopdata geoparr[eopsize];
double gjdeopstart;
double gjdspwstart;

double  sgn(double x)
{
    if (x < 0.0)
    {
        return -1.0;
    }
    else
    {
        return 1.0;
    }

}  // end sgn

/* -----------------------------------------------------------------------------
*
*                           function initspw
*
*  this function initializes the space weather data from the cssi files.
*
*  author        : david vallado                  719-573-2600   2 nov 2005
*
*  revisions
*
*  inputs          description                    range / units
*
*  outputs       :
*    spwarr      - array of spw data records
*    jdspwstart  - julian date of the start of the spwarr data
*
*  locals        :
*                -
*
*  coupling      :
*
*  references    :
*
*  -------------------------------------------------------------------------- */
void initspw(const char *file, spwdata spwarr[spwsize], double* jdspwstart)
{
    FILE *infile;
    char longstr[140];
    char str[9], blk[30],buff[1024];
    int numrecsobs=0, numrecspred, i, year, mon, day, tmp;

    // ---- open files
    infile = fopen(file, "r");
    if (infile==NULL) return;

    while (fgets(buff,sizeof(buff),infile)) {
        if (strstr(buff,"NUM_OBSERVED_POINTS")) {
            sscanf(buff,"%s %i",blk,&numrecsobs);
            break;
        }
    }
    if (numrecsobs<=0) return;

    // ---- find epoch date
    fgets( longstr,140,infile);

    // ---- process observed records
    for (i = 0; i < numrecsobs; i++)
    {
        // use d format for integers with leading 0's
        fscanf(infile,"%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4lf %2i %4i %6lf %2i %6lf %6lf %6lf %6lf %6lf \n",
               &spwarr[i].year, &spwarr[i].mon, &spwarr[i].day, &spwarr[i].bsrn, &spwarr[i].nd,
               &spwarr[i].kparr[0],&spwarr[i].kparr[1],&spwarr[i].kparr[2],&spwarr[i].kparr[3],
               &spwarr[i].kparr[4],&spwarr[i].kparr[5],&spwarr[i].kparr[6], &spwarr[i].kparr[7], &spwarr[i].sumkp,
               &spwarr[i].aparr[0],&spwarr[i].aparr[1],&spwarr[i].aparr[2],&spwarr[i].aparr[3],
               &spwarr[i].aparr[4],&spwarr[i].aparr[5],&spwarr[i].aparr[6], &spwarr[i].aparr[7], &spwarr[i].avgap,
               &spwarr[i].cp, &spwarr[i].c9, &spwarr[i].isn, &spwarr[i].adjf10, &spwarr[i].q,
               &spwarr[i].adjctrf81, &spwarr[i].adjlstf81, &spwarr[i].obsf10, &spwarr[i].obsctrf81, &spwarr[i].obslstf81 );
    }

    jday( spwarr[0].year,spwarr[0].mon,spwarr[0].day, 0,0,0.0, jdspwstart);
    jdspwstart = jdspwstart;

    fgets( longstr,140,infile);
    fgets( longstr,140,infile);
    fscanf(infile,"%s %i ", blk,&numrecspred);
    fgets( longstr,140,infile);

    // ---- process predicted records
    for (i = numrecsobs; i < numrecsobs + numrecspred; i++)
    {
        // use d format for integers with leading 0's
        fscanf(infile,"%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4lf %2i %4i %6lf %6lf %6lf %6lf %6lf %6lf \n",
               &spwarr[i].year, &spwarr[i].mon, &spwarr[i].day, &spwarr[i].bsrn, &spwarr[i].nd,
               &spwarr[i].kparr[0],&spwarr[i].kparr[1],&spwarr[i].kparr[2],&spwarr[i].kparr[3],
               &spwarr[i].kparr[4],&spwarr[i].kparr[5],&spwarr[i].kparr[6], &spwarr[i].kparr[7], &spwarr[i].sumkp,
               &spwarr[i].aparr[0],&spwarr[i].aparr[1],&spwarr[i].aparr[2],&spwarr[i].aparr[3],
               &spwarr[i].aparr[4],&spwarr[i].aparr[5],&spwarr[i].aparr[6], &spwarr[i].aparr[7], &spwarr[i].avgap,
               &spwarr[i].cp, &spwarr[i].c9, &spwarr[i].isn, &spwarr[i].adjf10,
               &spwarr[i].adjctrf81, &spwarr[i].adjlstf81, &spwarr[i].obsf10, &spwarr[i].obsctrf81, &spwarr[i].obslstf81 );
        spwarr[i].q = 0;
    }
    fclose( infile );
}

/* -----------------------------------------------------------------------------
*
*                           function initeop
*
*  this function initializes the earth orientation parameter data. the input
*  data files are from celestrak and the eoppn.txt file contains the nutation
*  daily values used for optimizing the speed of operation. note these nutation
*  values do not have corrections applied (dx/dy/ddpsi/ddeps) so a single
*  routine can process past and future data with these corrections. the
*  corrections are not known in the future. the files could be combined, but it
*  may make more sense to keep them separate because the values can be calculated
*  long into the future.
*
*  author        : david vallado                  719-573-2600   2 nov 2005
*
*  revisions
*
*  inputs          description                    range / units
*
*  outputs       :
*    eoparr      - array of eop data records
*    jdeopstart  - julian date of the start of the eoparr data
*
*  locals        :
*                -
*
*  coupling      :
*
*  references    :
*
*  -------------------------------------------------------------------------- */
void initeop
        (
                const char *file,
                eopdata eoparr[eopsize],
                double* jdeopstart
        )
{
    FILE *infile, *infile1;
    char longstr[140];
    char str[9], blk[19];
    double mjd;
    int numrecsobs, numrecspred, year, mon, day;
    long i;

    // ---- open files select compatible files!!
    infile  = fopen(file, "r");

    // ---- read number of data points
    fscanf(infile,"%19c %5i ", blk,&numrecsobs );

    // ---- find epoch date
    fgets( longstr,140,infile);

    // ---- process observed records
    for (i = 0; i < numrecsobs; i++)
    {
        // use d format for integers with leading 0's
        fscanf(infile,"%4i %3d %3d %6i %10lf %10lf %11lf %11lf %10lf %10lf %10lf %10lf %4i \n",
               &eoparr[i].year, &eoparr[i].mon, &eoparr[i].day, &eoparr[i].mjd,
               &eoparr[i].xp, &eoparr[i].yp, &eoparr[i].dut1, &eoparr[i].lod,
               &eoparr[i].ddpsi, &eoparr[i].ddeps, &eoparr[i].dx, &eoparr[i].dy, &eoparr[i].dat );

    }

    jday( eoparr[0].year, eoparr[0].mon, eoparr[0].day, 0,0,0.0, jdeopstart);
    jdeopstart = jdeopstart;

    fgets( longstr,140,infile);
    fgets( longstr,140,infile);
    fscanf(infile,"%20c %5i ", blk,&numrecspred );
    fgets( longstr,140,infile);

    // ---- process predicted records
    for (i = numrecsobs; i < numrecsobs + numrecspred; i++)
    {
        // use d format for integers with leading 0's
        fscanf(infile,"%4i %3d %3d %6i %10lf %10lf %11lf %11lf %10lf %10lf %10lf %10lf %4i \n",
               &eoparr[i].year, &eoparr[i].mon, &eoparr[i].day, &eoparr[i].mjd,
               &eoparr[i].xp, &eoparr[i].yp, &eoparr[i].dut1, &eoparr[i].lod,
               &eoparr[i].ddpsi, &eoparr[i].ddeps, &eoparr[i].dx, &eoparr[i].dy, &eoparr[i].dat );

    }

    fclose( infile );
}   // procedure initeop


/* -----------------------------------------------------------------------------
*
*                           function findeopparam
*
*  this routine finds the eop parameters for a given time. several types of
*  interpolation are available. the cio and iau76 nutation parameters are also
*  read for optimizing the speeed of calculations.
*
*  author        : david vallado                      719-573-2600   12 dec 2005
*
*  inputs          description                          range / units
*    jde         - julian date of epoch (0 hrs utc)
*    mfme        - minutes from midnight epoch
*    interp      - interpolation                        n-none, l-linear, s-spline
*    eoparr      - array of eop data records
*    jdeopstart  - julian date of the start of the eoparr data (set in initeop)
*
*  outputs       :
*    dut1        - delta ut1 (ut1-utc)                  sec
*    dat         - number of leap seconds               sec
*    lod         - excess length of day                 sec
*    xp          - x component of polar motion          rad
*    yp          - y component of polar motion          rad
*    ddpsi       - correction to delta psi (iau-76 theory) rad
*    ddeps       - correction to delta eps (iau-76 theory) rad
*    dx          - correction to x (cio theory)         rad
*    dy          - correction to y (cio theory)         rad
*    x           - x component of cio                   rad
*    y           - y component of cio                   rad
*    s           -                                      rad
*    deltapsi    - nutation longitude angle             rad
*    deltaeps    - obliquity of the ecliptic correction rad
*
*  locals        :
*                -
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2004,
* --------------------------------------------------------------------------- */
void findeopparam
        (
                double  jd,       double mfme,     char interp,
                eopdata eoparr[eopsize],           double jdeopstart,
                double* dut1,     int* dat,
                double* lod,      double* xp,      double* yp,
                double* ddpsi,    double* ddeps,   double* dx,   double* dy,
                double* x,        double* y,       double* s,
                double* deltapsi, double* deltaeps
        )
{
    long i, recnum;
    int  year, mon, day, idx, off1, off2;
    eopdata eoprec, lasteoprec, nexteoprec,  tempeoprec;
    double  fixf,   fixa, mjdt, recnumstart, jdeopstarto;

    // ---- read data for day of interest
    jd          = jd + mfme/1440.0;
    jdeopstarto = floor( jd - jdeopstart);
    recnum      = (int)(jdeopstarto);

    // check for out of bound values
    if ((recnum >= 1) && (recnum <= eopsize))
    {
        eoprec   = eoparr[recnum-1];

        // ---- set non-interpolated values
        *dut1     = eoprec.dut1;
        *dat      = eoprec.dat;
        *lod      = eoprec.lod;
        *xp       = eoprec.xp;
        *yp       = eoprec.yp;
        *ddpsi    = eoprec.ddpsi;
        *ddeps    = eoprec.ddeps;
        *dx       = eoprec.dx;
        *dy       = eoprec.dy;

        // ---- find nutation parameters for use in optimizing speed
        *x        = eoprec.x;
        *y        = eoprec.y;
        *s        = eoprec.s;
        *deltapsi = eoprec.deltapsi;
        *deltaeps = eoprec.deltaeps;

        // ---- do linear interpolation
        if (interp == 'l')
        {
            nexteoprec = eoparr[recnum];
            fixf = mfme / 1440.0;

            *dut1     = eoprec.dut1     + (nexteoprec.dut1     - eoprec.dut1    ) * fixf;
            *dat      = eoprec.dat      + (nexteoprec.dat      - eoprec.dat     ) * fixf;
            *lod      = eoprec.lod      + (nexteoprec.lod      - eoprec.lod     ) * fixf;
            *xp       = eoprec.xp       + (nexteoprec.xp       - eoprec.xp      ) * fixf;
            *yp       = eoprec.yp       + (nexteoprec.yp       - eoprec.yp      ) * fixf;
            *ddpsi    = eoprec.ddpsi    + (nexteoprec.ddpsi    - eoprec.ddpsi   ) * fixf;
            *ddeps    = eoprec.ddeps    + (nexteoprec.ddeps    - eoprec.ddeps   ) * fixf;
            *dx       = eoprec.dx       + (nexteoprec.dx       - eoprec.dx      ) * fixf;
            *dy       = eoprec.dy       + (nexteoprec.dy       - eoprec.dy      ) * fixf;
            *x        = eoprec.x        + (nexteoprec.x        - eoprec.x       ) * fixf;
            *y        = eoprec.y        + (nexteoprec.y        - eoprec.y       ) * fixf;
            *s        = eoprec.s        + (nexteoprec.s        - eoprec.s       ) * fixf;
            *deltapsi = eoprec.deltapsi + (nexteoprec.deltapsi - eoprec.deltapsi) * fixf;
            *deltaeps = eoprec.deltaeps + (nexteoprec.deltaeps - eoprec.deltaeps) * fixf;
        }

        // ---- do spline interpolations
        if (interp == 's')
        {
            off1 = 10;   // every 5 days data...
            off2 = 5;
            *dut1 = cubicinterp ( eoparr[recnum-off1].dut1, eoparr[recnum-off2].dut1, eoparr[recnum].dut1, eoparr[recnum+off2].dut1,
                                  eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                  mfme );
            *dat = cubicinterp ( eoparr[recnum-off1].dat, eoparr[recnum-off2].dat, eoparr[recnum].dat, eoparr[recnum+off2].dat,
                                 eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                 mfme );
            *lod = cubicinterp ( eoparr[recnum-off1].lod, eoparr[recnum-off2].lod, eoparr[recnum].lod, eoparr[recnum+off2].lod,
                                 eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                 mfme );
            *xp = cubicinterp ( eoparr[recnum-off1].xp, eoparr[recnum-off2].xp, eoparr[recnum].xp, eoparr[recnum+off2].xp,
                                eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                mfme );
            *yp = cubicinterp ( eoparr[recnum-off1].yp, eoparr[recnum-off2].yp, eoparr[recnum].yp, eoparr[recnum+off2].yp,
                                eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                mfme );
            *ddpsi = cubicinterp ( eoparr[recnum-off1].ddpsi, eoparr[recnum-off2].ddpsi, eoparr[recnum].ddpsi, eoparr[recnum+off2].ddpsi,
                                   eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                   mfme );
            *ddeps = cubicinterp ( eoparr[recnum-off1].ddeps, eoparr[recnum-off2].ddeps, eoparr[recnum].ddeps, eoparr[recnum+off2].ddeps,
                                   eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                   mfme );
            *dx = cubicinterp ( eoparr[recnum-off1].dx, eoparr[recnum-off2].dx, eoparr[recnum].dx, eoparr[recnum+off2].dx,
                                eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                mfme );
            *dy = cubicinterp ( eoparr[recnum-off1].dy, eoparr[recnum-off2].dy, eoparr[recnum].dy, eoparr[recnum+off2].dy,
                                eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                mfme );
            *x = cubicinterp ( eoparr[recnum-off1].x, eoparr[recnum-off2].x, eoparr[recnum].x, eoparr[recnum+off2].x,
                               eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                               mfme );
            *y = cubicinterp ( eoparr[recnum-off1].y, eoparr[recnum-off2].y, eoparr[recnum].y, eoparr[recnum+off2].y,
                               eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                               mfme );
            *s = cubicinterp ( eoparr[recnum-off1].s, eoparr[recnum-off2].s, eoparr[recnum].s, eoparr[recnum+off2].s,
                               eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                               mfme );
            *deltapsi = cubicinterp ( eoparr[recnum-off1].deltapsi, eoparr[recnum-off2].deltapsi, eoparr[recnum].deltapsi, eoparr[recnum+off2].deltapsi,
                                      eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                      mfme );
            *deltaeps = cubicinterp ( eoparr[recnum-off1].deltaeps, eoparr[recnum-off2].deltaeps, eoparr[recnum].deltaeps, eoparr[recnum+off2].deltaeps,
                                      eoparr[recnum-off1].mjd, eoparr[recnum-off2].mjd, eoparr[recnum].mjd, eoparr[recnum+off2].mjd,
                                      mfme );
        }
    }
        // set default values
    else
    {
        *dut1     = 0.0;
        *dat      = 33;
        *lod      = 0.0;
        *xp       = 0.0;
        *yp       = 0.0;
        *ddpsi    = 0.0;
        *ddeps    = 0.0;
        *dx       = 0.0;
        *dy       = 0.0;

        // ---- find nutation parameters for use in optimizing speed
// these could be set here, or in the program calling this...
//           x        = eoprec.x;
//           y        = eoprec.y;
//           s        = eoprec.s;
//           deltapsi = eoprec.deltapsi;
//           deltaeps = eoprec.deltaeps;
    }
}  // procedure findeopparam

/* -----------------------------------------------------------------------------
*
*                           function findatmosparam
*
*  this routine finds the atmospheric parameters for a given time.
*    ap/kp 3 hourly data is valid at 0000, 0300 hrs, etc
*    apavg and kpsum are valid at 1200 hrs
*    f107 and f107bar values are valid at 1700/2000 hrs depending on the date
*    ap arrays go 0-7, but msisarr goes 1-8 to match msis code and fortran source
*
*  author        : david vallado                      719-573-2600    2 dec 2005
*
*  inputs          description                          range / units
*    jd          - julian date of epoch (0 hrs utc)     days from 4713 bc
*    mfme        - minutes from midnight epoch          mins
*    interp      - interpolation           n-none, a-ap only, f-f10.7 only, b-both
*    fluxtype    - flux type               a-adjusted, o-observed
*    f81type     - flux 81-day avg type    l-last, c-centered
*    inputtype   - input type              a-actual, u - user   c - constant
*    spwarr      - array of space weather data
*    jdspwstart  - julian date of the start of the spwarr data (set in initspw)
*
*  outputs       :
*    f107        - f10.7 value (current day)
*    f107bar     - f10.7 81-day avg value
*    ap          - planetary effect array
*    avgap       - daily average ap planetary value
*    aparr       - last 8 values of 3-hourly ap
*    kp          - planetary effect array
*    sumkp       - daily kp sum of planetary values
*    kparr       - last 8 values of 3-hourly kp
*
*  locals        :
*    fluxtime    - minutes from midnight where f107 is valid (1020 or 1200)
*
*  coupling      :
*    none        -
*
*  references    :
*
*  -------------------------------------------------------------------------- */

void findatmosparam
        (
                double jd, double mfme, char interp, char fluxtype, char f81type, char inputtype,
                spwdata spwarr[spwsize], double jdspwstart,
                double* f107, double* f107bar,
                double* ap, double* avgap, double aparr[8],
                double* kp, double* sumkp, double kparr[8]
        )
{
    int     i, recnum, year, mon, day, idx, j;
    double  tf107,  tf107bar, tavgap;
    char    ftype,  fctrtype;
    spwdata spwrec, lastspwrec, nextspwrec, tempspwrec;
    double  fixf,   fixa, mjdt, recnumstart, fluxtime, jdspwstarto;

    // --------------------  implementation   ----------------------
    // ---- set flux time based on when measurments were taken
    // ---- before may 31, 1991, use 1700 hrs (1020 minutes)
    if (jd > 2448407.5)
        fluxtime = 1200.0;
    else
        fluxtime = 1020.0;

    // ---- process actual data
    if (inputtype == 'a')
    {
        // ---- read data for day of interest
        jd          = jd + mfme/1440.0;
        jdspwstarto = floor(jd - jdspwstart);
        recnum      = (int)(jdspwstarto);

        if (recnum < 1)
        {
            printf("%14.5lf before %14.5lf date in file, hit ctrl-c \n", jd, jdspwstart);
            scanf( "%lf", &jd );
        }
        // ---- set non-interpolated values
        spwrec     = spwarr[recnum];
        lastspwrec = spwarr[recnum-1];
        nextspwrec = spwarr[recnum+1];

        if (fluxtype == 'a')
        {
            *f107 = spwrec.adjf10;
            if (f81type == 'l')
                *f107bar = spwrec.adjlstf81;
            else
                *f107bar = spwrec.adjctrf81;
        }
        else
        {
            *f107 = spwrec.obsf10;
            if (f81type == 'l')
                *f107bar = spwrec.obslstf81;
            else
                *f107bar = spwrec.obsctrf81;
        }
        *avgap = spwrec.avgap;
        *sumkp = spwrec.sumkp;

        // ---- get last ap/kp array value from the current time value
        idx = floor(mfme/180.0); // values change at 0, 3, 6, ... hrs
        idx = (int)(idx);
        if (idx < 0) idx = 0;
        if (idx > 7) idx = 7;

        j = idx;
        for (i = 1; i <= 8; i++)
        {
            if (j >= 0)
            {
                aparr[8-i] = spwrec.aparr[j];
                kparr[8-i] = spwrec.kparr[j];
            }
            else
            {
                aparr[8-i] = lastspwrec.aparr[8+j];
                kparr[8-i] = lastspwrec.kparr[8+j];
            }
            j = j - 1;
        }
        *ap = spwrec.aparr[ idx ];
        *kp = spwrec.kparr[ idx ]*0.1;

        // ------------------------ do interpolation ------------------------
        if (interp != 'n')
        {
            if ((interp == 'f') | (interp == 'b'))
            {
                if (mfme > fluxtime-720.0) // go 12 hrs before...
                {
                    if (mfme > fluxtime)
                        tempspwrec = nextspwrec;
                    else
                        tempspwrec = lastspwrec;
                    fixf = (fluxtime - mfme) / 1440.0;
                }
                else
                {
                    tempspwrec = lastspwrec;
                    fixf = (mfme + (1440 - fluxtime)) / 1440.0;
                }
                if (fluxtype == 'a') // adjusted or observed values
                {
                    tf107 = tempspwrec.adjf10;
                    if (f81type == 'l')
                        tf107bar = tempspwrec.adjlstf81;
                    else
                        tf107bar = tempspwrec.adjctrf81;
                }
                else
                {
                    tf107 = tempspwrec.obsf10;
                    if (f81type == 'l')
                        tf107bar = tempspwrec.obslstf81;
                    else
                        tf107bar = tempspwrec.obsctrf81;
                }
                // ---- perform simple linear interpolation
                if (mfme <= fluxtime)
                {
                    if (mfme > fluxtime-720.0)
                    {
                        *f107    = *f107    + (tf107 - *f107) * fixf;
                        *f107bar = *f107bar + (tf107bar - *f107bar) * fixf;
                    }
                    else
                    {
                        *f107    = tf107    + (*f107 - tf107) * fixf;
                        *f107bar = tf107bar + (*f107bar - tf107bar) * fixf;
                    }
                }
                else
                {
                    *f107    = *f107    + (tf107 - *f107) * fixf;
                    *f107bar = *f107bar + (tf107bar - *f107bar) * fixf;
                }
                // ---- perform cubic splining
//                   f107 = cubicinterp ( lastspwrec.f107, f107, nextspwrec.f107, nextspwrec+1.f107,
//                                        1??, 2, 3, 4, mfme+?? );
//                   f107bar = cubicinterp ( lastspwrec.f107bar, f107bar, nextspwrec.f107bar, nextspwrec+1.f107bar,
//                                        1??, 2, 3, 4, mfme+?? );
//                   f107bar = cubicinterp ( bkp[idx-2], bkp[idx-1], bkp[idx], bkp[idx+1],
//                                        bap[idx-2], bap[idx-1], bap[idx], bap[idx+1],
//                                        apin );
//                   f107bar = cubicinterp ( bkp[idx-2], bkp[idx-1], bkp[idx], bkp[idx+1],
//                                        bap[idx-2], bap[idx-1], bap[idx], bap[idx+1],
//                                        apin );
            }

            // avgap is effective each 1200, sumkp is effective at 2400 hrs
            // sumkp is probably not right - ie, centered on 0hrs instead of 24. really need next day??
            if ((interp == 'a') | (interp == 'b'))
            {
                fixa = (720 - mfme) / 1440.0;
                if (mfme > 720)
                {
                    *avgap = *avgap - (nextspwrec.avgap - *avgap) * fixa;
                    *sumkp = nextspwrec.sumkp - (nextspwrec.sumkp - *sumkp) * (mfme/1440.0);
                }
                else
                {
                    *avgap = *avgap - (*avgap - lastspwrec.avgap) * fixa;
                    *sumkp = *sumkp - (*sumkp - lastspwrec.sumkp) * (1440.0-mfme)/1440.0;
                }

                // this fraction is the same for the remainder of calculations
                fixa = (fmod(mfme,180)) / 180.0;
                if (idx+1 < 8 )
                {
                    *ap = spwrec.aparr[idx] + (spwrec.aparr[idx+1]-spwrec.aparr[idx]) * fixa;
                    *kp = (spwrec.kparr[idx] + (spwrec.aparr[idx+1]-spwrec.aparr[idx]) * fixa)*0.1;
                }
                else
                {
                    *ap = spwrec.aparr[idx] + (nextspwrec.aparr[0]-spwrec.aparr[idx]) * fixa;
                    *kp = (spwrec.kparr[idx] + (nextspwrec.aparr[0]-spwrec.aparr[idx]) * fixa)*0.1;
                }

                // step down from idx through the 8 points
                j = idx;
                fixa = fmod(mfme, 90.0) / 180.0;
                for (i = 1; i <= 8; i++)
                {
                    if (j >= 0) // j = 0 .. 6
                    {
                        if (j+1 < 8)
                        {
                            aparr[8-i] = spwrec.aparr[j] + (spwrec.aparr[j+1] - spwrec.aparr[j]) * fixa;
                            kparr[8-i] = spwrec.kparr[j] + (spwrec.kparr[j+1] - spwrec.kparr[j]) * fixa;
                        }
                        else // j = 0
                        {
                            aparr[8-i] = spwrec.aparr[j] + (nextspwrec.aparr[0] - spwrec.aparr[j]) * fixa;
                            kparr[8-i] = spwrec.kparr[j] + (nextspwrec.kparr[0] - spwrec.kparr[j]) * fixa;
                        }
                    }
                    else
                    {    // j = -2 .. -7
                        if (j < -1)
                        {
                            aparr[8-i] = lastspwrec.aparr[8+j] + (lastspwrec.aparr[9+j] - lastspwrec.aparr[8+j]) * fixa;
                            kparr[8-i] = lastspwrec.kparr[8+j] + (lastspwrec.kparr[9+j] - lastspwrec.kparr[8+j]) * fixa;
                        }
                        else // j = -1
                        {
                            aparr[8-i] = lastspwrec.aparr[7] + (spwrec.aparr[0] - lastspwrec.aparr[7]) * fixa;
                            kparr[8-i] = lastspwrec.kparr[7] + (spwrec.kparr[0] - lastspwrec.kparr[7]) * fixa;
                        }
                    }

                    j = j - 1;
                }  // for
            }
        } // if interp != n
    }  // if interp = a
    else
        // ---- user input data
    if (inputtype == 'u')
    {
        // this is for data that may be simulated, or otherwise different from
        // the current noaa data

        // there could also be the interpolation stuff from above

    }
    else
        // ---- constant data
    if (inputtype == 'c')
    {
        // this data is the same all the time
        // leave the same as when it enters
    }

}


/* -----------------------------------------------------------------------------
*
*                           function kp2ap
*
*  this function converts kp to ap using cubic splines. notice the arrays go
*  beyond the range of values to permit endpoint evaluations without additional
*  logic. the arrays have an extra element so they will start at 1.
*
*  author        : david vallado                  719-573-2600   7 aug  2005
*
*  revisions
*
*  inputs          description                    range / units
*    kpin        - kp
*
*  outputs       :
*    kp2ap       - ap
*
*  locals        :
*    idx         - index of function value above the input value so the input
*                  value is between the 2nd and 3rd point
*
*  coupling      :
*    cubicspl    - perform the splining operation given 4 points
*
*  references    :
*    vallado       2004, 899-901
* --------------------------------------------------------------------------- */

double kp2ap
        (
                double kpin
        )
{
    static double bapt[32] =
            { 0, -0.00001, -0.001,
              0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32,
              39, 48, 56, 67, 80, 94, 111, 132, 154, 179, 207, 236, 300, 400, 900
            };
    static double bkpt[32] =
            { 0, -0.66666667,  -0.33333,
              0, 0.33333, 0.66667, 1, 1.33333, 1.66667, 2, 2.33333,
              2.66667, 3, 3.33333, 3.66667, 4, 4.33333,
              4.66667, 5, 5.33333, 5.66667, 6.0, 6.33333, 6.66667, 7,
              7.33333, 7.66667, 8, 8.33333, 8.666667, 9, 9.33333
            };

    double bap[32], bkp[32];
    int idx;

    memcpy( bap,  bapt,  32*sizeof(double) );
    memcpy( bkp,  bkpt,  32*sizeof(double) );

    idx = 1;
    while ((idx < 33) && (kpin > bkp[idx]))
    {
        idx = idx + 1;
    }

    if (idx > 2)
    {
        return cubicinterp ( bap[idx-2], bap[idx-1], bap[idx], bap[idx+1],
                             bkp[idx-2], bkp[idx-1], bkp[idx], bkp[idx+1],
                             kpin );
    } // if idx > 3
    else
        return 0.0;
}

/* -----------------------------------------------------------------------------
*
*                           function ap2kp
*
*  this function converts ap to kp using cubic splines. notice the values go
*  beyond the range of values to permit endpoint evaluations without additional
*  logic. the arrays have an extra element so they will start at 1.
*
*  author        : david vallado                  719-573-2600   7 aug  2005
*
*  revisions
*
*  inputs          description                    range / units
*    kpin        - kp
*
*  outputs       :
*    ap2kp       - ap
*
*  locals        :
*    idx         - index of function value above the input value so the input
*                  value is between the 2nd and 3rd point
*
*  coupling      :
*
*  references    :
*    vallado       2004, 899-901
* --------------------------------------------------------------------------- */

double ap2kp
        (
                double apin
        )
{
    static double bapt[32] =
            { 0, -0.00001, -0.001 ,
              0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32,
              39, 48, 56, 67, 80, 94, 111, 132, 154, 179, 207, 236, 300, 400, 900
            };
    static double bkpt[32] =
            { 0, -0.66666667,  -0.33333  ,
              0, 0.33333, 0.66667, 1, 1.33333, 1.66667, 2, 2.33333,
              2.66667, 3, 3.33333, 3.66667, 4, 4.33333,
              4.66667, 5, 5.33333, 5.66667, 6.0, 6.33333, 6.66667, 7,
              7.33333, 7.66667, 8, 8.33333, 8.666667, 9, 9.33333
            };

    double kp;
    double bap[32], bkp[32];
    int  idx;

    memcpy( bap,  bapt,  32*sizeof(double) );
    memcpy( bkp,  bkpt,  32*sizeof(double) );

    idx = 1;
    while ((idx < 33) && (apin > bap[idx]))
    {
        idx = idx + 1;
    }

    if (idx > 2)
    {
        return cubicinterp ( bkp[idx-2], bkp[idx-1], bkp[idx], bkp[idx+1],
                             bap[idx-2], bap[idx-1], bap[idx], bap[idx+1],
                             apin );
    } // if idxs > 3
    else
        return 0.0;
}

/* -----------------------------------------------------------------------------
*
*                           procedure jday
*
*  this procedure finds the julian date given the year, month, day, and time.
*    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
*
*  algorithm     : calculate the answer in one step for efficiency
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    year        - year                           1900 .. 2100
*    mon         - month                          1 .. 12
*    day         - day                            1 .. 28,29,30,31
*    hr          - universal time hour            0 .. 23
*    min         - universal time min             0 .. 59
*    sec         - universal time sec             0.0 .. 59.999
*
*  outputs       :
*    jd          - julian date                    days from 4713 bc
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
*
*  references    :
*    vallado       2007, 189, alg 14, ex 3-14
*
* --------------------------------------------------------------------------- */
void jday(int year, int mon, int day, int hr, int minute, double sec, double* jd)
{
    *jd = 367.0 * year -
          floor((7 * (year + floor((mon + 9) / 12.0))) * 0.25) +
          floor( 275 * mon / 9.0 ) +
          day + 1721013.5 +
          ((sec / 60.0 + minute) / 60.0 + hr) / 24.0;  // ut in days
    // - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
}

/* -----------------------------------------------------------------------------
*
*                           function cubicinterp
*
*  this function performs a cubic spline. four points are needed.
*
*  author        : david vallado                  719-573-2600   1 dec  2005
*
*  revisions
*
*  inputs          description                    range / units
*    valuein     - kp
*
*  outputs       :
*    out         - ap
*
*  locals        :
*                -
*
*  coupling      :
*    cubicspl
*
*  references    :
*    vallado       2007, 556
* --------------------------------------------------------------------------- */
double  cubicinterp(double p1a, double p1b, double p1c, double p1d, double p2a,
                    double p2b, double p2c, double p2d, double valuein)
{
    double kc0, kc1, kc2, kc3, ac0, ac1, ac2, ac3,
            r1r, r1i, r2r, r2i, r3r, r3i, value;

    // -------- assign function points ---------
    cubicspl(p1a, p1b, p1c, p1d, &ac0, &ac1, &ac2, &ac3);
    cubicspl(p2a, p2b, p2c, p2d, &kc0, &kc1, &kc2, &kc3);

    // recover the original function values
    // use the normalized time first, but at an arbitrary interval
    cubic(kc3, kc2, kc1, kc0-valuein, 'R', &r1r, &r1i, &r2r, &r2i, &r3r, &r3i);

    if ((r1r >= -0.000001) && (r1r <= 1.001))
    {
        value = r1r;
    }
    else
    {
        if ((r2r >= -0.000001) && (r2r <= 1.001))
        {
            value = r2r;
        }
        else
        {
            if ((r3r >= -0.000001) && (r3r <= 1.001))
            {
                value = r3r;
            }
            else
            {
                value = 0.0;
                printf("error in cubicinterp root %17.14f %11.7f %11.7f %11.7f \n",
                       valuein,r1r,r2r,r3r);
            }
        }
    }
    return (ac3*pow(value,3)  + ac2*value*value  + ac1*value + ac0);

} // cubicinterp

/* -----------------------------------------------------------------------------
*
*                           function cubicspl
*
*  this function performs cubic splining of an input zero crossing
*  function in order to find function values.
*
*  author        : david vallado                  719-573-2600     7 aug 2005
*
*  revisions
*                -
*  inputs          description                    range / units
*    p0,p1,p2,p3 - function values used for splining
*    t0,t1,t2,t3 - time values used for splining
*
*  outputs       :
*    acu0..acu3  - splined polynomial coefficients. acu3 t^3, etc
*
*  locals        : none
*
*  coupling      :
*    none
*
*  references    :
*    vallado       2007, 556
* --------------------------------------------------------------------------- */
void cubicspl(double p1, double p2, double p3, double p4, double* acu0,
              double* acu1, double* acu2, double* acu3)
{
    *acu0 = p2;
    *acu1 = -p1/3.0 - 0.5*p2 + p3 -p4/6.0;
    *acu2 = 0.5*p1 - p2 + 0.5*p3;
    *acu3 = -p1/6.0 + 0.5*p2 - 0.5*p3 + p4/6.0;
}

/* -----------------------------------------------------------------------------
*
*                           function cubic
*
*  this function solves for the three roots of a cubic equation.  there are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  the general form is y = a3x3 + b2x2 + c1x + d0.  note
*    that r1i will always be zero since there is always at least one real root.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  revisions
*    vallado     - convert to matlab              719-573-2600   18 dec 2002
*
*  inputs          description                    range / units
*    a3          - coefficient of x cubed term
*    b2          - coefficient of x squared term
*    c1          - coefficient of x term
*    d0          - constant
*    opt         - option for output              I all roots including imaginary
*                                                 R only real roots
*                                                 U only unique real roots (no repeated)
*
*  outputs       :
*    r1r         - real portion of root 1
*    r1i         - imaginary portion of root 1
*    r2r         - real portion of root 2
*    r2i         - imaginary portion of root 2
*    r3r         - real portion of root 3
*    r3i         - imaginary portion of root 3
*
*  locals        :
*    temp1       - temporary value
*    temp2       - temporary value
*    p           - coefficient of x squared term where x cubed term is 1.0
*    q           - coefficient of x term where x cubed term is 1.0
*    r           - coefficient of constant term where x cubed term is 1.0
*    delta       - discriminator for use with cardans formula
*    e0          - angle holder for trigonometric solution
*    phi         - angle used in trigonometric solution
*    cosphi      - cosine of phi
*    sinphi      - sine of phi
*
*  coupling      :
*    quadric     - quadratic roots
*
*  references    :
*    vallado       2007, 975
*
* --------------------------------------------------------------------------- */
void cubic(double a3, double b2, double c1, double d0, char opt, double* r1r,
           double* r1i, double* r2r, double* r2i, double* r3r, double* r3i)
{
    const double rad      = 57.29577951308230;
    const double onethird = 1.0/3.0;
    const double small    = 0.00000001;
    double temp1, temp2, p, q, r, delta, e0, cosphi, sinphi, phi;
    // ------------------------  implementation   --------------------------
    *r1r  = 0.0;
    *r1i  = 0.0;
    *r2r  = 0.0;
    *r2i  = 0.0;
    *r3r  = 0.0;
    *r3i  = 0.0;

    if (fabs(a3) > small)
    {
        // ------------- force coefficients into std form -------------------
        p= b2/a3;
        q= c1/a3;
        r= d0/a3;

        a3= onethird*( 3.0 *q - p*p );
        b2= (1.0 /27.0 )*( 2.0 *p*p*p - 9.0 *p*q + 27.0 *r );

        delta= (a3*a3*a3/27.0 ) + (b2*b2*0.25 );

        // -------------------- use cardans formula ------------------------
        if ( delta > small )
        {
            temp1= (-b2*0.5 )+sqrt(delta);
            temp2= (-b2*0.5 )-sqrt(delta);
            temp1= sgn(temp1)*pow( fabs(temp1),onethird );
            temp2= sgn(temp2)*pow( fabs(temp2),onethird );
            *r1r= temp1 + temp2 - p*onethird;

            if (opt=='I')
            {
                *r2r= -0.5 *(temp1 + temp2) - p*onethird;
                *r2i= -0.5 *sqrt( 3.0  )*(temp1 - temp2);
                *r3r= -0.5 *(temp1 + temp2) - p*onethird;
                *r3i= -*r2i;
            }
            else
            {
                *r2r= 99999.9;
                *r3r= 99999.9;
            }
        }
        else
        {
            // -------------------- evaluate zero point ------------------------
            if ( fabs( delta ) < small  )
            {
                *r1r= -2.0*sgn(b2)*pow(fabs(b2*0.5),onethird) - p*onethird;
                *r2r=      sgn(b2)*pow(fabs(b2*0.5),onethird) - p*onethird;
                if (opt=='U')
                    *r3r= 99999.9;
                else
                    *r3r= *r2r;
            }
            else
            {
                // --------------- use trigonometric identities ----------------
                e0     = 2.0 *sqrt(-a3*onethird);
                cosphi = (-b2/(2.0 *sqrt(-a3*a3*a3/27.0 )) );
                sinphi = sqrt( 1.0 -cosphi*cosphi );
                phi    = atan2( sinphi,cosphi );
                if (phi < 0.0)
                    phi = phi + 2.0*M_PI;
                *r1r= e0*cos( phi*onethird ) - p*onethird;
                *r2r= e0*cos( phi*onethird + 120.0 /rad ) - p*onethird;
                *r3r= e0*cos( phi*onethird + 240.0 /rad ) - p*onethird;
            } // if fabs(delta)...
        }  // if delta > small
    }  // if fabs > small
    else
    {
        quadric( b2, c1, d0, opt, r1r, r1i, r2r, r2i );
        *r3r  = 99999.9;
        *r3i  = 99999.9;
    }
}

/* -----------------------------------------------------------------------------
*
*                           function quadric
*
*  this function solves for the two roots of a quadric equation.  there are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  the general form is y = ax2 + bx + c.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  revisions
*    vallado     - convert to matlab              719-573-2600    3 dec 2002
*
*  inputs          description                    range / units
*    a           - coefficient of x squared term
*    b           - coefficient of x term
*    c           - constant
*    opt         - option for output              I all roots including imaginary
*                                                 R only real roots
*                                                 U only unique real roots (no repeated)
*
*  outputs       :
*    r1r         - real portion of root 1
*    r1i         - imaginary portion of root 1
*    r2r         - real portion of root 2
*    r2i         - imaginary portion of root 2
*
*  locals        :
*    discrim     - discriminate b2 - 4ac
*
*  coupling      :
*    none.
*
*  references    :
*    vallado       2007, 974
*
* ----------------------------------------------------------------------------*/
void quadric(double a, double b, double c, char opt, double* r1r, double* r1i,
             double* r2r, double* r2i)
{
    const double small    = 0.0000001;
    double delta, discrim;
    // --------------------  implementation   ----------------------
    *r1r = 0.0;
    *r1i = 0.0;
    *r2r = 0.0;
    *r2i = 0.0;

    discrim = b*b - 4.0 *a*c;

    // ---------------------  real roots  --------------------------
    if ( fabs(discrim) < small  )
    {
        *r1r = -b / ( 2.0 *a );
        *r2r = *r1r;
        if (opt=='U')
            *r2r = 99999.9;
    }
    else
    {
        if (fabs(a) < small)
            *r1r = -c/b;
        else
        {
            if ( discrim > 0.0  )
            {
                *r1r = ( -b + sqrt(discrim) ) / ( 2.0 *a );
                *r2r = ( -b - sqrt(discrim) ) / ( 2.0 *a );
            }
            else
            {
                // ------------------ complex roots --------------------
                if (opt=='I')
                {
                    *r1r = -b / ( 2.0 *a );
                    *r2r = *r1r;
                    *r1i =  sqrt(-discrim) / ( 2.0 *a );
                    *r2i = -sqrt(-discrim) / ( 2.0 *a );
                }
                else
                {
                    *r1r = 99999.9;
                    *r2r = 99999.9;
                }
            }
        }
    }
}

/* -----------------------------------------------------------------------------
*
*                           procedure finddays
*
*  this procedure finds the fractional days through a year given the year,
*    month, day, hour, minute and second.
*
*  algorithm     : set up array for the number of days per month
*                  find leap year - use 1900 because 2000 is a leap year
*                  check for a leap year
*                  loop to find the elapsed days in the year
*
*  author        : david vallado                  719-573-2600
*
*  inputs          description                    range / units
*    year        - year                           1900 .. 2100
*    mon         - month                          1 .. 12
*    day         - day                            1 .. 28,29,30,31
*    hr          - hour                           0 .. 23
*    min         - minute                         0 .. 59
*    sec         - second                         0.0 .. 59.999
*
*  outputs       :
*    days        - day of year plus fraction of a
*                    day                          days
*
*  locals        :
*    lmonth      - length of months of year
*    i           - index
*
*  coupling      :
*    none.
*
*  references    :
*    vallado       2007, 207, ex 3-12
*
* --------------------------------------------------------------------------- */
void finddays(int year, int month, int day, int hr, int minute, double sec, double* days)
{
    int lmonth[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int i;

    if (((year - 1900) % 4) == 0)
        lmonth[1] = 29;

    i    = 1;
    *days = 0.0;
    while ((i < month) && (i < 12))
    {
        *days = *days + lmonth[i-1];
        i = i + 1;
    }
    *days = *days + day + hr / 24.0 + minute / 1440.0 + sec / 86400.0;
}

