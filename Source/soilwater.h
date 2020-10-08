/*****************************************************************************
**
**  FILE:    soilwater.h
**
**  AUTHOR:  Melannie Hartman  9/7/93 - 8/21/96
**
*****************************************************************************/

#include <stdio.h>
#include "swconst.h"

#define FILEPATH_MAX 60   /* Max # of chars in a file path-name */

#define INPTSTRLEN 120    /* Max length of input file line */
#define OTPTSTRLEN 80     /* Max length of output file string */

typedef struct {
  int verbose;
  int debug;
} FLAG_S, *FLAG_SPT;

typedef struct
{
  char wbname[FILEPATH_MAX];  /* Water Balance File */
  FILE *fp_outf;
  FILE *fp_wb;
  FILE *fp_soiln;
  FILE *fp_co2;
  FILE *fp_swc;
  FILE *fp_wfps;
  FILE *fp_soiltavg;
  FILE *fp_soiltmin;
  FILE *fp_soiltmax;
  FILE *fp_stempdx;
  FILE *fp_bio;
  FILE *fp_wflux;
  FILE *fp_resp;
  FILE *fp_yearsum;
  FILE *fp_livec;
  FILE *fp_deadc;
  FILE *fp_soilc;
  FILE *fp_sysc;
  FILE *fp_tgmonth;
  FILE *fp_dN2lyr;
  FILE *fp_dN2Olyr;
  FILE *fp_dels;
  FILE *fp_dcsip;
  FILE *fp_harv;
  FILE *fp_cflows;
  FILE *fp_yrcflows;
  FILE *fp_daily;
  FILE *fp_nflux;
  FILE *fp_summary;
  FILE *fp_methane;
  FILE *fp_psyn;
/*!!  FILE *fp_snow; */
  int  write_wb;
  int  write_co2;
  int  write_soiln;
  int  write_swc;
  int  write_wfps;
  int  write_soiltavg;
  int  write_soiltmin;
  int  write_soiltmax;
  int  write_stempdx;
  int  write_bio;
  int  write_wflux;
  int  write_resp;
  int  write_yearsum;
  int  write_livec;
  int  write_deadc;
  int  write_soilc;
  int  write_sysc;
  int  write_tgmonth;
  int  write_dN2lyr;
  int  write_dN2Olyr;
  int  write_dels;
  int  write_dcsip;
  int  write_harvest;
  int  write_cflows;
  int  write_yrcflows;
  int  write_daily;
  int  write_nflux;
  int  write_summary;
  int  write_methane;
  int  write_psyn;
} FILES_S, *FILES_SPT;

typedef struct
{
  float  width[MAXLYR];
  float  depth[MAXLYR];
  float  dpthmn[MAXLYR];
  float  dpthmx[MAXLYR];
  float  bulkd[MAXLYR];
  float  fieldc[MAXLYR];
  float  wiltpt[MAXLYR];
  float  ecoeff[MAXLYR];
  float  tcoeff[MAXLYR];
  float  sandfrac[MAXLYR];
  float  clayfrac[MAXLYR];
  float  orgfrac[MAXLYR];
  float  swclimit[MAXLYR];
  float  satcond[MAXLYR];
  float  pH[MAXLYR];
  int    numlyrs;
  float  swcfc[MAXLYR];
  float  swcwp[MAXLYR];
  double swc[MAXLYR];
  double swcmin[MAXLYR];
  float  minpot[MAXLYR];
  float  wfps[MAXLYR];
  float  sumecoeff;
/*  float  sumtcoeff[NTDEPTHS];
  int    ntlyrs;
  int    lyrmin[NTDEPTHS];
  int    lyrmax[NTDEPTHS]; */
  int    nelyrs;
  float  thetas[MAXLYR];
  float  thetas_bd[MAXLYR];
  float  psis[MAXLYR];
  float  b[MAXLYR];
  float  binverse[MAXLYR];
  int    ubnd[CENTMAXLYR];
  int    lbnd[CENTMAXLYR];
} LAYERPAR_S, *LAYERPAR_SPT;

typedef struct
{
  int    usexdrvrs;
  float  sublimscale;
  float  reflec;
  float  albedo;
  float  fswcinit;
  float  dmpflux;
  float  hours_rain;
  float  hpotdeep;
  float  ksatdeep;
  float  rlatitude;
  float  cldcov[NMONTH+1];
  int    texture;
  float  tbotmn;
  float  tbotmx;
  float  dmp;
  float  timlag;
  float  Ncoeff;
  int    drainlag;
  int    jdayStart;
  int    jdayEnd;
  float  NO3_N2O_x;
  float  NO3_N2O_y;
  float  NO3_N2O_step;
  float  NO3_N2O_slope;
  float  MaxNitAmt;
  int    SnowFlag;
  float  netmn_to_no3;
  float  wfpsdnitadj;
  float  n2n2oadj;
  float  th2ocoef1;
  float  th2ocoef2;
  float  flood_N2toN2O;
  float  CO2_to_CH4;
  float  C6H12O6_to_CH4;
  float  frac_to_exudates;
  float  Aeh;
  float  Deh;
  float  Beh_flood;
  float  Beh_drain;
  float  zero_root_frac;
  float  elevation;
  float  slope;
  float  aspect;
  float  ehoriz;
  float  whoriz;
} SITEPAR_S, *SITEPAR_SPT;

typedef struct
{
  float stmtemp[MAXSTLYR];
  float soiltavg[MAXLYR];
  float soiltmin[MAXLYR];
  float soiltmax[MAXLYR];
} SOIL_S, *SOIL_SPT;

float c_shwave(int month, float rlatitude, int jday);

/*void calcdefac(int *texture, float *stemp, float *tfunc, float *bgwfunc,
               float *agdefac, float *bgdefac, float *avgwfps, float teff[],
               float *rprpet, int *idef);

void calcpet(int *jday, int *month, float *tmin, float *tmax, float *tavg,
             float *solrad, float *rhumid, float *windsp, float *snow,
             int *petswitch, float fwloss[4], float *sitlat, float tmn2m[],
             float tmx2m[], float *petdly);

void daylen(int *curday, float *sitlat, float *dylngth);*/

void fracbslos(float *fbse, float *fbst, float blivelai);

void getdiff(float *valsubt, double valtochk, double valmin);

void h2oflux(int jday, int numlyrs, double swc[MAXLYR], double swcmin[MAXLYR],
             float minpot[MAXLYR], float depth[MAXLYR], float width[MAXLYR],
             float satcond[MAXLYR],  LAYERPAR_SPT layers,
             float soiltavg[MAXLYR], float watrinput, float bserate,
             float *soilEvap, float hours_rain_param, float dmpflux,
             float *aet, float *outflow, float wfluxout[MAXLYR],
             float snowpack, float *runoffdly, float *basef, float *baseflow,
             int watertable, float hpotdeep, float ksatdeep);

void hwdrain(double swc[MAXLYR], double *drain_out, int numlyrs,
             float swcfc[MAXLYR], float wfluxout[MAXLYR], int watertable,
             float thetas_bd[MAXLYR], float width[MAXLYR]);

void initdaily(int month, float biolive, float biodead, float blitter,
               float *biomass, float *blivelai, float *vegcov, float *totagb,
               float aggreenc, float *stcrlai, LAYERPAR_SPT layers);

void initlyrs(char *soilname, LAYERPAR_SPT layers, FLAG_SPT flags,
              SITEPAR_SPT sitepar);

void initphotosyn(char *paramFile);

void initsite(char *sitename, SITEPAR_SPT sitepar, LAYERPAR_SPT layers,
              FLAG_SPT flags, float sradadj[NMONTH], float *tminslope,
              float *tminintercept, float *maxphoto, float *bioabsorp);

void initsrad(double elevation, double latitude, double aspect, double slope,
              double ehoriz, double whoriz, float daylength[NDAY]);

/*void initsw(float *sitlat, float swcinit[MAXLYR], int *usexdrvrs,
              int *numlyrs, double *frac_nh4_fert,  double *frac_no3_fert,
              int *texture);*/

void litstcr_evap(float *cwlit, float *cwstcr, float *petleft, float *aet,
                  float totlit, float totstcr);

float petrad(int jday, int month, float cldcov, float windsp, float rhumid,
             float rlatitude, float avgtemp, float reflec, float solrad,
             float albedo, float snowpack);

float pevapdly(float tmin, float tmax, float sitlat, float tmn2m[],
               float tmx2m[], int jday, int month);

void potbse(float *bserate, float totagb, float fbse, float petday,
            LAYERPAR_SPT layers, float aggreenc);

void potbst(float *bstrate, float h2ogef, float biolive, float biodead,
            float fbst, float petday, float co2val);

void pteevap(float *bserate, float *bstrate, float petleft);

void rainflux(float watrinput, int infil_time, double swc[MAXLYR],
              float thetas_bd[MAXLYR], float width[MAXLYR],
              float satcond[MAXLYR], float swcfc[MAXLYR], float incap,
              float soiltavg[MAXLYR], int numlyrs, double *runoff,
              double *cumtime, double *cuminfl, double *drain,
              float wfluxout[MAXLYR], int *impedence);

void setamov(float amovdly[CENTMAXLYR], int nlayer, float wfluxout[MAXLYR],
             int numlyrs, int lbnd[CENTMAXLYR]);

void setasmos(float asmos[CENTMAXLYR], int *nlayer, float swcinit[MAXLYR],
              int *numlyrs, float avh2o[3], float rwcf[CENTMAXLYR]);

void showlyrs(double swc[MAXLYR], int numlyrs);

void soiltransp(float transp[MAXLYR], int numlyrs, float bstrate,
                LAYERPAR_SPT layers, float *aet);

void snowmodel(int jday, float tmax, float tmin, float tavg, float tmpcrit,
               float ppt, float *pptadj, float netrad, float relhum,
               float windsp, float cldcov, float rlatitude, float albedo,
               float *snow, float *Emelt, float *accum, float *Esublim,
               float *petlocal, float sublimscale, float petday,
               float tmelt[2], float *snlq, double srad, float daylength);

void snowCent(float tmelt[2], float tave, float pptactual, float *pptsoil,
              float *snow, float *snlq, float *petleft, float *melt,
              float *accum, float *sublim, float tmin, float tmax,
              double srad, float daylength);

float svapor(float atemp);

float swpotentl(double swc, int ilyr, LAYERPAR_SPT layers, char callname[]);

float tanfunc(float z, float a, float b, float c, float d);

float trwtavg(int *nlaypg, LAYERPAR_SPT layers, FLAG_SPT flags, 
              SITEPAR_SPT sitepar);

float watrate(float swp, float petday, float a, float b);

/*void watrbal(int *jdy, float *time, float *ppt, float *accum, float *melt,
             float *wbswc1, float *wbswc2, float *evap, float *transp,
             float *sublim, float *intrcpt, float *outflow, float *snlq1,
             float *snlq2, float *snow, float *runoffdly);*/

void watreqn(float sand, float clay, float *thetas, float *psis, float *b,
             float *binverse);

/*void watrflow(int *jday, int *month, int *nlayer, int *nlaypg,
              float *avgtemp, float *tempmin, float *tempmax, float *solrad,
              float *rhumid, float *windsp, float *pptactual, float *biolive,
              float *blitter, float *biodead, float rwcf[CENTMAXLYR],
              float avh2o[3], float asmos[CENTMAXLYR], float *snowpack,
              float *snlq, float amovdly[CENTMAXLYR], float *pet,
              float *evaptot, float *trantot, float *stream1, float *basef,
              float *pottransp, float *baseflow, float *accum, float *melt,
              float *intrcpt, float *outflow, float tmelt[], float *sublim,
              float wfluxout[], float *time, float *strplt, float *co2val,
              float *tmns, float *tmxs, float *runoffdly, float *trandep,
              float *soiltavewk);*/

void watrlit(float watrinput, float *pptsoil, float *wintlit, float blitter);

void watrstcr(float *pptleft, float *wintstcr, float ppt, float vegcov);

void wfps(LAYERPAR_SPT layers);

float wfunc_pulse(float *ppt, float *snow);

/*void wrtsoiln(float *time, int *jday, double *ammonium, double nitrate[]);

void wrtco2(float *time, int *jday, float co2PPM[]);*/

void wrtswc(FILE *fp, float time, int jday, double swc[], float width[],
            int numlyrs);

void wrtstemp(FILE *fp, float time, int jday, float soiltavg[], int numlyrs);

void wrtwfps(FILE *fp, float time, int jday, float wfps[], int numlyrs,
             float width[]);

void soiltemp(int jday, float tmin, float tmax, float depth[MAXLYR],
              float width[MAXLYR], float fieldc[MAXLYR], float sand[MAXLYR],
              float clay[MAXLYR], float org[MAXLYR], float bulkd[MAXLYR],
              double swc[MAXLYR], int numlyrs, float soiltavg[MAXLYR],
              float soiltmin[MAXLYR], float soiltmax[MAXLYR],
              float stmtemp[MAXSTLYR], float tmns, float tmxs,
              float *soiltavewk, float srfctemp, float diurnal_range);

void surftemp(float elitst, float pmxtmp, float pmntmp, float pmxbio,
              float tempmax, float tempmin, float *tmxs, float *tmns,
              float *srfctemp, float daylength, float biolive, float biodead,
              float blitter, float woodb, float stsys, float ststart,
              float stamt, float time, float snowpack, float *diurnal_range,
              float litrcrbn, int SnowFlag);

/*void wrtbio(float *time, int *curday, float *aglivc, float *bglivc,
            float *aglivn, float *bglivn, float *rleavc, float *frootc,
            float *fbrchc, float *rlwodc, float *crootc);*/
