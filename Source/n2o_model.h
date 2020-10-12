/*****************************************************************************
**
**  FILE:    n2o_model.h
**
**  AUTHOR:  Melannie D. Hartman
**           Natural Resource Ecology Laboratory
**           Colorado State University
**           11/10/95, 2/13/97
**
*****************************************************************************/

#define COARSE 1           /* texture type for sandy soil */
#define MEDIUM 2           /* texture type for medium (loamy) soil */
#define FINE 3             /* texture type for fine soil */
#define VERYFINE 4         /* texture type for volcanic soil */

#define MIN_NH4_CONC 0.05  /* minimum NH4 concentration (ppm) */
#define MIN_NO3_CONC 0.05  /* minimum NO3 concentration (ppm) */

void nitrify(double *ammonium, double *nh4_2_no3, float *maxt, float *nreduce,
             float *pHscale);

void denitrify(float *newCO2, double *newNO3, double nitrate[],
               float tfluxout[], float *critflow, float frlechd[],
               float stream[], float *basef, float* stormf, double *inorglch,
               double *Dn2oflux, double *Dn2flux, float stdfieldc,
               float stdbulkd, float *efscltef, float co2PPM[],
               double dN2lyr[], double dN2Olyr[], int *jday,
               int is_saturated);

float nox_pulse(float *ppt, float *snow);

void getsoilprop(float *asand, float *asilt, float *aclay, float *bulkden,
                 float *fieldcap, int *texture);

float diffusiv(float *A, float *bulkden, float *wfps);

void leachdly(float tfluxout[], int numlyrs, double nitrate[], float critflow,
              float frlechd[], float stream[], float basef, float stormf,
              double *inorglch);

void methane_emission(float *aglivc, float *tmxbio, float *bglivc,
                      float *avgst_10cm, float *CH4_prod, float *CH4_Ep,
                      float *CH4_Ebl, float zero_root_frac);

void methane_oxidation(double *CH4_oxid, int *isdecid, int *isagri);

void methane_production(float *prev_bgprd, float *Com, float *avgst_10cm,
                        float *TI, float *SI, float *Cr, float *Eh,
                        float *Feh, float *CH4_prod, int *watertable,
                        int *watrflag);

/*void update_npool(int *clyr, float *amt, double *frac_nh4, double *frac_no3,
                  double *ammonium, double nitrate[], char subname[]);

void bal_npool(int *nlayer, float minerl[], double *ammonium,
               double nitrate[], double *inorglch);*/

float f_allometric(float x, float A[]);

float f_arctangent(float x, float A[]);

float f_exponential(float x, float A[]);

float f_gen_gompertz(float x, float A[]);

float f_logistic(float x, float A[]);

float f_gen_poisson_density(float x, float A[]);
