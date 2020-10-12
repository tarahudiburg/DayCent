
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:     calcPhotosyn.h
**
**  PURPOSE:  Header file for calcPhotosyn.c.
**
**  AUTHOR:   Adapted from code extracted from the Sipnet model
**            (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**            Cindy Keough - 04/2009
**
*****************************************************************************/

/***********************
**                    **
**  GLOBAL CONSTANTS  **
**                    **
***********************/
#define B            2.1f          /* coefficient that controls temperature */
                                   /* decrease at night */
                                   /* "Stream Temperature Investigations: */
                                   /* Field and Analytic Methods, Instream */
                                   /* Flow Information Paper No. 13" by */
                                   /* John M. Bartholow, p. 41. suggests */
                                   /* using a value of of 2.1 */
#define C_WEIGHT     12.0          /* molecular weight of carbon */
#define NUM_LAYERS   6             /* number of canopy layers used in the */
                                   /* calculation of light absorbed */
                                   /* NOTE:  6 layers gives approximately */
                                   /*        the same result as 100 layers */
#ifndef SEC_PER_DAY
  #define SEC_PER_DAY  86400.0     /* number of seconds in 24 hours */
#endif
#define TEN_9        1000000000.0  /* for conversions from nano  */
#define TIME_LAG_MAX  1.8f         /* time lag in maximum temperature after */
                                   /* noon (hours) */
                                   /* "A Model for Diurnal Variation in */
                                   /* Soil and Air Temperature" by William */
                                   /* J. Parton and Jesse A. Logan suggests */
                                   /* using a value of 1.8 */
#define TIME_LAG_MIN  0.0f         /* time lag for the minimum temperature */
                                   /* after sunrise (hours) */
#define TINY         0.000001      /* to avoid divide-by-zero errors */


/****************************
**                         **
**  STRUCTURE DEFINITIONS  **
**                         **
****************************/
typedef struct {
  double length;   /* length of this timestep (in days) (number of */
                   /* daylight hours divided by 24) */
  double par;      /* average par for this time step */
                   /* (Einsteins * m^-2 ground area * day^-1) */
                   /* NOTE: input is in Einsteins * m^-2 ground area, */
                   /* summed over entire time step */
  double tair;     /* avg. air temp for this time step (degrees C) */
  double vpd;      /* average vapor pressure deficit (kPa) */
  double minTemp;  /* minimum temperature for interval (degrees C) */
  double maxTemp;  /* maximum temperature for interval (degrees C) */
} CLIMATE_S, *CLIMATE_SPT;

/* parameter values read from crop or tree option */
typedef struct {
  /* photosynthesis */
  double aMax;                 /* maximum photosynthesis */
                               /* (nmol CO2 * g^-1 leaf * sec^-1) */
                               /* assuming maximum possible par, all */
                               /* intercepted, no temperature, water or vpd */
                               /* stress */
  double aMaxFrac;             /* avg. daily params.aMax as fraction of */
                               /* instantaneous */
  double baseFolRespFrac;      /* basal foliage resp. rate, as % of maximum */
                               /* net photosynthesis rate */
  double psnTMin;              /* minimum tempature at which net */
                               /* photosynthesis occurs (degrees C) */
  double psnTOpt;              /* optimal temperature at which net */
                               /* photosynthesis occurs (degrees C) */
  double dVpdSlope;            /* slope value in dVpd equation */
  double dVpdExp;              /* exponential value in dVpd equation */
  double halfSatPar;           /* par at which photosynthesis occurs at 1/2 */
                               /* theoretical maximum */
                               /* (Einsteins * m^-2 ground area * day^-1) */
  double attenuation;          /* light attenuation coefficient */
  double leafCSpWt;            /* grams of carbon in a square meter of leaf */
                               /* area */
  double cFracLeaf;            /* factor for converting leaf biomass to */
                               /* carbon */
                               /* (leaf biomass * cFracLeaf = leaf carbon) */
  double psnTMax;              /* maximum temperature at which net */
                               /* photosynthesis occurs, assumed */
                               /* symmetrical around psnTOpt (degrees C) */
  double growthDays1;          /* number of days after germination to start */
                               /* using aMaxScalar1 */
  double growthDays2;          /* number of days after germination to start */
                               /* using aMaxScalar2 */
  double growthDays3;          /* number of days after germination to start */
                               /* using aMaxScalar3 */
  double growthDays4;          /* number of days after germination to start */
                               /* using aMaxScalar4 */
  double aMaxScalar1;          /* multiplier used to adjust aMax based on */
                               /* growthDays1 days since germination */
  double aMaxScalar2;          /* multiplier used to adjust aMax based on */
                               /* growthDays2 days since germination */
  double aMaxScalar3;          /* multiplier used to adjust aMax based on */
                               /* growthDays3 days since germination */
  double aMaxScalar4;          /* multiplier used to adjust aMax based on */
                               /* growthDays4 days since germination */
} PSPARAMS_S, *PSPARAMS_SPT;


/***************************
**                        ** 
**  FUNCTION DEFINITIONS  **
**                        **
***************************/
void calc_vpd(double *average_vpd, double *average_temp, double maxTemp,
              double minTemp, float dayLength, float tminslope,
              float tminintercept);

void calcLightEff3 (double *lightEff, double lai, double par);

void calcMoistureEff(double etrans, double *dWater, double potETrans);

void calcTemperatureEff(double tair, double *dTemp);

void calcVpdEff(double vpd, double *dVpd);

double hourTemp(double *sunrise, double *sunset, int hour, double maxTemp,
                double minTemp, float dayLength);

FILE *openFile(char *name, char *mode);

void potPsn(double *potGrossPsn, double *lai, double par, double *dTemp,
            double *dVpd, double *lightEff, int *grwdys);

float svapor(float);