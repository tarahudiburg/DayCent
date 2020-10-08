
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:     calcSrad.h
**
**  PURPOSE:  Header file for calcSrad.c.
**
**  AUTHOR:   Adapted from Peter Thorton's code extracted from the Sipnet
**            model (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**            Cindy Keough - 04/2009
**
*****************************************************************************/

/***********************
**                    **
**  GLOBAL CONSTANTS  **
**                    **
***********************/
#define ABASE     -6.1e-5     /* vapor pressure effect on transmittance */
                              /* (1/Pa) */
#define B0          0.013     /* radiation parameter (dim) */
#define B1          0.201     /* radiation parameter (dim) */
#define B2          0.185     /* radiation parameter (dim) */
#define C           1.5       /* radiation parameter (dim) */
#define CP       1010.0       /* specific heat of air (J kg-1 K-1) */
#define DAYSOFF    11.25      /* julian day offset of winter solstice */
#define DIF_ALB     0.6       /* diffuse albedo for horizon correction */
                              /* (dim) */
#define EPS      0.62196351   /* unitless ratio of molec weights (MW/MA) */
#define G_STD    9.80665      /* standard gravitational accel. (m s-2) */ 
#define LR_STD   0.0065       /* standard temperature lapse rate (-K m-1) */
#define MA       28.9644e-3   /* molecular weight of air (kg mol-1) */
#define MINDECL -0.4092797    /* minimum declination (radians) */
#define P_STD    101325.0     /* standard pressure at 0.0 m elevation (Pa) */
#define PI       3.14159265   /* pi */
#define R        8.3143       /* gas law constant (m3 Pa mol-1 K-1) */
#define RADPERDAY 0.017214    /* radians of Earth orbit per julian day */
#define RADPERDEG 0.01745329  /* radians per degree */
#define RAIN_SCALAR  0.75     /* correction to trans. for rain day (dim) */
#define SECPERRAD 13750.9871  /* seconds per radian of hour angle */
#define SNOW_TCRIT   -6.0     /* critical temperature for snowmelt (deg C) */   
#define SNOW_TRATE  0.042     /* snowmelt rate (cm/degC/day) */
#define SRADDT 600.0          /* timestep for radiation routine (seconds) */
#define T_STD    288.15       /* standard temp at 0.0 m elevation (K) */ 
#define TBASE       0.870     /* max inst. trans., 0m, nadir, dry atm (dim) */
#define TDAYCOEF     0.45     /* daylight air temperature coefficient (dim) */


/****************************
**                         **
**  STRUCTURE DEFINITIONS  **
**                         **
****************************/
typedef struct {
  double site_lat;       /* site latitude, dec. degrees (- for south) */
  double site_elev;      /* site elevation, meters */
  double site_slp;       /* site slope, degrees */
  double site_asp;       /* site aspect, degrees */
  double site_ehoriz;    /* site east horizon, degrees */
  double site_whoriz;    /* site west horizon, degrees */
} PARAMETER_S, *PARAMETER_SPT;

typedef struct {
  double s_tdew;         /* site dewpoint temperature value (degrees C) */
  double s_tmax;         /* site tmax value */
  double s_tmin;         /* site tmin value */
  double s_tday;         /* site daylight temperature value */
  double s_prcp;         /* site prcp value */
  double s_srad;         /* site shortwave radiation value */
  double s_dayl;         /* site daylength value */
  double s_swe;          /* site snow water equivalent value (cm) */
} DATA_S, *DATA_SPT;

typedef struct {
  double ttmax0[366];
  double flat_potrad[366];
  double slope_potrad[366];
  double daylength[366];
} YEARDATA_S, *YEARDATA_SPT;

typedef struct {
  int indewpt;           /* input dewpoint temperature flag (0=NO, 1=YES) */
} CONTROL_S, *CONTROL_SPT;


/***************************
**                        **
**  FUNCTION DEFINITIONS  **
**                        **
***************************/
double atm_pres(double elev);

double calc_pet(double rad, double ta, double pa, double dayl);

void calc_srad_humidity(CONTROL_SPT ctrl, PARAMETER_SPT params, 
                        DATA_SPT sitedata, YEARDATA_SPT yeardata, int yday);

void calc_srad_humidity_iterative(CONTROL_SPT ctrl, PARAMETER_SPT params,
                                  DATA_SPT sitedata, YEARDATA_SPT yeardata,
                                  int yday);

void snowpack(DATA_SPT sitedata);
