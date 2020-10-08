/*  swconst.h */

#define SEC_PER_DAY 86400
#define SEC_PER_HOUR 3600
#define HOURS_PER_DAY 24
#define NDAY 366               /* dimension for daily arrays */
#define NMONTH 12              /* # of months in a year */
#define MAXLYR 21              /* Max # of soil layers (soil water model) */
#define CENTMAXLYR 10          /* Max # of century soil layers */
#define MAXSTLYR 200           /* Max # of soil temperature layers */
#define NTDEPTHS 4             /* Max # of soil regions */
#define VALMISSING 99          /* Used for value missing from histweth file */
#define BAR2CM 1024            /* 1 bar = 1024 cm H2O */
#define PARTDENS 2.65          /* Particle Density (g/cm3) */
#define PI 3.14159265
#define FRZSOIL -0.5           /* temperature at which soil is considered */
                               /* frozen (deg C) */
#define MIN_FRZN_COND 0.00005  /* minimum hydraulic conductivity in frozen */
                               /* soil (cm/sec) */

#define max(a,b)        (((a) > (b)) ? (a) : (b))
#define min(a,b)        (((a) < (b)) ? (a) : (b))


#define FNSOIL "soils.in"
#define FNSITE "sitepar.in"

#define CONVLAI 100    /* biomass needed to produce an LAI of 1 (g/m**2) */
