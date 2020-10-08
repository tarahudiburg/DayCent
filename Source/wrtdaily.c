/*****************************************************************************
**
**  FILE:      wrtdaily.c
**
**  FUNCTION:  void wrtdaily()
**
**  PURPOSE:   This function writes daily values for evapotranspiration,
**             defac, soil temperature, snow water, thermal units, and solar
**             radiation.
**
**  INPUTS:
**    agdefac    - surface decomposition factor based on temperature and
**                 moisture
**    aglivc     - above ground live carbon (g/m2)
**    aggreenc   - the amount of photosynthetic active carbon (g/m2)
**    bgdefac    - soil decomposition factor based on temperature and moisture
**    curday     - the day of the year (1..366)
**    hwstress   - water stress term used to determine if full maturity has
**                 been reached (0.0-1.0)
**    petdly     - potential evapotranspiration rate for the day (cm H2O)
**    scenfrac   - multiplier used to indicate the fraction of the aboveground
**                 live carbon that is photosynthetic active carbon (0.0-1.0)
**                   1.0 = no senescence has occurred, 100% photosynthetic
**                         active carbon
**                   0.0 = full senescence, 0% photosynthetic active carbon
**    snlq       - liquid snow water content (cm H2O)
**    snow       - snowpack water content (cm H2O)
**    srad       - total incoming shortwave radiation (W/m2)
**    stemp      - average soil temperature near the soil surface (degrees C)
**    thermunits - accumulator of thermal units for growing degree day
**                 implementation
**    time       - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_daily    - file pointer to daily.out output file
**    files->write_daily - flag to indicate if daily.out output file should
**                         be created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    None
**
**  CALLED BY:
**    dailymoist()
**
*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "soilwater.h"

    void wrtdaily(float *time, int *curday, float *petdly, float *agdefac,
                  float *bgdefac, float *stemp, float *snow, float *snlq,
                  float *thermunits, float *aglivc, float *aggreenc,
                  float *hwstress, float *scenfrac, double *srad)
    {

      extern FILES_SPT files;

      if (!files->write_daily) {
        return;
      }

      fprintf(files->fp_daily, "%8.2f  %4d", *time, *curday);
      fprintf(files->fp_daily, "%12.4f  ", *petdly);
      fprintf(files->fp_daily, "%12.4f  ", *agdefac);
      fprintf(files->fp_daily, "%12.4f  ", *bgdefac);
      fprintf(files->fp_daily, "%12.4f  ", *stemp);
      fprintf(files->fp_daily, "%12.4f  ", *snow);
      fprintf(files->fp_daily, "%12.4f  ", *snlq);
      fprintf(files->fp_daily, "%12.4f  ", *thermunits);
      fprintf(files->fp_daily, "%12.4f  ", *aglivc);
      fprintf(files->fp_daily, "%12.4f  ", *aggreenc);
      fprintf(files->fp_daily, "%12.4f  ", *hwstress);
      fprintf(files->fp_daily, "%12.4f  ", *scenfrac);
      fprintf(files->fp_daily, "%12.4f  ", *srad);
      fprintf(files->fp_daily, "\n");

      return;
    }
