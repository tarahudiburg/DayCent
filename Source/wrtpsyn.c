
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtpsyn.c
**
**  FUNCTION:  void wrtpsyn()
**
**  PURPOSE:   Write out the photosynthesis values. 
**
**  AUTHOR:    Cindy Keough  01/2009
** 
**  INPUTS:
**    aetdly         - actual evapotranspiration (cm H2O)
**    annPrecip      - average annual precipitation for site (cm)
**    average_temp   - average temperature for daylight hours (degrees C)
**    average_vpd    - average vapor pressure deficit (kPa)
**    curday         - the day of the year (1..366)
**    dailyPrecip    - precipitation for current day (cm)
**    daytime        - fraction of day that has sunlight (0.0 - 1.0)
**    crpGrossPsn    - gross photosynthesis, with water stress, for grass/crop
**                     system (g C * m^-2 ground area * day^-1)
**    crpLAI         - grass/crop system leaf area index
**    crpLightEff    - decrease in photosynthesis due to amount of light
**                     absorbed for grass/crop system
**    crpPotGrossPsn - potential photosynthesis, without water stress, for
**                     grass/crop system (g C * m^-2 ground area * day^-1)
**    crpdTemp       - decrease in photosynthesis due to temperature for
**                     grass/crop system
**    crpdVpd        - decrease in photosynthesis due to vapor pressure
**                     deficit for grass/crop system
**    crpdWater      - effect of water stress on photosynthesis for grass/crop
**                     system
**    forGrossPsn    - gross photosynthesis, with water stress, for forest
**                     system (g C * m^-2 ground area * day^-1)
**    forLAI         - forest system leaf area index
**    forLightEff    - decrease in photosynthesis due to amount of light
**                     absorbed for forest system
**    forPotGrossPsn - potential photosynthesis, without water stress, for
**                     forest system (g C * m^-2 ground area * day^-1)
**    fordTemp       - decrease in photosynthesis due to temperature for
**                     forest system
**    fordVpd        - decrease in photosynthesis due to vapor pressure
**                     deficit for forest system
**    fordWater      - effect of water stress on photosynthesis for forest
**                     system
**    maxTemp        - maximum temperature for day (degrees C)
**    minTemp        - minimum temperature for day (degrees C)
**    petdly         - potential evapotranspiration rate for day (cm H2O)
**    srad           - shortwave radiation value for day (W/m^2)
**    time           - simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files             - structure containing information about output files
**    files->fp_psyn    - file pointer to psyndy.out output file
**    files->write_psyn - flag to indicate if psyndy.out output file should
**                        be created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**     None
**
**  CALLED BY:
**     simsom()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    void wrtpsyn(float *time, int *curday, double *minTemp, double *maxTemp,
                 double *annPrecip, double *dailyPrecip, float *aetdly,
                 float *petdly, double *daytime, double *srad,
                 double *average_temp, double *average_vpd, double *crpLAI,
                 double *crpdTemp, double *crpdVpd, float *crpdWater,
                 double *crpLightEff, double *crpPotGrossPsn,
                 double *crpGrossPsn, double *forLAI, double *fordTemp,
                 double *fordVpd, float *fordWater, double *forLightEff,
                 double *forPotGrossPsn, double *forGrossPsn)
    {

      extern FILES_SPT files;

      if (!files->write_psyn) {
        goto ex;
      }

      fprintf(files->fp_psyn, "%6.2f  %8d  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *time, *curday, *minTemp, *maxTemp, *annPrecip, *dailyPrecip);
      fprintf(files->fp_psyn, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *aetdly, *petdly, *daytime, *srad, *average_temp);
      fprintf(files->fp_psyn, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *average_vpd, *crpLAI, *crpdTemp, *crpdVpd, *crpdWater);
      fprintf(files->fp_psyn, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *crpLightEff, *crpPotGrossPsn, *crpGrossPsn, *forLAI,
              *fordTemp);
      fprintf(files->fp_psyn, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n", 
              *fordVpd, *fordWater, *forLightEff, *forPotGrossPsn,
              *forGrossPsn);

ex:   return;
    }
