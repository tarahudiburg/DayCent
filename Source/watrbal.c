
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      watrbal.c
**
**  FUNCTION:  void watrbal()
**
**  PURPOSE:   Report all components of the water balance at the end of each
**             day.
**
**  AUTHOR:    Melannie Hartman
**             10/18/96
**
**  INPUTS:
**    accum     - the amount of snow added to the snowpack (cm H2O)
**    evap      - water evaporated from soil (cm H2O)
**    intrcpt   - amount of precipitation intercepted by the standing crop and
**                litter (cm H2O)
**    jdy       - current day of the year (1..366)
**    melt      - the amount of snow melted from the snowpack, if 
**                daily air temperature is warm enough (cm H2O)
**    outflow   - water that runs off, or drains out of the profile (cm H2O)
**    ppt       - precipitation for the day (cm)
**    runoffdly - amount of water (rain or snowmelt) which did not infiltrate
**                soil profile (cm)
**    snlq1     - the liquid water in the snowpack at the beginning of the day
**                (cm H2O)
**    snlq2     - the liquid water in the snowpack at the end of the day
**                (cm H2O)
**    snow      - current snowpack (equiv. cm H2O)
**    sublim    - amount of water sublimated from the snowpack (cm H2O)
**    transp    - water transpired from soil (cm H2O)
**    time      - simulation time (years)
**    wbswc1    - amount of water in the soil at the beginning of the day
**                (cm H2O)
**    wbswc2    - amount of water in the soil at the end of the day (cm H2O)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files           - structure containing information about output files
**    files->fp_wb    - file pointer to watrbal.out output file 
**    files->write_wb - flag to indicate if watrbal.out output file should be
**                      created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    balance - water balance (should equal zero)
**    dsnlq   - difference between the liquid water in the snowpack at the end
**              of the day and the liquid water in the snowpack at the
**              beginning of the day (cm H2O)
**
**  OUTPUTS:
**    None
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    void watrbal(int *jdy, float *time, float *ppt, float *accum, float *melt,
                 float *wbswc1, float *wbswc2, float *evap, float *transp,
                 float *sublim, float *intrcpt, float *outflow, float *snlq1,
                 float *snlq2, float *snow, float *runoffdly)
    {
      extern FILES_SPT files;
      float balance, dsnlq;

      if (!files->write_wb) {
        return;
      }

      /* swc2 = swc1 + ppt + melt - snow_accum - intrcpt - evap - transp - */
      /*        outflow */
      /* 0    = swc1 - swc2 + melt - snow_accum - intrcpt - evap - transp - */
      /*        outflow */

      dsnlq = *snlq2 - *snlq1;
      /* corrected the snow storage term in the balance equation
         Old equation:
         balance = *ppt + *melt - *accum - *intrcpt - *evap - *transp -
                *outflow + *wbswc1 - *wbswc2; 
         to properly account for rain-on-snow events. The snow storage term
         snowstor = snow1 - snow2 + snlq1 - snlq2  =  *ppt - *melt - sublim
         when snow is present. When snow is *NOT* present, melt = 0 doubling ppt.
         accum can NOT be used since it is snow but and does not reflect rain
         stored in snlq.
         balance =  ppt - snowstor - sublim + dswc -intrcpt-evap-transp-outflow
          =  ppt - sublim + dswc -intrcpt-evap-transp-outflow (NO snow)
          =  ppt - (ppt - *melt - sublim) - sublim + dswc -
             intrcpt-evap-transp-outflow (snow)
          =  *melt + dswc  -intrcpt-evap-transp-outflow (snow)
          The best indicator of snow is sublim <> 0
      */
      if(*sublim == 0.0) {balance = *ppt;} else {balance = *melt;}
      balance += (*wbswc1 - *wbswc2) - *intrcpt - *evap - *transp - *outflow;

      fprintf(files->fp_wb,"%7.2f %4d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f ",
              *time, *jdy, *ppt, -(*accum), dsnlq, *melt, -(*intrcpt),
              -(*evap));
      fprintf(files->fp_wb,"%7.3f %7.3f %8.5f %7.3f %9.5f ", -(*transp),
              -(*sublim), (*wbswc1 - *wbswc2), -(*outflow), balance);
      fprintf(files->fp_wb,"%7.3f %7.3f %7.3f \n", *snow, *snlq2, *runoffdly);

      return;
    }
