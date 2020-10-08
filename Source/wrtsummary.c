/*****************************************************************************
**
**  FILE:      wrtsummary.c
**
**  FUNCTION:  void wrtsummary()
**
**  PURPOSE:   This function writes daily values for daily climate, trace gas
**             fluxes, and CO respiration.
**
**  INPUTS:
**    CH4_oxid  - Methane oxidation (gCH4/ha)
**    CO2resp   - Heterotrophic CO2 respiration for the day (gCO2/ha)
**    curday    - the day of the year (1..366)
**    N2Oflux   - Nitrous oxide flux (gN/ha)
**    nit_amt   - Gross nitrification (gN/ha)
**    NOflux    - Nitric oxide flux (gN/ha)
**    ntCO2resp - "No-till" heterotrophic CO2 respiration removing cultivation
**                effects (gCO2/ha)
**    ppt       - Precipitation for day (cm)
**    tempmax   - Maximum temperature for day (degrees C)
**    tempmin   - Minimum temperature for day (degrees C)
**    time      - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_summary    - file pointer to summary.out output file
**    files->write_summary - flag to indicate if summary.out output file should
**                           be created, 0 = do not create, 1 = create
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

    void wrtsummary(float *time, int *curday, float *tempmax, float *tempmin,
                    float *ppt, double *N2Oflux, double *NOflux,
                    double *CH4_oxid, double *nit_amt, double *CO2resp,
                    float *ntCO2resp)
    {

      extern FILES_SPT files;

      if (!files->write_summary) {
        return;
      }

      fprintf(files->fp_summary, "%8.2f  %4d", *time, *curday);
      fprintf(files->fp_summary, "%12.4f  ", *tempmax);
      fprintf(files->fp_summary, "%12.4f  ", *tempmin);
      fprintf(files->fp_summary, "%12.4f  ", *ppt);
      fprintf(files->fp_summary, "%12.4f  ", *N2Oflux);
      fprintf(files->fp_summary, "%12.4f  ", *NOflux);
      fprintf(files->fp_summary, "%12.4f  ", *CH4_oxid);
      fprintf(files->fp_summary, "%12.4f  ", *nit_amt);
      fprintf(files->fp_summary, "%12.4f  ", *CO2resp);
      fprintf(files->fp_summary, "%12.4f  ", *ntCO2resp);
      fprintf(files->fp_summary, "\n");

      return;
    }
