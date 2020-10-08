
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtyrcflows.c
**
**  FUNCTION:  void wrtyrcflows()
**
**  PURPOSE:   This function writes daily values for tracking carbon flows
**             from decomposition.
**
**  INPUTS:
**    ametc1tosom11  - annual accumulator for carbon flow from surface
**                     metabolic pool to fast surface organic matter
**                     pool (g/m^2)
**    ametc2tosom12  - annual accumulator for  carbon flow from soil
**                     metabolic pool to fast soil organic matter pool
**                     (g/m^2)
**    asom11tosom21  - annual accumulator of carbon flow from fast surface
**                     organic matter pool to intermediate surface organic
**                     matter pool (g/m^2)
**    asom12tosom22  - annual accumulator of carbon flow from fast soil
**                     organic matter pool to intermediate soil organic matter
**                     pool (g/m^2)
**    asom12tosom3   - annual accumulator of carbon flow from fast soil
**                     organic matter pool to slow soil organic matter pool
**                     (g/m^2)
**    asom21tosom11  - annual accumulator of carbon flow from intermediate
**                     surface organic matter pool to fast surface organic
**                     matter pool (g/m^2)
**    asom21tosom22  - annual accumulator of carbon flow from intermediate
**                     surface organic matter pool to intermediate soil
**                     organic matter pool (g/m^2)
**    asom22tosom12  - annual accumulator of carbon flow from intermediate
**                     soil organic matter pool to fast soil organic matter
**                     pool (g/m^2)
**    asom22tosom3   - annual accumulator of carbon flow from intermediate
**                     soil organic matter pool to slow soil organic matter
**                     pool (g/m^2)
**    asom3tosom12   - annual accumulator of carbon flow from slow soil
**                     organic matter pool to fast soil organic matter pool
**                     (g/m^2)
**    astruc1tosom11 - annual accumulator for carbon flow from surface
**                     structural pool to fast surface organic matter
**                     pool (g/m^2)
**    astruc1tosom21 - annual accumulator for carbon flow from surface
**                     structural pool to intermediate surface organic
**                     matter pool (g/m^2)
**    astruc2tosom12 - annual accumulator for carbon flow from soil
**                     structural pool to fast soil organic matter pool
**                     (g/m^2)
**    astruc2tosom22 - annual accumulator for carbon flow from soil
**                     structural pool to intermediate soil organic
**                     matter pool (g/m^2)
**    awood1tosom11  - annual accumulator for carbon flow from dead
**                     fine branch pool to fast surface organic matter pool
**                     (g/m^2)
**    awood1tosom21  - annual accumulator for carbon flow from dead fine
**                     branch pool pool to intermediate surface organic matter
**                     pool (g/m^2)
**    awood2tosom11  - annual accumulator for carbon flow from dead large wood
**                     pool to fast surface organic matter pool (g/m^2)
**    awood2tosom21  - annual accumulator for carbon flow from dead large wood
**                     pool to intermediate surface organic matter pool
**                     (g/m^2)
**    awood3tosom12  - annual accumulator for carbon flow from dead coarse
**                     root pool to fast soil organic matter pool (g/m^2)
**    awood3tosom22  - annual accumulator for carbon flow from dead coarse
**                     root pool to intermediate soil organic matter pool
**                     (g/m^2)
**    time           - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_yrcflows    - file pointer to year_cflows.out output file
**    files->write_yrcflows - flag to indicate if yearcflows.out output file
**                            should be created, 0 = do not create, 1 = create
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

    void wrtyrcflows(float *time, float *asom11tosom21, float *asom12tosom22,
                     float *asom12tosom3, float *asom21tosom11,
                     float *asom21tosom22, float *asom22tosom12,
                     float *asom22tosom3, float *asom3tosom12,
                     float *ametc1tosom11, float *ametc2tosom12,
                     float *astruc1tosom11, float *astruc1tosom21,
                     float *astruc2tosom12, float *astruc2tosom22,
                     float *awood1tosom11, float *awood1tosom21,
                     float *awood2tosom11, float *awood2tosom21,
                     float *awood3tosom12, float *awood3tosom22)
    {

      extern FILES_SPT files;

      if (!files->write_yrcflows) {
        return;
      }

      fprintf(files->fp_yrcflows, "%8.2f", *time);
      fprintf(files->fp_yrcflows, "%12.6f  ", *asom11tosom21);
      fprintf(files->fp_yrcflows, "%12.6f  ", *asom12tosom22);
      fprintf(files->fp_yrcflows, "%12.6f  ", *asom12tosom3);
      fprintf(files->fp_yrcflows, "%12.6f  ", *asom21tosom11);
      fprintf(files->fp_yrcflows, "%12.6f  ", *asom21tosom22);
      fprintf(files->fp_yrcflows, "%12.6f  ", *asom22tosom12);
      fprintf(files->fp_yrcflows, "%12.6f  ", *asom22tosom3);
      fprintf(files->fp_yrcflows, "%12.6f  ", *asom3tosom12);
      fprintf(files->fp_yrcflows, "%12.6f  ", *ametc1tosom11);
      fprintf(files->fp_yrcflows, "%12.6f  ", *ametc2tosom12);
      fprintf(files->fp_yrcflows, "%12.6f  ", *astruc1tosom11);
      fprintf(files->fp_yrcflows, "%12.6f  ", *astruc1tosom21);
      fprintf(files->fp_yrcflows, "%12.6f  ", *astruc2tosom12);
      fprintf(files->fp_yrcflows, "%12.6f  ", *astruc2tosom22);
      fprintf(files->fp_yrcflows, "%12.6f  ", *awood1tosom11);
      fprintf(files->fp_yrcflows, "%12.6f  ", *awood1tosom21);
      fprintf(files->fp_yrcflows, "%12.6f  ", *awood2tosom11);
      fprintf(files->fp_yrcflows, "%12.6f  ", *awood2tosom21);
      fprintf(files->fp_yrcflows, "%12.6f  ", *awood3tosom12);
      fprintf(files->fp_yrcflows, "%12.6f  ", *awood3tosom22);
      fprintf(files->fp_yrcflows, "\n");

      return;
    }
