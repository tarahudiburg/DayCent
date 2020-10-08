/*****************************************************************************
**
**  FILE:      wrtcflows.c
**
**  FUNCTION:  void wrtcflows()
**
**  PURPOSE:   This function writes daily values for tracking carbon flows
**             from decomposition.
**
**  INPUTS:
**    curday         - the day of the year (1..366)
**    metc1tosom11   - carbon flow from surface metabolic pool to fast
**                     surface organic matter pool (g/m^2)
**    metc2tosom12   - carbon flow from soil metabolic pool to fast soil
**                     organic matter pool (g/m^2)
**    som11tosom21   - carbon flow from fast surface organic matter pool
**                     to intermediate surface organic matter pool (g/m^2)
**    som12tosom22   - carbon flow from fast soil organic matter pool
**                     to intermediate soil organic matter pool (g/m^2)
**    som12tosom3    - carbon flow from fast soil organic matter pool to
**                     slow soil organic matter pool (g/m^2)
**    som21tosom11   - carbon flow from intermediate surface organic
**                     matter pool to fast surface organic matter pool (g/m^2)
**    som21tosom22   - carbon flow from intermediate surface organic
**                     matter pool to intermediate soil organic matter
**                     pool (g/m^2)
**    som22tosom12   - carbon flow from intermediate soil organic matter
**                     pool to fast soil organic matter pool (g/m^2)
**    som22tosom3    - carbon flow from intermediate soil organic matter
**                     pool to slow soil organic matter pool (g/m^2)
**    som3tosom12    - carbon flow from slow soil organic matter pool to
**                     fast soil organic matter pool (g/m^2)
**    struc1tosom11  - carbon flow from surface structural pool to fast
**                    surface organic matter pool (g/m^2)
**    struc1tosom21  - carbon flow from surface structural pool to
**                     intermediate surface organic matter pool (g/m^2)
**    struc2tosom12  - carbon flow from soil structural pool to fast
**                     soil organic matter pool (g/m^2)
**    struc2tosom22  - carbon flow from soil structural pool to
**                     intermediate soil organic matter pool (g/m^2)
**    time           - current simulation time (years)
**    wood1tosom11   - carbon flow from dead fine branch pool to fast
**                     surface organic matter pool (g/m^2)
**    wood1tosom21   - carbon flow from dead fine branch pool to
**                     intermediate surface organic matter pool (g/m^2)
**    wood2tosom11   - carbon flow from dead large wood pool to fast
**                     surface organic matter pool (g/m^2)
**    wood2tosom21   - carbon flow from dead large wood pool pool to
**                     intermediate surface organic matter pool (g/m^2)
**    wood3tosom12   - carbon flow from dead coarse root pool to fast
**                     soil organic matter pool (g/m^2)
**    wood3tosom22   - carbon flow from dead coarse root pool to
**                     intermediate soil organic matter pool (g/m^2)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_cflows    - file pointer to cflows.out output file
**    files->write_cflows - flag to indicate if cflows.out output file should
**                          be created, 0 = do not create, 1 = create
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

    void wrtcflows(float *time, int *curday, float *som11tosom21,
                   float *som12tosom22, float *som12tosom3,
                   float *som21tosom11, float *som21tosom22,
                   float *som22tosom12, float *som22tosom3,
                   float *som3tosom12, float *metc1tosom11,
                   float *metc2tosom12, float *struc1tosom11,
                   float *struc1tosom21, float *struc2tosom12,
                   float *struc2tosom22, float *wood1tosom11,
                   float *wood1tosom21, float *wood2tosom11,
                   float *wood2tosom21, float *wood3tosom12,
                   float *wood3tosom22)
    {

      extern FILES_SPT files;

      if (!files->write_cflows) {
        return;
      }

      fprintf(files->fp_cflows, "%8.2f  %4d", *time, *curday);
      fprintf(files->fp_cflows, "%12.6f  ", *som11tosom21);
      fprintf(files->fp_cflows, "%12.6f  ", *som12tosom22);
      fprintf(files->fp_cflows, "%12.6f  ", *som12tosom3);
      fprintf(files->fp_cflows, "%12.6f  ", *som21tosom11);
      fprintf(files->fp_cflows, "%12.6f  ", *som21tosom22);
      fprintf(files->fp_cflows, "%12.6f  ", *som22tosom12);
      fprintf(files->fp_cflows, "%12.6f  ", *som22tosom3);
      fprintf(files->fp_cflows, "%12.6f  ", *som3tosom12);
      fprintf(files->fp_cflows, "%12.6f  ", *metc1tosom11);
      fprintf(files->fp_cflows, "%12.6f  ", *metc2tosom12);
      fprintf(files->fp_cflows, "%12.6f  ", *struc1tosom11);
      fprintf(files->fp_cflows, "%12.6f  ", *struc1tosom21);
      fprintf(files->fp_cflows, "%12.6f  ", *struc2tosom12);
      fprintf(files->fp_cflows, "%12.6f  ", *struc2tosom22);
      fprintf(files->fp_cflows, "%12.6f  ", *wood1tosom11);
      fprintf(files->fp_cflows, "%12.6f  ", *wood1tosom21);
      fprintf(files->fp_cflows, "%12.6f  ", *wood2tosom11);
      fprintf(files->fp_cflows, "%12.6f  ", *wood2tosom21);
      fprintf(files->fp_cflows, "%12.6f  ", *wood3tosom12);
      fprintf(files->fp_cflows, "%12.6f  ", *wood3tosom22);
      fprintf(files->fp_cflows, "\n");

      return;
    }
