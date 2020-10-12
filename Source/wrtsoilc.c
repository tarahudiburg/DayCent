
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtsoilc.c
**
**  FUNCTION:  void wrtsoilc()
**
**  PURPOSE:   Write out the soil carbon values. 
**
**  AUTHOR:    Cindy Keough  10/02
** 
**  INPUTS:
**    curday    - the day of the year (1..366)
**    metabc(2) - carbon in metabolic component of soil litter (g/m2)
**    som1c1    - carbon in surface active soil organic matter (g/m2)
**    som1c2    - carbon in soil active soil organic matter (g/m2)
**    som2c1    - carbon in surface slow soil organic matter (g/m2)
**    som2c2    - carbon in soil slow soil organic matter (g/m2)
**    som3c     - carbon in passive soil organic matter (g/m2)
**    strucc(2) - carbon in structural component of soil litter (g/m2)
**    time      - simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files              - structure containing information about output files
**    files->fp_soilc    - file pointer to soilc.out output file
**    files->write_soilc - flag to indicate if soilc.out output file should
**                         be created, 0 = do not create, 1 = create
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

    void wrtsoilc(float *time, int *curday, float *metabc2, float *strucc2,
                  float *som1c1, float *som1c2, float *som2c1, float *som2c2,
                  float *som3c)
    {
      extern FILES_SPT files;

      if (!files->write_soilc) {
        goto ex;
      }

      fprintf(files->fp_soilc, "%6.2f  %2d  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *time, *curday, *metabc2, *strucc2, *som1c1, *som1c2);
      fprintf(files->fp_soilc, "%10.4f  %10.4f %10.4f\n", 
              *som2c1, *som2c2, *som3c);

ex:   return;
    }
