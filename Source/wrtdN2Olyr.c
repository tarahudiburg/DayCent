/*****************************************************************************
**
**  FILE:      wrtdN2Olyr.c
**
**  FUNCTION:  void wrtdN2Olyr()
**
**  PURPOSE:   This function writes out daily N2O flux from denitrification
**             by layer (gN/m^2).
**
**  INPUTS:
**    dN2Olyr[] - N2O flux by layers (gN/m^2)
**    jday     - current day of the year (1..366)
**    time     - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_dN2Olyr    - file pointer to dN2Olyr.out output file
**    files->write_dN2Olyr - flag to indicate if dN2Olyr.out output file
**                           should be created, 0 = do not create, 1 = create
**    layers->numlyrs      - total number of layers in the soil water model
**                           soil profile
**
**  LOCAL VARIABLES:
**    ilyr - current layer in the soil profile
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

    void wrtdn2olyr(float *time, int *jday, double dN2Olyr[])
    {

      int    ilyr;
      extern LAYERPAR_SPT layers;
      extern FILES_SPT files;

      if (!files->write_dN2Olyr) {
        return;
      }

      fprintf(files->fp_dN2Olyr, "%8.2f  %4d  ", *time, *jday);

      for (ilyr=0; ilyr < layers->numlyrs; ilyr ++) {

        fprintf(files->fp_dN2Olyr, "%12.6f  ", dN2Olyr[ilyr]);
      }
      fprintf(files->fp_dN2Olyr, "\n");

      return;
    }
