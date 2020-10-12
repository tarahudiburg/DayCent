/*****************************************************************************
**
**  FILE:      wrtdN2lyr.c
**
**  FUNCTION:  void wrtdN2lyr()
**
**  PURPOSE:   This function writes out daily N2 flux from denitrification
**             by layer (gN/m^2).
**
**  INPUTS:
**    dN2lyr[] - N2 flux by layers (gN/m^2)
**    jday     - current day of the year (1..366)
**    time     - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_dN2lyr    - file pointer to dN2lyr.out output file
**    files->write_dN2lyr - flag to indicate if dN2lyr.out output file should
**                          be created, 0 = do not create, 1 = create
**    layers->numlyrs     - total number of layers in the soil water model
**                          soil profile
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

    void wrtdn2lyr(float *time, int *jday, double dN2lyr[])
    {

      int    ilyr;
      extern LAYERPAR_SPT layers;
      extern FILES_SPT files;

      if (!files->write_dN2lyr) {
        return;
      }

      fprintf(files->fp_dN2lyr, "%8.2f  %4d  ", *time, *jday);

      for (ilyr=0; ilyr < layers->numlyrs; ilyr ++) {

        fprintf(files->fp_dN2lyr, "%12.6f  ", dN2lyr[ilyr]);
      }
      fprintf(files->fp_dN2lyr, "\n");

      return;
    }
