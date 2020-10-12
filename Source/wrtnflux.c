/*****************************************************************************
**
**  FILE:      wrtnflux.c
**
**  FUNCTION:  void wrtnflux()
**
**  PURPOSE:   This function writes daily values for trace gas fluxes.
**
**  INPUTS:
**    curday     - the day of the year (1..366)
**    Dn2flux    - elemental inert nitrogen gas denitrification (gN/ha)
**    Dn2oflux   - nitrous oxide denitrification (gN/ha)
**    nflux_sum  - annual accumulator for nitrous oxide (gN/ha)
**    nflux_sum2 - annual accumulator for nitric oxide (gN/ha)
**    Nn2oflux   - nitrous oxide nitrification (gN/ha)
**    NOflux     - nitric oxide (gN/ha)
**    time       - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_nflux    - file pointer to nflux.out output file
**    files->write_nflux - flag to indicate if nflux.out output file should
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

    void wrtnflux(float *time, int *curday, double *Nn2oflux,
                  double *Dn2oflux, double *Dn2flux, double *NOflux,
                  double *nflux_sum, double *nflux_sum2)
    {

      extern FILES_SPT files;

      if (!files->write_nflux) {
        return;
      }

      fprintf(files->fp_nflux, "%8.2f  %4d", *time, *curday);
      fprintf(files->fp_nflux, "%12.4f  ", *Nn2oflux);
      fprintf(files->fp_nflux, "%12.4f  ", *Dn2oflux);
      fprintf(files->fp_nflux, "%12.4f  ", *Dn2flux);
      fprintf(files->fp_nflux, "%12.4f  ", *NOflux);
      fprintf(files->fp_nflux, "%12.4f  ", *nflux_sum);
      fprintf(files->fp_nflux, "%12.4f  ", *nflux_sum2);
      fprintf(files->fp_nflux, "\n");

      return;
    }
