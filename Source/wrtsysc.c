
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtsysc.c
**
**  FUNCTION:  void wrtsysc()
**
**  PURPOSE:   Write out the system carbon values. 
**
**  AUTHOR:    Cindy Keough  10/02
** 
**  INPUTS:
**    CO2resp - heterotrophic CO2 respiration (g/m2)
**    curday  - the day of the year (1..366)
**    deadc   - system dead carbon,
**              stdedc + metabc(1) + strucc(1) + wood1c + wood2c + wood3c
**              (g/m2)
**    livec   - system live carbon,
**              aglivc + bglivcj + bglivcm + rleavc + frootcj + frootcm +
**              fbrchc + rlwodc + crootc
**              (g/m2)
**    soilc   - system soil carbon,
**              metabc(2) + strucc(2) + som1c(1) + som1c(2) + som2c(1) +
**              som2c(2) + som3c
**              (g/m2)
**    sysc    - system carbon, livec + deadc + soilc (g/m2)
**    time    - simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files             - structure containing information about output files
**    files->fp_sysc    - file pointer to sysc.out output file
**    files->write_sysc - flag to indicate if sysc.out output file should
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

    void wrtsysc(float *time, int *curday, float *livec, float *deadc,
                 float *soilc, float *sysc, float *CO2resp)
    {
      extern FILES_SPT files;

      if (!files->write_sysc) {
        goto ex;
      }

      fprintf(files->fp_sysc, "%6.2f  %2d  %10.4f  %10.4f  %10.4f  %10.4f",
              *time, *curday, *livec, *deadc, *soilc, *sysc);
      fprintf(files->fp_sysc, " %10.4f\n", *CO2resp);

ex:   return;
    }
