
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtlivec.c
**
**  FUNCTION:  void wrtlivec()
**
**  PURPOSE:   Write out the live carbon values. 
**
**  AUTHOR:    Cindy Keough  10/02
** 
**  INPUTS:
**    aglivc  - above ground live carbon (g/m2)
**    bglivcj - below ground juvenile fine root live carbon (g/m2)
**    bglivcm - below ground mature fine root live carbon (g/m2)
**    crootc  - coarse root live carbon (g/m2)
**    curday  - the day of the year (1..366)
**    fbrchc  - fine branch live carbon (g/m2)
**    frootcj - forest juvenile fine root live carbon (g/m2)
**    frootcm - forest mature fine root live carbon (g/m2)
**    rleavc  - leaf live carbon (g/m2)
**    rlowdc  - large wood live carbon (g/m2)
**    time    - simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files              - structure containing information about output files
**    files->fp_livec    - file pointer to livec.out output file
**    files->write_livec - flag to indicate if livec.out output file should
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

    void wrtlivec(float *time, int *curday, float *aglivc, float *bglivcj,
                  float *bglivcm, float *rleavc, float *frootcj,
                  float *frootcm, float *fbrchc, float *rlwodc, float *crootc)
    {
      extern FILES_SPT files;

      if (!files->write_livec) {
        goto ex;
      }

      fprintf(files->fp_livec, "%6.2f  %2d  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *time, *curday, *aglivc, *bglivcj, *bglivcm, *rleavc);
      fprintf(files->fp_livec, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n", 
              *frootcj, *frootcm, *fbrchc, *rlwodc, *crootc);

ex:   return;
    }
