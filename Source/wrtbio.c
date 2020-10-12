
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtbio.c
**
**  FUNCTION:  void wrtbio()
**
**  PURPOSE:   Write out the biomass values. 
**
**  AUTHOR:    Melannie Hartman  6/93
** 
**  INPUTS:
**    aglivc  - above ground live carbon (g/m2)
**    aglivn  - amount of nitrogen in above ground live (g/m2)
**    bglivcj - below ground juvenile fine root live carbon (g/m2)
**    bglivcm - below ground mature fine root live carbon (g/m2)
**    bglivnj - amount of nitrogen in juvenile fine root live (g/m2)
**    bglivnm - amount of nitrogen in mature file root live (g/m2)
**    crootc  - forest coarse root carbon (g/m2)
**    curday  - the day of the year (1..366)
**    fbrchc  - forest fine branch carbon (g/m2)
**    frootcj - forest juvenile fine root carbon (g/m2)
**    frootcm - forest mature fine root carbon (g/m2)
**    h2ogef1 - water effect on crop/grass production
**    h2ogef2 - water effect on forest production
**    rleavc  - forest live leaf carbon (g/m2)
**    rlwodc  - forest large wood carbon (g/m2)
**    time    - simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files            - structure containing information about output files
**    files->fp_bio    - file pointer to bio.out output file
**    files->write_bio - flag to indicate if bio.out output file should be
**                       created, 0 = do not create, 1 = create
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

    void wrtbio(float *time, int *curday, float *aglivc, float *bglivcj,
                float *bglivcm, float *aglivn, float *bglivnj, float *bglivnm,
                float *rleavc, float *frootcj, float *frootcm, float *fbrchc,
                float *rlwodc, float *crootc, float *h2ogef1, float *h2ogef2)
    {
      extern FILES_SPT files;

      if (!files->write_bio) {
        return;
      }

      fprintf(files->fp_bio, "%6.2f  %2d  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *time, *curday, *aglivc, *bglivcj, *bglivcm, *aglivn);
      fprintf(files->fp_bio, "%10.4f  %10.4f  %10.4f  ",
               *bglivnj, *bglivnm, *rleavc);
      fprintf(files->fp_bio, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f", 
                *frootcj, *frootcm, *fbrchc, *rlwodc, *crootc);
      fprintf(files->fp_bio, "%10.6f  %10.6f\n", 
              *h2ogef1, *h2ogef2);

      return;
    }
