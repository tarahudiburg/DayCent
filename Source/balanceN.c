
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      balanceN.c
**
**  FUNCTION:  void bal_npool()
**
**  PURPOSE:   This routine checks to total nitrogen balance between
**             minerl[][N] and ammonium and nitrate; it updates the 
**             distribution of minerl[][N] when ammonium and nitrate 
**             have been updated.
**
**  INPUTS:
**    ammonium  - total ammonium in soil mineral pool (gN/m2)
**    inorglch  - N from organic leaching of stream flow (base flow + storm 
**                flow) (g/m2)
**    minerl[]  - soil minerl pool, 1..nlayer (g/m2)
**                Minerl is a 2-D array in FORTRAN.  FORTRAN stores multi-
**                dimensioned arrays in column-major order.  In this routine
**                it is a 1-D array, where minerl(clyr+1,iel+1) in FORTRAN is
**                equivalent to minerl[(iel)*CMXLYR+(clyr)] in this
**                routine.
**    nitrate[] - total nitrate in soil mineral pool, distributed among soil
**                water model layers (gN/m2)
**    nlayer    - number of layers in Century soil profile
**
**  GLOBAL VARIABLES:
**    CMXLYR     - maximum number of Century soil layers (10)
**    MXSWLYR    - maximum number of soil water model layers (21)
**
**  EXTERNAL VARIABLES:
**    layers          - soil water soil layer structure
**    layers->lbnd[]  - the index of the lower soil water model layer which
**                      corresponds to clyr in Century
**    layers->numlyrs - total number of layers in the soil water model soil
**                      profile
**    layers->ubnd[]  - the index of the upper soil water model layer which 
**                      corresponds to layer clyr in Century
**
**  LOCAL VARIABLES:
**    clyr       - Century soil mineral layer (1..nlayer)
**    iel        - current element, set to 0 for nitrogen
**    ilyr       - current layer in the soil profile
**    minerl_sum - minerl[], summed over soil layers
**    npool_sum  - nitrate[], summed over soil layers
**
**  OUTPUTS:
**    minerl[] - soil minerl pool, 1..nlayer (g/m2)
**               Minerl is a 2-D array in FORTRAN.  FORTRAN stores multi-
**               dimensioned arrays in column-major order.  In this routine
**               it is a 1-D array, where minerl(clyr+1,iel+1) in FORTRAN is
**               equivalent to minerl[(iel)*CMXLYR+(clyr)] in this
**               routine.
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "soilwater.h"

    void bal_npool(int *nlayer, float minerl[], double *ammonium,
                   double nitrate[], double *inorglch)
    {
      int    debug = 0;
      int    clyr, ilyr;
      int    iel;           /* Nitrogen */
      double npool_sum;
      double minerl_sum;

      extern LAYERPAR_SPT layers;

      /* Initializations */
      iel = 0;
      npool_sum = 0.0;
      minerl_sum = 0.0;

      /* Sum ammonium and nitrate pools and compare to sum of minerl pool */

      npool_sum = *ammonium;
      for (ilyr=0; ilyr < MXSWLYR; ilyr ++) {
        npool_sum += nitrate[ilyr];
      }

      for (clyr=0; clyr < CMXLYR; clyr++) {
        minerl_sum += minerl[(iel)*CMXLYR + (clyr)];
      }
      minerl_sum -= *inorglch;

      if (fabs(npool_sum - minerl_sum) > 1.0E-5) {
        if (debug) {
          fprintf(stdout, "Nitrogen soil pools are unbalanced.\n");
          fprintf(stdout, "npool_sum = %13.10f   minerl_sum = %13.10f  ",
                  npool_sum, minerl_sum);
          fprintf(stdout, "npool-minrl = %13.10f\n ", npool_sum - minerl_sum);
          fprintf(stderr, "Nitrogen soil pools are unbalanced.\n");
          fprintf(stderr, "npool_sum = %13.10f   minerl_sum = %13.10f  ",
                  npool_sum, minerl_sum);
          fprintf(stderr, "npool-minrl = %13.10f\n ", npool_sum - minerl_sum);
          for (clyr=0; clyr <= *nlayer; clyr++) {
            fprintf(stdout, "minerl[%1d][N] = %13.10f\n", clyr+1,  
                    minerl[(iel)*CMXLYR + (clyr)]);
            fprintf(stderr, "minerl[%1d][N] = %13.10f\n", clyr+1,  
                    minerl[(iel)*CMXLYR + (clyr)]);
          }
          for (ilyr=0; ilyr <= layers->numlyrs; ilyr ++) {
            fprintf(stdout, "nitrate[%1d] = %13.10f\n", ilyr, nitrate[ilyr]);
            fprintf(stderr, "nitrate[%1d] = %13.10f\n", ilyr, nitrate[ilyr]);
          }
          fprintf(stdout, "ammonium = %13.10f\n", *ammonium);
          fprintf(stderr, "ammonium = %13.10f\n", *ammonium);
        }
/*        if (fabs(npool_sum - minerl_sum) > 1.0E-3) { */
/*        if (fabs(npool_sum - minerl_sum) > 0.01) { */
/*        if (fabs(npool_sum - minerl_sum) > 0.05) { */
        if (fabs(npool_sum - minerl_sum) > 0.1) {
          fprintf(stdout, "ERROR!  N balance boundary check exceeded.\n");
          fprintf(stdout, "npool_sum = %13.10f   minerl_sum = %13.10f  ",
                  npool_sum, minerl_sum);
          fprintf(stdout, "npool-minrl = %13.10f\n ", npool_sum - minerl_sum);
          for (clyr=0; clyr <= *nlayer; clyr++) {
            fprintf(stdout, "minerl[%1d][N] = %13.10f\n", clyr+1,  
                    minerl[(iel)*CMXLYR + (clyr)]);
          }
          for (ilyr=0; ilyr <= layers->numlyrs; ilyr ++) {
            fprintf(stdout, "nitrate[%1d] = %13.10f\n", ilyr, nitrate[ilyr]);
          }
          fprintf(stdout, "ammonium = %13.10f\n", *ammonium);
/*          exit(1); */
        }
      }

      for (clyr=0; clyr < CMXLYR; clyr++) {
        minerl[(iel)*CMXLYR + (clyr)] = 0.0f;
      }

      /* minerl(1,N) = ammonium */
      clyr = 0;
      minerl[(iel)*CMXLYR+(clyr)] = (float)*ammonium;

      for (clyr=0; clyr < *nlayer; clyr++) {
        for(ilyr = layers->ubnd[clyr]; ilyr <= layers->lbnd[clyr]; ilyr++) {
          minerl[(iel)*CMXLYR+(clyr)] += (float)nitrate[ilyr];
        }
      }
      minerl[(iel)*CMXLYR+(*nlayer)] = (float)nitrate[layers->numlyrs];

      return;
    }
