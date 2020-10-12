
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      soiltransp.c
**
**  FUNCTION:  void soiltransp()
**
**  PURPOSE:   Transpire water from soil
**
**  REWRITE:   Melannie Hartman  9/23/93 - 9/28/93
**
**  HISTORY:
**    04/30/92 (SLC)
**    08/13/92 (SLC) Changed call to function which checks lower bound on soil
**                   water content.  Replaced call to "chkzero" with the
**                   function "getdiff".
**    11/27/95 (MDH) Added aet update to this routine.
**
**  INPUTS:
**    bstrate          - transpiration rate (cm H2O/day)
**    layers           - soil water soil layer structure
**    layers->swc[]    - soil water content by layer (cm H2O)
**    layers->swcmin[] - lower bound on soil water content by layer (cm H2O)
**                       swc will not be allowed to drop below this minimum
**    layers->tcoeff[] - transpiration water absoption coefficients by layer
**                       (ND)
**    numlyrs          - number of layers in the soil water model soil profile
**                       that correspond to the top nlaypg Century layers
**
**  GLOBAL VARIABLES:
**    MXSWLYR   - maximum number of soil water model layers (21)
**
**  LOCAL VARIABLES:
**    callname  - call name for subroutine
**    crootwp   - root water potential (not used, perhaps later though?)
**    ftransp   - weighted average of soil water potential
**    ilyr      - current layer in the soil profile
**    sumswpf   - sum of swpfrac over all layers
**    swp[]     - soil water potential at the layer
**    swpfrac[] - ratio of transpiration coefficient to swp per layer
**
**  OUTPUTS:
**    aet           - the actual evapotranspiration so far on this day
**                    (cm H2O)
**    layers->swc[] - soil water content by layer (cm H2O)
**    transp[]      - water transpired from the soil by layer (cm H2O).
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    getdiff()   - check that soil water content does not fall below
**                  the limit, if so adjust the parameter to be subtracted
**    swpotentl() - given its soil water content calculate the soil water
**                  potential of a soil layer
** 
*****************************************************************************/

#include "soilwater.h"

    void soiltransp(float transp[MXSWLYR], int numlyrs, float bstrate,
                    LAYERPAR_SPT layers, float *aet)
    {
      float swp[MXSWLYR];
      float sumswpf;
      float swpfrac[MXSWLYR];
      float ftransp;
      int   ilyr;
      static char *callname = "soiltransp";

/*      printf("Entering soiltransp...\n"); */

      sumswpf = 0.0f;

      for (ilyr=0; ilyr < numlyrs; ilyr++) {
        swp[ilyr] = swpotentl(layers->swc[ilyr], ilyr, layers, callname);
        swpfrac[ilyr] = layers->tcoeff[ilyr] / swp[ilyr];
        sumswpf = sumswpf + swpfrac[ilyr];
      }

      for (ilyr=0; ilyr < numlyrs; ilyr++) {
        ftransp = swpfrac[ilyr] / sumswpf;
/*        crootwp = crootwp+ftransp*swp[ilyr] */

        transp[ilyr] = ftransp * bstrate;
        if (transp[ilyr] < 0.0) {
          printf("transp[%1d] = %5.2f\n", ilyr, transp[ilyr]);
          printf("ftransp = %5.2f\n", ftransp);
          printf("bstrate = %5.2f\n", bstrate);
        }

        /* 8/13/92 (SLC) before subtracting, check if difference is going */
        /* to fall below minimum soil water content, don't allow more than */
        /* this to be subtracted. Adjust transp[ilyr] if necessary. */

        getdiff(&transp[ilyr], layers->swc[ilyr], layers->swcmin[ilyr]);
        if (transp[ilyr] < 0.0) {
/*          printf("After getdiff: transp[%1d] = %5.2f\n", ilyr,
                 transp[ilyr]); */
          transp[ilyr] = 0.0f;
        }

        layers->swc[ilyr] = layers->swc[ilyr] - transp[ilyr];
        *aet += transp[ilyr];

      }

/*      printf("Exitting soiltransp...\n"); */

      return;
    }
