
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      setamov.c
**
**  FUNCTION:  void setamov()
**
**  PURPOSE:   To set amovdly (passed to Century) using the value wfluxout
**             (daily soil water model).
**
**  INPUTS:
**    lbnd[]   - the index of the lower soil water model layer which
**               corresponds to clyr in Century
**    nlayer   - number of layers in Century soil profile
**    numlyrs  - total number of layers in the soil water model soil profile
**    wfluxout - the amount of water moving through the bottom of a soil layer
**               (cm H2O) (positive is downward, negative is upward)
**
**  GLOBAL VARIABLES:
**    CMXLYR   - maximum number of Century soil layers (10)
**    MXSWLYR  - maximum number of soil water model layers (21)
**
**  LOCAL VARIALBES:
**    clyr -  Century soil layer (1..nlayer)
**
**  OUTPUTS:
**    amovdly - the amount of water moving through the bottom of a soil
**              layer (cm H2O) (positive is downward, negative is upward)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "swconst.h"

    void setamov(float amovdly[CMXLYR], int nlayer,
                 float wfluxout[MXSWLYR], int numlyrs, int lbnd[CMXLYR])
    {
      int clyr;

      for(clyr = 0; clyr < nlayer; clyr++) {
        amovdly[clyr] = wfluxout[lbnd[clyr]];
      }

      for(clyr=nlayer; clyr<CMXLYR; clyr++) {
        amovdly[clyr] = 0.0f;
      }

      return;
    }
