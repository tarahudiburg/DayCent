
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      potbst.c
**
**  FUNCTION:  void potbst()
**
**  PURPOSE:   Calculate potential transpiration rate.
**             See 2.11 in ELM doc.
**
**  REWRITE:   Melannie Hartman  9/22/93 - 9/28/93
**
**  HISTORY:
**    04/30/92 (SLC)
**    07/01/92 (SLC) Fixed bug.  Equation for bstrate not right.
**    09/01/92 (SLC) Allow transpiration, even if biodead is zero.  This
**                   was in the original model, we will compute shade if biodead
**                   is greater than deadmax.  Otherwise, shadeaf = 1.0
**    11/09/95 (MDH) Added PILPS modifications.
**    11/16/01 (CAK) Soil water potential computed at field capacity 
**    07/25/11 (CAK) Using 0.0 to 1.0 value for swpavg to replace the call to
**                   watrate.
**    01/22/14 (CAK) Using a calculated water effect on transpiration
**                   replacing swpavg with h2ogef.
**
**  INPUTS:
**    biodead - dead above-ground biomass (g/m2)
**    biolive - live above-ground biomass (g/m2)
**    co2val  - CO2 effect on transpiration, added 8/14/98 -mdh
**    fbst    - fraction of water loss from transpiration
**    h2ogef  - 0.0 to 1.0 value returned from trwtavg for wettest soil layer
**              within the plant rooting zone
**    petday  - potential evapotranspiration for the day (cm H2O)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    deadmax - maximum biomass of dead (g/m2), before shade has any affect
**    par1    - intermediate parameter for shade effect calculation
**    par2    - intermediate parameter for shade effect calculation
**    scale1  - scale for shade affect
**    shadeaf - shade affect on transpiration rate
**    trpar1  - input paramter to watrate
**
**  OUTPUTS:
**    bstrate - bare-soil transpiration loss rate (cm/day)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    tanfunc() - tangent function
**    watrate() - compute transpiration rate
**
*****************************************************************************/

#include "soilwater.h"

    void potbst(float *bstrate, float h2ogef, float biolive, float biodead,
                float fbst, float petday, float co2val)
    {
      float scale1 = 0.3f;
/*      float trpar1 = 28.0f; */
      float trpar1 = 12.0f; 
      float deadmax = 150.0f;
      float par1, par2;
      float shadeaf;

      if (biolive <= 0.0) {
        *bstrate=0.0f;
      } else {
        if (biodead >= deadmax) {
          par1 = tanfunc(biolive, 300.0f, 12.0f, 34.0f, 0.002f);
          par2 = tanfunc(biodead, 300.0f, 12.0f, 34.0f, 0.002f);

          /* calculate affect of shading by dead */

          shadeaf = (par1/par2)*(1-scale1) + scale1;
          if (shadeaf > 1) {
            shadeaf = 1.0f;
          }
        } else {
          shadeaf = 1.0f;
        }

/*        *bstrate = watrate(swpavg, petday, trpar1, 0.07f) *
                   shadeaf * petday * fbst; */
/*        *bstrate = watrate(swpavg, petday, trpar1, 0.10f) *
                   shadeaf * petday * fbst * co2val; */
/*        *bstrate = swpavg * shadeaf * petday * fbst * co2val; */
        *bstrate = h2ogef * shadeaf * petday * fbst * co2val;
      }

      return;
    }
