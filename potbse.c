
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      potbse.c
**
**  FUNCTION:  void potbse()
**
**  PURPOSE:   Calculate potential bare soil evaporation rate.
**             See 2.11 in ELM doc.
**
**  REWRITE:   Melannie Hartman  9/22/93 - 9/28/93
**
**  HISTORY:
**    04/30/92 (SLC)
**    08/27/92 (SLC) Put in a check so that bserate cannot become
**                   negative.  If total aboveground biomass (i.e.
**                   litter+biomass) is > ABG_LIMIT, bserate=0.
**    11/09/95 (MDH) Added PILPS modifications (avswp adjustment).
**    05/18/01 (CAK) ecoeff[] no longer being used in this calculation
**    06/12/14       use photosynthetic active carbon rather rather than
**                   aboveground live carbon
**
**  INPUTS:
**    aggreenc          - the amount of photosynthetic active carbon (g/m2)
**    fbse              - fraction of water loss from bare soil evaporation
**    layers            - soil water soil layer structure
**    layers->ecoeff[]  - bare-soil evaporation water absorption coefficients
**                        by layer
**    layers->nelyrs    - number of layers to consider in evaporation
**    layers->sumecoeff - sum of evaporation coefficients
**    layers->swc[]     - the current day's soil water content by layer
**                        (cm H2O)
**    layers->width[]   - the thickness of soil water model layers (cm)
**    petday            - potential evapotranspiration rate for the day
**                        (cm H2O)
**    totagb            - sum of aboveground biomass and litter (g/m2)
**
**   GLOBAL VARIABLES:
**
**   LOCAL VARIABLES:
**    ABG_LIMIT    - upper limit on total aboveground biomass (g/m2)
**                   when biomass values are above this value assume the soil
**                   surface is completely covered with litter
**    avswp        - average soil water potential over all layers
**    bserate_min  - minimum value for bare soil evaporation rate (cm H2O/day)
**    callname     - call name for subroutine
**    evpar1       - input parameter to watrate
**    ilyr         - current layer in the soil profile
**    max_pet_rate - maximum potential evapotranspiration rate (cm H2O/day) 
**    swpsum       - soil water potential summed for all layers considered in
**                   evaporation 
**
**  OUTPUTS:
**    bserate - bare soil evaporation loss rate (cm H2O/day)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    swpotentl() - given its soil water content calculate the soil water
**                  potential of a soil layer
**
*****************************************************************************/

#include "soilwater.h"

/* #define ABG_LIMIT 300 */
#define ABG_LIMIT 1000

    void potbse(float *bserate, float totagb, float fbse, float petday,
                LAYERPAR_SPT layers, float aggreenc)
    {
/*      float evpar1 = 15.0f; */
/*      float evpar1 = 8.0f; */
/*      float avswp; */
/*      float swpsum = 0.0f; */
      float bserate_min = 0.05f;
      float max_pet_rate = 0.70f;
/*      int ilyr; */
/*      static char *callname = "potbse"; */

/*      for(ilyr=0; ilyr < layers->nelyrs; ilyr++) {
        swpsum += layers->width[ilyr]*layers->ecoeff[ilyr] *
                  swpotentl(layers->swc[ilyr], ilyr, layers,callname);
      }

      avswp = swpsum / layers->sumecoeff; */

      /* Change from PILPS runs. Bill Parton. 12/2/94 */
        
/*      avswp = min(avswp, swpotentl(layers->swc[1], 1, layers, callname)); */

      /* 8/27/92 (SLC) if totagb > ABG_LIMIT, assume soil surface is */
      /* completely covered with litter and that bare soil */
      /* evaporation is inhibited. */

      if (totagb >= ABG_LIMIT) {
        *bserate = bserate_min;
      } else {

/*        *bserate = petday * watrate(avswp, petday, evpar1, 0.06)*
                   (1 - (totagb/ABG_LIMIT))*fbse; */
/*        *bserate = petday * watrate(avswp, petday, evpar1, 0.10)*
                   (1 - (totagb/ABG_LIMIT))*fbse; */
/*        *bserate = max(bserate_min, (petday*(1-(totagb/ABG_LIMIT))*fbse)); */
/*        *bserate = max(bserate_min,
                      (petday*(1-((totagb-biolive)/ABG_LIMIT))*fbse)); */
        *bserate = max(bserate_min,
                      (petday*(1-((totagb-(aggreenc*2.5f))/ABG_LIMIT))*fbse));

        if (*bserate > max_pet_rate * petday) {
          *bserate = max_pet_rate * petday;
        }
      }

      return;
    }
