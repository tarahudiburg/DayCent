
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      methane_emission.c
**
**  FUNCTION:  void methane_emission()
**
**  PURPOSE:   Calculate methane emission via plants and bubbles.
**
**  DESCRIPTION:
**    Methane emission is a function of aboveground plant size and root
**    biomass.
**
**  REFERENCE:
**    A semi-empirical model of methane emission from flooded rice paddy
**    soils.  Y. Huang, R.L. Sass, and F.M. Fisher, Jr.  1998.  Global Change
**    Biology 4: 247-268.
**   
**    Modeling methane emission from rice paddies with various agricultrual
**    pratices.  Y. Huang, W. Zhang, X. Zheng, J. Li, and Y. Yu.  Journal of
**    Geophysical Research.  2004.  19: DO8113
**
**  INPUTS:
**    aglivc         - aboveground live carbon (g/m^2)
**    avgst_10cm     - average soil temperature in top 10 cm of soil profile
**                     (degrees C)
**    bglivc         - belowground live carbon (g/m^2)
**    CH4_Ep         - methane emitted via plants (g/m^2)
**    CH4_Ebl        - methane emitted via bubbles (g/m^2)
**    CH4_prod       - methane production (gC/m2/day)
**    tmxbio         - theoretical maximum biomass where no more methane is
**                     emitted by the plant (g biomass/m^2)
**    zero_root_frac - fraction of CH4 emitted via bubbles when there is no
**                     root biomass (0.0 - 1.0)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    agbiomass       - aboveground live biomass (g biomass/m^2)
**    bgbiomass       - belowground live biomass (g biomass/m^2)
**    CH4_prod_remain - methane production minus methane emitted via plants
**                      (g/m^2)
**    Fp              - fraction of produced methane emitted via plants
**                      (0.0 - 0.55)
**    P0              - criterion when bubbles occur (g/m^2)
**   
**  OUTPUTS:
**    CH4_Ep  - methane emitted via plants (g/m^2)
**    CH4_Ebl - methane emitted via bubbles (g/m^2)
**
**  CALLED BY:
**    trace_gas_model()
**
**  CALLS:
**    None
**  
*****************************************************************************/

#include "soilwater.h"
#include "n2o_model.h"
#include <math.h>
#include <stdio.h>

    void methane_emission(float *aglivc, float *tmxbio, float *bglivc,
                          float *avgst_10cm, float *CH4_prod, float *CH4_Ep,
                          float *CH4_Ebl, float zero_root_frac)
    {
      float P0 = 0.002f;
      float agbiomass, bgbiomass;
      float Fp;
      float CH4_prod_remain;

      /* Methane emission via plants */
      if (*aglivc > 0.0) {
        agbiomass = *aglivc * 2.5f;
        Fp = 0.55f * (float)pow(max(1.0 - (agbiomass / *tmxbio),0.0), 0.25);
        *CH4_Ep = Fp * *CH4_prod;
      } else {
        *CH4_Ep = 0.0f;
      }

      /* Methane emission via bubbles */
      CH4_prod_remain = max(*CH4_prod - *CH4_Ep - P0, 0.0f);
      *CH4_Ebl = 0.0f;
      if (*avgst_10cm >= 1.0) {
        if (CH4_prod_remain > 0.0) {
          if (*bglivc > 1.0) {
            bgbiomass = *bglivc * 2.5f;
            *CH4_Ebl = 0.7f * (CH4_prod_remain) *
                       (float)log(*avgst_10cm)/bgbiomass;
          } else {
            /* CH4 emitted via bubbles when there is no root biomass, */
            /* user defined */
            *CH4_Ebl = zero_root_frac * (CH4_prod_remain) *
                       (float)log(*avgst_10cm);
          }
        }
      }
      if (*CH4_Ebl > CH4_prod_remain) {
        *CH4_Ebl = CH4_prod_remain;
      }

      return;
    }