
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      methane_production.c
**
**  FUNCTION:  void methane_production()
**
**  PURPOSE:   Calculate methane production.
**
**  DESCRIPTION:
**    Methane production is a function of soil texture, soil temperature,
**    litter and SOM decomposition rates, aboveground live carbon, and soil
**    redox potential.
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
**    Modeling methane emissions from rice paddies.  M. Cao, J.B. Dent, and
**    O.W. Heal.  1995.  Global Biogeochemical Cycles.  Vol 9, No. 2, pages
**    183-195.
**    
**  INPUTS:
**    Com        - sum of CO2 losses from heterotrophic decomposition of
**                 metabc(1), metabc(2), strucc(1), strucc(2), som1c(1),
**                 som1c(2), som2c(1), som2c(2), and som3c (g/m^2)
**    Eh         - effect of water management on soil redox potential (Eh)
**                 (-250.0 mv - 300.0 mv)
**    prev_bgprd - previous day's fine root carbon production (g/m^2)
**    watertable - 1 = simulate water table (flooding stage),
**                 0 = no water table (drainage stage)
**    watrflag   - 1 = water added to the system automatically to bring the
**                     soil water content to saturation
**                 0 = rain and irrigation events used to add water to the
**                     soil potentially bringing it to saturation
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    layers                  - soil water soil layer structure
**    layers->sandfrac[]      - sand fraction in soil layer, 0.0 - 1.0
**    layers->width[]         - the thickness of soil water model layers (cm)
**    sitepar                 - site specific parameters structure for soil
**                              water model
**    sitepar->Aeh            - differential coefficient for calculation of
**                              soil Eh
**    sitepar->Beh_drain      - upper-limit value of Eh during drainage course
**                              (mv)
**    sitepar->Beh_flood      - low-limit value for Eh during flooding course
**                              (mv)
**    sitepar->C6H12O6_to_CH4 - reaction of anaerobic carbohydrate
**                              fermentation with methanogenesis (mole weight)
**    sitepar->Deh            - differential coefficient for calculation of
**                              soil Eh
**    soil                    - soil temperature structure
**    soil->soiltavg[]        - average soil temperature of layer (degrees C)
**
**  LOCAL VARIABLES:
**    aglivc_biomass - aboveground live carbon converted from grams of carbon
**                     to grams of biomass
**    avgsand_10cm   - average sand content in top 3 soil layers as defined by
**                     the soils.in file (percentage 1.0 - 100.0) 
**    Q10            - temperature coefficient for methane emission
**    VI             - variety index, relative difference in plant methane
**                     production
**
**  OUTPUTS:
**    avgst_10cm - average soil temperature in top 10 cm of soil profile
**                 (degrees C)
**    CH4_prod   - methane production (gC/m2/day)
**    Cr         - carbohydrates derived from rice plants (g/m^2)
**    Eh         - effect of water management on soil redox potential (Eh)
**                 (-250.0 mv - 300.0 mv)
**    Feh        - reduction factor of effect of soil redox potential (Eh) on
**                 methane production (0.0 - 1.0)
**    SI         - soil texture index for methane production (0.0 - 1.0)
**    TI         - soil temperature index for methane production (0.0 - 1.0)
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

    void methane_production(float *prev_bgprd, float *Com, float *avgst_10cm,
                            float *TI, float *SI, float *Cr, float *Eh,
                            float *Feh, float *CH4_prod, int *watertable,
                            int *watrflag)
    {
      float VI = 1.0f;
      float Q10 = 3.0f;
      float avgsand_10cm;

      extern LAYERPAR_SPT layers;
      extern SITEPAR_SPT sitepar;
      extern SOIL_SPT soil;

      /* Average sand content in top 3 soils.in soil layers (10 cm) */
      avgsand_10cm = (layers->sandfrac[0] * layers->width[0] + 
                      layers->sandfrac[1] * layers->width[1] + 
                      layers->sandfrac[2] * layers->width[2]) /
                     (layers->width[0] + layers->width[1] + layers->width[2]);
      avgsand_10cm *= 100.0f;
      /* Soil texture index for methane production */
      *SI = 0.325f + 0.0225f * avgsand_10cm;

      /* Average soil temperature in top 3 soils.in soil layers (10 cm) */
      *avgst_10cm = (soil->soiltavg[0] * layers->width[0] + 
                     soil->soiltavg[1] * layers->width[1] + 
                     soil->soiltavg[2] * layers->width[2]) /
                    (layers->width[0] + layers->width[1] + layers->width[2]);
      /* Soil temperature index for methane production */
      *TI = min(*avgst_10cm, 30);
      *TI = (float)pow(Q10, ((*TI - 30.0) / 10.0));

      /* Use the previous day's fine root production to calculate */
      /* carbohydrates derived from rice plants */
      *Cr = sitepar->frac_to_exudates * (*prev_bgprd);

      /* Soil redox potential */
      if (*watertable == 1) {
        /* Flooding stage */
        if (*watrflag == 1) {
          /* Water added automatically to keep soil at saturation */
          *Eh = *Eh - sitepar->Deh * (sitepar->Aeh + min(1.0f, *Com)) *
                (*Eh - sitepar->Beh_flood);
        } else {
          /* Water added via rain and irrigation events */
          *Eh = -20.0f;
        }
      } else {
        /* Drainage stage */
        *Eh = *Eh - sitepar->Deh * (sitepar->Aeh + 0.7f) *
              (*Eh - sitepar->Beh_drain);
      }

      /* Reduction factor of effect of soil redox potential on methane */
      /* production */
      if (*Eh < -150.0) {
        *Feh = (float)exp(-1.7 * ((150.0 + (-150.0)) / 150.0));
      } else {
        *Feh = (float)exp(-1.7 * ((150.0 + *Eh) / 150.0));
      }

      /* Methane production */
      *CH4_prod = sitepar->C6H12O6_to_CH4 * (*Feh) * (*Com + (*TI * (*Cr)));

      return;
    }