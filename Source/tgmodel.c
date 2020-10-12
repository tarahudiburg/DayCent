
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      tgmodel.c
**
**  FUNCTION:  void trace_gas_model()
**
**  PURPOSE:   Trace Gas Model - calculates daily N2, N2O, NO fluxes from the
**             soil.
**
**  HISTORY:   New module routine, pulled from dailymoist. -mdh 6/20/00
**
**  INPUTS:
**    afiel      - field capacity in the top 10 cm (vol. frac)
**    aglivc     - aboveground live carbon (g/m^2)
**    ammonium   - ammonium concentration (gN/m^2)
**    avgst_10cm - average soil temperature in top 10 cm of soil profile
**                 (degrees C)
**    avgwfps    - avg wfps in top 3 soil layers (~top 10 cm) 0-1
**    basef      - fraction of base flow
**    bglivc     - belowground live carbon (g/m^2)
**    btolai     - biomass to LAI conversion factor for trees
**    bulkd      - bulk density (g/cm^3)
**    CH4_Ep     - methane emitted via plants (gCH4-C/m^2/day)
**    CH4_Ebl    - methane emitted via bubbles (gCH4-C/m^2/day)
**    CH4_oxid   - methane oxidation (gCH4-C/ha/day)
**    CH4_prod   - methane production (gCH4-C/m^2/day)
**    clay       - fraction of clay (0.0 - 1.0)
**    co2PPM[]   - CO2 concentration by layer (ppm)
**    Com        - sum of CO2 losses from heterotrophic decomposition of
**                 metabc(1), metabc(2), strucc(1), strucc(2), som1c(1),
**                 som1c(2), som2c(1), som2c(2), and som3c (g/m^2)
**    Cr         - carbohydrates derived from rice plants (g/m^2)
**    critflow   - amount of flow between lyrs to leach (cm H2O)
**    crpstore   - crop storage, NO absorped by canopy is added to this pool
**    Dn2flux    - denitrification N2 flux (gN/m^2/day)
**    Dn2oflux   - denitrification N2O flux (gN/m^2/day)
**    dN2lyr[]   - N2 flux from Denitrification by layer (gN/m2/day)
**    dN2Olyr[]  - N2O flux from Denitrification by layer (gN/m2/day)
**    efscltef   - effective soil cultivation effect.
**                 This is the C flow weighted clteff for correction of
**                 the CO2 concentration when used as a surrogate for O2
**                 availability. This reduces the CO2 produced by tillage
**                 reflecting enhanced gas permeability from the tillage.
**                 Modification recommended S Del Grosso KLK 14Aug2012
**    Eh         - effect of water management on soil redox potential (Eh)
**                 (-250.0 mv - 300.0 mv)
**    Feh        - reduction factor of effect of soil redox potential (Eh) on
**                 methane production (0.0 - 1.0)
**    forstore   - forest storage, NO absorped by canopy is added to this pool
**    frlechd[]  - leaching fractions (?)
**    inorglch   - inorganic N leaching (gN/m^2) (?)
**    isdecid    - flag, set to 1 for deciduous forest, 0 otherwise
**    isagri     - flag, set to 1 for agricultural system or grassland that
**                 has been fertilized or cultivated, 0 otherwise
**    jday       - current day of year (1-366)
**    maxt       - max of montly max temperatures (deg C)
**    newCO2     - amount of soil respiration (gC/m^2/day)
**    newminrl   - net mineralization (gN/m^2/day)
**    nit_amt    - gross nitrification (gN/m^2/day)
**    nitrate[]  - nitrate concentration by lyr (gN/m^2)
**    Nn2oflux   - nitrification N2O flux (gN/m^2/day)
**    NOflux     - NO flux (gN/m^2/day)
**    nreduce    - reduction factor on nitrification rates due to
**                 nitrification inhibitors added with fertilizer (0-1)
**    pHscale    - optional scalar on soil pH
**    ppt        - daily precip (cm)
**    prev_bgprd - previous day's fine root carbon production (g/m^2)
**    rleavc     - leaf carbon (g/m^2)
**    sand       - fraction of sand (0.0 - 1.0)
**    SI         - soil texture index for methane production (0.0 - 1.0)
**    silt       - fraction of silt (0.0 - 1.0)
**    snow       - snow cover (cm SWE)
**    stormf     - fraction of storm flow
**    stream[]   - stream flow (?)
**    texture    - soil texture index (see swconst.h)
**    TI         - soil temperature index for methane production (0.0 - 1.0)
**    time       - current time in years
**    tmxbio     - theoretical maximum biomass where no more methane is
**                 emitted by the plant (g biomass/m^2)
**    watertable - flag, 1 = simulate water table, 0 = no water table
**    watrflag   - 1 = water added to the system automatically to bring the
**                     soil water content to saturation
**                 0 = rain and irrigation events used to add water to the
**                     soil potentially bringing it to saturation
**    wfluxout[] - amount of flow between lyrs (cm H2O)
**
**  GLOBAL VARIABLES:
**    MXSWLYR - maximum number of soil water model layers (21)
**    PI      - pi (3.14159265)
**
**  EXTERNAL VARIABLES:
**    sitepar                - site parameter structure
**    sitepar->NO3_N2O_slope - slope, controls how fast the percentage of NO3
**                             from mineralization that is lost as N2O changes
**                             as a function of water filled pore space
**    sitepar->NO3_N2O_step  - step between minimum and maximum percentage of
**                             NO3 from mineralization that is lost as N2O
**    sitepar->NO3_N2O_x     - x-inflection, water filled pore space where 1/2
**                             of maximum percentage of NO3 from
**                             mineralization is lost as N2O
**    sitepar->NO3_N2O_y     - y-inflection, percentage of NO3 from
**                             mineralization that is lost as N2O at
**                             x-inflection water filled pore space
**    sitepar->netmn_to_no3  - fraction of new net mineralization that goes to
**                             NO3 (0.0-1.0)
**
**  LOCAL VARIABLES:
**    A[]              - parameters to Parton-Innis functions
**    avg_fc           - average field capacity in the top 3 soil layers
**                       (volumetric)
**    avg_vswc         - average soil water content in top 3 soil layers 
**                       (volumetric)
**    canopy_reduction - reduction factor for NO absorbed by canopy (0-1)
**    dDO              - normalized diffusivity in aggregate soil media (0-1)
**    debug            - flag to set debugging mode, 0 = off, 1 = on
**    grass_lai        - amount of LAI in grass/crop leaf canopy (m^2/m^2)
**    ilyr             - current layer in the soil profile
**    is_saturated     - flag, 0 = soil is not saturated,
**                       1 = soil is saturated
**    krainNO          - increase of NO due to moisture and rain >=1.0
**    newNH4           - new NH4 produced (gN/m^2/day)
**    newNO3           - new NO3 produced (gN/m^2/day)
**    NH4_to_NO        - ammonium converted to NO to reach potential NO flux
**                       (gN/m^2/day)
**    nh4_to_no3       - amount of NH4 converted to NO3 by nitrification
**                       (gN/m^2/day)
**    NO_N2O_ratio     - NO/N2O ratio <= 1.0
**    NOabsorp         - NO absorbed by canopy (gN/m^2/day)
**    NOsoil           - total NO flux emitted by the soil (gN/m^2/day)
**    npool_sum        - sum of N in nitrate and ammonium pools (gN/m^2)
**    potential_NOflux - maximum possible NO flux based on NO/N2O and krain
**                       (gN/m^2/day)
**    stdbulkd         - standard bulk density based on soil texture (g/cm^3)
**    stdfieldc        - standard field capacity based on soil texture
**                       (vol.frac)
**    total_lai        - total LAI at site (m^2/m^2)
**    tree_lai         - amount of LAI in tree leaf canopy (m^2/m^2)
**    turnovfrac       - fraction of new NO3 that goes to N2O
**
**  OUTPUTS:
**    ammonium   - ammonium concentration (gN/m^2)
**    avgst_10cm - average soil temperature in top 10 cm of soil profile
**                 (degrees C)
**    CH4_Ep     - methane emitted via plants (gCH4-C/m^2/day)
**    CH4_Ebl    - methane emitted via bubbles (gCH4-C/m^2/day)
**    CH4_oxid   - methane oxidation (gCH4-C/ha/day)
**    CH4_prod   - methane production (gCH4-C/m^2/day)
**    co2PPM[]   - CO2 concentration by layer (ppm)
**    Cr         - carbohydrates derived from rice plants (g/m^2)
**    crpstore   - crop storage, NO absorped by canopy is added to this pool
**    dN2lyr[]   - N2 flux from Denitrification by layer (gN/m2/day)
**    dN2Olyr[]  - N2O flux from Denitrification by layer (gN/m2/day)
**    Dn2flux    - denitrification N2 flux (gN/m^2/day)
**    Dn2oflux   - denitrification N2O flux (gN/m^2/day)
**    Eh         - effect of water management on soil redox potential (Eh)
**                 (-250.0 mv - 300.0 mv)
**    Feh        - reduction factor of effect of soil redox potential (Eh) on
**                 methane production (0.0 - 1.0)
**    forstore   - forest storage, NO absorped by canopy is added to this pool
**    inorglch   - inorganic N leaching (gN/m^2) (?)
**    nit_amt    - gross nitrification (gN/m^2/day)
**    nitrate[]  - nitrate concentration by lyr (gN/m^2)
**    Nn2oflux   - nitrification N2O flux (gN/m^2/day)
**    NOflux     - NO flux (gN/m^2/day)
**    SI         - soil texture index for methane production (0.0 - 1.0)
**    stream[]   - stream flow (?)
**    TI         - soil temperature index for methane production (0.0 - 1.0)
**
**  CALLED BY:
**    dailymoist() - FORTRAN routine
**
**  CALLS:
**    denitrify()   - calculate N2O flux and N2:N2O ratio due to
**                    denitrification
**    diffusiv()    - estimate normalized diffusivity in soils
**    getsoilprop() - determine soil texture class based on percent of sand,
**                    silt, clay
**    nitrify()     - nitrification, produces N2O gas
**    nox_pulse()   - increase NO due to moisture and rain
**
*****************************************************************************/

#include "soilwater.h"
#include "n2o_model.h"
#include "swconst.h"
#include <math.h>

    void trace_gas_model(double *newminrl, double *ammonium, double nitrate[],
                         int *texture, float *sand, float *silt, float *clay,
                         float *afiel, float *bulkd, float *maxt, float *ppt,
                         float *snow, float *avgwfps, float *stormf,
                         float *basef, float frlechd[], float stream[],
                         double *inorglch, float *critflow, float wfluxout[],
                         float *newCO2, float co2PPM[], float *efscltef,
                         float *time, double *NOflux, double *Nn2oflux,
                         double *Dn2oflux, double *Dn2flux, double *CH4_oxid,
                         int *isdecid, int *isagri, float *aglivc,
                         float *rleavc, float *btolai, float *crpstore,
                         float *forstore, double *nit_amt, float *nreduce,
                         int *jday, float *pHscale, double dN2lyr[],
                         double dN2Olyr[], float *prev_bgprd, float *Com,
                         float *avgst_10cm, float *TI, float *SI, float *Cr,
                         float *Eh, float *Feh, float *CH4_prod,
                         float *CH4_Ep, float *CH4_Ebl, int *watertable,
                         int *watrflag, float *bglivc, float *tmxbio)
    {

      /* Local Variables */

      int    debug = 0;
      int    ilyr;
      int    is_saturated;
      double turnovfrac;
      double newNH4;
      double newNO3;
      double nh4_to_no3;
      double krainNO;
      double potential_NOflux;
      double dDO;
      float  stdbulkd;
      float  stdfieldc;
      double NO_N2O_ratio;
      double NH4_to_NO;
      double npool_sum;
      float canopy_reduction;
      float grass_lai;
      double NOabsorp;
      float total_lai;
      float tree_lai;
      float  A[4];
      float avg_fc, avg_vswc;

      extern LAYERPAR_SPT layers;
      extern SITEPAR_SPT sitepar;

      *Nn2oflux = 0.0;
      *NOflux = 0.0;
      *Dn2oflux = 0.0;
      *Dn2flux = 0.0;

      /* Compute fraction of new mineralization that is converted to NH4 */
      /* and NO3 */

      if (debug) {
        printf("newminrl = %6.4f\n", *newminrl);
      }

      if (*newminrl <= 0.0) {

        /* Immobilization */
        /* Distribute N loss proportionally between ammonium and nitrate   */
        /* layers.  There is no check that these N pools won't go negative */
        /* once immobilization is accounted for.  It is assumed that the   */
        /* immobilization calculation by the decomp model is moderated by  */
        /* the supply of minerl N.                                         */

        npool_sum = (*ammonium > 0.0) ? *ammonium : 0.0;
        for (ilyr=0; ilyr < MXSWLYR; ilyr ++) {
          npool_sum += (nitrate[ilyr] > 0.0) ? nitrate[ilyr] : 0.0;
        }
        if (*ammonium > 0.0) {
          *ammonium += *newminrl * (*ammonium / npool_sum);
        }
        for (ilyr=0; ilyr < MXSWLYR; ilyr ++) {
          if (nitrate[ilyr] > 0.0) {
            nitrate[ilyr] += *newminrl * (nitrate[ilyr] / npool_sum);
          }
        }
        newNH4 = 0.0;
        newNO3 = 0.0;
      } else {
        /* Mineralization */
        newNH4 = *newminrl * (1.0 - sitepar->netmn_to_no3);
        newNO3 = *newminrl * sitepar->netmn_to_no3;
      }

      if (debug) {
        printf("newNH4 = %6.4f\n", newNH4);
        printf("newNO3 = %6.4f\n", newNO3);
      }

      *ammonium += newNH4;

      /* Compute the amount of NH4 that is converted to NO3 due to */
      /* nitrification */

      nitrify(ammonium, &nh4_to_no3, maxt, nreduce, pHscale);
      *nit_amt = nh4_to_no3;

      if (debug) {
        printf("texture = %1d\n", *texture);
        printf("nh4_to_no3 = %6.4f\n", nh4_to_no3);
        printf("maxt = %6.4f\n", *maxt);
      }

      /* Compute fraction of new NO3 that is converted to N2O and NO */

      krainNO = nox_pulse(ppt, snow);

      getsoilprop(sand, silt, clay, &stdbulkd, &stdfieldc, texture);

      /* Use standard field capacity and bulk density according */
      /* to the soil class in the texture triangle -mdh 10/26/99 */
/*      dDO = diffusiv(afiel(1), bulkd, *avgwfps) */
      /* No, change back to soils.in field capacity and bulk density. */
      /* -mdh 6/20/00 */
/*      dDO = diffusiv(&stdfieldc, &stdbulkd, avgwfps); */
      dDO = diffusiv(afiel, bulkd, avgwfps);

      newNO3 += nh4_to_no3;

      if (newNO3 > 1.0E-30) {
        /* Use an arctangent function to control the percentage of NO3 that */
        /* is lost as N2O, cak - 08/28/2014 */
        A[0] = sitepar->NO3_N2O_x;
        A[1] = sitepar->NO3_N2O_y;
        A[2] = sitepar->NO3_N2O_step;
        A[3] = sitepar->NO3_N2O_slope;
        turnovfrac = f_arctangent(*avgwfps, A);
        /* We need a fraction here, rather than a percentage */
        turnovfrac /= 100.0;
        *Nn2oflux = newNO3 * turnovfrac;
        newNO3 -= *Nn2oflux;

        /* Another update to NO flux calculation -mdh 10/26/99 */

/*        NO_N2O_ratio = 15.23 + (35.45*atan(0.676*PI*(10*dDO-1.86)))/PI; */
        NO_N2O_ratio = 8.0 + (18.0*atan(0.75*PI*(10*dDO-1.86)))/PI;
        /* If this is an agricultural system adjust the NO to N2O ratio */
        /* cak - 01/28/03 */
        if (*isagri) {
/*          NO_N2O_ratio *= 0.2; */
          NO_N2O_ratio *= 0.5;
        }
        potential_NOflux = NO_N2O_ratio * *Nn2oflux * krainNO;

        if (potential_NOflux <= newNO3) {
          *NOflux = potential_NOflux;
          newNO3 -= *NOflux;
        } else {
          /* take N out of ammonimum to get max NOflux possible */
          NH4_to_NO = min(*ammonium, (potential_NOflux-newNO3));
          *NOflux = newNO3 + NH4_to_NO;
          *ammonium -= NH4_to_NO;
          newNO3 = 0;
        }

        if (*NOflux < 1.0E-30) {
          *NOflux = 0.0;
        }

      } else {
        NO_N2O_ratio = 0.0;
      }

      /* If the volumetric soil water content in the top 3 soil layers is */
      /* greater than field capacity consider the soil to be saturated */
      avg_vswc = (float)(layers->swc[0] / layers->width[0] +
                         layers->swc[1] / layers->width[1] +
                         layers->swc[2] / layers->width[2]) / 3.0f;
      avg_fc = (layers->fieldc[0] * layers->width[0]  +
                layers->fieldc[1] * layers->width[1] +
                layers->fieldc[2] * layers->width[2]) /
               (layers->width[0] + layers->width[1] + layers->width[2]);
      if (avg_vswc <= avg_fc) {
        is_saturated = 0;
      } else {
        is_saturated = 1;
      }

      /* Compute the N2O flux (Dn2oflux) and N2 flux (Dn2flux) due to */
      /* denitrification */
      denitrify(newCO2, &newNO3, nitrate, wfluxout, critflow, frlechd,
                stream, basef, stormf, inorglch, Dn2oflux, Dn2flux,
                stdfieldc, stdbulkd, efscltef, co2PPM, dN2lyr, dN2Olyr, jday,
                is_saturated);

      /* Now compute NOflux from denitrification (new calculation */
      /* -mdh 6/1/00 */
/*      potential_NOflux = NO_N2O_ratio * *Dn2oflux * krainNO; */
      /* For denitrification, krainNO is >= 1.0 -mdh 6/22/00 */

      potential_NOflux = NO_N2O_ratio * *Dn2oflux * min(1.0, krainNO);

      if (potential_NOflux <= *ammonium) {
        /* Take all N out of ammonimum pool */
        *NOflux += potential_NOflux;
        *ammonium -= potential_NOflux;
      } else {
        /* Take N out of available ammonium, then convert some Dn2oflux to */
        /* NOflux */
        *NOflux += *ammonium;
        potential_NOflux -= *ammonium;
        *ammonium = 0.0;
        if (potential_NOflux <= *Dn2oflux) {
          *NOflux += potential_NOflux;
          *Dn2oflux -= potential_NOflux;
        }
      }

      /* Compute the amount of the soil NO flux that is absorped by the */
      /* canopy, cak - 09/23/03 */
      grass_lai = (*aglivc * 2.5f) / CONVLAI;
      tree_lai = (*rleavc * 2.5f) * *btolai;
      total_lai = grass_lai + tree_lai;
      if (total_lai > 0.0) {
        /* canopy_reduction appears to be the reabsorbed fraction.
           This equation is a parabola with the minimum about 8.28. It reduces
           absorption for LAI above that so limit LAI to 8.0. The problem seems
           to be the simple biomass to LAI ratio. 200 bu/acre corn is
           aglivC > 1100 and an LAI > 34. Obviously its not all leaves. */
          canopy_reduction = (total_lai > 8.0)? 0.4428f:
                (float)(0.0077 * pow(total_lai,2) + -0.13 * total_lai + 0.99);

        /* We need to retain the soil flux value
        NOsoil = *NOflux;
         *NOflux *= canopy_reduction;  NOabsorp = NOsoil - *NOflux; */
        /*   KLK 2Nov12  algebraically reduce the above steps to: */
        NOabsorp = *NOflux * (1.0f - canopy_reduction);

        /* NO absorped by canopy goes to crop storage and forest storage */
        *crpstore += (float)(NOabsorp * (grass_lai / total_lai));
        *forstore += (float)(NOabsorp * (tree_lai /  total_lai));

        *NOflux -= NOabsorp; /* reduce the NOflux by absorption*/
      }

      if (*NOflux < 1.0E-30) {
        *NOflux = 0.0;
      }
      if (*Nn2oflux < 1.0E-30) {
        *Nn2oflux = 0.0;
      }
      if (*Dn2oflux < 1.0E-30) {
        *Dn2oflux = 0.0;
      }
      if (*Dn2flux < 1.0E-30) {
        *Dn2flux = 0.0;
      }

      *CH4_oxid = 0.0f;
      *CH4_prod = 0.0f;
      *CH4_Ep = 0.0f;
      *CH4_Ebl = 0.0f;
/*      if (is_saturated == 0) { */
        /* Calculate methane oxidation */
        methane_oxidation(CH4_oxid, isdecid, isagri);
/*      } else { */
        /* Calculate methane production */
        *Com = (1.0f / sitepar->C6H12O6_to_CH4) * sitepar->CO2_to_CH4 * (*Com);
        methane_production(prev_bgprd, Com, avgst_10cm, TI, SI, Cr, Eh, Feh,
                           CH4_prod, watertable, watrflag);
        /* Calculate methane emission via plants and bubbles */
        if (*CH4_prod > 0.0) {
          methane_emission(aglivc, tmxbio, bglivc, avgst_10cm, CH4_prod, CH4_Ep,
                           CH4_Ebl, sitepar->zero_root_frac);
        }
/*      } */

      return;
    }
