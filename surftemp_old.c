
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      surftemp.c
**
**  FUNCTION:  void surftemp()
**
**  PURPOSE:   Compute the minimum, maximum, and average soil surface
**             temperatures
**
**  HISTORY:
**    Rewritten from Fortran to C to enable this routine to be called from
**    the watrflow subroutine.  Moved the code that is calculating the
**    insulation effects of snow on surface soil temperature from the soiltemp
**    routine for inclusion in this subroutine.
**    Cindy Keough
**    01/2010
**
**  INPUTS:
**    aglivb        - above ground live biomass (g/m2)
**    daylength     - amount of daylight hours (1..24) 
**    diurnal_range - the difference between the maximum and minimum surface
**                    soil temperatures (deg C)
**    elitst        - effect of litter on soil temperature relative
**                    to live and standing dead biomass (fix.100)
**    litrcrbn      - amount of carbon in the surface litter (sum of som1c(1),
**                    som2c(1), strucc(1), and metabc(1)
**    pmntmp        - effect of biomass on minimum surface temperature
**                    (fix.100)
**    pmxbio        - maximum biomass for soil temperature calculations
**                    (fix.100)
**    pmxtmp        - effect of biomass on maximum surface temperature
**                    (fix.100)
**    sfclit        - above ground litter biomass (g/m2)
**    SnowFlag      - effect of snow on soil surface temperature calc
**                    1 = insulating effect, 0 = no insulating effect
**    snowpack      - current snowpack (equiv. cm H2O)
**    srfctemp      - average soil surface temperature (degrees C)
**    stamt         - amount of warming applied to soil surface temperature
**                    (degrees C)
**    stdead        - above ground dead biomass (g/m2)
**    ststart       - start time for soil warming simulation
**    stsys         - flag to indicate if soil warming will be simulated
**                      0 = no, 1 = yes
**    time          - current simulation time
**    tmax          - maximum air temperature (deg C)
**    tmin          - minimum air temperature (deg C)
**    tmns          - minimum soil surface temperature (degrees C)
**    tmxs          - maximum soil surface temperature (degrees C)
**    woodb         - wood biomass, fine branch + large wood (g/m^2)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    basetemp  - impact of surface litter carbon on mineral soil surface
**                temperature
**    bio       - leaf biomass (g/m^2)
**    prevstemp - previously calculated value of srfctemp
**    snowmult  - effect of snowdepth on surface temperature.  It is
**                smallest with deep snow and multiples average air
**                temperature.  The more snow, the closer soil surface
**                temperature is to freezing (frac).
**    tmns_leaf - minimum temperature with leaf shading (degrees C)
**    tmns_mlt  - fraction used to compute weighted mimimum surface soil
**                temperature (0.5 - 0.95)
**    tmns_wood - minimum temperature with wood shading (degrees C)
**    tmxs_leaf - maximum temperature with leaf shading (degrees C)
**    tmxs_mlt  - fraction used to compute weighted maximum surface soil
**                temperature (1.0 - tmns_mlt)
**    tmxs_wood - maximum temperature with wood shading (degrees C)
**    woodbio   - wood biomass (g/m^2)
**
**  OUTPUTS:
**    diurnal_range - the difference between the maximum and minimum surface
**                    soil temperatures (deg C)
**    srfctemp      - average soil surface temperature (degrees C)
**    tmns          - minimum soil surface temperature (degrees C)
**    tmxs          - maximum soil surface temperature (degrees C)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <math.h>
#include <stdlib.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

    void surftemp(float elitst, float pmxtmp, float pmntmp, float pmxbio,
                  float tmax, float tmin, float *tmxs, float *tmns,
                  float *srfctemp, float daylength, float aglivb,
                  float stdead, float sfclit, float woodb, float stsys,
                  float ststart, float stamt, float time, float snowpack,
                  float *diurnal_range, float litrcrbn, int SnowFlag)
    {

      float bio, woodbio;
      float tmns_mlt, tmxs_mlt;
      float tmns_leaf, tmns_wood, tmxs_leaf, tmxs_wood;
      float snowmult;
      float basetemp;
      static float prevstemp = 100.0f;

      /* Leaf biomass */
      bio = aglivb + stdead + elitst * sfclit;
      bio = min(bio, pmxbio);

      /* Wood biomass */
      woodbio = min(woodb, 5000.0f);

      /* Maximum temperature with leaf shading */
      tmxs_leaf = tmax+(25.4f/(1.0f+18.0f*(float)exp(-0.20*tmax))) *
                  (float)(exp(pmxtmp*bio)-0.13);
      /* Minimum temperature with leaf shading */
      tmns_leaf = tmin+pmntmp*bio-1.78f;

      /* Maximum temperature with wood shading */
      tmxs_wood = tmax+(25.4f/(1.0f+18.0f*(float)exp(-0.20*tmax)))*
                  (float)(exp(pmxtmp*0.1f*woodbio)-0.13f);
      /* Minimum temperature with wood shading */
      tmns_wood = tmin+pmntmp*0.1f*woodbio-1.78f;

      *tmxs = min(tmxs_leaf, tmxs_wood);
      *tmns = max(tmns_leaf, tmns_wood);

      /* Average surface temperature */
      /* Use day length to compute the values for the multipliers on */
      /* minimum and maximum soil surface temperature, cak - 04/25/03 */
      /* Weigh surface temperature more towards tmns in the winter between */
      /* the fall and spring equinox when nights are longer than days. */
      if (daylength < 12.0f) {
        tmns_mlt = ((12.0f - daylength) * 3.0f + 12.0f) / 24.0f;
      } else {
        tmns_mlt = ((12.0f - daylength) * 1.2f + 12.0f) / 24.0f;
      }
      tmns_mlt = min(0.95f, tmns_mlt);
      tmns_mlt = max(0.05f, tmns_mlt);
      tmxs_mlt = 1.0f - tmns_mlt;

      /* Calculate the impact of surface litter carbon on mineral soil */
      /* surface temperature */
      basetemp = min((2.0f - (-2.0f)) / (3000.0f - 0.0f) *
                     (litrcrbn - 3000.0f) + 2.0f, 0.0f);

      /* Compute insulation effects of snow on surface soil temperature */
      /* 11/30/95 (Bill Parton). */
      if ((snowpack <= 0.00000001) || (SnowFlag == 0)) {
        *srfctemp = tmxs_mlt* *tmxs + tmns_mlt * *tmns;
        *diurnal_range = *tmxs - *tmns;
      } else {
        if ((tmin + tmax) / 2.0 >= 0.0) {
          /* if there is snow, and average air temperature gets above */
          /* freezing, average soil surface temperature stays at freezing */
          *srfctemp = basetemp;
          *diurnal_range = 0.3f * (tmax - tmin) / 2.0f;
          if (*diurnal_range / 2.0 + *srfctemp > 0.0) {
            *diurnal_range = basetemp * *srfctemp;
          }
        } else {
          /* average air temperature is below freezing */
          snowmult = -0.15f * snowpack + 1.0f;
          if (snowmult < 0.0) {
            snowmult = 0.0f;
          }
          *srfctemp = basetemp + (0.3f * (tmin + tmax) / 2.0f) * snowmult;
          *diurnal_range = 0.3f * (tmax - tmin) * snowmult;
          if (*diurnal_range / 2.0 + *srfctemp > 0) {
            *diurnal_range = basetemp * *srfctemp;
          }
        }
        /* If the soil has frozen prior to snow accumulating on the site */
        /* the insulating effect of the snow will hold the low frozen soil */
        /* temperature until the snow melts, cak 09/29/2011 */
        if (prevstemp < *srfctemp) {
          *srfctemp = prevstemp;
        }
        prevstemp = *srfctemp;
      }

      /* Added code to allow warming of the soil surface temperature, */
      /* cak - 07/02/03 */
      if (stsys > 0 && time >= ststart) {
        *srfctemp = *srfctemp + stamt;
      }

      return;
    }
