
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
**                    -1 = minimum temperature calculated using snow temp buffering at
**                         snow base during run.
**                         for comparison ONLY; run level retention is a bug.
**                     0 = no snow insulation.
**                     1 = snow temp buffering from Parton paper; no
**                         tave > 0 basetemp = -2 + littrC/750
**                         tave < 0 basetemp + 0.13 * tave * (1-3*snowpack/20)
**                     2 = snow insulation around basetemp; may exhibit freeze thaw
**                     3 = low temperature snow insulation; based around minimum of
**                         snowfall or basetemp; may exhibit freeze thaw
**                     4 = constant at basetemp
**                     5 = minimum temperature calculated using snow temp buffering at
**                         snow base under current snowpack.
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
**                temperature under snow
**    bio       - leaf biomass (g/m^2)
**    sfaltmp   - surface temp under snow, for latching models this is the minimum;
**                large value flags no snow
**    basetemp  - temp at base of snow layer
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
#include <stdio.h>
/*void abortmssg(char buffer[256]);*/

#ifndef min
  #define min(a,b) (((a) < (b)) ? (a) : (b))
  #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

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
      static float sfaltmp = 256.0f;

      if(sfaltmp == 256.0f) {
        sfaltmp = 128.0f;
        if     (SnowFlag ==-1) {printf("SnowFlag = -1; run level minimum insulated temp\n");}
        else if(SnowFlag == 0) {printf("SnowFlag = 0;  no snow insulation\n");}
        else if(SnowFlag == 1) {printf("SnowFlag = 1;  snow insulation with basetemp max; Parton paper\n");}
        else if(SnowFlag == 2) {printf("SnowFlag = 2;  snow insulation, around basetemp\n");}
        else if(SnowFlag == 3) {printf("SnowFlag = 3;  low temp snow insulation, around the minimum of basetemp or snowfall temperature\n");}
        else if(SnowFlag == 4) {printf("SnowFlag = 4;  constant; basetemp\n");}
        else if(SnowFlag == 5) {printf("SnowFlag = 5;  snow cover minimum insulated temp\n");}
        else {
          fprintf(stderr,"unknown SnowFlag in 'surftemp'. Expecting 0-5. SnowFlag = %d\n", SnowFlag);
          exit(1);
        }
      }

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
      tmxs_mlt =  0.5f - tmns_mlt + 0.5f;
      *srfctemp = tmxs_mlt * *tmxs + tmns_mlt * *tmns; /* temperature assuming no snow */

      /* Surface temperature without snow */
      if (snowpack <= 0.00000001 || SnowFlag == 0) {
        *diurnal_range = *tmxs - *tmns;
        if(SnowFlag > 1) {sfaltmp = 128.0f;} /* under water, no snow fall temp */
      }

      /* Surface temperature with insulation effects of snow 11/30/95 (Bill Parton). */
      else {

        /* Snow base temperature
           impact of surface litter carbon on mineral soil surface temperature
           On Jan 17, 2013, at 2:56 PM, Keough,Cynthia wrote:
           Modify the constants in the equation that is calculating effect of surface
           litter on snow soil surface temperature.
           Code before:
           basetemp = min((2.0f - (-2.0f)) / (1000.0f - 0.0f) *
                         (litrcrbn - 1000.0f) + 2.0f, 0.0f);
                    = min(4.0f * (litrcrbn/ 1000.0f - 1.0f) + 2.0f, 0.0f);

           Revised equation
           basetemp = min(4.0f * (litrcrbn/3000.0f - 1.0f) + 2.0f, 0.0f);
           NOTE basetemp ~ -2.0C for normal litter litrcrbn KLK 6/2016
        */
        basetemp = min(-2.0f + litrcrbn/750.0f, 0.0f); /* Simplified algebra */

        if(sfaltmp == 128.0f) {sfaltmp = *srfctemp;} /* save temperature at snowfall */

        if(SnowFlag == 3) {
          basetemp = min(basetemp, sfaltmp); /* low temperature model base temperature */
        }

        /* these are the pieces of partons snow temperature fit.
           The If splits these for constant and reactive the sub models, 3 and 4 */
        if ((SnowFlag != 4 && tmin + tmax < 0.0) || SnowFlag == 2 || SnowFlag == 3) {
          /* average air temperature is below freezing */
          snowmult = max(1.0f - 0.15f * snowpack , 0.0f);
          *srfctemp = basetemp + (0.3f * (tmin + tmax) / 2.0f) * snowmult;
          /* *srfctemp = min(*srfctemp, 0.0); */
          *diurnal_range = 0.3f * (tmax - tmin) * snowmult;
          if (*diurnal_range / 2.0 + *srfctemp > 0.0) {
            *diurnal_range = basetemp * *srfctemp;
          }
        }

        else if (SnowFlag == 4 || tmin + tmax >= 0.0) {
          /* if there is snow, and average air temperature gets above
             freezing, average soil surface temperature stays at freezing
             NOTE: applies to SnowFlag -1,1,5;  0,2,3 have been handled so omit them from if */
          *srfctemp = basetemp;
          *diurnal_range = 0.3f * (tmax - tmin) / 2.0f;
          if (*diurnal_range / 2.0 + *srfctemp > 0.0) {
            *diurnal_range = basetemp * *srfctemp;
          }
        }

        /* code to hole srfctemp less than a season or run minimum temperature
           ratchets down the minimum temperature seen under the snow layer
           The original coding (SnowFlag == -1) didn't reset the minimum and it continued
           to decrease over THE ENTIRE RUN!
           Obviously A BUG, snow insulation doesn't carry past snowmelt. KLK 6/2016
        */
        if(SnowFlag == -1  || SnowFlag == 5) {
          /* If the soil has frozen prior to snow accumulating on the site
             the insulating effect of the snow will hold the low frozen soil
             temperature until the snow melts, cak 09/29/2011 */
          if (sfaltmp < *srfctemp) {*srfctemp = sfaltmp;}
          sfaltmp = *srfctemp; /* ratchet the previous down */
        }
        /* printf("surftemp = %f\n", *srfctemp); */
      }

      /* Added code to allow warming of the soil surface temperature, */
      /* cak - 07/02/03 */
      if (stsys > 0 && time >= ststart) {
        *srfctemp = *srfctemp + stamt;
      }

      return;
    }
