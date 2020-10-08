
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:  watrflow.c
**
**  FUNCTION:  void watrflow()
**
**  PURPOSE: Water-flow submodel.  This submodel is a rewrite of a
**           model originally written by William Parton.  It simulates
**           the flow of water through the plant canopy and soil.
**           See "Abiotic Section of ELM", as a reference.
** 
**  AUTHOR:  Susan Chaffee    4/30/92 
** 
**  REWRITE:  Melannie Hartman  9/20/93 - 8/21/96
**
**  HISTORY:
**    11/16/01 (CAK) Soil water potential computed at field capacity in
**                   the trwtavg subroutine.
**    05/19/08 (CAK) Soil water potential returned from trwtavg subroutine
**                   is calculated using the wettest soil layer in the plant
**                   rooting zone
**    01/22/14 (CAK) Calculate the water effect of transpiration using the
**                   same equation as used when calculating the water effect
**                   on potential production, the coefficients for this
**                   calculation are read from the sitepar.in file.
**
**  INPUTS:
**    accum      - the amount of snow added to the snowpack (cm H2O)
**    aetdly     - actual evapotranspiration (cm H2O)
**    aggreenc   - the amount of photosynthetic active carbon (g/m2)
**    amovdly    - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**    asmos[]    - soil water content by layer (cm H2O)
**    avgtemp    - average air temperature for the day (deg C)
**    avh2o[]    - water available for plant growth (avh2o[0]), plant survival
**                 (avh2o[1]), and in the first two Century soil layers 
**                 (avh2o[2])
**    basef      - the fraction of soil water content of the soil layer below
**                 the bottom of the soil profile which is lost via base flow
**    baseflow   - soil water content water lost to base flow from the soil
**                 layer directly below the bottom layer of the soil profile 
**    biodead    - above ground dead biomass (g/m2)
**    biolive    - above ground live biomass (g/m2)
**    blitter    - above ground litter (g/m2)
**    co2val     - CO2 effect on transpiration.   Added 8/14/98 -mdh
**    daylength  - amount of daylight hours (1..24) 
**    elitst     - effect of litter on soil surface temperature
**    evaptot    - total amount of water evaporated from all soil layers
**                 (cm H2O)
**    intrcpt    - amount of precipitation intercepted by the standing crop
**                 and litter (cm H2O)
**    jday       - current day of the year (1..366)
**    litrcrbn   - amount of carbon in the surface litter (sum of som1c(1),
**                 som2c(1), strucc(1), and metabc(1)
**    melt       - the amount of water added to the soil when melted from the
**                 snowpack, if daily air temperature is warm enough (cm H2O)
**    month      - current month of the year (1..12)
**    nlayer     - number of layers in Century soil profile
**    nlaypg     - number of Century soil layers used for plant growth and
**                 root death
**    outflow    - water that runs off, or drains out of the profile (cm H2O)
**    pet        - daily potential evapotranspiration rates (cm H2O)
**    pmxbio     - maximum dead biomass level for soil surface temperature
**                 calculation
**    pmntmp     - effect of biomass on minimum soil surface temperature
**    pmxtmp     - effect of biomass on maximum soil surface temperature
**    pottransp  - bare-soil transpiration loss rate (cm H2O/day)
**    pptactual  - the current day's precipitation (cm H2O)
**    rhumid     - average relative humidity for the day (% 1..100)
**    runoffdly  - amount of water (rain or snowmelt) which did not infiltrate
**                 soil profile (cm)
**    rwcf[]     - relative water content by layer
**    snlq       - the liquid water in the snowpack (cm H2O)
**    snowpack   - current snowpack (equiv. cm H2O)
**    soiltavewk - average soil temperature in the second soils.in soil layer
**                 over the previous 7 days (deg C)
**    solrad     - total incoming shortwave radiation (langleys/day)
**    srad       - shortwave radiation value for day (W/m^2)
**    srfctemp   - soil surface temperature (deg C)
**    stamt      - amount of warming applied to soil surface temperature
**                 (degrees C)
**    stcrlai    - standing crop leaf area index
**    stream1    - runoff plus baseflow (runoff has replaced stormflow)
**    strplt     - year in simulation to start output
**    ststart    - start time for soil warming simulation
**    stsys      - flag to indicate if soil warming will be simulated
**                 0 = no, 1 = yes
**    sublim     - amount of water sublimated from the snowpack (cm H2O)
**    tadep      - depth in the soil profile above which transpiration occurs
**                 (cm)
**    tempmax    - maximum air temperature for the day (deg C)
**    tempmin    - minimum air temperature for the day (deg C)
**    time       - simulation time (years)
**    tmelt[]    - tmelt[0] = melting temperature (deg C), if temperature is
**                 >= this value, snow is allowed to melt and is added to the
**                 precipitation
**                 tmelt[1] = the slope of the melting equation
**                 (cm snow / deg C)
**    tmns       - minimum soil surface temperature (deg C)
**    tmxs       - maximum soil surface temperature (deg C)
**    tranlayers - number of soil layers used in transpiration calculations
**    trantot    - total amount of water transpired from all soil layers
**                 (cm H2O)
**    watertable - 1 = simulate water table, 0 = no water table
**    watr2sat   - amount of water automatically added to the system to bring
**                 the soil water content in the full soil profile to
**                 saturation
**    watrflag   - 1 = water added to the system automatically to bring the
**                     soil water content to saturation
**                 0 = rain and irrigation events used to add water to the
**                     soil potentially bringing it to saturation
**    wfluxout[] - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**    windsp     - average daily windspeed at 2 meters (mph)
**    woodb      - wood biomass, fine branch + large wood (g/m^2)
**
**  GLOBAL VARIABLES:
**    CENTMAXLYR - maximum number of Century soil layers (10)
**    MAXLYR     - maximum number of soil water model layers (21)
**
**  EXTERNAL VARIABLES:
**    files                 - structure containing information about output
**                            files
**    files->fp_soiltavg    - file pointer to soiltavg.out output file
**    files->fp_soiltmax    - file pointer to soiltmax.out output file
**    files->fp_soiltmin    - file pointer to soiltmin.out output file
**    files->fp_stempdx     - file pointer to stemp_dx.out output file
**    files->fp_swc         - file pointer to vswc.out output file
**    files->fp_wfps        - file pointer to wfps.out output file
**    files->write_soiltavg - flag to indicate if soiltavg.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmax - flag to indicate if soiltmax.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmin - flag to indicate if soiltmin.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_stempdx  - flag to indicate if stemp_dx.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_swc      - flag to indicate if vswc.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_wfps     - flag to indicate if wfps.out output file should
**                            be created, 0 = do not create, 1 = create
**    flags                 - structure containing debugging flags
**    flags->debug          - flag to set debugging mode, 0 = off, 1 = on
**    layers                - soil water soil layer structure
**    layers->bulkd[]       - bulk density by layer (g/cm3)
**    layers->clayfrac[]    - clay fraction in soil layer, 0.0 - 1.0
**    layers->depth[]       - the distance from the surface to the middle of
**                            the soil layer (cm)
**    layers->ecoeff[]      - bare-soil evaporation water absorption
**                            coefficients by layer
**    layers->fieldc[]      - volumetric water content at field capacity for
**                            layer (cm H2O/cm of soil)
**    layers->lbnd[]        - the index of the lower soil water model layer
**                            which corresponds to given soil layer in Century
**    layers->minpot[]      - minimum matric potential by layer based on
**                            swcmin (-cm)
**    layers->nelyrs        - number of layers to consider in evaporation
**    layers->numlyrs       - total number of layers in the soil water model
**                            soil profile
**    layers->orgfrac[]     - organic matter in soil layer, fraction 0.0 - 1.0
**    layers->sandfrac[]    - sand fraction in soil layer, 0.0 - 1.0
**    layers->satcond[]     - saturated hydraulic conductivity by layer
**                            (cm/sec)
**    layers->sumecoeff     - sum of evaporation coefficients
**    layers->swc[]         - soil water content by layer (cm H2O)
**    layers->swcfc[]       - volumetric soil water content at field capacity
**                            for layer (cm H2O/cm of soil)
**    layers->swcmin[]      - lower bound on soil water content by layer
**                            (cm H2O) soil water content will not be allowed
**                            to drop below this minimum
**    layers->tcoeff[]      - transpiration water absoption coefficients by
**                            layer (ND)
**    layers->width[]       - the thickness of soil water model layers (cm)
**    layers->wfps[]        - water-filled pore space by layer
**                            (fraction of a porespace that is filled with
**                            water, 0.0-1.0)
**    sitepar               - site specific parameters structure for soil
**                            water model
**    sitepar->albedo       - fraction of light reflected by snow
**    sitepar->cldcov[]     - average cloud cover for the month (%, 1..100)
**    sitepar->dmpflux      - damping factor for soil water flux (in h2oflux)
**    sitepar->hours_rain   - the duration of the rainfall/snowmelt event
**                            (hours)
**    sitepar->hpotdeep     - hydraulic water potential of deep storage
**                            layer, the more negative the number the dryer
**                            the soil layer (units?)
**    sitepar->ksatdeep     - saturated hydraulic conductivity of deep
**                            storage layer (cm/sec)
**    sitepar->rlatitude    - latitude of the site (in radians)
**    sitepar->sublimscale  - multiplier to scale sublimation
**    sitepar->usexdrvrs    - 0 = use air temperature to drive PET rates
**                            1 = use extra weather drivers (solrad, rhumid,
**                                windsp) for PET calculation, 
**                            2 = use extra weather drivers (srad, vpd)
**                                for photosynthesis calculation
**                            3 = use extra drivers for both PET and
**                                photosynthesis calculations
**    sitepar->watertable   - 1 = simulate water table, 0 = no water table
**    soil                  - soil temperature structure
**    soil->soiltavg[]      - average soil temperature of layer (degrees C)
**    soil->soiltmax[]      - maximum soil temperature by layer (degrees C)
**    soil->soiltmin[]      - minimum soil temperature by layer (degrees C)
**    soil->stmtemp[]       - soil surface temperature (degrees Celsius)
**
**  LOCAL VARIABLES:
**    blivelai      - live biomass leaf area index
**    biomass       - above ground biomass (g/m2)
**    bserate       - bare soil evaporation rates (cm/day)
**    bstrate       - bare soil transpiration rates (cm/day)
**    cwlit         - cummulative water in litter to date (cm H2O)
**    cwstcr        - cummulative water in standing crops to date (cm H2O)
**    diurnal_range - the difference between the maximum and minimum surface
**                    soil temperatures (deg C)
**    fbse          - fraction of bare soil water loss by evaporation
**    fbst          - fraction of bare soil water loss by transpiration
**    h2ogef        - water effect on transpiration
**    ii            - index
**    ilyr          - current layer in the soil profile
**    intrcpt_limit - limit on interception, total interception of the
**                    precipitation by the standing crop and litter will not
**                    be allowed to exceed 30% of the PET
**    petleft       - water left to be evaporated from soil layers after,
**                    water has been evaporated/transpired from standing crop
**                    and litter
**    pptsoil       - remaining precipitation after interception (cm H2O)
**    soilEvap      - water evaporated from soil surface (cm H2O)
**    stlyrs        - number of soil temperature layers to output
**    sumintrcpt    - sum of water intercepted by standing crop and litter
**                    for the day (cm H2O)
**    swcsat        - the soil water content of a layer at saturation (cm H2O)
**    swctemp[]     - temporary storage location for layers->swc[] array
**                    values
**    totagb        - total monthly above ground biomass (g/m2)
**    totlit        - total water intercepted by litter, to date (cm H2O)
**    totstcr       - total water intercepted by standing crop, to date
**                    (cm H2O)
**    transp[]      - water transpired from each layer (cm H2O)
**    vegcov        - vegetation cover (units?)
**    watrinput     - precipitation left after interception by standing crop
**                    (cm H2O)
**    wintlit       - water intercepted by litter for the day (cm H2O)
**    wintstcr      - water intercepted by standing crops for the day (cm H2O)
**    wtosat        - the amount of water to be added to a layer to bring it
**                    to the desired level of saturation (cm H2O)
**
**  OUTPUTS:
**    aetdly     - actual evapotranspiration (cm H2O)
**    accum      - the amount of snow added to the snowpack (cm H2O)
**    amovdly    - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**    asmos[]    - soil water content by layer (cm H2O)
**    avh2o[]    - water available for plant growth (avh2o[0]), plant survival
**                 (avh2o[1]), and in the first two Century soil layers 
**                 (avh2o[2])
**    baseflow   - soil water content water lost to base flow from the soil
**                 layer directly below the bottom layer of the soil profile 
**    evaptot    - total amount of water evaporated from all soil layers
**                 (cm H2O)
**    intrcpt    - amount of precipitation intercepted by the standing crop
**                 and litter (cm H2O)
**    melt       - the amount of water added to the soil when melted from the
**                 snowpack, if daily air temperature is warm enough (cm H2O)
**    outflow    - water that runs off, or drains out of the profile (cm H2O)
**    pottransp  - bare-soil transpiration loss rate (cm H2O/day)
**    runoffdly  - amount of water (rain or snowmelt) which did not infiltrate
**                 soil profile (cm)
**    rwcf[]     - relative water content by layer
**    snlq       - the liquid water in the snowpack (cm H2O)
**    snowpack   - current snowpack (equiv. cm H2O)
**    soiltavewk - average soil temperature in the second soils.in soil layer
**                 over the previous 7 days (deg C)
**    srfctemp   - soil surface temperature (deg C)
**    stcrlai    - standing crop leaf area index
**    stream1    - runoff plus baseflow (runoff has replaced stormflow)
**    sublim     - amount of water sublimated from the snowpack (cm H2O)
**    trantot    - total amount of water transpired from all soil layers
**                 (cm H2O)
**    watr2sat   - amount of water automatically added to the system to bring
**                 the soil water content in the full soil profile to
**                 saturation
**    wfluxout[] - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    fracbslos()    - calculate fraction of water loss from bare soil 
**                     evaporation and transpiration
**    h2oflux()      - water-flow submodel
**    initdaily()    - initialize daily values for the soil water model
**    litstcr_evap() - evaporate water from litter and standing crop
**    pteevap()      - adjust bare-soil evaporation and transpiration rates so
**                     that the day's total evaporation/transpiration does not
**                     exceed the day's PET, also increase the day's AET by
**                     the bare-soil evaporation/transpiration rates
**    potbse()       - calculate potential bare soil evaporation rate
**    potbst()       - calculate potential transpiration rate
**    setamov()      - set amovdly (passed to Century) using the value
**                     wfluxout (daily soil water model)
**    setasmos()     - set asmos, avh2o and rfwc (Century variables) from swc,
**                     the soil water content in the daily soil water model
**    showlyrs()     - print the soil water content by layer
**    snowCent()     - accumulate a snow pack if temperatures are low enough,
**                     calculate melt and sublimation
**    snowmodel()    - compute snow melt and sublimation (used when using
**                     extra weather drivers)
**    soiltemp()     - calculate the daily average, maximum, and minimum soil
**                     temperature for a specified number of soil depths
**    soiltransp()   - transpire water from soil
**    surftemp()     - compute the minimum, maximum, and average soil surface
**                     temperatures
**    trwtavg()      - compute soil water potential of wettest soil layer in
**                     plant rooting zone to be used for transpiration
**                     calculations
**    watrlit()      - calculate water intercepted by litter
**    watrstcr()     - calculate the water intercepted by standing crop
**    wfps()         - compute the daily water-filled pore space by layer
**    wrtstemp()     - write out the soil temperature by layer
**    wrtswc()       - write out the soil water content by layer
**    wrtwfps()      - write out the water-filled pore space by layer
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "soilwater.h"

    void watrflow(int *jday, int *month, int *nlayer, int *nlaypg,
                  int *watertable, int *watrflag, float *avgtemp,
                  float *tempmin, float *tempmax, float *solrad,
                  float *rhumid, float *windsp, float *pptactual,
                  float *biolive, float *blitter, float *biodead,
                  float rwcf[CENTMAXLYR], float avh2o[3],
                  float asmos[CENTMAXLYR], float *snowpack, float *snlq,
                  float amovdly[CENTMAXLYR], float *pet, float *evaptot,
                  float *trantot, float *stream1, float *basef,
                  float *pottransp, float *baseflow, float *accum,
                  float *melt, float *intrcpt, float *outflow, float tmelt[],
                  float *sublim, float wfluxout[], float *time, float *strplt,
                  float *co2val, float *tmns, float *tmxs, float *runoffdly,
                  float *tadep, float *soiltavewk, float *daylength,
                  float *woodb, float *elitst, float *pmxtmp, float *pmntmp,
                  float *pmxbio, float *srfctemp, float *stsys,
                  float *ststart, float *stamt, float *litrcrbn,
                  float *aggreenc, float *watr2sat, float *aetdly,
                  float *h2ogef, float *stcrlai, double *srad)
    {
      extern FLAG_SPT flags;
      extern LAYERPAR_SPT layers;
      extern SITEPAR_SPT sitepar;
      extern SOIL_SPT soil;
      extern FILES_SPT files;

      int ilyr;
      int stlyrs = 63;
      float fbse;
      float fbst;
      float totlit;
      float totstcr;
      float petleft;
      float sumintrcpt;
      float intrcpt_limit;
      float watrinput;
      float pptsoil;
      float bserate;
      float bstrate;
      float wintstcr;
      float wintlit;
      static float cwstcr = 0.0f;
      static float cwlit = 0.0f; 
      float soilEvap;
      float transp[MAXLYR];
      float biomass, blivelai, vegcov, totagb;
      float swctemp[MAXLYR];
      float diurnal_range;
      int   ii, tranlayers;
      float swcsat;
      float wtosat;

      *sublim = 0.0f;
      *pottransp = 0.0f;
      *h2ogef = 0.0f;
      fbse = 0.0f;
      fbst = 0.0f;
      totlit = 0.0f;
      totstcr = 0.0f;
      petleft = 0.0f;
      watrinput = 0.0f;
      pptsoil = 0.0f;
      bserate = 0.0f;
      bstrate = 0.0f;
      wintstcr = 0.0f;
      wintlit = 0.0f;
      *aetdly = 0.0f;

      if (flags->debug > 1) {
        printf("Entering function watrflow\n");
      }

      for(ilyr=0; ilyr < layers->numlyrs; ilyr++) {
        transp[ilyr] = 0.0f;
      }
         
      soilEvap = 0.0f;
      *evaptot = 0.0f;
      *trantot = 0.0f;
      *melt = 0.0f;
      *accum = 0.0f;
      *outflow = 0.0f;

      /* Call daily initialization routine */
      initdaily(*month, *biolive, *biodead, *blitter, &biomass, &blivelai,
                &vegcov, &totagb, *aggreenc, stcrlai, layers);

      for(ilyr=0; ilyr < layers->numlyrs; ilyr++) {
        if (layers->swc[ilyr] + 0.0000001 < layers->swcmin[ilyr]) {
          printf("h2oflux:  swcmin[%1d] > swc[%1d]\n", ilyr, ilyr);
          printf("swcmin[%1d] = %1f, swc[%1d] = %1f\n", ilyr,
                 layers->swcmin[ilyr], ilyr, layers->swc[ilyr]);
/*          exit(1); */
        }
      }
/*      printf("pet in watrflow = %8.4f\n", *pet); */
      petleft = *pet;

      /* Calculate soil water potential for wettest soil layer in the */
      /* plant rooting zone, cak - 05/19/08 */
      *h2ogef = trwtavg(nlaypg, layers, flags, sitepar);

      /* if snow is to be accumulated, call the snow routine */

/*      printf("tempmax = %6.2f\n", *tempmax);
      printf("avgtemp = %6.2f\n", *avgtemp);
      printf("pptactual = %6.2f\n", *pptactual);
      printf("solrad = %6.2f\n", *solrad);
      printf("rhumid = %6.2f\n", *rhumid);
      printf("windsp = %6.2f\n", *windsp);
      printf("petleft = %8.4f\n", petleft); */

/*!!
fprintf(files->fp_snow, "%7.2f %4d %7.3f %7.3f %7.3f %7.3f %7.3f",
        *time, *jday, *avgtemp, *pptactual, *pet, *snowpack, *snlq);
!!*/
      if (sitepar->usexdrvrs == 1 || sitepar->usexdrvrs == 3) {
        snowmodel(*jday, *tempmax, *tempmin, *avgtemp, -1.0f, *pptactual,
                  &pptsoil, *solrad, *rhumid, *windsp,
                  sitepar->cldcov[*month], sitepar->rlatitude,
                  sitepar->albedo, snowpack, melt, accum, sublim, &petleft,
                  sitepar->sublimscale, petleft, tmelt, snlq, *srad,
                  *daylength);
      } else {
        snowCent(tmelt, *avgtemp, *pptactual, &pptsoil, snowpack, snlq,
                 &petleft, melt, accum, sublim, *tempmin, *tempmax, *srad,
                 *daylength);
      }

      *aetdly += *sublim;

      if (petleft < 0.0) {
        fprintf(stderr, "ERROR in PET/AET balance.  petleft = %12.10f\n",
                petleft);
        exit(1);
      }

      /* calculate the water intercepted by the standing crop and litter.
         For each, return the amount intercepted, along with the precip. left over
         Exclude interception is no precip, and no snowfall 
         -- Exclude interception for rain-on-snow where the snowpack (snowpack>o)
            absorbs the precip (*melt or pptsoil ==0)  to prevent double counting KLK
         -- reversing the sense of the if made it easier to include new conditions
         KLK  12May2012
      */
      if (*pptactual == 0.0  || *accum != 0.0  || 
          (*snowpack > 0.0 && *melt == 0.0)) {
        wintstcr = 0.0f;
        wintlit = 0.0f;
      } else {
        watrstcr(&pptsoil, &wintstcr, *pptactual, vegcov);
        watrinput = *pptactual - wintstcr;
        watrlit(watrinput, &pptsoil, &wintlit, *blitter);
      }

      /* Limit total water interception to 30% PET */
      /*      evaporation is petleft limited... why use PET for this limit?  KLK */
      intrcpt_limit = 0.30f * (*pet);
      sumintrcpt = wintstcr + wintlit;
      if (sumintrcpt > intrcpt_limit) {
        /*RARE condition -- limit interception to unallocated water in rain on snow.
          This is probably sufficient to solve the fully absorbed case, pptsoil == 0
          but the if above saves a lot of 0 calculations. KLK 12May2011
         */
        if(intrcpt_limit > *melt  &&  *snowpack > 0.0) {intrcpt_limit = *melt;} 
        wintstcr *= intrcpt_limit/sumintrcpt; 
        wintlit *= intrcpt_limit/sumintrcpt; 
           
        /* More water is available for infiltration now */
        pptsoil += sumintrcpt - wintstcr - wintlit;
      }

      *intrcpt = wintstcr + wintlit;

      /* Calculate the fraction of water loss from bare soil evaporation */
      /* and transpiration. */

      fracbslos(&fbse, &fbst, blivelai);

      /* Calculate the potential bare soil evaporation rate            */
      /* If there is a snowpack, no bare soil evaporation will occur,  */
      /* therefore set bserate to zero (MDH) - 1/14/94                 */

      if (*snowpack > 0) {
        bserate = 0.0f;
      } else { 
        potbse(&bserate, totagb, fbse, petleft, layers, *aggreenc);
      }

      /* Calculate the potential bare soil transpiration rate */
      /* Pass the soil water potential computed at field capacity to the */
      /* potbst subroutine, see call to trwtavg above, cak - 11/16/01 */
      /* Pass the soil water potential of the wettest soil layer in the */
      /* plant rooting zone, see call to trwtavg above, cak - 05/19/08 */

      potbst(&bstrate, *h2ogef, *biolive, *biodead, fbst, petleft, *co2val);

      /* Sum cumulative water intercepted by litter */

      cwlit += wintlit;
      totlit = cwlit;

      /* Sum cumulative water intercepted by standing crop */

      cwstcr += wintstcr;
      totstcr = cwstcr;

      /*  Evaporate as much water as possible first from the standing */
      /*  crop and then litter.  Total evaporation/transporation will */
      /*  not exceed the PET for the day */

      litstcr_evap(&cwlit, &cwstcr, &petleft, aetdly, totlit, totstcr);

      /*  Reduce the potential bare-soil evapotranspiration rates */
      /*  if necessary to prevent total evapotranspiration from */
      /*  exceeding the day's PET.  Adjust AET also. */

      pteevap(&bserate, &bstrate, petleft);

      *pottransp = bstrate;
      /* Limit transpiration as a result of root resistance to water */
      /* flow, cak - 10/13/2009 */
      if (*pottransp > 0.36) {
        *pottransp = (0.42f - 0.36f) / (6.0f - 0.36f) * (petleft - 6.0f) +
                     0.42f;
      }

      /* At this point, soil water contents are equal to yesterday's values */

      /* Calculate soil surface temperature */
      surftemp(*elitst, *pmxtmp, *pmntmp, *pmxbio, *tempmax, *tempmin, tmxs,
               tmns, srfctemp, *daylength, *biolive, *biodead, *blitter,
               *woodb, *stsys, *ststart, *stamt, *time, *snowpack,
               &diurnal_range, *litrcrbn, sitepar->SnowFlag);

      /* Calculate soil temperature for all soil layers */
      soiltemp(*jday, *tempmin, *tempmax, layers->depth, layers->width,
               layers->fieldc, layers->sandfrac, layers->clayfrac,
               layers->orgfrac, layers->bulkd, layers->swc, layers->numlyrs,
               soil->soiltavg, soil->soiltmin, soil->soiltmax, soil->stmtemp,
               *tmns, *tmxs, soiltavewk, *srfctemp, diurnal_range);

      if (flags->debug > 1) {
        printf("Before h2oflux: ");
        showlyrs(layers->swc, layers->numlyrs);
      }

      h2oflux(*jday, layers->numlyrs, layers->swc, layers->swcmin,
              layers->minpot, layers->depth, layers->width, layers->satcond,
              layers, soil->soiltavg, pptsoil, bserate, &soilEvap,
              sitepar->hours_rain, sitepar->dmpflux, aetdly, outflow,
              wfluxout, *snowpack, runoffdly, basef, baseflow,
              *watertable, sitepar->hpotdeep, sitepar->ksatdeep);

/*!!
fprintf(files->fp_snow, "%7.3f\n", outflow);
!!*/

      if (flags->debug > 1) {
        printf("After h2oflux: ");
        showlyrs(layers->swc, layers->numlyrs);
      }

      /* If petleft > 0, evaporate/transpire water from the soil... */
      /* the total water already evaporated/transpired from the  */
      /* standing crop and litter is less than the day's PET. */

/*      printf("petleft = %8.4f\n", petleft); */

      if (petleft > 0) { 

        /* Transpiration occurs in only down to the calculated tadep depth */
        /* Compute the number of DayCent layers that are in the calculated */
        /* tadep depth.  Pass this value to the soiltransp subroutine */
        /* in place of layers->numlyrs.  cak - 01/29/03 */
        ii = 0;
        while ((layers->dpthmx[ii] <= *tadep) && (ii < layers->numlyrs)) {
          ii++;
        }
        tranlayers = ii;

        soiltransp(transp, tranlayers, bstrate, layers, aetdly);
      } else {
        bserate=0.0f;
        bstrate=0.0f;
      }

      *evaptot += soilEvap;
      for(ilyr=0; ilyr < layers->numlyrs; ilyr++) {
        *trantot += transp[ilyr];
      }
      /* stormflow has been replaced by runoff */
/*      *stormflow = wfluxout[layers->numlyrs-1] * (*stormf); */

      *stream1 = *outflow;

      /* When simulating a water table, if desired, add enough water to */
      /* keep the soil at saturation, cak - 02/02/2012 */
      *watr2sat = 0.0f;
      if ((*watertable == 1) && (*watrflag == 1)) {
        for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
          swcsat = (1.0f - layers->bulkd[ilyr] / (float)PARTDENS) *
                   layers->width[ilyr];
          wtosat = swcsat - (float)layers->swc[ilyr];
          if (wtosat > 0.0) {
            layers->swc[ilyr] += wtosat;
            *watr2sat += wtosat;
          }
        }
      }

      for(ilyr=0; ilyr <= layers->numlyrs; ilyr++) {
        swctemp[ilyr] = (float)(layers->swc[ilyr]);
      }
      setasmos(asmos, nlayer, swctemp, &layers->numlyrs, avh2o, rwcf);
      setamov(amovdly, *nlayer, wfluxout, layers->numlyrs, layers->lbnd);

      if (*time >= *strplt) {

        wfps(layers);

        if (files->write_swc) {
          wrtswc(files->fp_swc, *time, *jday, layers->swc, layers->width,
                 layers->numlyrs);
        }

        if (files->write_wfps) {
          wrtwfps(files->fp_wfps, *time, *jday, layers->wfps, layers->numlyrs,
                  layers->width);
        }

        if (files->write_soiltavg) {
          wrtstemp(files->fp_soiltavg, *time, *jday, soil->soiltavg,
                   layers->numlyrs);
        }

        if (files->write_soiltmin) {
          wrtstemp(files->fp_soiltmin, *time, *jday, soil->soiltmin,
                   layers->numlyrs);
        }

        if (files->write_soiltmax) {
          wrtstemp(files->fp_soiltmax, *time, *jday, soil->soiltmax,
                   layers->numlyrs);
        }

        if (files->write_stempdx) {
          wrtstemp(files->fp_stempdx, *time, *jday, soil->stmtemp, stlyrs);
        }

      }

      if (flags->debug > 1) {
        printf("Exiting function watrflow\n");
      }

      return;
    }
