
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      calcSrad.c
**
**  FUNCTION:  void calc_srad()
**
**  PURPOSE:   Calculate short wave radiation for current day.
**
**  AUTHOR:    Adapted from Peter Thorton's code extracted from the Sipnet
**             model (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    curday      - current day of the year (1..366)
**    dailyPrecip - precipitation for current day (cm)
**    maxTemp     - maximum temperature for current day (degrees C)
**    minTemp     - minimum temperature for current day (degrees C)
**    tdew        - dewpoint temperature for current day (degrees C)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    ctrl              - program control variables
**    ctrl->indewpt     - dewpoint temperature input flag, 0=NO, 1=YES
**    params            - site parameters
**    sitedata          - site data structure
**    sitedata->s_tmax  - site maximum temperature value (degrees C)
**    sitedata->s_tmin  - site minimum temperature value (degrees C)
**    sitedata->s_prcp  - site daily precipitation values (cm)
**    sitedata->s_tdew  - site dewpoint temperature value (degrees C)
**    yeardata          - structure containing year long data ararys of
**                        ttmax0, potential radiation, and day length
**
**  LOCAL VARIABLES:
**    humid - site humidity value, VPD or VP (Pa)
**    yday  - index for accessing yeardata array values (0..365)
**
**  OUTPUTS:
**    srad - shortwave radiation value for the current day (W/m2)
**
**  CALLED BY:
**    simsom()
**
**  CALLS:
**    calc_prcp
**    calc_srad_humidity
**    calc_srad_humidity_iterative
**    snowpack
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calcSrad.h"

    void calc_srad(double *maxTemp, double *minTemp, double *dailyPrecip,
                   int *curday, double *tdew, double *srad)
    {

      extern CONTROL_SPT ctrl;
      extern DATA_SPT sitedata;
      extern YEARDATA_SPT yeardata;
      extern PARAMETER_SPT params;

      int    yday;
      double tmean;

      /* Initialize site parameters */

      /* assign meteorological data passed by the DayCent model file into */
      /* sitedata structure */
      sitedata->s_tmax = *maxTemp;
      sitedata->s_tmin = *minTemp;
      sitedata->s_tdew = *tdew;
      sitedata->s_prcp = *dailyPrecip;
      tmean = (sitedata->s_tmax + sitedata->s_tmin)/2.0;
      sitedata->s_tday = ((sitedata->s_tmax - tmean)*TDAYCOEF) + tmean;

      yday = *curday - 1;

      /* estimate daily snowpack */
      snowpack(sitedata);

      /* test for the presence of Tdew observations, and branch to the */
      /* appropriate srad and humidity algorithms */
      if (ctrl->indewpt) {
        /* estimate srad and humidity using real Tdew data */
        calc_srad_humidity(ctrl, params, sitedata, yeardata, yday);
      } else {  /* no dewpoint temperature data */
        /* estimate srad and humidity with iterative algorithm */
        calc_srad_humidity_iterative(ctrl, params, sitedata, yeardata, yday);
      }

      *srad = sitedata->s_srad;

      return;
    }  /* end of calc_srad() */


/*****************************************************************************
**
**  FUNCTION:  void snowpack()
**
**  PURPOSE:   Estimates the accumulation and melt of snow for radiation
**             algorithm corrections.
**
**  AUTHOR:    Adapted from Peter Thorton's code extracted from the Sipnet
**             model (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    sitedata - site data structure
**
**  GLOBAL VARIABLES:
**    SNOW_TCRIT - critical temperature for snowmelt, degrees C (-6.0)
**    SNOW_TRATE - snowmelt rate, cm/degC/day (0.042)
**
**  EXTERNAL VARIABLES:
**    sitedata         - site data structure
**    sitedata->s_prcp - site precipitation value (cm)
**    sitedata->s_swe  - site snow water equivalent value (cm)
**    sitedata->s_tmin - site minimum temperature value (degrees C)
**
**  LOCAL VARIABLES:
**    newsnow  - new snow water equivalent (cm)
**    snowmelt - snow melt water equivalent (cm)
**    snowpack - snowpack water equivalent (cm)
**
**  OUTPUTS:
**    sitedata->s_swe - site snow water equivalent value (cm)
**
**  CALLED BY:
**    calc_srad()
**
**  CALLS:
**    None
**
*****************************************************************************/

    void snowpack(DATA_SPT sitedata)
    {

      double snowpack,newsnow,snowmelt;

      /* first pass to initialize SWE */
      snowpack = 0.0;
      newsnow = 0.0;
      snowmelt = 0.0;
      if (sitedata->s_tmin <= SNOW_TCRIT) {
        newsnow = sitedata->s_prcp;
      } else {
        snowmelt = SNOW_TRATE * (sitedata->s_tmin - SNOW_TCRIT);
      }
      snowpack += newsnow - snowmelt;
      if (snowpack < 0.0) {
        snowpack = 0.0;
      }
      sitedata->s_swe = snowpack;

      return;
   }  /* end of snowpack() */


/*****************************************************************************
**
**  FUNCTION:  void calc_srad_humidity()
**
**  PURPOSE:   Calculate solar radiation and humidity for current day.  When
**             dewpoint temperature observations are available radiation and
**             humidity can be estimated directly.
**
**  AUTHOR:    Adapted from Peter Thorton's code extracted from the Sipnet
**             model (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    ctrl     - program control variables
**    params   - site parameters
**    sitedata - site data structure
**    yday     - current day of year, index for yeardata arrays
**    yeardata - structure containing year long data ararys of ttmax0,
**               potential radiation, and day length
**
**  GLOBAL VARIABLES:
**    ABASE       - vapor pressure effect on transmittance, 1/Pa (-6.1e-5)
**    B0          - radiation parameter, dimensionless (0.013)
**    B1          - radiation parameter, dimensionless (0.201)
**    B2          - radiation parameter, dimensionless (0.185)
**    C           - radiation parameter, dimensionless (1.5)
**    DIF_ALB     - diffuse albedo for horizon correction, dimensionless (0.6)
**    RADPERDEG   - radians per degree (0.01745329)
**    RAIN_SCALAR - correction to transmittance for rain day, dimensionless
**                  (0.75)
**
**  EXTERNAL VARIABLES:
**    params                   - site parameters
**    params->site_ehoriz      - site east horizon, degrees
**    params->site_slp         - site slope, degrees
**    params->site_whoriz      - site west horizon, degrees
**    sitedata                 - site data structure
**    sitedata->s_tmax         - site maximum temperature value (degrees C)
**    sitedata->s_tmin         - site minimum temperature value (degrees C)
**    sitedata->s_prcp         - site daily precipitation values (cm)
**    sitedata->s_dayl         - site daylength value (seconds)
**    sitedata->s_srad         - site shortwave radiation value (W/m2)
**    sitedata->s_swe          - site snow water equivalent value (cm)
**    sitedata->s_tday         - site daylight average temperature value
**                               (degrees C)
**    sitedata->s_tdew         - site dewpoint temperature value (degrees C)
**    yday                     - current day of year, index for yeardata
**                               arrays
**    yeardata                 - structure containing year long data ararys of
**                               ttmax0, potential radiation, and day length
**    yeardata->daylength[]    - length of day light (seconds)
**    yeardata->flat_potrad[]  - daylight average flux density for a flat
**                               surface
**    yeardata->slope_potrad[] - daylight average flux density for the slope
**    yeardata->ttmax0[]       - maximum daily total transmittance
**
**  LOCAL VARIABLES:
**    avg_horizon    - average horizon angle of sky proportion for diffuse
**                     radiation calculation
**    b              - effect of diurnal temperature range on daily total
**                     transmittance
**    dtr            - diurnal temperature range (degrees C)
**    horizon_scalar - horizon scalar of sky proportion for diffuse
**                     radiation calculation
**    pdif           - fraction of radiation that is diffuse
**    pdir           - fraction of radiation that is direct
**    pva            - vapor pressure
**    sc             - snow correction
**    sky_prop       - sky proportion for diffuse radiation
**    slope_excess   - excess slope of sky proportion for diffuse
**                     radiation calculation
**    slope_scalar   - slope scalar of sky proportion for diffuse
**                     radiation calculation
**    srad1          - direct radiation
**    srad2          - diffuse radiation
**    t_final        - final daily total transmittance
**    t_fmax         - proportion of daily maximum transmittance
**    t_tmax         - maximum transmittance corrected for vapor pressure
**    tmax           - maximum temperature for day
**    tmin           - minimum temperature for day
**
**  OUTPUTS:
**    sitedata->s_dayl - site daylength value (seconds)
**    sitedata->s_srad - site shortwave radiation value (W/m2)
**
**  CALLED BY:
**    calc_srad()
**
**  CALLS:
**    None
**
*****************************************************************************/

    void calc_srad_humidity(CONTROL_SPT ctrl, PARAMETER_SPT params,
                            DATA_SPT sitedata, YEARDATA_SPT yeardata,
                            int yday)
    {

      double pva;
      double dtr;
      double tmax,tmin;
      double t_tmax,b,t_fmax;
      double t_final,pdif,pdir,srad1,srad2; 
      double sc;
      double sky_prop;
      double avg_horizon, slope_excess;
      double horizon_scalar, slope_scalar;

      /* estimate radiation using Tdew observations */
      /* calculate diurnal temperature range for transmittance calculations */
      tmax = sitedata->s_tmax;
      tmin = sitedata->s_tmin;
      if (tmax < tmin) {
        tmax = tmin;
      }
      dtr = tmax-tmin;

      /*****************************************
       *                                       *
       * start of the main radiation algorithm *
       *                                       *
       *****************************************/

      /* STEP (1) calculate pressure ratio (site/reference) = f(elevation) */
      /* done in initsrad subroutine */

      /* STEP (2) correct initial transmittance for elevation */
      /* done in initsrad subroutine */

      /* STEP (3) build 366-day array of ttmax0, potential radiation, and */
      /*  daylength */
      /* done in initsrad subroutine */

      /* STEP (4)  calculate the sky proportion for diffuse radiation */
      /* uses the product of spherical cap defined by average horizon angle */
      /* and the great-circle truncation of a hemisphere. this factor does */
      /* not vary by yearday. */
      avg_horizon = (params->site_ehoriz + params->site_whoriz)/2.0;
      horizon_scalar = 1.0 - sin(avg_horizon * RADPERDEG);
      if (params->site_slp > avg_horizon) {
        slope_excess = params->site_slp - avg_horizon;
      } else {
        slope_excess = 0.0;
      }
      if (2.0*avg_horizon > 180.0) {
        slope_scalar = 0.0;
      } else {
        slope_scalar = 1.0 - (slope_excess/(180.0 - 2.0*avg_horizon));
        if (slope_scalar < 0.0) {
          slope_scalar = 0.0;
        }
      }
      sky_prop = horizon_scalar * slope_scalar;

      /* STEP (5)  final calculation of daily total radiation */
      /* correct this day's maximum transmittance for vapor pressure */
      pva = 610.7 * exp(17.38 * sitedata->s_tdew /
                        (239.0 + sitedata->s_tdew));
      t_tmax = yeardata->ttmax0[yday] + ABASE * pva;

      /* b parameter from DTR */
      b = B0 + B1 * exp(-B2 * dtr);

      /* proportion of daily maximum transmittance */
      t_fmax = 1.0 - 0.9 * exp(-b * pow(dtr,C));

      /* correct for precipitation if this is a rain day */
      if (sitedata->s_prcp) {
        t_fmax *= RAIN_SCALAR;
      }

      /* final daily total transmittance */
      t_final = t_tmax * t_fmax;

      /* estimate fraction of radiation that is diffuse, on an */
      /* instantaneous basis, from relationship with daily total */
      /* transmittance in Jones (Plants and Microclimate, 1992) Fig 2.8, */
      /* p. 25, and Gates (Biophysical Ecology, 1980) Fig 6.14, p. 122. */
      pdif = -1.25*t_final + 1.25;
      if (pdif > 1.0) {
        pdif = 1.0;
      }
      if (pdif < 0.0) {
        pdif = 0.0;
      }

      /* estimate fraction of radiation that is direct, on an instantaneous */
      /* basis */
      pdir = 1.0 - pdif;

      /* the daily total radiation is estimated as the sum of the following */
      /* two components: */
      /* 1. The direct radiation arriving during the part of the day when */
      /*    there is direct beam on the slope. */
      /* 2. The diffuse radiation arriving over the entire daylength (when */
      /*    sun is above ideal horizon). */

      /* component 1 (direct) */
      srad1 = yeardata->slope_potrad[yday] * t_final * pdir;

      /* component 2 (diffuse) */
      /* includes the effect of surface albedo in raising the diffuse */
      /* radiation for obstructed horizons */
      srad2 = yeardata->flat_potrad[yday] * t_final * pdif *
              (sky_prop + DIF_ALB*(1.0-sky_prop));

      /* snow pack influence on radiation */
      if (sitedata->s_swe > 0.0) {
         /* snow correction in J/m2/day */
        sc = (1.32 + 0.096 * sitedata->s_swe) * 1e6;
        /* convert to W/m2 and check for zero daylength */
        if (yeardata->daylength[yday] > 0.0) {
          sc /= yeardata->daylength[yday];
        } else {
          sc = 0.0;
        }
        /* set a maximum correction of 100 W/m2 */
        if (sc > 100.0) {
          sc = 100.0;
        }
      } else {
        sc = 0.0;
      }

      /* save daily radiation and daylength */
      sitedata->s_srad = srad1 + srad2 + sc;
      sitedata->s_dayl = yeardata->daylength[yday];

      return;
    } /* end of calc_srad_humidity() */


/*****************************************************************************
**
**  FUNCTION:  void calc_srad_humidity_iterative()
**
**  PURPOSE:   Calculate solar radiation and humidity for current day.  In the
**             of dew point temperature observations an iterative estimation
**             of shortwave radiation and humidity is required.
**
**  AUTHOR:    Adapted from Peter Thorton's code extracted from the Sipnet
**             model (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    ctrl     - program control variables
**    params   - site parameters
**    sitedata - site data structure
**    yday     - current day of year, index for yeardata arrays
**    yeardata - structure containing year long data ararys of ttmax0,
**               potential radiation, and day length
**
**  GLOBAL VARIABLES:
**    ABASE       - vapor pressure effect on transmittance, 1/Pa (-6.1e-5)
**    B0          - radiation parameter, dimensionless (0.013)
**    B1          - radiation parameter, dimensionless (0.201)
**    B2          - radiation parameter, dimensionless (0.185)
**    C           - radiation parameter, dimensionless (1.5)
**    DIF_ALB     - diffuse albedo for horizon correction, dimensionless (0.6)
**    RADPERDEG   - radians per degree (0.01745329)
**    RAIN_SCALAR - correction to transmittance for rain day, dimensionless
**                  (0.75)
**
**  EXTERNAL VARIABLES:
**    params                   - site parameters
**    params->site_ehoriz      - site east horizon, degrees
**    params->site_elev        - site elevation (meters)
**    params->site_slp         - site slope, degrees
**    params->site_whoriz      - site west horizon, degrees
**    sitedata                 - site data structure
**    sitedata->s_tmax         - site maximum temperature value (degrees C)
**    sitedata->s_tmin         - site minimum temperature value (degrees C)
**    sitedata->s_prcp         - site daily precipitation values (cm)
**    sitedata->s_dayl         - site daylength value (seconds)
**    sitedata->s_srad         - site shortwave radiation value (W/m2)
**    sitedata->s_swe          - site snow water equivalent value (cm)
**    sitedata->s_tday         - site daylight average temperature value
**    sitedata->s_tmin         - site minimum temperature value (degrees C)
**                               (degrees C)
**    sitedata->s_tdew         - site dewpoint temperature value (degrees C)
**    yday                     - current day of year, index for yeardata
**                               arrays
**    yeardata                 - structure containing year long data ararys of
**                               ttmax0, potential radiation, and day length
**    yeardata->daylength[]    - length of day light (seconds)
**    yeardata->flat_potrad[]  - daylight average flux density for a flat
**                               surface
**    yeardata->slope_potrad[] - daylight average flux density for the slope
**    yeardata->ttmax0[]       - maximum daily total transmittance
**
**  LOCAL VARIABLES:
**    ann_pet        - annual potential evapotranspiration
**    ann_prcp       - annual total precipitation
**    avg_horizon    - average horizon angle of sky proportion for diffuse
**                     radiation calculation
**    b              - effect of diurnal temperature range on daily total
**                     transmittance
**    dtr            - diurnal temperature range (degrees C)
**    effann_prcp    - effective annual precipitation
**    horizon_scalar - horizon scalar of sky proportion for diffuse
**                     radiation calculation
**    pa             - air pressure at site
**    pdif           - fraction of radiation that is diffuse
**    pdir           - fraction of radiation that is direct
**    pet            - potential evapotranspiration
**    pva            - vapor pressure
**    ratio          - ratio of PET/effann_prcp
**    ratio2         - ratio squared
**    ratio3         - ratio cubed
**    sc             - snow correction
**    sky_prop       - sky proportion for diffuse radiation
**    slope_excess   - excess slope of sky proportion for diffuse
**                     radiation calculation
**    slope_scalar   - slope scalar of sky proportion for diffuse
**                     radiation calculation
**    srad1          - direct radiation
**    srad2          - diffuse radiation
**    sum_pet        - sum of potential evapotranspiration for time step
**    t_final        - final daily total transmittance
**    t_fmax         - proportion of daily maximum transmittance
**    t_tmax         - maximum transmittance corrected for vapor pressure
**    tdew           - dew point temperature (degrees C)
**    tdewk          - dew point correction factor
**    tmax           - maximum temperature for day (degrees C)
**    tmin           - minimum temperature for day (degrees C)
**    tmink          - minimum temperature correction factor for dew point
**                     calculations
**
**  OUTPUTS:
**    sitedata->s_dayl - site daylength value (seconds)
**    sitedata->s_srad - site shortwave radiation value (W/m2)
**
**  CALLED BY:
**    calc_srad()
**
**  CALLS:
**    atm_pres()
**    calc_pet()
**
*****************************************************************************/

    void calc_srad_humidity_iterative(CONTROL_SPT ctrl, PARAMETER_SPT params,
                                      DATA_SPT sitedata,
                                      YEARDATA_SPT yeardata, int yday)
    {

      double dtr;
      double t_fmax, tdew;
      double ann_prcp,effann_prcp;
      double sum_pet,ann_pet;
      double tmax,tmin;
      double pva,t_tmax,b;
      double tmink,pet,ratio,ratio2,ratio3,tdewk;
      double t_final,pdif,pdir,srad1,srad2; 
      double pa;
      double sc;
      double sky_prop;
      double avg_horizon, slope_excess;
      double horizon_scalar, slope_scalar;

      /* calculate diurnal temperature range for transmittance calculations */
      tmax = sitedata->s_tmax;
      tmin = sitedata->s_tmin;
      if (tmax < tmin) {
        tmax = tmin;
      }
      dtr = tmax-tmin;

      /* calculate the annual total precip for decision between simple and */
      /* arid-corrected humidity algorithm */
      ann_prcp = sitedata->s_prcp * 365.25;
      if (ann_prcp == 0.0) {
        ann_prcp = 1.0;
      }

      /* Since we are receiving one day of precipitation data from DayCent */
      /* use a simple total scaled to effective annual precip */
      effann_prcp = sitedata->s_prcp * 365.25;
      /* if the effective annual precip for this period is less than 8 cm, */
      /* set the effective annual precip to 8 cm to reflect an arid */
      /* condition, while avoiding possible division-by-zero errors and */
      /* very large ratios (PET/Pann) */
      if (effann_prcp < 8.0) {
        effann_prcp = 8.0;
      }

      /*****************************************
       *                                       *
       * start of the main radiation algorithm *
       *                                       *
       *****************************************/

      /* before starting the iterative algorithm between humidity and */
      /* radiation, calculate all the variables that don't depend on */
      /* humidity so they only get done once. */

      /* STEP (1) calculate pressure ratio (site/reference) = f(elevation) */
      /* done in initsrad subroutine */

      /* STEP (2) correct initial transmittance for elevation */
      /* done in initsrad subroutine */

      /* STEP (3) build 366-day array of ttmax0, potential radidation, and */
      /* daylength */
      /* done in initsrad subroutine */

      /* STEP (4)  calculate the sky proportion for diffuse radiation */
      /* uses the product of spherical cap defined by average horizon angle */
      /* and the great-circle truncation of a hemisphere. this factor does */
      /* not vary by yearday. */
      avg_horizon = (params->site_ehoriz + params->site_whoriz)/2.0;
      horizon_scalar = 1.0 - sin(avg_horizon * RADPERDEG);
      if (params->site_slp > avg_horizon) {
        slope_excess = params->site_slp - avg_horizon;
      } else {
        slope_excess = 0.0;
      }
      if (2.0*avg_horizon > 180.0) {
        slope_scalar = 0.0;
      } else {
        slope_scalar = 1.0 - (slope_excess/(180.0 - 2.0*avg_horizon));
        if (slope_scalar < 0.0) {
          slope_scalar = 0.0;
        }
      }
      sky_prop = horizon_scalar * slope_scalar;

      /* b parameter, and t_fmax not varying with Tdew, so these can be */
      /* calculated once, outside the iteration between radiation and */
      /* humidity estimates. */
      /* b parameter from DTR */
      b = B0 + B1 * exp(-B2 * dtr);

      /* proportion of daily maximum transmittance */
      t_fmax = 1.0 - 0.9 * exp(-b * pow(dtr,C));

      /* correct for precipitation if this is a rain day */
      if (sitedata->s_prcp) {
        t_fmax *= RAIN_SCALAR;
      }

      /* As a first approximation, calculate radiation assuming that */
      /* Tdew = Tmin */
      tdew = sitedata->s_tmin;
      pva = 610.7 * exp(17.38 * tdew / (239.0 + tdew));
      t_tmax = yeardata->ttmax0[yday] + ABASE * pva;

      /* final daily total transmittance */
      t_final = t_tmax * t_fmax;

      /* estimate fraction of radiation that is diffuse, on an */
      /* instantaneous basis, from relationship with daily total */
      /* transmittance in Jones (Plants and Microclimate, 1992) Fig 2.8, */
      /* p. 25, and Gates (Biophysical Ecology, 1980) Fig 6.14, p. 122. */
      pdif = -1.25*t_final + 1.25;
      if (pdif > 1.0) {
        pdif = 1.0;
      }
      if (pdif < 0.0) {
        pdif = 0.0;
      }

      /* estimate fraction of radiation that is direct, on an instantaneous */
      /* basis */
      pdir = 1.0 - pdif;

      /* the daily total radiation is estimated as the sum of the following */
      /* two components: */
      /* 1. The direct radiation arriving during the part of the day when */
      /*    there is direct beam on the slope. */
      /* 2. The diffuse radiation arriving over the entire daylength (when */
      /*    sun is above ideal horizon). */

      /* component 1 */
      srad1 = yeardata->slope_potrad[yday] * t_final * pdir;

      /* component 2 (diffuse) */
      /* includes the effect of surface albedo in raising the diffuse */
      /* radiation for obstructed horizons */
      srad2 = yeardata->flat_potrad[yday] * t_final * pdif *
              (sky_prop + DIF_ALB*(1.0-sky_prop)); 

      /* snow pack influence on radiation */
      if (sitedata->s_swe > 0.0) {
        /* snow correction in J/m2/day */
        sc = (1.32 + 0.096 * sitedata->s_swe) * 1e6;
        /* convert to W/m2 and check for zero daylength */
        if (yeardata->daylength[yday] > 0.0) {
          sc /= yeardata->daylength[yday];
        } else {
          sc = 0.0;
        }
        /* set a maximum correction of 100 W/m2 */
        if (sc > 100.0) {
          sc = 100.0;
        }
      } else {
        sc = 0.0;
      }

      /* save daily radiation and daylength */
      sitedata->s_srad = srad1 + srad2 + sc;
      sitedata->s_dayl = yeardata->daylength[yday];

      /* estimate annual PET first, to decide which humidity algorithm */
      /* should be used */
      /* estimate air pressure at site */
      pa = atm_pres(params->site_elev);
      sum_pet =
        calc_pet(sitedata->s_srad,sitedata->s_tday,pa,sitedata->s_dayl);
      ann_pet = sum_pet * 365.25;

      /* humidity algorithm decision: */
      /* PET/prcp >= 2.5 -> arid correction */
      /* PET/prcp <  2.5 -> use tdew-tmin, which is already finished */
      if (ann_pet/ann_prcp >= 2.5) {
        /* Estimate Tdew using the initial estimate of radiation for PET */
        tmink = sitedata->s_tmin + 273.15;
        pet = sum_pet;

        /* calculate ratio (PET/effann_prcp) and correct the dewpoint */
        ratio = pet/effann_prcp;
        ratio2 = ratio*ratio;
        ratio3 = ratio2*ratio;
        tdewk = tmink*(-0.127 + 1.121*
          (1.003 - 1.444*ratio + 12.312*ratio2 - 32.766*ratio3) + 0.0006*
          (dtr));
        tdew = tdewk - 273.15;

        /* Revise estimate of radiation using new Tdew */
        pva = 610.7 * exp(17.38 * tdew / (239.0 + tdew));
        t_tmax = yeardata->ttmax0[yday] + ABASE * pva;

        /* final daily total transmittance */
        t_final = t_tmax * t_fmax;

        /* estimate fraction of radiation that is diffuse, on an */
        /* instantaneous basis, from relationship with daily total */
        /* transmittance in Jones (Plants and Microclimate, 1992) Fig 2.8, */
        /* p. 25, and Gates (Biophysical Ecology, 1980) Fig 6.14, p. 122. */
        pdif = -1.25*t_final + 1.25;
        if (pdif > 1.0) {
          pdif = 1.0;
        }
        if (pdif < 0.0) {
          pdif = 0.0;
        }

        /* estimate fraction of radiation that is direct, on an */
        /* instantaneous basis */
        pdir = 1.0 - pdif;

        /* the daily total radiation is estimated as the sum of the */
        /* following two components: */
        /* 1. The direct radiation arriving during the part of the day when */
        /*    there is direct beam on the slope. */
        /* 2. The diffuse radiation arriving over the entire daylength */
        /*    (when sun is above ideal horizon). */

        /* component 1 */
        srad1 = yeardata->slope_potrad[yday] * t_final * pdir;

        /* component 2 (diffuse) */
        /* includes the effect of surface albedo in raising the diffuse */
        /* radiation for obstructed horizons */
        srad2 = yeardata->flat_potrad[yday] * t_final * pdif *
                (sky_prop + DIF_ALB*(1.0-sky_prop)); 

        /* snow pack influence on radiation */
        if (sitedata->s_swe > 0.0) {
          /* snow correction in J/m2/day */
          sc = (1.32 + 0.096 * sitedata->s_swe) * 1e6;
          /* convert to W/m2 and check for zero daylength */
          if (yeardata->daylength[yday] > 0.0) {
            sc /= yeardata->daylength[yday];
          } else {
            sc = 0.0;
          }
          /* set a maximum correction of 100 W/m2 */
          if (sc > 100.0) {
            sc = 100.0;
          }
        } else {
          sc = 0.0;
        }

        /* save daily radiation */
        sitedata->s_srad = srad1 + srad2 + sc;

        /* Revise estimate of Tdew using new radiation */
        tmink = sitedata->s_tmin + 273.15;
        pet = calc_pet(sitedata->s_srad,sitedata->s_tday,pa,sitedata->s_dayl);

        /* calculate ratio (PET/effann_prcp) and correct the dewpoint */
        ratio = pet/effann_prcp;
        ratio2 = ratio*ratio;
        ratio3 = ratio2*ratio;
        tdewk = tmink*(-0.127 + 1.121*
          (1.003 - 1.444*ratio + 12.312*ratio2 - 32.766*ratio3) + 0.0006*
          (dtr));
        tdew = tdewk - 273.15;
      }  /* end of arid-correction humidity and radiation estimation */

      return;
    } /* end of calc_srad_humidity_iterative() */


/*****************************************************************************
**
**  FUNCTION:  double atm_pres()
**
**  PURPOSE:   Calculate the atmospheric pressure (Pa) as a function of
**             elevation (m).  From the discussion on atmospheric statics in:
**             Iribane, J.V., and W.L. Godson, 1981. Atmospheric
**             Thermodynamics,  2nd Edition.  D. Reidel Publishing Company,
**             Dordrecht, The Netherlands.  (p. 168)
**
**  AUTHOR:    Adapted from Peter Thorton's code extracted from the Sipnet
**             model (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    elev - elevation (meters)
**
**  GLOBAL VARIABLES:
**    G_STD  - standard gravitational accel, m s-2, (9.80665)
**    LR_STD - standard temperature lapse rate, -K m-1, (0.0065)
**    MA     - molecular weight of air, kg mol-1, (28.9644e-3)
**    P_STD  - standard pressure at 0.0 m elevation, Pa, (101325.0)
**    R      - gas law constant, m3 Pa mol-1 K-1 (8.3143)
**    T_STD  - standard temp at 0.0 m elevation, K, (288.15)
**
**  EXTERNAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    pa - atmospheric pressure (Pa)
**    t1 - first term in calculation
**    t2 - second term in calculation
**
**  OUTPUTS:
**    pa - atmospheric pressure (Pa)
**
**  CALLED BY:
**    calc_srad_humidity_iterative()
**
**  CALLS:
**    None
**
*****************************************************************************/

    double atm_pres(double elev)
    {

      double t1,t2;
      double pa;

      t1 = 1.0 - (LR_STD * elev)/T_STD;
      t2 = G_STD / (LR_STD * (R / MA));
      pa = P_STD * pow(t1,t2);

      return(pa);
    }  /* end of atm_pres() */


/*****************************************************************************
**
**  FUNCTION:  double calc_pet()
**
**  PURPOSE:   Calculate the potential evapotranspiration for aridity
**             corrections in calc_vpd(), according to Kimball et al., 1997
**
**  AUTHOR:    Adapted from Peter Thorton's code extracted from the Sipnet
**             model (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    dayl - daylength (seconds)
**    pa   - air pressure (Pa)
**    rad  - daylight average incident shortwave radiation (W/m2)
**    ta   - daylight average air temperature (degrees C)
**
**  GLOBAL VARIABLES:
**    EPS - ratio of molec weights, unitless MW/MA (0.62196351)
**
**  EXTERNAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    dt    - offset for saturation vapor pressure calculation
**    gamma - psychrometer parameter (Pa K-1)
**    lhvap - latent heat of vaporization of water (J kg-1)
**    pet   - potential evapotranspiration (kg m-2 day-1)
**    pvs1  - saturated vapor pressures (Pa)
**    pvs2  - saturated vapor pressures (Pa)
**    rnet  - absorbed shortwave radiation avail. for ET (W m-2)
**    s     - slope of saturated vapor pressure curve (Pa K-1)
**    t1    - air temperature offset (deg C)
**    t2    - air temperature offset (deg C)
**
**  OUTPUTS:
**    pet - potential evapotranspiration (centimeters/day) 
**
**  CALLED BY:
**    calc_srad_humidity_iterative()
**
**  CALLS:
**    None
**
*****************************************************************************/

    double calc_pet(double rad, double ta, double pa, double dayl)
    {

      double rnet;
      double lhvap;
      double gamma;
      double dt = 0.2;
      double t1, t2;
      double pvs1, pvs2;
      double pet;
      double s;

      /* calculate absorbed radiation, assuming albedo = 0.2 and ground */
      /* heat flux = 10% of absorbed radiation during daylight */
      rnet = rad * 0.72;

      /* calculate latent heat of vaporization as a function of ta */
      lhvap = 2.5023e6 - 2430.54 * ta;

      /* calculate the psychrometer parameter:                              */
      /*   gamma = (cp pa)/(lhvap epsilon)                                  */
      /* where:                                                             */
      /*   cp       specific heat of air (J/kg K)                           */
      /*   epsilon - ratio of molecular weights of water and air (unitless) */
      gamma = CP * pa / (lhvap * EPS);

      /* estimate the slope of the saturation vapor pressure curve at ta */
      /* temperature offsets for slope estimate */
      t1 = ta+dt;
      t2 = ta-dt;

      /* calculate saturation vapor pressures at t1 and t2, using formula */
      /* from Abbott, P.F., and R.C. Tabony, 1985. The estimation of */
      /* humidity parameters.  Meteorol. Mag., 114:49-56. */
      pvs1 = 610.7 * exp(17.38 * t1 / (239.0 + t1));
      pvs2 = 610.7 * exp(17.38 * t2 / (239.0 + t2));

      /* calculate slope of pvs vs. T curve near ta */
      s = (pvs1-pvs2) / (t1-t2);

      /* calculate PET using Priestly-Taylor approximation, with */
      /* coefficient set at 1.26.  Units of result are kg/m^2/day, */
      /* equivalent to mm water/day */
      pet = (1.26 * (s/(s+gamma)) * rnet * dayl)/lhvap;

      /* return a value in centimeters/day, because this value is used in a */
      /* ratio to annual total precip, and precip units are centimeters */
      return (pet/10.0);
    }  /* end of calc_pet() */
