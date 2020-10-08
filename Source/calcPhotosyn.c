

/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      calcPhotosyn.c
**
**  FUNCTION:  void calcPsnEffects()
**
**  PURPOSE:   Calculate the temperature and VPD effects on photosynthesis.
**             These calculations will be the same for grass/crops or trees so
**             there is no need to calculate them separately for each system.
**
**  AUTHOR:    Adapted from code extracted from the Sipnet model
**             (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    average_temp  - average temperature for daylight hours (degrees C)
**    average_vpd   - average vapor pressure deficit (kPa)
**    dayLength     - length of day in hours
**    dTemp         - decrease in photosynthesis due to temperature
**    dVpd          - decrease in photosynthesis due to vapor pressure deficit
**    dVpdExp       - exponential value in dVpd equation
**    dVpdSlope     - slope value in dVpd equation
**    dWater        - effect of water stress on photosynthesis
**    etrans        - actual evapotranspiration for the day (cm H2O) or
**                    transpiration rate as a function of potential
**                    evapotranspiration and soil water potential
**    h2ogef3       - effect of water stress on photosynthesis
**    maxTemp       - maximum temperature for day (degrees C)
**    minTemp       - minimum temperature for day (degrees C)
**    obs_vpd       - observed vapor pressure deficit value for the day as
**                    read from the weather data file (kPa)
**    potETrans     - potential evapotranspiration for the day (cm H2O) or 1.0
**    psnTMax       - maximum temperature at which net photosynthesis
**                    occurs, assumed symmetrical around psnTOpt
**                    (degrees C)
**    psnTMin       - minimum tempature at which net photosynthesis occurs
**                    (degrees C)
**    tminintercept - intercept term used to adjust minimum temperature for
**                    calculating VPD at dewpoint
**    tminslope     - slope term used to adjust minimum temperature for
**                    calculating VPD at dewpoint
**
**  GLOBAL VARIABLES:
**    TINY - to avoid divide-by-zero errors (0.000001)
**
**  EXTERNAL VARIABLES:
**    climate             - climate parameters
**    climate->maxTemp    - maximum temperature for day (degrees C)
**    climate->minTemp    - minimum temperature for day (degrees C)
**    climate->tair       - average air temperature for this time step
**                          (degrees C)
**    climate->vpd        - average vapor pressure deficit for this time step
**                          (kPa)
**    psparams            - parameter values read from the photosyn.in file
**    psparams->dVpdExp   - exponential value in dVpd equation
**    psparams->dVpdSlope - slope value in dVpd equation
**    psparams->psnTMax   - maximum temperature at which net photosynthesis
**                          occurs, assumed symmetrical around psnTOpt
**                          (degrees C)
**    psparams->psnTMin   - minimum tempature at which net photosynthesis
**                          occurs (degrees C)
**    psparams->psnTOpt   - optimal temperature at which net photosynthesis
**                          occurs (degrees C)
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    average_temp - average temperature for daylight hours (degrees C)
**    average_vpd  - average vapor pressure deficit (kPa)
**    dTemp        - decrease in photosynthesis due to temperature
**    dVpd         - decrease in photosynthesis due to vapor pressure deficit
**    dWater       - effect of water stress on photosynthesis
**
**  CALLED BY:
**    potprod()
**
**  CALLS:
**    calc_vpd()
**    calcMoistureEff()
**    calcTemperatureEff()
**    calcVpdEff()
**
*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "calcPhotosyn.h"

    void calcPsnEffects(double *maxTemp, double *minTemp, double *etrans,
                        double *potETrans, double *average_temp,
                        double *average_vpd, double *dTemp, double *dVpd,
                        float *dWater, float *dayLength, double *obs_vpd,
                        double *dVpdExp, double *dVpdSlope, double *psnTMin,
                        double *psnTOpt, float *tminslope,
                        float *tminintercept, float *h2ogef3)
    {

      extern CLIMATE_SPT climate;
      extern PSPARAMS_SPT psparams;

      climate->maxTemp = *maxTemp;
      climate->minTemp = *minTemp;

      psparams->dVpdExp = *dVpdExp;
      psparams->dVpdSlope = *dVpdSlope;
      psparams->psnTMin = *psnTMin;
      psparams->psnTOpt = *psnTOpt;
      /* calculate psnTMax, assumed symmetrical */
      psparams->psnTMax = psparams->psnTOpt +
                          (psparams->psnTOpt - psparams->psnTMin);

      /* Calculate average temperature and vapor pressure deficit for */
      /* daylight hours */
      calc_vpd(average_vpd, average_temp, climate->maxTemp, climate->minTemp,
               *dayLength, *tminslope, *tminintercept);
      /* If we have an observed vpd value for the current day replace the */
      /* calculated value with the observed vpd */
      if (*obs_vpd <= -99.0) {
        /* do nothing */
      } else {
        *average_vpd = *obs_vpd;
      }
      climate->tair = *average_temp;
      climate->vpd = *average_vpd;
      /* avoid divide by zero */
      if (climate->vpd < TINY) {
         climate->vpd = TINY;
      }

      /* calculate decrease in photosynthesis due to temperature */
      calcTemperatureEff(climate->tair, dTemp);

      /* calculate decrease in photosynthesis due to vapor pressure deficit */
      calcVpdEff(climate->vpd, dVpd);

      /* factor for effect of water stress on photosynthesis */
/*      calcMoistureEff(*etrans, dWater, *potETrans); */
      /* the factor for the effect of water stress on photosynthesis is now */
      /* being calculated in the watrflow subroutine based on the soil */
      /* content of the wettest soil layer in the plant rooting zone */
      /* CAK - 01/22/2014 */
      *dWater = *h2ogef3;

    } /* end of calcPsnEffects() */


/*****************************************************************************
**
**  FUNCTION:  void calc_vpd()
**
**  PURPOSE:   Calculate vapor pressure deficit for daylight hours.
**
**  AUTHOR:    Craig Maxwell - 02/2009
**
**  INPUTS:
**    average_temp  - average temperature for daylight hours (degrees C)
**    average_vpd   - average vapor pressure deficit (kPa)
**    dayLength     - length of day in hours
**    maxTemp       - maximum temperature for day (degrees C)
**    minTemp       - minimum temperature for day (degrees C)
**    tminintercept - intercept term used to adjust minimum temperature for
**                    calculating VPD at dewpoint
**    tminslope     - slope term used to adjust minimum temperature for
**                    calculating VPD at dewpoint
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**      hour               - index, current hour of day
**      hour_counter       - index, number of daylight hours in day
**      min_vapor_pressure - vapor pressure based on adjusted minimum
**                           temperature value for the day (mmHg)
**      rounded_sunrise    - time of sunrise rounded to nearest integer
**                           (24 hour clock)
**      rounded_sunset     - time of sunset rounded to nearest integer
**                           (24 hour clock)
**      sunrise            - time of sunrise (24 hour clock)
**      sunset             - time of sunset (24 hour clock)
**      temp_sum           - sum of temperature values for each hour, used in
**                           average_vpd calculation (degrees C)
**      temperature        - calculated temperature for a given hour
**                           (degrees C)
**      vapor_pressure     - vapor pressure for the current hour based on
**                           calculated hourly temperature value (mmHg)
**      vpd                - vapor pressure deficit for current hour (mmHg)
**      vpd_sum            - sum of vpd values for each hour, used in
**                           average_vpd calculation (mmHg)
**      newMinTemp         - adjusted minimum temperature value (degrees C)
**
**  OUTPUTS:
**    average_temp - average temperature for daylight hours (degrees C)
**    average_vpd  - average vapor pressure deficit (kPa)
**
**  CALLED BY:
**    calcPsnEffects()
**
**  CALLS:
**    hourTemp()
**    svapor()
**
*****************************************************************************/
    void calc_vpd(double *average_vpd, double *average_temp, double maxTemp,
                  double minTemp, float dayLength, float tminslope,
                  float tminintercept)
    {

      /* Using 12 for the first call to hourTemp since we're just looking */
      /* for the sunrise and sunset times. */
      int hour = 12;
      int hour_counter = 0;
      double min_vapor_pressure;
      int rounded_sunrise, rounded_sunset;
      double sunrise, sunset;
      double temp_sum = 0;
      double temperature;
      double vapor_pressure;
      double vpd;
      double vpd_sum = 0;
      double newMinTemp;

      /* Call hourTemp to get sunrise and sunset. */
      temperature = hourTemp(&sunrise, &sunset, hour, maxTemp, minTemp,
                    dayLength);

      rounded_sunrise = (int)(sunrise + 0.5);
      rounded_sunset = (int)(sunset + 0.5);

      /* Adjust the minimum temperature used to calculate the VPD at */
      /* dewpoint.  This is necessary for dry sites where the dewpoint */
      /* can be less than the actual minimum temperature for the day */
      newMinTemp = tminslope * minTemp + tminintercept;
      min_vapor_pressure = svapor((float)newMinTemp);

      for (hour=rounded_sunrise; hour<=rounded_sunset; ++hour) {
        ++hour_counter;

        /* Calculate the temperature for this hour. */
        temperature = hourTemp(&sunrise, &sunset, hour, maxTemp, minTemp,
                               dayLength);

        vapor_pressure = svapor((float)temperature);

        vpd = vapor_pressure - min_vapor_pressure;

        temp_sum += temperature;
        vpd_sum += vpd;
      }

      *average_temp = temp_sum / hour_counter;
      *average_vpd = vpd_sum / hour_counter;

      /* Convert mmHg to kPa */
      *average_vpd *= (101.325 / 760.0);
    }  /* end of calc_vpd() */
   
   
/*****************************************************************************
**
**  FUNCTION:  double hourTemp()
**
**  PURPOSE:   Calculate the temperature for a given hour.
**
**  AUTHOR:    Craig Maxwell - 02/2009
**
**  INPUTS:
**    dayLength - length of day in hours
**    hour      - hour of day (24 hour clock)
**    maxTemp   - maximum temperature for day (degrees C)
**    minTemp   - minimum temperature for day (degrees C)
**    sunrise   - time of sunrise (24 hour clock)
**    sunset    - time of sunset (24 hour clock)
**
**  GLOBAL VARIABLES:
**    B                          - coefficient that controls temperature
**                                 decrease at night (2.1)
**                                 "Stream Temperature Investigations: Field
**                                 and Analytic Methods, Instream Flow
**                                 Information Paper No. 13" by John M.
**                                 Bartholow, p. 41. suggests using a value of
**                                 of 2.1
**    TIME_LAG_MAX               - time lag in maximum temperature after noon
**                                 (hours) (1.8)
**                                 "A Model for Diurnal Variation in Soil and
**                                 Air Temperature" by William J. Parton and
**                                 Jesse A. Logan suggests using a value of
**                                 1.8
**    TIME_LAG_MIN               - time lag for the minimum temperature after
**                                 sunrise (hours) (0)
**
**  EXTERNAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ddy                        - daylength adjusted by the minimum time lag
**    hrsAfterMinTempUntilSunset - number of hours after the minimum
**                                 temperature until sunset
**    hrsAfterSunsetUntilMinTemp - number of hours after sunset until the time
**                                 of the minimum temperature
**    nightLength                - length of night in hours
**    sunsetTemp                 - temperature at sunset (degrees C)
**    temp                       - temperature for the given hour (degrees C)
**
**  OUTPUTS:
**    temp - temperature for the given hour (degrees C)
**
**  CALLED BY:
**    calc_vpd()
**
**  CALLS:
**    None
**
*****************************************************************************/
    double hourTemp(double *sunrise, double *sunset, int hour, double maxTemp,
                    double minTemp, float dayLength)
    {

      float hrsAfterMinTempUntilSunset;
      float hrsAfterSunsetUntilMinTemp;
      float nightLength;
      float sunsetTemp;
      float temp;
      float ddy;

      nightLength = 24 - dayLength;

      /* Determine if the hour is during the day or night */
      *sunrise = 12 - (dayLength / 2) + TIME_LAG_MIN;
      *sunset = 12 + (dayLength / 2);
      if ((hour > *sunrise) && (hour < *sunset)) {
        /* Calculate temperature for a day time hour. */
        hrsAfterMinTempUntilSunset = hour - (float)*sunrise;
        temp = (float)((maxTemp - minTemp) *
                        sin((3.14*hrsAfterMinTempUntilSunset) /
                       (dayLength+(2*TIME_LAG_MAX))) + minTemp);
      } else {
        /* Calculate temperature for a nght time hour. */
        if (hour > (float)*sunset) {
          hrsAfterSunsetUntilMinTemp = hour - (float)*sunset;
        }
        if (hour < (float)*sunrise) {
          hrsAfterSunsetUntilMinTemp = (24 - (float)*sunset) + hour;
        }
        ddy = dayLength - TIME_LAG_MIN;
        sunsetTemp = (float)((maxTemp - minTemp) * sin((3.14*ddy) /
                             (dayLength+(2*TIME_LAG_MAX))) + minTemp);
        temp = (float)(minTemp + (sunsetTemp - minTemp) *
                       exp(-(B * hrsAfterSunsetUntilMinTemp) / nightLength));
      }

      return temp;
    }  /* end of hourTemp() */


/*****************************************************************************
**
**  FUNCTION:  void calcTemperatureEff()
**
**  PURPOSE:   Calculate decrease in photosynthesis due to temperature.
**
**  AUTHOR:    Adapted from code extracted from the Sipnet model
**             (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    dTemp - decrease in photosynthesis due to temperature
**    tair  - average air temperature (degrees C)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    psparams          - parameter values read from the photosyn.in file
**    psparams->psnTMax - maximum temperature at which net photosynthesis
**                        occurs, assumed symmetrical around psnTOpt
**                        (degrees C)
**    psparams->psnTMin - minimum tempature at which net photosynthesis occurs
**                        (degrees C)
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    dTemp - decrease in photosynthesis due to temperature
**
**  CALLED BY:
**    calcPsnEffects()
**
**  CALLS:
**    None
**
*****************************************************************************/
    void calcTemperatureEff(double tair, double *dTemp)
    {

      extern PSPARAMS_SPT psparams;

      *dTemp = (psparams->psnTMax - tair)*(tair - psparams->psnTMin)/
              pow((psparams->psnTMax - psparams->psnTMin)/2.0, 2);
      if (*dTemp < 0) {
        *dTemp = 0.0;
      }
    }


/*****************************************************************************
**
**  FUNCTION:  void calcVpdEff()
**
**  PURPOSE:   Calculate decrease in photosynthesis due to vapor pressure
**             deficit.
**
**  AUTHOR:    Adapted from code extracted from the Sipnet model
**             (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    dVpd - decrease in photosynthesis due to vapor pressure deficit
**    vpd  - average vapor pressure deficit (kPa)
**
**  GLOBAL VARIABLES:
**    TINY - to avoid divide-by-zero errors (0.000001)
**
**  EXTERNAL VARIABLES:
**    psparams            - parameter values read from the photosyn.in file
**    psparams->dVpdSlope - slope value in dVpd equation
**    psparams->dVpdExp   - exponential value in dVpd equation
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    dVpd - decrease in photosynthesis due to vapor pressure deficit
**
**  CALLED BY:
**    calcPsnEffects()
**
**  CALLS:
**    None
**
*****************************************************************************/
    void calcVpdEff(double vpd, double *dVpd)
    {

      extern PSPARAMS_SPT psparams;

      /*   dVpd = 1.0 - psparams->dVpdSlope * pow(vpd, psparams->dVpdExp); */
      *dVpd = psparams->dVpdSlope * exp(psparams->dVpdExp * vpd);
      if (*dVpd < TINY) {
        *dVpd = 0;
      } else {
        *dVpd = *dVpd / psparams->dVpdSlope;
      }
    }


/*****************************************************************************
**
**  FUNCTION:  void calcMoistureEff()
**
**  PURPOSE:   Calculate factor for effect of water stress on photosynthesis
**             (factor between 0 and 1).
**
**  AUTHOR:    Adapted from code extracted from the Sipnet model
**             (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**   dWater      - effect of water stress on photosynthesis
**   etrans      - actual evapotranspiration (cm H2O)  or
**                 transpiration rate as a function of potential
**                 evapotranspiration and soil water potential 
**   potETrans   - potential evapotranspiration (cm H2O) or 1.0
**
**  GLOBAL VARIABLES:
**    TINY - to avoid divide-by-zero errors (0.000001)
**
**  EXTERNAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**   dWater - effect of water stress on photosynthesis
**
**  CALLED BY:
**    calcPsnEffects()
**
**  CALLS:
**    None
**
*****************************************************************************/
/*    void calcMoistureEff(double etrans, double *dWater, double potETrans) */
/*    {                                                                     */
/*                                                                          */
/*      if (potETrans > TINY) {                                             */
/*        if (etrans/potETrans > 1.0) {                                     */
/*          *dWater = 1.0;                                                  */
/*        } else {                                                          */
/*          *dWater = etrans/potETrans;                                     */
/*        }                                                                 */
/*      } else {                                                            */
/*        *dWater = 0.0;                                                    */
/*      }                                                                   */
/*    }  */ /* end of calcMoistureEff()                                     */


/*****************************************************************************
**
**  FUNCTION:  void calcPhotosyn()
**
**  PURPOSE:   Calculate photosynthesis value for current day.
**
**  AUTHOR:    Adapted from code extracted from the Sipnet model
**             (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    daylenght   - the length of daylight  (in days) (number of daylight
**                  hours divided by 24)
**    dTemp       - decrease in photosynthesis due to temperature
**    dVpd        - decrease in photosynthesis due to vapor pressure deficit
**    dWater      - effect of water stress on photosynthesis
**    GrossPsn    - gross photosynthesis, with water stress
**                  (g C * m^-2 ground area * day^-1)
**    grwdys      - number of days photosynthesis has occurred during the
**                  current growing season
**    lai         - leaf area
**    lightEff    - decrease in photosynthesis due to amount of light absorbed
**    potGrossPsn - potential photosynthesis, without water stress
**                  (g C * m^-2 ground area * day^-1)
**    srad        - shortwave radiation value for day (W/m^2)
**
**  GLOBAL VARIABLES:
**    TINY - to avoid divide-by-zero errors (0.000001)
**
**  EXTERNAL VARIABLES:
**    climate                   - climate parameters
**    climate->length           - length of this timestep (in days) (number of
**                                daylight hours divided by 24)
**    climate->par              - average par for this time step
**                                (Einsteins * m^-2 ground area * day^-1)
**    psparams                  - parameter values read from the photosyn.in
**                                file
**    psparams->aMax            - maximum photosynthesis
**                                (nmol CO2 * g^-1 leaf * sec^-1)
**    psparams->aMaxFrac        - average daily psparams.aMax as fraction of
**                                instantaneous
**    psparams->aMaxScalar1     - multiplier used to adjust aMax based on
**                                growthDays1 days since germination
**    psparams->aMaxScalar2     - multiplier used to adjust aMax based on
**                                growthDays2 days since germination
**    psparams->aMaxScalar3     - multiplier used to adjust aMax based on
**                                growthDays3 days since germination
**    psparams->aMaxScalar4     - multiplier used to adjust aMax based on
**                                growthDays4 days since germination
**    psparams->attenuation     - light attenuation coefficient
**    psparams->baseFolRespFrac - basal foliage respiration rate, as % of
**                                of maximum net photosynthesis rate
**    psparams->cFracLeaf       - factor for converting leaf biomass to carbon
**                                (leaf biomass * cFracLeaf = leaf carbon)
**    psparams->growthDays1     - number of days after germination to start
**                                using aMaxScalar1
**    psparams->growthDays2     - number of days after germination to start
**                                using aMaxScalar2
**    psparams->growthDays3     - number of days after germination to start
**                                using aMaxScalar3
**    psparams->growthDays4     - number of days after germination to start
**                                using aMaxScalar4
**    psparams->halfSatPar      - par at which photosynthesis occurs at 1/2
**                                theoretical maximum
**                                (Einsteins * m^-2 ground area * day^-1)
**    psparams->leafCSpWt       - grams of carbon in a square meter of leaf
**                                area
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    dWater      - effect of water stress on photosynthesis
**    GrossPsn    - gross photosynthesis, with water stress
**                  (g C * m^-2 ground area * day^-1)
**    lightEff    - decrease in photosynthesis due to amount of light absorbed
**    potGrossPsn - potential photosynthesis, without water stress
**                  (g C * m^-2 ground area * day^-1)
**
**  CALLED BY:
**    potprod()
**
**  CALLS:
**    potPsn()
**
*****************************************************************************/
    void calcPhotosyn(double *srad, double *daylength, double *lai,
                      double *potGrossPsn, double *GrossPsn, double *dTemp,
                      double *dVpd, float *dWater, double *lightEff,
                      int *grwdys, double *aMaxScalar1, double *aMaxScalar2,
                      double *aMaxScalar3, double *aMaxScalar4, double *aMax,
                      double *aMaxFrac, double *attenuation,
                      double *baseFolRespFrac, double *cFracLeaf,
                      double *halfSatPar, double *leafCSpWt,
                      double *growthDays1, double *growthDays2,
                      double *growthDays3, double *growthDays4)
    {

      extern CLIMATE_SPT climate;
      extern PSPARAMS_SPT psparams;

      climate->length = *daylength;

      psparams->aMaxScalar1 = *aMaxScalar1;
      psparams->aMaxScalar2 = *aMaxScalar2;
      psparams->aMaxScalar3 = *aMaxScalar3;
      psparams->aMaxScalar4 = *aMaxScalar4;
      psparams->aMax = *aMax;
      psparams->aMaxFrac = *aMaxFrac;
      psparams->attenuation = *attenuation;
      psparams->baseFolRespFrac = *baseFolRespFrac;
      psparams->cFracLeaf = *cFracLeaf;
      psparams->halfSatPar = *halfSatPar;
      psparams->leafCSpWt = *leafCSpWt;
      psparams->growthDays1 = *growthDays1;
      psparams->growthDays2 = *growthDays2;
      psparams->growthDays3 = *growthDays3;
      psparams->growthDays4 = *growthDays4;

      /* convert srad (W/m^2) to par (W/m^2) */
      climate->par = *srad * 0.5;
      /* 1 W/m^2 = 1.8 micromole/m^2/sec */
      climate->par *= 1.8;
      /* convert micromoles to moles */
      climate->par /= 1E6;
      /* convert from Einsteins/m^2/sec to Einsteins/m^2/day */ 
      climate->par *= 86400;
      climate->par *= (1.0/climate->length);
   
      /* gross photosynthesis without water effect */
      potPsn(potGrossPsn, lai, climate->par, dTemp, dVpd, lightEff, grwdys);

      if (*potGrossPsn < TINY) {
        /* we don't have any photosynthesis */
        *GrossPsn = *potGrossPsn;  
      } else {
        /* gross phtosynthesis with water stress */
        *GrossPsn = *potGrossPsn * *dWater;
      }

    } /* end of calcPhotosyn() */


/*****************************************************************************
**
**  FUNCTION:  void potPsn()
**
**  PURPOSE:   Calculate gross photosynthesis without water effect
**             (g C * m^-2 ground area * day^-1).
**
**  AUTHOR:    Adapted from code extracted from the Sipnet model
**             (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    dTemp       - decrease in photosynthesis due to temperature
**    dVpd        - decrease in photosynthesis due to vapor pressure deficit
**    grwdys      - number of days photosynthesis has occurred during the
**                  current growing season
**    lai         - leaf area
**    lightEff    - decrease in photosynthesis due to amount of light absorbed
**    par         - average par for this time step
**                  (Einsteins * m^-2 ground area * day^-1)
**    potGrossPsn - potential photosynthesis, without water stress
**
**  GLOBAL VARIABLES:
**    C_WEIGHT    - molecular weight of carbon (12.0)
**    SEC_PER_DAY - number of seconds in 24 hours (86400.0)
**    TEN_9       - for conversions from nano (1000000000.0)
**
**  EXTERNAL VARIABLES:
**    psparams                  - parameter values read from the photosyn.in
**                                file
**    psparams->aMax            - maximum photosynthesis
**                                (nmol CO2 * g^-1 leaf * sec^-1)
**    psparams->aMaxFrac        - average daily psparams.aMax as fraction of
**                                instantaneous
**    psparams->aMaxScalar1     - multiplier used to adjust aMax based on
**                                growthDays1 days since germination
**    psparams->aMaxScalar2     - multiplier used to adjust aMax based on
**                                growthDays2 days since germination
**    psparams->aMaxScalar3     - multiplier used to adjust aMax based on
**                                growthDays3 days since germination
**    psparams->aMaxScalar4     - multiplier used to adjust aMax based on
**                                growthDays4 days since germination
**    psparams->baseFolRespFrac - basal foliage respiration rate, as % of
**                                of maximum net photosynthesis rate
**    psparams->cFracLeaf       - factor for converting leaf biomass to carbon
**                                (leaf biomass * cFracLeaf = leaf carbon)
**    psparams->growthDays1     - number of days after germination to start
**                                using aMaxScalar1
**    psparams->growthDays2     - number of days after germination to start
**                                using aMaxScalar2
**    psparams->growthDays3     - number of days after germination to start
**                                using aMaxScalar3
**    psparams->growthDays4     - number of days after germination to start
**                                using aMaxScalar4
**    psparams->leafCSpWt       - grams of carbon in a square meter of leaf
**                                area
**
**  LOCAL VARIABLES:
**    conversion  - convert from (nmol CO2 * g^-1 leaf * sec^-1) to
**                  (g C * m^-2 ground area * day^-1)
**    grossAMax   - maximum possible gross respiration
**                  (nmol CO2 * g^-1 leaf * sec^-1)
**    respPerGram - base foliar respiration (nmol CO2 * g^-1 leaf * sec^-1)
**
**  OUTPUTS:
**    lightEff    - decrease in photosynthesis due to amount of light absorbed
**    potGrossPsn - potential photosynthesis, without water stress
**                  (g C * m^-2 ground area * day^-1)
**
**  CALLED BY:
**    calcPhotosyn()
**
**  CALLS:
**    calcLightEff3()
**
*****************************************************************************/
    void potPsn(double *potGrossPsn, double *lai, double par, double *dTemp,
                double *dVpd, double *lightEff, int *grwdys)
    {

      extern PSPARAMS_SPT psparams;

      double conversion;
      double grossAMax;
      double respPerGram;
      double scalar;

      /* Add a multiplier to aMax that varies depending on how many */
      /* days of growth there have been in the current growing season, */
      /* cak - 11/04/2009 */
      if (*grwdys < psparams->growthDays2) {
        scalar = (psparams->aMaxScalar2 - psparams->aMaxScalar1) /
                 (psparams->growthDays2 - psparams->growthDays1) *
                 (*grwdys - psparams->growthDays2) + psparams->aMaxScalar2;
      } else if (*grwdys >= psparams->growthDays2 &&
                 *grwdys < psparams->growthDays3) {
        scalar = (psparams->aMaxScalar3 - psparams->aMaxScalar2) /
                 (psparams->growthDays3 - psparams->growthDays2) *
                 (*grwdys - psparams->growthDays3) + psparams->aMaxScalar3;
      } else if (*grwdys >= psparams->growthDays3 &&
                 *grwdys < psparams->growthDays4) {
        scalar = (psparams->aMaxScalar4 - psparams->aMaxScalar3) /
                 (psparams->growthDays4 - psparams->growthDays3) *
                 (*grwdys - psparams->growthDays4) + psparams->aMaxScalar4;
      } else {
        scalar = psparams->aMaxScalar4;
      }

      respPerGram = psparams->baseFolRespFrac * (psparams->aMax * scalar);
      /* foliar respiration, unmodified by temp, etc. */
      grossAMax = (psparams->aMax * scalar) * psparams->aMaxFrac +
                  respPerGram;

      /* calculate decrease in photosynthesis due to amount of light */
      /* absorbed */
      calcLightEff3(lightEff, *lai, par);

      /* to convert units */
      conversion = C_WEIGHT * (1.0/TEN_9) *
                   (psparams->leafCSpWt/psparams->cFracLeaf) * *lai *
                   SEC_PER_DAY;
      *potGrossPsn = grossAMax * *dTemp * *dVpd * *lightEff * conversion;
    }  /* end of potPsn() */


/*****************************************************************************
**
**  FUNCTION:  void calcLightEff3()
**
**  PURPOSE:   Another method to calculate amount of light absorbed using
**             light attenuation (as in PnET).
**             NOTE:  Use Simpson's method to approximate the integral.  We
**                    are essentially integrating over the canopy, from top to
**                    bottom, but it's an ugly integral so we'll approximate
**                    it numerically.  Here we use Simpson's method to
**                    approximate the integral so we need an odd number of
**                    points.  This means that NUM_LAYERS must be EVEN because
**                    we loop from layer = 0 to layer = NUM_LAYERS.  As a
**                    reminder, Simpson's rule approximates the integral as:
**                      (h/3) * (y(0) + 4y(1) + 2y(2) + 4y(3) + ... +
**                      2y(n-2) + 4y(n-1) + y(n)),
**                    where h is the distance between each x value (here
**                    h = 1/NUM_LAYERS).  We keep track of the running sum in
**                    cumLightEff.  Information on the distribution of LAI
**                    with height is available as of March 2007.
**                    Contact:  Dr. Maggie Prater, Maggie.Prater@colorado.edu
**
**  AUTHOR:    Adapted from code extracted from the Sipnet model
**             (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    lai      - leaf area
**    lightEff - average light effect
**    par      - average par for this time step
**               (Einsteins * m^-2 ground area * day^-1)
**
**  GLOBAL VARIABLES:
**    NUM_LAYERS - number of canopy layers used in the calculation of light
**                 absorbed (6)
**                 NOTE:  6 layers gives approximately the same result as 100
**                        layers
**
**  EXTERNAL VARIABLES:
**    climate               - climate parameters
**    climate->par          - average par for this time step
**                            (Einsteins * m^-2 ground area * day^-1)
**    psparams              - parameter values read from the photosyn.in file
**    psparams->attenuation - light attenuation coefficient
**    psparams->halfSatPar  - par at which photosynthesis occurs at 1/2
**                            theoretical maximum
**                            (Einsteins * m^-2 ground area * day^-1)
**
**  LOCAL VARIABLES:
**    coeff          - current coefficient in Simpson's rule
**    cumLai         - lai from current layer up
**    cumLightEff    - light effect the running sum
**    currLightEff   - light effect for the current layer
**    layer          - counter
**    lightIntensity - light intensity for current layer
**
**  OUTPUTS:
**    lightEff - average light effect
**
**  CALLED BY:
**    potPsn()
**
**  CALLS:
**    None
**
*****************************************************************************/
    void calcLightEff3 (double *lightEff, double lai, double par)
    {

      extern CLIMATE_SPT climate;
      extern PSPARAMS_SPT psparams;

      int layer;
      double cumLai;
      double lightIntensity;
      double currLightEff, cumLightEff;
      int coeff;

      if (lai > 0 && par > 0) {
        /* must have at least some leaves and some light */
        cumLai = 0.0;
        cumLightEff = 0.0;
        layer = 0;
        coeff = 1;

        while (layer <= NUM_LAYERS) {
          /* lai from this layer up */
          cumLai = lai * ((double)layer / NUM_LAYERS);
          /* between 0 and par */
          lightIntensity = climate->par * exp(-1.0 * psparams->attenuation * cumLai);
           /* between 0 and 1 */
          currLightEff = (1 - pow(2, (-1.0 * lightIntensity/psparams->halfSatPar)));
          /* when lightIntensity = halfSatPar, currLightEff = 1/2 */
          cumLightEff += coeff * currLightEff;
 
          /* now move to the next layer */
          layer++;
          /* coeff. goes 1, 4, 2, 4, ..., 2, 4, 2 */
          coeff = 2*(1 + layer%2);
        }
        /* last value should have had a coefficient of 1, but actually had */
        /* a coefficient of 2, so subtract 1 */
        cumLightEff -= currLightEff;
 
        /* multiplying by (h/3) in Simpson's rule */
        *lightEff = cumLightEff/(3.0*NUM_LAYERS);
      } else {
        /* no leaves or no light! */
        *lightEff = 0;
      }
    }  /* end of calcLightEff3() */
