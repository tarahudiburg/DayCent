
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initsite_tg.c
**
**  FUNCTION:  void initsite()
**
**  PURPOSE:   Read in site specific parameters (from sitepar.in)
**
**  AUTHOR:    Susan Chaffee    March 10, 1992
**
**  REWRITE:   Melannie Hartman  9/9/93 - 9/23/93
**
**  HISTORY:
**    8/13/92 (SLC) - Change the way the parameter for minimum soil water
**                    content is used.  No longer a function of wilting point,
**                    now it simply gives the percent of water at each layer.
**    5/08/08 (CAK) - Add checks for end of file, if the end of the file is
**                    reached before all of the variables have been read and
**                    initialized this will result in a fatal error.
**
**  INPUTS:
**    flags          - structure containing debugging flags
**    flags->debug   - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose - flag to set verbose debugging mode, 0 = off, 1 = on
**    sitename       - data file name containing site parameters
**    sitepar        - site specific parameters structure for soil water model
**    layers         - soil water soil layer structure
**
**  GLOBAL VARIABLES:
**    INPTSTRLEN - maximum length of input file line (120)
**    NTDEPTHS   - maximum number of soil regions (4)
**
**  LOCAL VARIABLES:
**    errmsg[]   - string containing error message
**    fp_in      - pointer to input file
**    imo        - current month (1..12)
**    inptline[] - line read from input file
**    m          - month
**
**  OUTPUTS:
**    Read from file sitename:
**    bioabsorp              - litter biomass at full absorption of radiation
**                             (grams biomass)
**    maxphoto               - maximum carbon loss due to photodecomposition
**                             (ug C/KJ srad)
**    sitepar->albedo        - fraction of light reflected by snow
**    sitepar->aspect        - site aspect (degrees)
**    sitepar->Aeh           - differential coefficient for calculation of
**                             soil Eh
**    sitepar->Beh_drain     - upper-limit value of Eh during drainage
**                             course (mv)
**    sitepar->Beh_flood     - low-limit value for Eh during flooding
**                             course (mv)
**    sitepar->C6H12O6_to_CH4 - reaction of anaerobic carbohydrate
**                              fermentation with methanogenesis
**                              (mole weight)
**    sitepar->cldcov[]      - average cloud cover for the month (%, 1..100)
**    sitepar->CO2_to_CH4    - fraction of CO2 from heterotrophic soil
**                             respiration that is used to produce CH4
**                             (0.0 - 1.0)
**    sitepar->Deh           - differential coefficient for calculation of
**                             soil Eh
**    sitepar->dmp           - damping factor for calculating soil temperature
**                             by layer
**    sitepar->dmpflux       - damping factor for soil water flux (in h2oflux)
**    sitepar->drainlag      - number of days that soil drainage should lag
**                             behind rain/irrigation/melt events (1-5)
**    sitepar->ehoriz        - site east horizon, degrees
**    sitepar->elevation     - site elevation (meters)
**    sitepar->frac_to_exudates - fraction of fine root production that is
**                                root exudates
**    sitepar->fswcinit      - initial soil water content, fraction of field
**                             capacity (0.0 - 1.0)
**    sitepar->hours_rain    - the duration of the rainfall/snowmelt event
**                             (hours)
**    sitepar->hpotdeep      - hydraulic water potential of deep storage
**                             layer, the more negative the number the dryer
**                             the soil layer (units?)
**    sitepar->jdayEnd       - the day of the year to end the turning off of
**                             the restriction of the CO2 effect on
**                             denitrification
**    sitepar->jdayStart     - the day of the year to start the turning off of
**                             the restriction of the CO2 effect on
**                             denitrification
**    sitepar->ksatdeep      - saturated hydraulic conductivity of deep
**                             storage layer (cm/sec)
**    sitepar->MaxNitAmt     - maximum daily nitrification amount (gN/m^2)
**    sitepar->N2N2Oadj      - multiplier on N2/N2O ratio
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
**    sitepar->Ncoeff        - minimum water/temperature limitation
**                             coefficient for nitrification
**    sitepar->netmn_to_no3  - fraction of new net mineralization that goes to
**                             NO3 (0.0-1.0)
**    sitepar->reflec        - fraction of light reflected by vegetation
**    sitepar->slope         - site slope (degrees)
**    sitepar->SnowFlag      - effect of snow on soil surface temperature calc
**                             1 = insulating effect, 0 = no insulating effect
**    sitepar->sublimscale   - multiplier to scale sublimation
**    sitepar->tbotmn        - minimum emperature for bottom soil layer
**                             (degrees C)
**    sitepar->tbotmx        - maximum emperature for bottom soil layer
**                             (degrees C)
**    sitepar->texture       - texture classification for trace gas model
**                             (1 = coarse, 2 = medium, 3 = fine)
**    sitepar->th2ocoef1     - relative soil water content required for 50% of
**                             maximum transpiration
**    sitepar->th2ocoef2     - 4 times the slope at relative soil water
**                             content required for 50% of maximum
**                             transpiration
**    sitepar->timlag        - days from Jan 1 to coolest temp at bottom of
**                             soil (days)
**    sitepar->usexdrvrs     - 0 = use air temperature to drive PET rates
**                             1 = use extra weather drivers (solrad, rhumid,
**                                 windsp) for PET calculation
**                             2 = use extra weather drivers (srad, vpd)
**                                 for photosynthesis calculation
**                             3 = use extra drivers for both PET and
**                                 photosynthesis calculations
**    sitepar->wfpsdnitadj   - adjustment on inflection point for water filled
**                             pore space effect on denitrification curve
**    sitepar->whoriz        - site west horizon, degrees
**    sitepar->zero_root_frac - fraction of CH4 emitted via bubbles when there
**                              is zero root biomass (0.0-1.0)
**    sradadj[]              - solar radiation adjustment for cloud cover and
**                             transmission coeffient
**    tminintercept          - intercept used to adjust minimum temperature
**                             for calculating VPD at dewpoint
**    tminslope              - slope used to adjust minimum temperature
**                             for calculating VPD at dewpoint
**
**  EXAMPLE INPUT FILE:
**  0        / 0 = don't use extra weather drivers for PET 
**             1 = use extra weather drivers for PET (solrad, rhumid, windsp)
**             2 = use extra weather drivers (solrad, vpd) for photosynthesis
**             3 = use extra drivers for both PET and photosynthesis
**  1.0      / sublimscale
**  0.18     / reflec - vegetation reflectivity/albedo (frac)
**  0.65     / albedo - snow albedo (frac)
**  0.90     / fswcinit - initial swc, fraction of field capacity
**  0.000001 / dmpflux - in h2oflux routine (0.000001 = original value)
**  4        / hours_rain - duration of each rain event
**  0        / # of days between rainfall event and drainage of soil (-1=computed)
**  -200     / hpotdeep - hydraulic water potential of deep storage layer
**             (units?)
**  0.0002   / ksatdeep - saturated hydraulic conductivity of deep storage
**             layer (cm/sec)
**  1  58    / cldcov[month] - cloud cover (%)
**  2  58
**  3  58
**  4  58
**  5  58
**  6  58
**  7  58
**  8  58
**  9  58
**  10 58
**  11 58
**  12 58
**  5.0 16.4 / min and max temperature for bottom soil layer (degrees C)
**  0.003    / damping factor for calculating soil temperature by layer
**  30.0     / timlag, days from Jan 1 to coolest temp at bottom of soil (days)
**  0.03     / min water/temperature limitation coefficient for nitrify
**  0 0      / turn off resp restraint on denit between these days of the year
**  0.7      / x-inflection, WFPS where 1/2 of max % of NO3 is lost as N2O
**  2.0      / y-inflection, percentage of NO3 lost as N2O at x-inflection WFPS
**  4.0      / step between minimum and maximum percentage of NO3 lost as N2O
**  5.0      / slope, controls how fast the percent changes as a function of WFPS
**  0.4      / maximum daily nitrification amount (gN/m^2)
**  1        / snow effect on soil surface temp: 0 = not insulating, 1 = insulating
**  0.2      / fraction of new net mineralization that goes to NO3 (0.0-1.0)
**  1.0      / adjustment on inflection point for WFPS effect on denit
**  1.0      / N2/N2O ratio adjustment coefficient
**  0.378    / relative swc required for 50% of maximum transpiration
**  9.0      / 4 times the slope at relative swc required for 50% of max transp
**  2.0      / N2/N2O ratio for flooded state
**  0.5      / fraction of CO2 from soil respiration used to produce CH4
**  0.27     / reaction of anaerobic carbohydrate fermentation with methanogenesis
**             (mole weight C6H12O6 to CH4)
**  0.45     / fraction of root production that is root exudates
**  0.23     / differential coefficient (Aeh)
**  0.16     / differential coefficient (Deh)
**  -250.0   / low-limit value for Eh during flooding course (mv)
**  300.0    / upper-limit value of Eh during drainage course (mv)
**  7.0      / fraction CH4 emitted via bubbles when zero root biomass (0.0-1.0)
**  1525     / elevation, meters
**  0.0      / site slope, degrees
**  0.0      / site aspect, degrees
**  0.0      / site east horizon, degrees
**  0.0      / site west horizon, degrees
**  1   0.63 / sradadj[month] srad adjust for cloud cover & transmission coeff
**  2   0.63 
**  3   0.63 
**  4   0.63 
**  5   0.63 
**  6   0.63 
**  7   0.63 
**  8   0.63 
**  9   0.63 
**  10  0.63 
**  11  0.63
**  12  0.63
**  1.0      / slope for adjusting minimum temperature for VPD dewpoint calc
**  0.0      / intercept for adjusting minimum temperature for VPD dewpoint calc
**  1.16     / maximum carbon loss due to photodecomposition (ug C/KJ srad)
**  200.0    / litter biomass for full absorption of solar radiation (g biomass)
**
**  CALLED BY:
**    initsw()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"
#include <stdlib.h>
#include <stdio.h>

    void initsite(char *sitename, SITEPAR_SPT sitepar, LAYERPAR_SPT layers,
                  FLAG_SPT flags, float sradadj[NMONTH], float *tminslope,
                  float *tminintercept, float *maxphoto, float *bioabsorp)
    {

      int  imo, m;
      char inptline[INPTSTRLEN];
      char errmsg[INPTSTRLEN];
      FILE *fp_in;

      if (flags->debug) {
        printf("Entering function initsite\n");
      }

      if ((fp_in = fopen(sitename, "r")) == NULL) {
        sprintf(errmsg, "Cannot open file %s\n", sitename);
        perror(errmsg);
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%d", &sitepar->usexdrvrs);
        if (flags->debug) {
          printf("usexdrvrs: %1d\n", sitepar->usexdrvrs);
        }
        if (sitepar->usexdrvrs < 0 || sitepar->usexdrvrs > 3) {
          printf("ERROR:  Invalid value for usexdrvrs\n");
          printf("in sitepar.in file.  Ending simulation!\n");
          exit(1);
        }
      } else {
        printf("ERROR:  Missing extra weather drivers value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->sublimscale);
        if (flags->debug) {
          printf("sublimscale: %f\n", sitepar->sublimscale);
        }
      } else {
        printf("ERROR:  Missing sublimation scalar value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->reflec);
        if (flags->debug) {
          printf("reflec: %f\n", sitepar->reflec);
        }
      } else {
        printf("ERROR:  Missing fraction of light reflected by vegetation\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->albedo);
        if (flags->debug) {
          printf("albedo: %f\n", sitepar->albedo);
        }
      } else {
        printf("ERROR:  Missing fraction of light reflected by\n");
        printf("snow value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->fswcinit);
        if (flags->debug) {
          printf("fswcinit: %f\n", sitepar->fswcinit);
        }
      } else {
        printf("ERROR:  Missing initial soil water content value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->dmpflux);
        if (flags->debug) {
          printf("dmpflux: %f\n", sitepar->dmpflux);
        }
      } else {
        printf("ERROR:  Missing damping factor for soil water flux value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->hours_rain);
        if (flags->debug) {
          printf("hours_rain: %f\n", sitepar->hours_rain);
        }
      } else {
        printf("ERROR:  Missing duration of the rainfall/snowmelt event\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Allow user to set number of days between rainfall event and */
      /* drainage of soil profile.  If a value of -1 is entered set the */
      /* number of days to drainage based in the soil texture.  Constrain */
      /* the number of days to drainage to be <=5 to prevent numerical */
      /* instabilities in the h2oflux subroutine.  cak - 02/13/04 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%d", &sitepar->drainlag);
        if (sitepar->drainlag < 0) {
          sitepar->drainlag = sitepar->texture - 1;
        }
        if (sitepar->drainlag > 5) {
          printf("lag period for drainage too long, setting to max value\n");
          sitepar->drainlag = 5;
        }
        if (flags->debug) {
          printf("drainlag: %d\n", sitepar->drainlag);
        }
      } else {
        printf("ERROR:  Missing number of days that soil drainage should\n");
        printf("lag behind rain/irrigation/melt events value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->hpotdeep);
        if (flags->debug) {
          printf("hpotdeep: %f\n", sitepar->hpotdeep);
        }
      } else {
        printf("ERROR:  Missing hydraulic water potential of deep storage\n");
        printf("layer value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->ksatdeep);
        if (flags->debug) {
          printf("ksatdeep: %f\n", sitepar->ksatdeep);
          printf("Cloud cover:\n");
        }
      } else {
        printf("ERROR:  Missing saturated hydraulic conductivity of deep\n");
        printf("storage layer value in sitepar.in file.\n");
        printf("Ending simulation!\n");
        exit(1);
      }

      for(imo=1; imo<=12; imo++) {
        if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
          sscanf(inptline, "%d %f", &m, &sitepar->cldcov[imo]);
          if (flags->debug) {
            printf("cloud cover: %d  %6.2f\n", imo, sitepar->cldcov[imo]);
          }
        } else {
          printf("ERROR:  Missing average cloud cover value\n");
          printf("in sitepar.in file.  Ending simulation!\n");
          exit(1);
        }
      }

      /* The texture parameter is being set in the initlyrs subroutine */
      /* CAK - 05/31/01 */
      /* The texture parameter has been replaced by the minimum and */
      /* maximum temperature values for the bottom soil layer */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
/*      sscanf(inptline, "%d", &sitepar->texture); */
        sscanf(inptline, "%f %f", &sitepar->tbotmn, &sitepar->tbotmx);
        if (flags->debug) {
          printf("tbotmn: %f\n", sitepar->tbotmn);
          printf("tbotmx: %f\n", sitepar->tbotmx);
        }
        if (sitepar->tbotmn > sitepar->tbotmx) {
          fprintf(stderr, "Error in input file %s.\n", sitename);
          fprintf(stderr, "the minimum soil temperature at the bottom of\n");
          fprintf(stderr, "the soil is greater than the maximum soil\n");
          fprintf(stderr, "temperature at the bottom of the soil.\n");
          exit(1);
        }
      } else {
        printf("ERROR:  Missing minimum and maximum temperature for\n");
        printf("bottom soil layer values in sitepar.in file.\n");
        printf("Ending simulation!\n");
        exit(1);
      }

      /* Added dmp parameter to the sitepar.in file, read in the value */
      /* for this parameter, cak - 12/16/02 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->dmp);
        if (flags->debug) {
          printf("damping factor: %f\n", sitepar->dmp);
        }
      } else {
        printf("ERROR:  Missing damping factor for calculating soil\n");
        printf("temperature value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added timlag parameter to the sitepar.in file, read in the value */
      /* for this parameter, cak - 04/24/03 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->timlag);
        if (flags->debug) {
          printf("lag time: %f\n", sitepar->timlag);
        }
      } else {
        printf("ERROR:  Missing days from Jan 1 to coolest temp at bottom\n");
        printf("of soil value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added Ncoeff parameter to the sitepar.in file, read in the value */
      /* for this parameter, cak - 04/08/03 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->Ncoeff);
        if (flags->debug) {
          printf("water/temp coeff: %f\n", sitepar->Ncoeff);
        }
      } else {
        printf("ERROR:  Missing coefficient for nitrification value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%d %d", &sitepar->jdayStart, &sitepar->jdayEnd);
        if (flags->debug) {
          printf("jdayStart: %d\n", sitepar->jdayStart);
          printf("jdayEnd: %d\n", sitepar->jdayEnd);
        }
      } else {
        printf("ERROR:  Missing days of year to turn off respiration\n");
        printf("restraint on denitrification value in sitepar.in file.\n");
        printf("Ending simulation!\n");
        exit(1);
      }

      /* Added 4 coefficients to parameterize the arctangent function that */
      /* used to control the amount of NO3 from mineralization is lost as */
      /* N2O, cak - 08/28/2014 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->NO3_N2O_x);
        if (flags->debug) {
          printf("WFPS where 1/2 of max NO3 is lost as N2O: %f\n",
                 sitepar->NO3_N2O_x);
        }
      } else {
        printf("ERROR:  Missing WFPS where 1/2 of max NO3 is lost as\n");
        printf("N2O value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->NO3_N2O_y);
        if (flags->debug) {
          printf("percent NO3 lost as N2O at x-inflection WFPS: %f\n",
                 sitepar->NO3_N2O_y);
        }
      } else {
        printf("ERROR:  Missing percent NO3 lost as N2O at x-inflection\n");
        printf("WFPS in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->NO3_N2O_step);
        if (flags->debug) {
          printf("step between min and max NO3 lost as N2O: %f\n",
                 sitepar->NO3_N2O_step);
        }
      } else {
        printf("ERROR:  Missing step between min and max NO3 lost as\n");
        printf("N2O in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->NO3_N2O_slope);
        if (flags->debug) {
          printf("slope, how fast percent changes as function of WFPS: %f\n",
                 sitepar->NO3_N2O_slope);
        }
      } else {
        printf("ERROR:  Missing slope, how fast percent changes as\n");
        printf("function of WFPS in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added coefficient to control the maximum daily nitrification */
      /* amount - cak 09/21/2011 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->MaxNitAmt);
        if (flags->debug) {
          printf("Maximum daily nitrifcation amount: %f\n",
                 sitepar->MaxNitAmt);
        }
      } else {
        printf("ERROR:  Missing maximum daily nitrifcation amount\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added flag to turn off the insulating effect of snow on soil */
      /* surface temperature - cak 10/07/2011 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%d", &sitepar->SnowFlag);
        if (flags->debug) {
          printf("Snow insulating flag: %d\n",
                 sitepar->SnowFlag);
        }
      } else {
        printf("ERROR:  Missing snow insulating flag\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added coefficient to control the fraction of new net */
      /* mineralization that goes to NO3 - cak 01/09/2012 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->netmn_to_no3);
        if (flags->debug) {
          printf("Fraction of new net mineralization to NO3: %f\n",
                 sitepar->netmn_to_no3);
        }
      } else {
        printf("ERROR:  Missing fraction of new net mineralization to NO3\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added adjustment on inflection point for the water filled pore */
	  /* space effect on denitrification curve - cak 01/06/2013 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->wfpsdnitadj);
        if (flags->debug) {
          printf("Adjustment on WFPS inflection point for denit: %f\n",
                 sitepar->wfpsdnitadj);
        }
      } else {
        printf("ERROR:  Missing adjustment on WFPS inflection point for\n");
        printf("denitrification in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added N2/N2O ratio adjustment coefficient - cak 01/06/2013 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->n2n2oadj);
        if (flags->debug) {
          printf("N2/N2O ratio adjustment coefficient: %f\n",
                 sitepar->n2n2oadj);
        }
      } else {
        printf("ERROR:  Missing N2/N2O ratio adjustment coefficient\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added coefficients for relative water content on transpiration */
      /* equation, cak - 01/22/2014 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->th2ocoef1);
        if (flags->debug) {
          printf("rswc required for half of max transp: %f\n",
                 sitepar->th2ocoef1);
        }
      } else {
        printf("ERROR:  Missing rswc required for half of max transp\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->th2ocoef2);
        if (flags->debug) {
          printf("4 times slope at rwc for half of max transp: %f\n",
                 sitepar->th2ocoef2);
        }
      } else {
        printf("ERROR:  Missing 4 times slope at rwc for half of max\n");
        printf("transp value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added N2/N2O ratio that is to be used when the system is */
      /* flooded - cak 02/21/2012 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->flood_N2toN2O);
        if (flags->debug) {
          printf("N2/N2O ratio for flooded state: %f\n",
                 sitepar->flood_N2toN2O);
        }
      } else {
        printf("ERROR:  Missing N2/N2O ratio for flooded state\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added fraction of CO2 from soil respiration that is used to */
      /* produce methane */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->CO2_to_CH4);
        if (flags->debug) {
          printf(" Fraction of CO2 from heterotrophic soil respiration that");
          printf(" is used to produce CH4: %f\n", sitepar->CO2_to_CH4);
        }
      } else {
        printf("ERROR:  Missing fraction of CO2 from soil respiration\n");
        printf("that is used to produce CH4\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added coefficient to represent the reaction of anaerobic */
      /* carbohydrate fermentation with methanogenesis - cak 02/06/2012 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->C6H12O6_to_CH4);
        if (flags->debug) {
          printf(" Reaction of anaerobic carbohydrate fermentation with");
          printf(" methanogenesis: %f\n", sitepar->C6H12O6_to_CH4);
        }
      } else {
        printf("ERROR:  Missing reaction of anaerobic carbohydrate\n");
        printf("fermentation with methanogenesis\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added coefficient to represent the fraction of fine root */
      /* production that is root exudates - cak 03/22/2012 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->frac_to_exudates);
        if (flags->debug) {
          printf(" Fraction of fine root production that is root");
          printf(" exudates: %f\n", sitepar->frac_to_exudates);
        }
      } else {
        printf("ERROR:  Missing fraction of fine root production that is\n");
        printf("root exudates in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added differential coefficient for calculating Eh */
      /* Aeh - cak 02/06/2012 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->Aeh);
        if (flags->debug) {
          printf("Differential coefficient Aeh: %f\n",
                 sitepar->Aeh);
        }
      } else {
        printf("ERROR:  Missing differential coefficient Aeh\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added differential coefficient for calculating Eh */
      /* Deh - cak 02/06/2012 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->Deh);
        if (flags->debug) {
          printf("Differential coefficient Deh: %f\n",
                 sitepar->Deh);
        }
      } else {
        printf("ERROR:  Missing differential coefficient Deh\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added low-limit value for Eh during flooding course */
      /* Beh_flood - cak 02/06/2012 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->Beh_flood);
        if (flags->debug) {
          printf("low-limit value for Eh during flooding: %f\n",
                 sitepar->Beh_flood);
        }
      } else {
        printf("ERROR:  Missing low-limit value for Eh during flooding\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added upper-limit value for Eh during drainage course */
      /* Beh_drain - cak 02/06/2012 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->Beh_drain);
        if (flags->debug) {
          printf("upper-limit value for Eh during drainage: %f\n",
                 sitepar->Beh_drain);
        }
      } else {
        printf("ERROR:  Missing low-limit value for Eh during drainage\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added fraction CH4 emitted via bubbles when zero root biomass */
      /* cak 08/14/2012 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->zero_root_frac);
        if (flags->debug) {
          printf("frac CH4 emitted via bubbles when no root biomass: %f\n",
                 sitepar->zero_root_frac);
        }
      } else {
        printf("ERROR:  Missing frac CH4 emitted via bubbles when no root\n");
        printf("biomass in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added elevation (meters) to sitepar.in file, */
      /* cak - 04/15/2009 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->elevation);
        if (flags->debug) {
          printf("Elevation (meters): %f\n", sitepar->elevation);
        }
      } else {
        printf("ERROR:  Missing elevation (meters) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      /* Added site slope, aspect, east horizon, and west horizon */
      /* to sitepar.in file, */
      /* cak - 05/08/2009 */
      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->slope);
        if (flags->debug) {
          printf("Site slope (degrees): %f\n", sitepar->slope);
        }
      } else {
        printf("ERROR:  Missing site slope (degrees) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->aspect);
        if (flags->debug) {
          printf("Site aspect (degrees): %f\n", sitepar->aspect);
        }
      } else {
        printf("ERROR:  Missing site aspect (degrees) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->ehoriz);
        if (flags->debug) {
          printf("Site east horizon (degrees): %f\n", sitepar->ehoriz);
        }
      } else {
        printf("ERROR:  Missing site east horizon (degrees) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", &sitepar->whoriz);
        if (flags->debug) {
          printf("Site west horizon (degrees): %f\n", sitepar->whoriz);
        }
      } else {
        printf("ERROR:  Missing site west horizon (degrees) value\n");
        printf("in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      for(imo=1; imo<=12; imo++) {
        if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
          sscanf(inptline, "%d  %f", &m, &sradadj[imo-1]);
          if (flags->debug) {
            printf("Solar radiation adjustment factor: %d %f\n", imo,
                   sradadj[imo-1]);
          }
        } else {
          printf("ERROR:  Missing solar radiation adjustment factor\n");
          printf("value in sitepar.in file.  Ending simulation!\n");
          exit(1);
        }
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", tminslope);
        if (flags->debug) {
          printf("Slope for VPD min temp adjustment: %f\n", *tminslope);
        }
      } else {
        printf("ERROR:  Slope for VPD minimum temperature adjustment\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", tminintercept);
        if (flags->debug) {
          printf("Intercept for VPD min temp adjust: %f\n", *tminintercept);
        }
      } else {
        printf("ERROR:  Intercept for VPD minimum temperature adjustment\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", maxphoto);
        /* Convert ug C to grams C */
        *maxphoto *= 0.000001f;
        if (flags->debug) {
          printf("Maximum carbon loss due to photodecomp: %f\n", *maxphoto);
        }
      } else {
        printf("ERROR:  Missing maximum carbon loss due to photodecomp\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        sscanf(inptline, "%f", bioabsorp);
        if (flags->debug) {
          printf("Litter biomass at full absorption of radiation: %f\n",
                 *bioabsorp);
        }
      } else {
        printf("ERROR:  Missing litter biomass at full absorption of srad\n");
        printf("value in sitepar.in file.  Ending simulation!\n");
        exit(1);
      }

      if (fgets(inptline, INPTSTRLEN, fp_in) != NULL) {
        printf("ERROR:  Incorrect number of lines in sitepar.in file.\n");
        printf("Ending simulation!\n");
        exit(1);
      }

      fclose(fp_in);

      if (flags->debug) {
        printf("Exiting function initsite\n");
      }

      return;
    }
