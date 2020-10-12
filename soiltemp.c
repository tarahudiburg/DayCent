
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      soiltemp.c
**
**  FUNCTION:  void soiltemp()
**
**  PURPOSE:   This subroutine calculates the daily average, maximum, and
**             minimum soil temperature (t[ii][j]) for a specified number of
**             depths (nd).  The inputs include the maximum and minimum air
**             temperature (tmax and tmin), the standing plant biomass
**             (biomass), the soil temperature at the bottom of the soil
**             profile, and the thermal diffusivity of the soil. The model is
**             described in a paper by Parton(1981).
**
**  HISTORY:
**    Modified for use with the Trace Gas Model
**    Melannie Hartman
**    4/97
**
**  INPUTS:
**    bulkd[]       - bulk density by layer (g/cm3)
**    clay[]        - the fraction of clay in soil layer
**    depth[]       - the distance from the surface to the middle of the soil
**                    layer (cm)
**    diurnal_range - the difference between the maximum and minimum surface
**                    soil temperatures (deg C)
**    fieldc[]      - volumetric water content at field capacity for
**                    layer (cm H2O/cm of soil)
**    jday          - current day of the year (1..366)
**    numlyrs       - total number of layers in the soil water model soil
**                    profile
**    org[]         - the fraction of organic matter in soil layer
**    sand[]        - the fraction of sand in soil layer
**    soiltavewk    - average soil temperature in the second soils.in soil
**                    layer over the previous 7 days (deg C)
**    soiltavg[]    - average soil temperature by layer (deg C)
**    soiltmax[]    - maximum soil temperature by layer (deg C)
**    soiltmin[]    - minimum soil temperature by layer (deg C)
**    srfctemp      - soil surface temperature (deg C)
**    stmtemp[]     - the average soil temperature of the soil temperature
**                    model layers
**    swc[]         - soil water content of the soil layer (cm H2O)
**                    at time when tbotmn was observed
**    tmax          - maximum air temperature for the day (deg C - 2m)
**    tmin          - minimum air temperature for the day (deg C - 2m)
**    tmns          - minimum soil surface temperature (deg C)
**    tmxs          - maximum soil surface temperature (deg C)
**    width[]       - the thickness of soil water model layers (cm)
**
**  GLOBAL VARIABLES:
**    MXSWLYR     - maximum number of soil water model layers
**    MAXSTLYR    - maximum number of 5 centimeter layers for the soil
**                  temperature model (200)
**    PI          - pi (3.14159265)
**    SEC_PER_DAY - number of seconds in a day (86400)
**
**  EXTERNAL VARIABLES:
**    sitepar->dmp    - time step correction factor, relates to how fast the
**                      heat gets into/out of the soil
**    sitepar->tbotmn - minimum soil temperature at bottom layer (dbot) for
**                      year (degrees C)
**    sitepar->tbotmx - maximum soil temperature at bottom layer (dbot) for
**                      year (degrees C)
**    sitepar->timlag - time lag in days from the beginning of the year
**                      to the occurrence of the coolest temperature at the
**                      bottom of the soil profile
**
**  LOCAL VARIABLES:
**    a, b, c, d     - intermediate variables for calculations
**    arrayindx      - position in the circular avgtemparray array
**    asilt[]        - the fraction of silt in soil layer (0.0 - 1.0)
**    avgtemparray[] - array used to store the average soil temperatures in
**                     the 2nd soil layer over the previous 7 days
**    avtd           - the average thermal diffusivity used to calculate
**                     maximum and minimum soil temperature as a function of
**                     depth
**    dbot           - depth at bottom of the soil profile
**                     depth(numlyrs) + 5 cm (cm)
**    deltat         - the change in soil temperature between today and
**                     yesterday (deg C)
**    diff1          - intermediate variable for calculations
**    differ         - the difference between the maximum and minimum surface
**                     soil temperatures (deg C)
**    dtemp[]        - the daily change of temperature (deg c/day) by the soil
**                     temperature model soil layer
**    dummy1, dummy2 - intermediate variables for calculations
**    dx             - the depth interval for the soil temperature model soil
**                     layers (cm)
**    ii, kk, ll     - loop control variables
**    ierror         - error condition flag
**    k1, k2         - intermediate variables for calculations
**    maxdepth       - the depth of the bottom layer in the soil profile (cm)
**    m1, m2         - intermediate variables for calculations
**    nd             - the number of soil temperature layers (calculated by
**                     model)
**    ndd            - nd + 1
**    prevstavg[]    - soiltavg[] from the previous day (deg C)
**    soilenrgy      - total energy absorbed/released from the soil (cal)
**                     Negative soilenergy represents heat going into the
**                     soil.  Positive soilenergy represents heat going to
**                     warm the atmosphere.
**    startofrun     - flag to indicate the start of the run
**    t[][]          - the average ([][0]),maximum ([][1]),and minimum ([][2])
**                     soil temperature for the soil layer (deg C)
**    tdif[]         - thermal diffusivity of the soil by soil temperature
**                     model layer
**    tem1           - intermediate variable for calculations
**    vclay          - fraction of soil volume made up of clay
**    vh2oc          - fraction of soil volume made up of water
**    vmuck          - fraction of soil volume made up of organic matter
**    volheat[]      - volume heat of soil (cal/cm3/deg)*(cm3*deg)=cal
**    vsand          - fraction of soil volume made up of sand
**    vsilt          - fraction of soil volume made up of silt
**
**  OUTPUTS:
**   soiltavewk - average soil temperature in the second soils.in soil layer
**                over the previous 7 days (deg C)
**   soiltavg[] - average soil temperature by layer (deg C)
**   soiltmax[] - maximum soil temperature by layer (deg C)
**   soiltmin[] - minimum soil temperature by layer (deg C)
**   stmtemp[]  - the average soil temperature of the soil temperature model
**                layers
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    therm() - calculate thermal diffusivity
**
*****************************************************************************/

#include <math.h>
#include "soilwater.h"
#include "stmtemp.h"

    void soiltemp(int jday, float tmin, float tmax, float depth[MXSWLYR],
                  float width[MXSWLYR], float fieldc[MXSWLYR],
                  float sand[MXSWLYR], float clay[MXSWLYR], float org[MXSWLYR],
                  float bulkd[MXSWLYR], double swc[MXSWLYR], int numlyrs,
                  float soiltavg[MXSWLYR], float soiltmin[MXSWLYR],
                  float soiltmax[MXSWLYR], float stmtemp[MAXSTLYR], float tmns,
                  float tmxs, float *soiltavewk, float srfctemp,
                  float diurnal_range)
    {
      extern SITEPAR_SPT sitepar;

      float maxdepth;
      float tdif[MAXSTLYR], dtemp[MAXSTLYR], t[MXSWLYR][3];
      int   ierror, ii, ll, kk, nd, ndd;
      float differ, diff1, dx, deltat;
      float a, b, c, d, dummy1, dummy2;
      float dbot, tem1, avtd;
      float vsand, vsilt, vclay, vh2oc, vmuck;
      float asilt[MXSWLYR];
      float volheat[MXSWLYR];
      float prevstavg[MXSWLYR];
      int   k1, k2;
      float m1, m2;
      float soilenrgy;
      static int   arrayindx;
      static float avgtemparray[7];
      static int   startofrun = 1;

      /* Intialize variables */

      soilenrgy = 0.0f;

      for (ii=0; ii<MXSWLYR; ii++) {
        prevstavg[ii] = soiltavg[ii];
      }

      ierror = 0 ;

      /* specify the width of each soil temperature layer  */
      dx = 2.0f;

      /* Calculate the number of soil layers nd */
      maxdepth=0.0f;
      for(ii=0; ii<numlyrs; ii++) {
        maxdepth=maxdepth+width[ii];
      }

      dbot=maxdepth+5.0f;
      nd=(int)(dbot/dx -1.0f);
      ndd=nd+1;

      /* Calculate thermal diffusivity */
      therm(numlyrs, width, depth, bulkd, fieldc, swc, nd, stmtemp, tdif,
            sand, clay, org, tmin, tmax, dx);
 
      /* Calculate change of temperature for each depth (dtemp[kk])         */
      /* the soil depth between layers is 5cm and the time step=86400 sec.  */
      /* 86400 is the number of seconds in a day and is equal to the period */
      /* of oscillation or time step.                                       */
      tem1=sitepar->dmp*tdif[0]*SEC_PER_DAY/(dx*dx);
      dummy1=stmtemp[0]-2.0f*stmtemp[1]+stmtemp[2];
      if ((dummy1 > -.3e-12) && (dummy1 < .3e-12)) {
        dummy1=0.0f;
      }
      dtemp[0]=tem1*dummy1;

      for (kk=1; kk<nd; kk++) {
        tem1=sitepar->dmp*tdif[kk]*SEC_PER_DAY/(dx*dx);
        dummy2=stmtemp[kk]+dtemp[kk-1]-2.0f*stmtemp[kk+1]+stmtemp[kk+2];
        if ((dummy2 > -.3e-12) && (dummy2 < .3e-12)) {
          dummy2=0.0f;
        }
        dtemp[kk]=tem1*dummy2;
      }

      /* The soil surface temperature has been calculated in surftemp */
      stmtemp[0] = srfctemp;

      /* Calculate the updated value for the average soil temperature */
      for (ll=1; ll<ndd; ll++) {
        stmtemp[ll]=stmtemp[ll]+dtemp[ll-1];
        if ((stmtemp[ll] > 50) || (stmtemp[ll] < -50)) {
          printf("Problem in soiltemp - invalid numbers\n");
          printf("day = %1d, stmtemp(%1d) = %7.2f\n", jday, ll, stmtemp[ll]);
          ierror = -1;
        }
      }

      if (ierror == -1) {
/*        write(*,61)'stmtemp ',(stmtemp[ll],ll=1,30);
        write(*,61)'tdif  ',(tdif[ll],ll=1,30); */
      }

      /* Specify the soil temperature at depth(numlyrs) + 5cm         */
      /* stmtemp[nd] is set to sin function, using max and min        */
      /* dtemp(tbotmx & tbotmn) at the specified layer, and lag time  */
      /* in days from the beginning of the year to the coldest        */
      /* time period(timlag). these parameters are based on           */
      /* Pawnee 1971-79 ave monthly soil temps at 183 cm.             */
      /* Use the values for tbotmn and tbotmx as read from the        */
      /* sitepar.in file rather than hard coded values - cak 12/12/02 */
      /* Use the value for timlag as read from the sitepar.in file,   */
      /* cak - 04/24/03 */
      a=(sitepar->tbotmx-sitepar->tbotmn)/2.0f;
      b=(2.0f*(float)PI)/365.0f;
      c=((365.0f*0.75f)-sitepar->timlag)*b;
      d=(sitepar->tbotmx+sitepar->tbotmn)/2.0f;
      stmtemp[nd]=a*(float)sin((double)(b*jday+c))+d;

      /* Calculate the average,maximum and minimum soil temperature(t[ii,j])*/
      /* for the iith soil depth.  The depth is defined by depth[ii].  The  */
      /* average soil temperature is calculated by linearly interpolating   */
      /* between soil temperatures calculated at 5cm intervals.  Based on   */
      /* the fourier heat transfer equation.                                */
      avtd=0.00215f;
      for(ii=0; ii<numlyrs; ii++) {
        if(depth[ii] == 0.) {
          t[ii][1]=tmxs;
          t[ii][2]=tmns;
          t[ii][0]=(tmxs+tmns)/2.0f;
        } else {
          kk = (int)(depth[ii]/dx);
          k1 = (int)max((depth[ii]-dx/2.0f)/dx, 0.0f);
          k2 = (int)((depth[ii]+dx/2.0f)/dx);
          m1 = (dx*k2-(depth[ii]-dx/2.0f))/dx;
          m2 = (depth[ii]+dx/2.0f - dx*k2)/dx;
          t[ii][0] = (m1*stmtemp[k1] + m2*stmtemp[k2]);
 
          /* Calculate the maximum and minimum soil temperature at depth    */
          /* depth[ii].use eqn presented by parton(1983),which is a function*/
          /* of the average thermal diffusivity in the top 15cm, and the    */
          /* soil depth(depth[ii]),and the diurnal variation at the soil    */
          /* surface.                                                       */
          differ = diurnal_range;
          diff1=-depth[ii]* (float)pow((0.00005/(double)avtd),0.5);
          if (diff1 < -60.0) {
            diff1 = -60.0f;
          }
          t[ii][1]=t[ii][0]+differ*(float)exp((double)diff1)/2.0f;
          t[ii][2]=t[ii][0]-differ*(float)exp((double)diff1)/2.0f;
        }
      }

      /* Save today's values of soil temperature for use tomorrow */
      for(ii=0; ii<numlyrs; ii++) {

        /* New volume heat calculation by Bill Parton 9/94  */
        /* The volume of soil is 100cm*100cm*width[ii]      */
        asilt[ii] = 1.0f-sand[ii]-clay[ii]-org[ii];
        vsilt = bulkd[ii]*asilt[ii]/2.65f;
        vsand=bulkd[ii]*sand[ii]/2.65f;
        vclay=bulkd[ii]*clay[ii]/2.65f;
        vmuck=bulkd[ii]*org[ii]/1.30f;
        vh2oc = (float)swc[ii]/width[ii];

        /* If the temperature drops in the soil, energy is      */
        /* given of back to the atmosphere, and is therefore    */
        /* positive.                                            */
        deltat = prevstavg[ii] - t[ii][0] ;

        volheat[ii] = (0.20f*(vsand+vclay+vsilt)*2.65f + 0.30f*vmuck*1.30f +
                       vh2oc) * deltat*(width[ii]*100*100);

        soiltavg[ii] = t[ii][0];
        soiltmax[ii] = t[ii][1];
        soiltmin[ii] = t[ii][2];

        /* Negative soilenergy represents heat going into the soil          */
        /* Positive soilenergy represents heat going to warm the atmosphere */
        soilenrgy = soilenrgy + volheat[ii];
      }

      /* If it is the start of the run initialize all of the values in the */
      /* avgtemparray array using the soil temperature value calculated on */
      /* the first day of the run, cak - 04/27/04 */
      if (startofrun) {
        for (arrayindx = 0; arrayindx < 7; arrayindx++) {
          avgtemparray[arrayindx] = soiltavg[1];
        }
        arrayindx = arrayindx % 7;
        startofrun = 0;
      }

      /* Code added to compute average soil temperature for the 2nd layer */
      /* over the previous 7 days.  Store the last 7 days worth of */
      /* soil temperatures in the 2nd layer in the circular array, */
      /* cak - 04/27/04 */
      avgtemparray[arrayindx] = soiltavg[1];
      arrayindx = arrayindx + 1;
      if (arrayindx > 6) {
          arrayindx = 0;
      }

      /* Code added to compute average soil temperature for the 2nd layer */
      /* over the previous 7 days, this value is used in the scheduling of */
      /* the start of plant growth when using the growing degree day */
      /* implementation, cak - 04/27/04 */
      *soiltavewk = 0.0f;
      for (ii = 0; ii < 7; ii++) {
        *soiltavewk = *soiltavewk + avgtemparray[ii];
      }
      *soiltavewk = *soiltavewk / 7.0f;

      return;
    }
