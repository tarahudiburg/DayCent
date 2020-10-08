
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initsrad.c
**
**  FUNCTION:  void initsrad()
**
**  PURPOSE:   Initialize the solar radiation submodel.  Build a 366-day
**             arrays of ttmax0, potential radiation and day length values for
**             use in the solar radiation calculations.
**
**  AUTHOR:    Adapted from Peter Thorton's code extracted from the Sipnet
**             model (Sipnet code written by Bill Sacks at wsacks@wisc.edu)
**             Cindy Keough - 04/2009
**
**  INPUTS:
**    aspect      - site aspect (degrees)
**    ehoriz      - site east horizon (degrees)
**    elevation   - site elevation (meters)
**    latitude    - site latitude (decimal degrees)
**    slope       - site slope (degrees)
**    whoriz      - site west horizon (degrees)
**    daylength[] - length of day light (hours)
**
**  GLOBAL VARIABLES:
**    DAYSOFF   - julian day offset of winter solstice (11.25)
**    MINDECL   - minimum declination, radians (-0.4092797)
**    PI        - pi (3.14159265)
**    RADPERDAY - radians of Earth orbit per julian day (0.017214)
**    RADPERDEG - radians per degree (0.01745329)
**    SECPERRAD - seconds per radian of hour angle (13750.9871)
**    SRADDT    - timestep for radiation routine, seconds (600.0)
**    TBASE     - max inst. trans., 0m, nadir, dry atm, dimensionless (0.870)
**
**  EXTERNAL VARIABLES:
**    ctrl                     - program control variables
**    ctrl->indewpt            - dewpoint temperature input flag, 0=NO, 1=YES
**    params                   - site parameters
**    params->base_elev        - site elevation, meters
**    params->site_asp         - site aspect, degrees, (0=N,90=E,180=S,270=W)
**    params->site_ehoriz      - site east horizon, degrees
**    params->site_elev        - site elevation, meters
**    params->site_lat         - site latitude, decimal degrees, (- for south)
**    params->site_slp         - site slope, degrees
**    params->site_whoriz      - site west horizon, degrees
**    yeardata                 - year long data ararys of ttmax0, potential
**                               radiation, and day length
**    yeardata->daylength[]    - length of day light (seconds)
**    yeardata->flat_potrad[]  - daylight average flux density for a flat
**                               surface
**    yeardata->slope_potrad[] - daylight average flux density for the slope
**    yeardata->ttmax0[]       - maximum daily total transmittance
**
**  LOCAL VARIABLES:
**    am               - optical air mass
**    ami              - index for accessing value from optical airmass by
**                       degrees array, optam
**    bsg1             - first term in beam-slope geometry calculation
**    bsg2             - second term in beam-slope geometry calculation
**    bsg3             - thrid term in beam-slope geometry calculation
**    cbsa             - cosine of beam-slope angle
**    cosasp           - cosine of site aspect
**    cosdecl          - cosine of declination
**    cosegeom         - cosine of site latitude * cosine of declination
**    cosh             - cosine of hour angle
**    coshss           - daylight hours, -1 = 24 hours of daylight
**                                        1 =  0 hours of daylight
**    coslat           - cosine of site latitude
**    cosslp           - cosine of site slope
**    coszeh           - cosine of zenith angle for east horizon
**    coszwh           - cosine of zenith angle for west horizon
**    cza              - cosine of solar zenith angle
**    decl             - declination
**    dh               - hour-angle increment
**    dir_beam_topa    - extraterrestrial radiation perpendicular to beam
**    dir_flat_topa    - potential flat surface radiation, top of atmosphere
**    dt               - timestep
**    hh               - hour-angle loop index
**    ii               - loop control variable
**    hss              - hour angle at sunset (radians)
**    lat              - Site latitude (decimal degrees)
**    optam[]          - optical airmass by degrees
**    pratio           - pressure ratio
**    sc               - solar constant as a function of yearday (W/m^2)
**    sinasp           - sine of site aspect
**    sindecl          - sine of declination
**    sinegeom         - sine of site latitude * sine of declination
**    sinh             - sine of hour angle
**    sinlat           - sine of site latitude
**    sinslp           - sine of site slope
**    sum_flat_potrad  - total potential radiation on a flat surface for ideal
**                       horizons
**    sum_slope_potrad - sun between east and west horizons, and direct on
**                       slope
**    sum_trans        - daily total transmittance
**    t1               - first term in pressure ratio calculation
**    t2               - second term in pressure ratio calculation
**    trans1           - initial transmittance corrected for elevation
**    trans2           - instantaneous transmittance corrected for optical air
**                       mass
**
**  OUTPUTS:
**    ctrl->indewpt            - dewpoint temperature input flag, 0=NO, 1=YES
**    daylenght[]              - length of day light (hours)
**    params->site_asp         - site aspect, degrees, (0=N,90=E,180=S,270=W)
**    params->site_ehoriz      - site east horizon, degrees
**    params->site_elev        - site elevation, meters
**    params->site_lat         - site latitude, decimal degrees, (- for south)
**    params->site_slp         - site slope, degrees
**    params->site_whoriz      - site west horizon, degrees
**    yeardata->daylength[]    - length of day light (seconds)
**    yeardata->flat_potrad[]  - daylight average flux density for a flat
**                               surface
**    yeardata->slope_potrad[] - daylight average flux density for the slope
**    yeardata->ttmax0[]       - maximum daily total transmittance
**
**  CALLED BY:
**    initsw()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calcSrad.h"
#include "swconst.h"

CONTROL_S Control;                   /* program control variables */
CONTROL_SPT ctrl = &Control;
PARAMETER_S Parameters;              /* site parameters */
PARAMETER_SPT params = &Parameters;
DATA_S Site_Data;                    /* site climate data */
DATA_SPT sitedata = &Site_Data;
YEARDATA_S Year_Data;                /* 366-day arrays of ttmax0, potential */
YEARDATA_SPT yeardata = &Year_Data;  /* radiation, and day length */

    void initsrad(double elevation, double latitude, double slope,
                  double aspect, double ehoriz, double whoriz,
                  float daylength[NDAY])
    {

      int    ami, ii;
      double t1,t2;
      double pratio;
      double lat,coslat,sinlat,dt,dh,hh;
      double cosslp,sinslp,cosasp,sinasp;
      double cza,cbsa,coszeh,coszwh;
      double decl,cosdecl,sindecl,cosegeom,sinegeom,coshss,hss;
      double bsg1,bsg2,bsg3;
      double sc,dir_beam_topa;
      double sum_flat_potrad, sum_slope_potrad, sum_trans;
      double cosh,sinh;
      double dir_flat_topa,am;
      double trans1,trans2;
      /* optical airmass by degrees */
      double optam[21] = {2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07, 4.37,
                          4.72, 5.12, 5.60, 6.18, 6.88, 7.77, 8.90, 10.39,
                          12.44, 15.36, 19.79, 26.96, 30.00};

      /* Initial the variables in the control structure */
      ctrl->indewpt = 0;  /* Dewpoint temperature input? (0=NO, 1=YES) */

      /* Initialize site parameters */
      params->site_lat = latitude;    /* Site latitude, degrees */
      params->site_elev = elevation;  /* Site elevation, meters */
      params->site_slp = slope;       /* Site slope, degrees */
      params->site_asp = aspect;      /* Site aspect, degrees */
      params->site_ehoriz = ehoriz;   /* Site east horizon, degrees */
      params->site_whoriz = whoriz;   /* Site west horizon, degrees */

      /* do basic error checking on parameters */
      /* check dewpoint input flag */
      if ((ctrl->indewpt < 0) || (ctrl->indewpt > 1)) {
        printf("WARNING: input dewpoint flag should be 0 or 1: assuming");
        printf(" dewpoint input ...\n");
        ctrl->indewpt = 1;
      }

      /* error checking for site parameters */
      if ((params->site_lat < -90.0) || (params->site_lat > 90.0)) {
        printf("ERROR: site latitude must be in the range -90.0 - 90.0");
        printf(" degrees ...\n");
        exit(1);
      }
      if (params->site_elev > 5000) {
        printf("WARNING: site elev = %.1lf m: be sure to use meters,",
               params->site_elev);
        printf(" not feet\n");
        exit(1);
      }
      if (params->site_slp > 60.0) {
        printf("WARNING: site slope = %.1lf deg: be sure to use deg,",
               params->site_slp);
        printf(" not %%\n");
        exit(1);
      }
      if (params->site_slp < 0) {
        printf("ERROR: site slope must be >= 0.0\n");
        exit(1);
      }
      if ((params->site_asp < 0.0) || (params->site_asp > 360.0)) {
        printf("ERROR: site aspect must be in the range 0.0 - 360.0");
        printf(" degrees ...\n");
        exit(1);
      }
      if ((params->site_ehoriz < 0.0) || (params->site_ehoriz > 180.0)) {
        printf("ERROR: site east horizon must be in the range 0.0 - 180.0");
        printf(" degrees ...\n");
        exit(1);
      }
      if ((params->site_whoriz < 0.0) || (params->site_whoriz > 180.0)) {
        printf("ERROR: site west horizon must be in the range 0.0 - 180.0");
        printf(" degrees ...\n");
        exit(1);
      }

      /* STEP (1) calculate pressure ratio (site/reference) = f(elevation) */
      t1 = 1.0 - (LR_STD * params->site_elev)/T_STD;
      t2 = G_STD / (LR_STD * (R/MA));
      pratio = pow(t1,t2);

      /* STEP (2) correct initial transmittance for elevation */ 
      trans1 = pow(TBASE,pratio);

      /* STEP (3) build 366-day array of ttmax0, potential rad, and */
      /* daylength */
      /* precalculate the transcendentals */
      lat = params->site_lat;
      /* check for (+/-) 90 degrees latitude, throws off daylength calc */
      lat *= RADPERDEG;
      if (lat > 1.5707) lat = 1.5707;
      if (lat < -1.5707) lat = -1.5707;
      coslat = cos(lat);
      sinlat = sin(lat);
      cosslp = cos(params->site_slp * RADPERDEG);
      sinslp = sin(params->site_slp * RADPERDEG);
      cosasp = cos(params->site_asp * RADPERDEG);
      sinasp = sin(params->site_asp * RADPERDEG);
      /* cosine of zenith angle for east and west horizons */
      coszeh = cos(1.570796 - (params->site_ehoriz * RADPERDEG));
      coszwh = cos(1.570796 - (params->site_whoriz * RADPERDEG));

      /* sub-daily time and angular increment information */
      dt = SRADDT;                /* set timestep */ 
      dh = dt / SECPERRAD;        /* calculate hour-angle step */

      /* begin loop through yeardays */
      for (ii=0 ; ii<365 ; ii++) {
        /* calculate cos and sin of declination */
        decl = MINDECL * cos(((double)ii + DAYSOFF) * RADPERDAY);
        cosdecl = cos(decl);
        sindecl = sin(decl);

        /* do some precalculations for beam-slope geometry (bsg) */
        bsg1 = -sinslp * sinasp * cosdecl;
        bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl;
        bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl;

        /* calculate daylength as a function of lat and decl */
        cosegeom = coslat * cosdecl;
        sinegeom = sinlat * sindecl;
        coshss = -(sinegeom) / cosegeom;
        if (coshss < -1.0) coshss = -1.0; /* 24-hr daylight */
        if (coshss > 1.0) coshss = 1.0;   /* 0-hr daylight */
        hss = acos(coshss);               /* hour angle at sunset (radians) */
        /* daylength (seconds) */
        yeardata->daylength[ii] = 2.0 * hss * SECPERRAD;
        daylength[ii] = (float)yeardata->daylength[ii] / SEC_PER_HOUR;

        /* solar constant as a function of yearday (W/m^2) */
        sc = 1368.0 + 45.5*sin((2.0*PI*(double)ii/365.25) + 1.7);
        /* extraterrestrial radiation perpendicular to beam, total over the */
        /* timestep */
        dir_beam_topa = sc * dt;

        sum_trans = 0.0;
        sum_flat_potrad = 0.0;
        sum_slope_potrad = 0.0;

        /* begin sub-daily hour-angle loop, from -hss to hss */
        for (hh=-hss ; hh<hss ; hh+=dh) {
          /* precalculate cos and sin of hour angle */
          cosh = cos(hh);
          sinh = sin(hh);

          /* calculate cosine of solar zenith angle */
          cza = cosegeom * cosh + sinegeom;

          /* calculate cosine of beam-slope angle */
          cbsa = sinh * bsg1 + cosh * bsg2 + bsg3;

          /* check if sun is above a flat horizon */
          if (cza > 0.0) {
            /* when sun is above the ideal (flat) horizon, do all the */
            /* flat-surface calculations to determine daily total */
            /* transmittance, and save flat-surface potential radiation for */
            /* later calculations of diffuse radiation */

            /* potential radiation for this time period, flat surface, top */
            /* of atmosphere */
            dir_flat_topa = dir_beam_topa * cza;

            /* determine optical air mass */
            am = 1.0/(cza + 0.0000001);
            if (am > 2.9) {
              ami = (int)(acos(cza)/RADPERDEG) - 69;
              if (ami < 0) ami = 0;
              if (ami > 20) ami = 20;
              am = optam[ami];
            }

            /* correct instantaneous transmittance for this optical air */
            /* mass */
            trans2 = pow(trans1,am);

            /* instantaneous transmittance is weighted by potential */
            /* radiation for flat surface at top of atmosphere to get daily */
            /* total transmittance */
            sum_trans += trans2 * dir_flat_topa;

            /* keep track of total potential radiation on a flat surface */
            /* for ideal horizons */
            sum_flat_potrad += dir_flat_topa;

            /* keep track of whether this time step contributes to */
            /* component 1 (direct on slope) */
            if ((hh<0.0 && cza>coszeh && cbsa>0.0) ||
                (hh>=0.0 && cza>coszwh && cbsa>0.0)) {
              /* sun between east and west horizons, and direct on slope. */
              /* this period contributes to component 1 */
              sum_slope_potrad += dir_beam_topa * cbsa;
            }

          } /* end if sun above ideal horizon */
        } /* end of sub-daily hour-angle loop */

        /* calculate maximum daily total transmittance and daylight average */
        /* flux density for a flat surface and the slope */
        if (yeardata->daylength[ii]) {
           yeardata->ttmax0[ii] = sum_trans / sum_flat_potrad;
           yeardata->flat_potrad[ii] =
             sum_flat_potrad / yeardata->daylength[ii];
           yeardata->slope_potrad[ii] =
             sum_slope_potrad / yeardata->daylength[ii];
        } else {
           yeardata->ttmax0[ii] = 0.0;
           yeardata->flat_potrad[ii] = 0.0;
           yeardata->slope_potrad[ii] = 0.0;
        }
      } /* end of ii=365 days loop */

      /* force yearday 366 = yearday 365 */
      yeardata->ttmax0[365] = yeardata->ttmax0[364];
      yeardata->flat_potrad[365] = yeardata->flat_potrad[364];
      yeardata->slope_potrad[365] = yeardata->slope_potrad[364];
      yeardata->daylength[365] = yeardata->daylength[364];
      daylength[365] = (float)yeardata->daylength[365] / SEC_PER_HOUR;
    }
