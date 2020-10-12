
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      trwtavg.c
**
**  FUNCTION:  float trwtavg()
**
**  PURPOSE:   Return the soil water effect on transpiration (h2ogef). This 
**             value is also the soil water effect on photosynthesis for the
**             photosynthesis version of DayCent that calls calcPsnEffects.
**
**  AUTHOR:    Susan Chaffee  4/30/92
**
**  REWRITE:   Melannie Hartman  9/20/93
**
**  HISTORY:
**    11/09/95 (MDH) Added PILPS modifications.
**    05/19/08 (CAK) Modify subroutine to remove the code that is using
**                   shallow, intermediate, deep, and very deep soil depths to
**                   calculate a weighted average value for calculating
**                   transpiration.  Replace this code with a routine that
**                   will return the soil water potential of the wettest soil
**                   layer within the plant rooting zone to be used for
**                   calculating transpiration.
**    07/25/11 (CAK) Return a 0.0 to 1.0 multiplier instead of the soil water
**                   potential of the wettest soil layer within the plant
**                   rooting zone.
**    01/22/14 (CAK) Calculate the water effect of transpiration using the
**                   same equation as used when calculating the water effect
**                   on potential production, the coefficients for this
**                   calculation are read from the sitepar.in file.
**    05/07/14 (MDH) Comment: the water effect on transpiration is a 
**                   multiplier 0-1, and is no longer a soil water potential. 
**                   I updated the comments in this file to reflect this change.
**
**  INPUTS:
**    flags              - structure containing debugging flags
**    flags->debug       - flag to set debugging mode, 0 = off, 1 = on
**    layers             - soil water soil layer structure
**    layers->fieldc[]   - volumetric water content at field capacity for
**                         layer (cm H2O/cm of soil)
**    layers->lbnd[]     - the index of the lower soil water model layer which
**                         corresponds clyr in Century
**    layers->swc[]      - soil water content by layer (cm H2O)
**    layers->width[]    - the thickness of soil water model layers (cm)
**    layers->wiltpt[]   - volumetric soil water content at wilting point for
**                         a layer (cm H2O/cm of soil)
**    nlaypg             - number of layers in plant rooting zone
**    sitepar->th2ocoef1 - relative soil water content required for 50% of
**                         maximum transpiration
**    sitepar->th2ocoef2 - 4 times the slope at relative soil water content
**                         required for 50% of maximum transpiration
**
**  GLOBAL VARIABLES:
**
**  LOCAL VARIABLES:
**    callname - call name for subroutine
**    h2ogef   - soil water effect on transpiration
**    ilyr     - current layer in the soil profile
**    maxrwcf  - relative water content of wettest layer in plant rooting zone
**    rwcf     - relative water content by layer 
**
**  OUTPUTS:
**    h2ogef - soil water effect on transpiration, 0.0 to 1.0 multiplier
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include "soilwater.h"

    float trwtavg(int *nlaypg, LAYERPAR_SPT layers, FLAG_SPT flags,
                  SITEPAR_SPT sitepar)
    {
      int ilyr;
      double maxrwcf;
      float h2ogef;
      double rwcf;
      /*static char *callname = "trwtavg"; */

      if (flags->debug > 2) 
      {
        printf("Entering function trwtavg\n");
      }

      /* Find the wettest layer in the plant rooting zone */
      maxrwcf = 0.0f;
      for(ilyr=0; ilyr <= layers->lbnd[*nlaypg-1]; ilyr++) 
      {
        rwcf = (layers->swc[ilyr]/layers->width[ilyr] - layers->wiltpt[ilyr])/
               (layers->fieldc[ilyr] - layers->wiltpt[ilyr]);
        rwcf = max(0.0, rwcf);
        if (rwcf > maxrwcf) 
        {
          maxrwcf = rwcf;
        }
      }
      maxrwcf = min(1.0, maxrwcf);

      /* Calculate the soil water effect on transpiration */
      h2ogef = (float)(1.0/(1.0+exp(sitepar->th2ocoef2*
                                    (sitepar->th2ocoef1-maxrwcf))));

      if (flags->debug > 2) 
      {
        printf("Exiting function trwtavg\n");
      }

      return(h2ogef);
    }
