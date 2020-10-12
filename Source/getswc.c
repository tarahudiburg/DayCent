
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**  
**  FILE:      getswc.c
**
**  FUNCTION:  void getswc()
**
**  PURPOSE:   Pass the soil water content by layer to the calling routine.
**
**  AUTHOR:    Cindy Keough  07/20/2011 
**
**  INPUTS:
**    swc[]   - soil water content by layer (cm H2O)
**    numlyrs - total number of layers in the soil water model soil profile
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ilyr - current layer in the soil profile
**
**  CALLED BY:
**    wrtextsite()
**
**  CALLS
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"


    void getswc(float swcextend[], int *numlyrs)
    {
      int ilyr;

      extern LAYERPAR_SPT layers;

      for(ilyr=0; ilyr < *numlyrs; ilyr++) {
        swcextend[ilyr] = (float)layers->swc[ilyr];
      }

      return;
    }
