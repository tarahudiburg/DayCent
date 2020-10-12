/*****************************************************************************
**
**  stmtemp.h
**
**  Soil temperature submodel output variables
**
**  MXSWLYR  - maximum number of soil water model layers
**  MAXSTLYR - maximum number of 2 centimeter layers for the soil
**             temperature model (200)
**
**  Melannie Hartman
**  11/95
**
*****************************************************************************/

#include "swconst.h"

void therm(int nlyrs, float width[MXSWLYR], float depth[MXSWLYR], 
           float bulkden[MXSWLYR], float fieldc[MXSWLYR], double swc[MXSWLYR], 
           int nd, float stmtemp[MAXSTLYR], float tdif[MAXSTLYR], 
           float asand[MXSWLYR], float aclay[MXSWLYR], float aorg[MXSWLYR],
           float tmin, float tmax, float dx);
