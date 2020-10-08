
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtmethane.c
**
**  FUNCTION:  void wrtmethane()
**
**  PURPOSE:   Write out the methane production, emission, and oxidation
**             values. 
**
**  AUTHOR:    Cindy Keough 02/2012
** 
**  INPUTS:
**    aglivc      - above ground live carbon (g/m^2)
**    avgst_10cm  - average soil temperature in top 10 cm of soil profile
**                  (degrees C)
**    bglivcj     - below ground juvenile fine root live carbon (g/m^2)
**    bglivcm     - below ground mature fine root live carbon (g/m^2)
**    CH4_Ep      - methane emitted via plants (g/m^2)
**    CH4_Ebl     - methane emitted via bubbles (g/m^2)
**    CH4_oxid    - methane oxidation (g/m^2)
**    CH4_prod    - methane production (g/m^2)
**    Com         - sum of CO2 losses from heterotrophic decomposition of
**                  metabc(1), metabc(2), strucc(1), strucc(2), som1c(1),
**                  som1c(2), som2c(1), som2c(2), and som3c (g/m^2)
**    Cr          - carbohydrates derived from rice plants (g/m^2)
**    curday      - the day of the year (1..366)
**    Eh          - effect of water management on soil redox potential (Eh)
**                  (-250.0 mv - 300.0 mv)
**    Feh         - reduction factor of effect of soil redox potential (Eh) on
**                  methane production (0.0 - 1.0)
**    irri        - irrigation for current day (cm/H2O)
**    ppt         - precipitation for current day (cm/H2O)
**    prev_mcprd1 - previous day's aboveground live carbon production (g/m^2)
**    prev_mcprd2 - previous day's juvenile fine root carbon production
**                  (g/m^2)
**    prev_mcprd3 - previous day's mature fine root carbon production (g/m^2)
**    SI          - soil texture index for methane production (0.0 - 1.0)
**    time        - simulation time (years)
**    TI          - soil temperature index for methane production (0.0 - 1.0)
**    watr2sat    - amount of water automatically added to bring soil profile
**                  to full saturation (cm/H2O)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files                - structure containing information about output
**                           files
**    files->fp_methane    - file pointer to methane.out output file
**    files->write_methane - flag to indicate if methane.out output file
**                           should be created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**     None
**
**  CALLED BY:
**     dailymoist()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    void wrtmethane(float *time, int *curday, float *aglivc, float *bglivcj,
                    float *bglivcm, float *prev_mcprd1, float *prev_mcprd2,
                    float *prev_mcprd3, float *Com, float *ppt, float *irri,
                    float *watr2sat, float *avgst_10cm, float *TI, float *SI,
                    float *Cr, float *Eh, float *Feh, float *CH4_prod,
                    float *CH4_Ep, float *CH4_Ebl, double *CH4_oxid)
    {
      extern FILES_SPT files;

      if (!files->write_methane) {
        return;
      }

      fprintf(files->fp_methane, "%6.2f  %2d  %12.4f  %12.4f  %12.4f  ",
              *time, *curday, *aglivc, *bglivcj, *bglivcm);
      fprintf(files->fp_methane, "%12.4f  %12.4f  %12.4f  ",
              *prev_mcprd1, *prev_mcprd2, *prev_mcprd3);
      fprintf(files->fp_methane, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ", 
              *Com, *ppt, *irri, *watr2sat, *avgst_10cm);
      fprintf(files->fp_methane, "%12.4f  %12.4f  %12.6f  %12.4f  %12.4f  ",
              *TI, *SI, *Cr, *Eh, *Feh);
      fprintf(files->fp_methane, "%12.6f  %12.6f  %12.6f  %12.6f\n",
              *CH4_prod, *CH4_Ep, *CH4_Ebl, *CH4_oxid);
      return;
    }
