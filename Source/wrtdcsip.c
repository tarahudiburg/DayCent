/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtdcsip.c
**
**  FUNCTION:  void wrtdcsip()
**
**  PURPOSE:   Write out the values for the SIPNET model. 
**
**  AUTHOR:    Cindy Keough  04/09
** 
**  INPUTS:
**    accum     - the amount of snow added to the snowpack (cm H2O)
**    aglivc    - above ground live carbon (g/m^2)
**    bglivcj   - juvenile fine root live carbon for crop/grass (g/m^2)
**    bglivcm   - mature fine root live carbon for crop/grass (g/m^2)
**    CO2resp   - heterotrophic CO2 respiration (g/m^2)
**    crootc    - coarse root live carbon for forest (g/m^2)
**    dayofyr   - current simulation day of the year (1..366)
**    evapdly   - water evaporated from soil (cm H2O)
**    fbrchc    - fine branch live carbon for forest (g/m^2)
**    frootcj   - juvenile fine root live carbon for forest (g/m^2)
**    frootcm   - mature fine root live carbon for forest (g/m^2)
**    intrcpt   - amount of precipitation intercepted by the standing crop and
**                litter (cm H2O)
**    mcprd1    - daily NPP for shoots for grass/crop system (gC/m^2)
**    mcprd2    - daily NPP for juvenile roots for grass/crop system (gC/m^2)
**    mcprd3    - daily NPP for mature roots for grass/crop system (gC/m^2)
**    melt      - the amount of snow melted from the snowpack, if 
**                daily air temperature is warm enough (cm H2O)
**    mfprd1    - daily NPP for live leaves for tree system (gC/m^2)
**    mfprd2    - daily NPP for live juvenile fine roots for tree system (gC/m^2)
**    mfprd6    - daily NPP for live mature fine roots for tree system (gC/m^2)
**    mfprd3    - daily NPP for live fine branches for tree system (gC/m^2)
**    mfprd4    - daily NPP for live large wood for tree system (gC/m^2)
**    mfprd5    - daily NPP for live coarse roots for tree system (gC/m^2)
**    metabc1   - carbon in metabolic component of surface litter (g/m^2)
**    metabc2   - carbon in metabolic component of soil litter (g/m^2)
**    nee       - net ecosystem exchange (npptot - CO2resp)
**    npptot    - summation of all production values (gC/m^2)
**    outflow   - water that runs off, or drains out of the profile (cm H2O)
**    petdly    - potential evapotranspiration rate for day (cm H2O)
**    ppt       - precipitation for the day (cm)
**    rleavc    - leaf live carbon for forest (g/m^2)
**    rlwodc    - large wood live carbon for forest (g/m^2)
**    runoffdly - amount of water (rain or snowmelt) which did not infiltrate
**                soil profile (cm H2O)
**    snlq      - the liquid water in the snowpack (cm H2O)
**    snow      - current snowpack (equiv. cm H2O)
**    som1c1    - carbon in surface active soil organic matter (g/m^2)
**    som1c2    - carbon in soil active soil organic matter (g/m^2)
**    som2c1    - carbon in surface slow soil organic matter (g/m^2)
**    som2c2    - carbon in soil slow soil organic matter (g/m^2)
**    som3c     - carbon in passive soil organic matter (g/m^2)
**    stdedc    - standing dead carbon (g/m^2)
**    stemp     - soil surface temperature (Celsius) 
**    strucc1   - carbon in structural component of surface litter (g/m^2)
**    strucc2   - carbon in structural component of soil litter (g/m^2)
**    sublim    - amount of water sublimated from the snowpack (cm H2O)
**    time      - simulation time (years)
**    tlai      - LAI of the tree leaves
**    totsysc   - total system carbon, summation of all live carbon, dead
**                carbon, and soil organic matter carbon pools (g/m^2)
**    trandly   - water transpired from soil (cm H2O)
**    wc_2cm    - water holding capacity of a 2 cm soil layer (cm H2O)
**    wc_3cm    - water holding capacity of a 3 cm soil layer (cm H2O)
**    wc_5cm    - water holding capacity of a 5 cm soil layer (cm H2O)
**    wc_10cm   - water holding capacity of a 10 cm soil layer (cm H2O)
**    wc_15cm   - water holding capacity of a 15 cm soil layer (cm H2O)
**    wc_30cm   - water holding capacity of a 30 cm soil layer (cm H2O)
**    wood1c    - dead fine branch carbon (g/m^2)
**    wood2c    - dead large wood carbon (g/m^2)
**    wood3c    - dead coarse root carbon (g/m^2)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files              - structure containing information about output files
**    files->fp_dcsip    - file pointer to dc_sip.out output file
**    files->write_dcsip - flag to indicate if dc_sip.out output file should
**                         be created, 0 = do not create, 1 = create
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

    void wrtdcsip(float *time, int *dayofyr, float *trandly, float *evapdly,
                  float *intrcpt, float *sublim, float *outflow,
                  float *runoffdly, float *ppt, float *accum, float *melt,
                  float *snow, float *snlq, float *petdly, float *stemp,
                  float *wc_2cm, float *wc_3cm, float *wc_5cm, float *wc_10cm,
                  float *wc_15cm, float *wc_30cm, float *CO2resp,
                  float *mcprd1, float *mcprd2, float *mcprd3, float *mfprd1,
                  float *mfprd2, float *mfprd6, float *mfprd3, float *mfprd4,
                  float *mfprd5, float *npptot, float *nee, float *aglivc,
                  float *bglivcj, float *bglivcm, float *rleavc,
                  float *frootcj, float *frootcm, float *fbrchc,
                  float *rlwodc, float *crootc, float *tlai, float *stdedc,
                  float *wood1c, float *wood2c, float *wood3c, float *strucc1,
                  float *metabc1, float *strucc2, float *metabc2,
                  float *som1c1, float *som1c2, float *som2c1, float *som2c2,
                  float *som3c, float *totsysc)
    {
      extern FILES_SPT files;

      if (!files->write_dcsip) {
        goto ex;
      }

      fprintf(files->fp_dcsip, "%.2f,%d,", *time, *dayofyr);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
              *trandly, *evapdly, *intrcpt, *sublim, *outflow, *runoffdly);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
              *ppt, *accum, *melt, *snow, *snlq, *petdly, *stemp);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,%.4f,",
              *wc_2cm,  *wc_3cm, *wc_5cm, *wc_10cm, *wc_15cm);
      fprintf(files->fp_dcsip, "%.4f,%.4f,", *wc_30cm, *CO2resp);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,%.4f,", 
              *mcprd1, *mcprd2, *mcprd3, *mfprd1, *mfprd2);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,", 
              *mfprd6, *mfprd3, *mfprd4, *mfprd5);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,%.4f,",
              *npptot, *nee, *aglivc, *bglivcj, *bglivcm);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,",
              *rleavc, *frootcj, *frootcm, *fbrchc);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
              *rlwodc, *crootc, *tlai, *stdedc, *wood1c, *wood2c);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,%.4f,",
              *wood3c, *strucc1, *metabc1, *strucc2, *metabc2);
      fprintf(files->fp_dcsip, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,\n",
              *som1c1, *som1c2, *som2c1, *som2c2, *som3c, *totsysc);

ex:   return;
    }
