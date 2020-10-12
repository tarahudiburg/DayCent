/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtharvest.c
**
**  FUNCTION:  void wrtharvest()
**
**  PURPOSE:   Write out the state of the system at a harvest event. 
**
**  AUTHOR:    Cindy Keough  08/10
** 
**  INPUTS:
**    accrst       - annual accumulator of carbon in straw removed during
**                   harvest (sum of crmvst) (g/m^2)
**    accrste1     - annual accumulator of nitrogen in straw removed during
**                   harvest (sum of ermvst(N)) (g/m^2)
**    accrste2     - annual accumulator of phosphorus in straw removed during
**                   harvest (sum of ermvst(P)) (g/m^2)
**    accrste3     - annual accumulator of sulfur in straw removed during
**                   harvest (sum of ermvst(P)) (g/m^2)
**    agcacc       - growing season accumulator for aboveground carbon
**                   production, reset to 0.0 on LAST event (g/m^2)
**    addsdc       - amount of carbon transferred from aboveground live carbon
**                   pool (aglivc) to standing dead carbon pool (stdedc) due
**                   to a grain harvest event (g/m^2)
**    addsde1      - amount of nitrogen transferred from aboveground live
**                   nitrogen pool (aglive(1)) to standing dead nitrogen pool
**                   (stdede(1)) due to a grain harvest event (g/m^2) 
**    addsde2      - amount of phosphorus transferred from aboveground live
**                   phosphorus pool (aglive(2)) to standing dead phosphorus
**                   pool (stdede(2)) due to a grain harvest event (g/m^2) 
**    addsde3      - amount of sulfur transferred from aboveground live sulfur
**                   pool (aglive(3)) to standing dead sulfur pool (stdede(3))
**                   due to a grain harvest event (g/m^2) 
**    bgcjacc      - growing season accumulator for juvenile fine root carbon
**                   production, reset to 0.0 on LAST event (g/m^2)
**    bgcmacc      - growing season accumulator for mature fine root carbon
**                   production, reset to 0.0 on LAST event (g/m^2)
**    cgracc       - annual accumulator for carbon in harvested grain and
**                   tubers (g/m^2)
**    cgrain       - amount of carbon in harvested grain and tubers (g/m^2)
**    crmvst       - amount of carbon removed as straw during harvest (sum of
**                   cstraw and stdstraw) (g/m^2)
**    crpval       - numerical representation of the current crop
**    cstraw       - amount of carbon removed from aboveground live pool as
**                   straw during harvest (g/m^2)
**    ctubesj      - amount of carbon removed from juvenile fine root carbon
**                   pool (bglivcj) as tubers (g/m^2)
**    ctubesm      - amount of carbon removed from mature fine root carbon
**                   pool (bglivcm) as tubers (g/m^2)
**    dayofyr      - current simulation day of the year (1..366)
**    egracc1      - annual accumulator of nitrogen in harvested grain and
**                   tubers (g/m^2)
**    egracc2      - annual accumulator of phosphorus in harvested grain and
**                   tubers (g/m^2)
**    egracc3      - annual accumulator of sulfur in harvested grain and
**                   tubers (g/m^2)
**    egrain1      - amount of nitrogen in harvested grain and tubers (g/m^2)
**    egrain2      - amount of phosphorus in harvested grain and tubers
**                   (g/m^2)
**    egrain3      - amount of sulfur in harvested grain and tubers (g/m^2)
**    ermvst1      - amount of nitrogen in straw removed during harvest (sum
**                   of estraw1 and estdstraw1) (g/m^2)
**    ermvst2      - amount of phosphorus in straw removed during harvest (sum
**                   of estraw2 and estdstraw2) (g/m^2)
**    ermvst3      - amount of sulfur in straw removed during harvest (sum of
**                   estraw3 and estdstraw3) (g/m^2)
**    esoillittrj1 - amount of dead juvenile fine root nitrogen transferred to
**                   soil litter pool (metabe(2,1) and struce(2,1)) due to a
**                   harvest event (g/m^2)
**    esoillittrj2 - amount of dead juvenile fine root phosphorus transferred
**                   to soil litter pool (metabe(2,2) and struce(2,2)) due to
**                   a harvest event (g/m^2)
**    esoillittrj3 - amount of dead juvenile fine root sulfur transferred to
**                   soil litter pool (metabe(2,3) and struce(2,3)) due to a
**                   harvest event (g/m^2)
**    esoillittrm1 - amount of dead mature fine root nitrogen transferred to
**                   soil litter pool (metabe(2,1) and struce(2,1)) due to a
**                   harvest event (g/m^2)
**    esoillittrm2 - amount of dead mature fine root phosphorus transferred to
**                   soil litter pool (metabe(2,2) and struce(2,2)) due to a
**                   harvest event (g/m^2)
**    esoillittrm3 - amount of dead mature fine root sulfur transferred to
**                   soil litter pool (metabe(2,3) and struce(2,3)) due to a
**                   harvest event (g/m^2)
**    esrfclittrj1 - amount of dead juvenile fine root nitrogen transferred to
**                   surface litter pool (metabe(1,1) and struce(1,1)) due to
**                   a harvest event (g/m^2)
**    esrfclittrj2 - amount of dead juvenile fine root phosphorus transferred
**                   to surface litter pool (metabe(1,2) and struce(1,2)) due
**                   to a harvest event (g/m^2)
**    esrfclittrj3 - amount of dead juvenile fine root sulfur transferred to
**                   surface litter pool (metabe(1,3) and struce(1,3)) due to
**                   a harvest event (g/m^2)
**    esrfclittrm1 - amount of dead mature fine root nitrogen transferred to
**                   surface litter pool (metabe(1,1) and struce(1,1)) due to
**                   a harvest event (g/m^2)
**    esrfclittrm2 - amount of dead mature fine root phosphorus transferred to
**                   surface litter pool (metabe(1,2) and struce(1,2)) due to
**                   a harvest event (g/m^2)
**    esrfclittrm3 - amount of dead mature fine root sulfur transferred to
**                   surface litter pool (metabe(1,3) and struce(1,3)) due to
**                   a harvest event (g/m^2)
**    estdstraw1   - amount of nitrogen in straw removed from standing dead
**                   pool during harvest (g/m^2)
**    estdstraw2   - amount of phosphorus in straw removed from standing dead
**                   pool during harvest (g/m^2)
**    estdstraw3   - amount of sulfur in straw removed from standing dead pool
**                   during harvest (g/m^2)
**    estraw1      - amount of nitrogen in straw removed from aboveground live
**                   pool during harvest (g/m^2)
**    estraw2      - amount of phosphorus in straw removed from aboveground
**                   live pool during harvest (g/m^2)
**    estraw3      - amount of sulfur in straw removed from aboveground live
**                   pool during harvest (g/m^2)
**    etubesj1     - amount of nitrogen removed from juvenile fine root
**                   nitrogen pool (bglivej(1)) as tubers (g/m^2)
**    etubesj2     - amount of phosphorus removed from juvenile fine root
**                   phosphorus pool (bglivej(2)) as tubers (g/m^2)
**    etubesj3     - amount of sulfur removed from juvenile fine root sulfur
**                   pool (bglivej(3)) as tubers (g/m^2)
**    etubesm1     - amount of nitrogen removed from mature fine root nitrogen
**                   pool (bglivem(1)) as tubers (g/m^2)
**    etubesm2     - amount of phosphorus removed from mature fine root
**                   phosphorus pool (bglivem(2)) as tubers (g/m^2)
**    etubesm3     - amount of sulfur removed from mature fine root sulfur
**                   pool (bglivem(3)) as tubers (g/m^2)
**    fertapp1     - amount of nitrogen fertilizer applied since pervious HARV
**                   event (g/m^2)
**    fertapp2     - amount of phosphorus fertilizer applied since pervious
**                   HARV event (g/m^2)
**    fertapp3     - amount of sulfur fertilizer applied since pervious HARV
**                   event (g/m^2)
**    irrapp       - amount of irrigation applied since pervious HARV event
**                   (cm H2O)
**    omadapp      - amount of carbon added to the system through organic
**                   matter addition events since the previous HARV event
**                   (g/m2)
**    omaeapp1     - amount of nitrogen added to the system through organic
**                   matter addition events since the previous HARV event
**                   (g/m2)
**    omaeapp2     - amount of for phosphorus added to the system through
**                   organic matter addition events since the previous HARV
**                   event (g/m2)
**    omaeapp3     - amount of for sulfur added to the system through organic
**                   matter addition events since the previous HARV event
**                   (g/m2)
**    resid        - amount of residue straw carbon added to the surface
**                   litter pool (metabc(1) and strucc(1)) during a grain
**                   harvest event (g/m^2)
**    reside1      - amount of residue straw nitrogen added to the surface
**                   litter pool (metabe(1,1) and struce(1,1)) during a grain
**                   harvest event (g/m^2)
**    reside2      - amount of residue straw phosphorus added to the surface
**                   litter pool (metabe(1,2) and struce(1,2)) during a grain
**                   harvest event (g/m^2)
**    reside3      - amount of residue straw sulfur added to the surface
**                   litter pool (metabe(1,3) and struce(1,3)) during a grain
**                   harvest event (g/m^2)
**    soillittrj   - amount of dead juvenile fine root carbon transferred to
**                   soil litter pool (metabc(2) and strucc(2)) due to a
**                   harvest event (g/m^2)
**    soillittrm   - amount of dead mature fine root carbon transferred to
**                   soil litter pool (metabc(2) and strucc(2)) due to a
**                   harvest event (g/m^2)
**    srfclittrj   - amount of dead juvenile fine root carbon transferred to
**                   surface litter pool (metabc(1) and strucc(1)) due to a
**                   harvest event (g/m^2)
**    srfclittrm   - amount of dead mature fine root carbon transferred to
**                   surface litter pool (metabc(1) and strucc(1)) due to a
**                   harvest event (g/m^2)
**    stdstraw     - amount of carbon removed from standing dead pool as straw
**                   during harvest (g/m^2)
**    strmac1      - annual accumulator for cm H2O in stream flow
**    strmac2      - annual accumulator for nitrogen from mineral leaching of
**                   stream flow (g/m^2)
**    strmac3      - annual accumulator for phosphorus from mineral leaching
**                   of stream flow (g/m^2)
**    strmac4      - annual accumulator for sulfur from mineral leaching of
**                   stream flow (g/m^2)
**    strmac5      - annual accumulator for carbon from organic leaching of
**                   stream flow (g/m^2)
**    strmac6      - annual accumulator for nitrogen from organic leaching of
**                   stream flow (g/m^2)
**    strmac7      - annual accumulator for phosphorus from organic leaching
**                   of stream flow (g/m^2)
**    strmac8      - annual accumulator for sulfur from organic leaching of
**                   stream flow (g/m^2)
**    time         - simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files                - structure containing information about output
**                           files
**    files->fp_harv       - file pointer to harvest.csv output file
**    files->write_harvest - flag to indicate if harvest.csv output file
**                           should be created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**     None
**
**  CALLED BY:
**     harvst()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    void wrtharvest(float *time, int *dayofyr, float *crpval, float *agcacc,
                    float *bgcjacc, float *bgcmacc, float *cgrain,
                    float *egrain1, float *egrain2, float *egrain3,
                    float *crmvst, float *ermvst1, float *ermvst2,
                    float *ermvst3, float *cstraw, float *estraw1,
                    float *estraw2, float *estraw3, float *stdstraw,
                    float *estdstraw1, float *estdstraw2, float *estdstraw3,
                    float *addsdc, float *addsde1, float *addsde2,
                    float *addsde3, float *resid, float *reside1,
                    float *reside2, float *reside3, float *irrapp,
                    float *fertapp1, float *fertapp2, float *fertapp3,
                    float *omadapp, float *omaeapp1, float *omaeapp2,
                    float *omaeapp3, float *strmac1, float *strmac2,
                    float *strmac3, float *strmac4, float *strmac5,
                    float *strmac6, float *strmac7, float *strmac8,
                    float *cgracc, float *egracc1, float *egracc2,
                    float *egracc3, float *accrst, float *accrste1,
                    float *accrste2, float *accrste3, float *ctubesj,
                    float *etubesj1, float *etubesj2, float *etubesj3,
                    float *ctubesm, float *etubesm1, float *etubesm2,
                    float *etubesm3, float *srfclittrj,
                    float *esrfclittrj1, float *esrfclittrj2,
                    float *esrfclittrj3, float *soillittrj,
                    float *esoillittrj1, float *esoillittrj2,
                    float *esoillittrj3, float *srfclittrm,
                    float *esrfclittrm1, float *esrfclittrm2,
                    float *esrfclittrm3, float *soillittrm,
                    float *esoillittrm1, float *esoillittrm2,
                    float *esoillittrm3)
    {
      extern FILES_SPT files;

      if (!files->write_harvest) {
        goto ex;
      }

      fprintf(files->fp_harv, "%.2f,%d,", *time, *dayofyr);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
              *crpval, *agcacc, *bgcjacc, *bgcmacc, *cgrain, *egrain1);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
              *egrain2, *egrain3, *crmvst, *ermvst1, *ermvst2, *ermvst3);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
              *cstraw,  *estraw1, *estraw2, *estraw3, *stdstraw, *estdstraw1);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,",
              *estdstraw2, *estdstraw3, *addsdc, *addsde1, *addsde2);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,", 
              *addsde3, *resid, *reside1, *reside2, *reside3, *irrapp);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,",
              *fertapp1, *fertapp2, *fertapp3, *omadapp, *omaeapp1);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
              *omaeapp2, *omaeapp3, *strmac1, *strmac2, *strmac3, *strmac4);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
              *strmac5, *strmac6, *strmac7, *strmac8, *cgracc, *egracc1);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,",
              *egracc2, *egracc3, *accrst, *accrste1);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,",
              *accrste2, *accrste3, *ctubesj);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,",
              *etubesj1, *etubesj2, *etubesj3, *ctubesm, *etubesm1);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,",
              *etubesm2, *etubesm3, *srfclittrj, *esrfclittrj1,
              *esrfclittrj2);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,",
              *esrfclittrj3, *soillittrj, *esoillittrj1, *esoillittrj2,
              *esoillittrj3);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f,%.4f,%.4f,",
              *srfclittrm, *esrfclittrm1, *esrfclittrm2, *esrfclittrm3, *soillittrm);
      fprintf(files->fp_harv, "%.4f,%.4f,%.4f\n",
              *esoillittrm1, *esoillittrm2, *esoillittrm3);

ex:   return;
    }
