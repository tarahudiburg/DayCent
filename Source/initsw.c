
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initsw.c
**
**  FUNCTION:  void initsw()
**
**  PURPOSE:   Initialize the soil water model
**
**  INPUTS:
**    sitlat - latitude (degrees)
**
**  GLOBAL VARIABLES:
**    BAR2CM   - conversion factor for bars to centimeters H2O (1024)
**               (1 bar = 1024 cm H2O)
**    FNSITE   - file name for site specific input parameters (sitepar.in)
**    FNSOIL   - file name for soil layer structure input file (soils.in)
**    MAXLYR   - maximum number of soil water model layers (21)
**    MAXSTLYR - maximum number of 5 centimeter layers for the soil
**               temperature model (200)
**    PI       - pi (3.14159265)
**
**  EXTERNAL VARIABLES:
**    files                  - structure containing information about output
**                             files
**    flags                  - structure containing debugging flags
**    layers                 - soil water soil layer structure
**    layers->numlyrs        - total number of layers in the soil water model
**                             soil profile
**    layers->swcfc[]        - volumetric soil water content at field capacity
**                             for layer (cm H2O/cm of soil)
**    layers->swclimit[]     - minimum volumetric soil water content of a
**                             layer, fraction 0.0 - 1.0
**    layers->swcwp[]        - volumetric soil water content at wilting point
**                             for layer (cm H2O)
**    layers->width[]        - the thickness of soil water model layers (cm)
**    sitepar                - site specific parameters structure for soil
**                             water model
**    sitepar->fswcinit      - initial soil water content, fraction of field
**                             capacity (0.0 - 1.0)
**    sitepar->usexdrvrs     - 0 = use air temperature to drive PET rates
**                             1 = use extra weather drivers (solrad, rhumid,
**                                 windsp) for PET calculation, 
**                             2 = use extra weather drivers (srad, vpd)
**                                 for photosynthesis calculation
**                             3 = use extra drivers for both PET and
**                                 photosynthesis calculations
**    soil                   - soil temperature structure
**
**  LOCAL VARIABLES:
**    callname - call name for subroutine
**    errmsg[] - string containing error message
**    ilyr     - current layer in the soil profile
**    latitude - site latitude (decimal degrees)
**    lcnt     - count of number of input file lines read
**    line[]   - buffer containing line read from input file
**    MAXL     - maximum length of line read from input file
**    wrt      - flag, 0 = do not write to output file,
**               1 = do not write to output file
**
**  OUTPUTS:
**    bioabsorp             - litter biomass at full absorption of radiation
**                            (grams biomass)
**    daylength[]           - length of day light (hours)
**    files->fp_bio         - file pointer to bio.out output file
**    files->fp_cflows      - file pointer to cflows.out output file
**    files->fp_co2         - file pointer to co2.out output file
**    files->fp_daily       - file pointer to daily.out output file
**    files->fp_dcsip       - file pointer to dc_sip.csv output file
**    files->fp_deadc       - file pointer to deadc.out output file
**    files->fp_dels        - file pointer to dels.out output file
**    files->fp_dN2lyr      - file pointer to dN2lyr.out output file
**    files->fp_dN2Olyr     - file pointer to the dN2Olyr.out output file
**    files->fp_harv        - file pointer to the harvest.csv output file
**    files->fp_livec       - file pointer to livec.out output file
**    files->fp_methane     - file pointer to methane.out output file
**    files->fp_nflux       - file pointer to nflux.out output file
**    files->fp_outf        - file pointer to outfiles.in input file
**    files->fp_psyn        - file pointer to psyn.out output file
**    files->fp_resp        - file pointer to resp.out output file
**    files->fp_soilc       - file pointer to soilc.out output file
**    files->fp_soiln       - file pointer to soiln.out output file
**    files->fp_soiltavg    - file pointer to soiltavg.out output file
**    files->fp_soiltmax    - file pointer to soiltmax.out output file
**    files->fp_soiltmin    - file pointer to soiltmin.out output file
**    files->fp_stempdx     - file pointer to stemp_dx.out output file
**    files->fp_summary     - file pointer to summary.out output file
**    files->fp_swc         - file pointer to vswc.out output file
**    files->fp_sysc        - file pointer to sysc.out output file
**    files->fp_tgmonth     - file pointer to the tgmonth.out output file
**    files->fp_wb          - file pointer to watrbal.out output file
**    files->fp_wflux       - file pointer to wflux.out output file
**    files->fp_wfps        - file pointer to wfps.out output file
**    files->fp_yearsum     - file pointer to year_summary.out output file
**    files->fp_yrcflows    - file pointer to year_cflows.out output file
**    files->write_bio      - flag to indicate if bio.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_cflows   - flag to indicate if cflows.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_co2      - flag to indicate if co2.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_daily    - flag to indicate if daily.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_dcsip    - flag to indicate if dc_sip.csv output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_deadc    - flag to indicate if deadc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_dels     - flag to indicate if dels.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_dN2lyr   - flag to indicate if dN2lyr.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_dN2Olyr  - flag to indicate if dN2Olyr.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_harvest  - flag to indicate if harvest.csv output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_livec    - flag to indicate if livec.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_methane  - flag to indicate if methane.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_nflux    - flag to indicate if nflux.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_psyn     - flag to indicate if psyn.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_resp     - flag to indicate if resp.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soilc    - flag to indicate if soilc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiln    - flag to indicate if  output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soiltavg - flag to indicate if soiltavg.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmax - flag to indicate if soiltmax.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmin - flag to indicate if soiltmin.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_stempdx  - flag to indicate if stemp_dx.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_summary  - flag to indicate if summary.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_swc      - flag to indicate if vswc.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_sysc     - flag to indicate if sysc.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_tgmonth  - flag to indicate if tgmonth.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wb       - flag to indicate if watrbal.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wflux    - flag to indicate if wflux.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_wfps     - flag to indicate if wfps.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_yearsum  - flag to indicate if year_summary.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_yrcflow  - flag to indicate if year_cflows.out output file
**                            should be created, 0 = do not create, 1 = create
**    flags->debug          - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose        - flag to set verbose debugging mode, 0 = off,
**                            1 = on
**    layers->minpot[]      - minimum matric potential by layer based on
**                            swcmin (-cm)
**    layers->swc[]         - soil water content by layer (cm H2O)
**    layers->swcmin[]      - lower bound on soil water content by layer
**                            (cm H2O) swc will not be allowed to drop below
**                            this minimum
**    maxphoto              - maximum carbon loss due to photodecomposition
**                            (ug C/KJ srad)
**    numlyrs               - total number of layers in the soil water model
**                            soil profile 
**    sitepar->rlatitude    - latitude of the site (in radians)
**    soil->soiltavg[]      - average soil temperature of layer (degrees C)
**    soil->soiltmax[]      - maximum soil temperature by layer (degrees C)
**    soil->soiltmin[]      - minimum soil temperature by layer (degrees C)
**    soil->stmtemp[]       - the average soil temperature of the soil
**                            temperature model layers (degrees C)
**    sradadj[]             - solar radiation adjustment for cloud cover and
**                            transmission coeffient
**    swcinit[]             - initial soil water content by layer (cm H2O)
**    texture               - texture classification for trace gas model
**                            (1 = coarse, 2 = medium, 3 = fine)
**    tminintercept         - intercept used to adjust minimum temperature
**                            for calculating VPD at dewpoint
**    tminslope             - slope used to adjust minimum temperature
**                            for calculating VPD at dewpoint
**    usexdrvrs             - 0 = use air temperature to drive PET rates
**                            1 = use extra weather drivers (solrad, rhumid,
**                                windsp) for PET calculation, 
**                            2 = use extra weather drivers (srad, vpd)
**                                for photosynthesis calculation
**                            3 = use extra drivers for both PET and
**                                photosynthesis calculations
** 
**  CALLED BY:
**    detiv
**
**  CALLS:
**    initlyrs()  - read in layers of the soil structure for the site
**                  (from soils.in)
**    initsite()  - read in site specific parameters (from sitepar.in)
**    initsrad()  - initialize the solar radiation submodel
**    swpotentl() - given its soil water content calculate the soil water
**                  potential of a soil layer
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "soilwater.h"
#include "calcPhotosyn.h"

#define MAXL 150

FLAG_S flagstruct;
FLAG_SPT flags = &flagstruct;
LAYERPAR_S lyrstruct;
LAYERPAR_SPT layers = &lyrstruct;
SITEPAR_S sitestruct;
SITEPAR_SPT sitepar = &sitestruct;
SOIL_S soilstruct;
SOIL_SPT soil = &soilstruct;
FILES_S filestruct;
FILES_SPT files = &filestruct;
PSPARAMS_S PSParameters;
PSPARAMS_SPT psparams = &PSParameters;
CLIMATE_S ClimateVars;
CLIMATE_SPT climate = &ClimateVars;

    void initsw(float *sitlat, float swcinit[MAXLYR], int *usexdrvrs,
                int *numlyrs, int *texture, float daylength[NDAY],
                float sradadj[NMONTH], float *tminslope, float *tminintercept,
                float *maxphoto, float *bioabsorp, int *ext_flag,
                float swcextend[MAXLYR])
    {

      int wrt, ilyr, lcnt;
      static char *callname = "initsw";
      char errmsg[100], line[MAXL];
      double latitude;

      flags->debug = 1;
      flags->verbose = 1;

      initlyrs(FNSOIL, layers, flags, sitepar);
      initsite(FNSITE, sitepar, layers, flags, sradadj, tminslope,
               tminintercept, maxphoto, bioabsorp);
      latitude = *sitlat;
      initsrad(sitepar->elevation, latitude, sitepar->slope, sitepar->aspect,
               sitepar->ehoriz, sitepar->whoriz, daylength);

      sitepar->rlatitude = *sitlat * (float)(PI/180.0);
      *texture = sitepar->texture;
   
      if (flags->verbose) {
        printf("sitlat = %6.2f\n", *sitlat);
        printf("rlatitude = %6.2f\n", sitepar->rlatitude);
      }
 
      for (ilyr=0; ilyr<layers->numlyrs; ilyr++) {

        if (*ext_flag == 0) {
          layers->swc[ilyr] = sitepar->fswcinit * layers->swcfc[ilyr];
        } else {
          layers->swc[ilyr] = swcextend[ilyr];
        }
        swcinit[ilyr] = (float)layers->swc[ilyr];
   
        /* Set the lower limit on soil water content and potential. */ 

        layers->swcmin[ilyr] = min(layers->swclimit[ilyr]*layers->width[ilyr],
                                   layers->swcwp[ilyr]);
        layers->minpot[ilyr] = -swpotentl(layers->swcmin[ilyr],ilyr,layers,
                                          callname)*BAR2CM;
     
        printf("%2s  %8s  %8s  %8s\n", "ly", "swcinit", "swcmin", "minpot");
        if (flags->verbose) {
          printf("%2d  %8.4f  %8.4f  %8.4f\n", ilyr, swcinit[ilyr],
                 layers->swcmin[ilyr], layers->minpot[ilyr]);
        }
      }

      for (ilyr=0; ilyr<MAXLYR; ilyr++) {
        soil->soiltavg[ilyr] = 0.0f;
        soil->soiltmin[ilyr] = 0.0f;
        soil->soiltmax[ilyr] = 0.0f;
      }

      for (ilyr=0; ilyr<MAXSTLYR; ilyr++) {
        soil->stmtemp[ilyr] = 0.0f;
      }
   
      layers->swc[layers->numlyrs] = 0.0;
      swcinit[layers->numlyrs] = (float)layers->swc[layers->numlyrs];
      *usexdrvrs = sitepar->usexdrvrs;
      *numlyrs = layers->numlyrs;

      if ((files->fp_outf = fopen("outfiles.in", "r")) == NULL) {
        sprintf(errmsg, "Cannot open file %s\n", "outfiles.in");
        perror(errmsg);
        exit(1);
      }

      files->write_bio = 0;
      files->write_soiln = 0;
      files->write_soiltavg = 0;
      files->write_soiltmax = 0;
      files->write_soiltmin = 0;
      files->write_stempdx = 0;
      files->write_swc = 0;
      files->write_wb = 0;
      files->write_wfps = 0;
      files->write_co2 = 0;
      files->write_wflux = 0;
      files->write_resp = 0;
      files->write_yearsum = 0;
      files->write_livec = 0;
      files->write_deadc = 0;
      files->write_soilc = 0;
      files->write_sysc = 0;
      files->write_tgmonth = 0;
      files->write_dN2lyr = 0;
      files->write_dN2Olyr = 0;
      files->write_dels = 0;
      files->write_dcsip = 0;
      files->write_harvest = 0;
      files->write_cflows = 0;
      files->write_yrcflows = 0;
      files->write_daily = 0;
      files->write_nflux = 0;
      files->write_summary = 0;
      files->write_methane = 0;
      files->write_psyn = 0;

      lcnt = 0;
      while( fgets(line, MAXL, files->fp_outf) != NULL) {
        printf("%s", line);
        lcnt++;
        if (lcnt > 1) {
          sscanf(line, "%d", &wrt);
          printf("wrt = %d\n", wrt);
        }
        if (lcnt == 2)  files->write_bio = wrt;
        if (lcnt == 3)  files->write_soiln = wrt;
        if (lcnt == 4)  files->write_soiltavg = wrt;
        if (lcnt == 5)  files->write_soiltmax = wrt;
        if (lcnt == 6)  files->write_soiltmin = wrt;
        if (lcnt == 7)  files->write_stempdx = wrt;
        if (lcnt == 8)  files->write_swc = wrt;
        if (lcnt == 9)  files->write_wb = wrt;
        if (lcnt == 10) files->write_wfps = wrt;
        if (lcnt == 11) files->write_co2 = wrt;
        if (lcnt == 12) files->write_wflux = wrt;
        if (lcnt == 13) files->write_resp = wrt;
        if (lcnt == 14) files->write_yearsum = wrt;
        if (lcnt == 15) files->write_livec = wrt;
        if (lcnt == 16) files->write_deadc = wrt;
        if (lcnt == 17) files->write_soilc = wrt;
        if (lcnt == 18) files->write_sysc = wrt;
        if (lcnt == 19) files->write_tgmonth = wrt;
        if (lcnt == 20) files->write_dN2lyr = wrt;
        if (lcnt == 21) files->write_dN2Olyr = wrt;
        if (lcnt == 22) files->write_dels = wrt;
        if (lcnt == 23) files->write_dcsip = wrt;
        if (lcnt == 24) files->write_harvest = wrt;
        if (lcnt == 25) files->write_cflows = wrt;
        if (lcnt == 26) files->write_yrcflows = wrt;
        if (lcnt == 27) files->write_daily = wrt;
        if (lcnt == 28) files->write_nflux = wrt;
        if (lcnt == 29) files->write_summary = wrt;
        if (lcnt == 30) files->write_methane = wrt;
        if (lcnt == 31) files->write_psyn = wrt;
      }

      if (files->write_soiltavg) {
        files->fp_soiltavg = fopen("soiltavg.out", "w"); 
      }
      if (files->write_soiltmax) {
        files->fp_soiltmax = fopen("soiltmax.out", "w"); 
      }
      if (files->write_soiltmin) {
        files->fp_soiltmin = fopen("soiltmin.out", "w"); 
      }
      if (files->write_stempdx) {
        files->fp_stempdx = fopen("stemp_dx.out", "w");  
      }
      if (files->write_swc) {
        files->fp_swc = fopen("vswc.out", "w"); 
      }
      if (files->write_wfps) {
        files->fp_wfps = fopen("wfps.out", "w"); 
      }
    
      if (files->write_wb) {
        files->fp_wb = fopen("watrbal.out", "w"); 
        fprintf(files->fp_wb, "0=ppt+dswc-intrcpt-evap-transp-outflow");
        fprintf(files->fp_wb, " (when sublim = 0)\n");
        fprintf(files->fp_wb, "0=melt+dswc-intrcpt-evap-transp-outflow");
        fprintf(files->fp_wb, " (when sublim > 0)\n");
        fprintf(files->fp_wb, "%8s %8s %7s %7s %7s %7s %7s %7s %7s %7s %9s",
                "time", "dayofyr", "ppt", "accum", "dsnlq", "melt", "intrcpt",
                "evap", "transp", "sublim", "dswc");
        fprintf(files->fp_wb, " %7s %9s", "outflow", "balance");
        fprintf(files->fp_wb, " %7s %7s %7s \n", "snow", "snlq", "runoff");
      }

      if (files->write_soiln) {
        printf("Open soiln.out\n");
        files->fp_soiln = fopen("soiln.out", "w"); 
        fprintf(files->fp_soiln, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "ammonium", "NO3_ppm[0]", "NO3_ppm[1]",
                "NO3_ppm[2]", "NO3_ppm[3]");
        fprintf(files->fp_soiln, "%12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
                "NO3_ppm[4]", "NO3_ppm[5]", "NO3_ppm[6]", "NO3_ppm[7]",
                "NO3_ppm[8]", "NO3_ppm[9]", "etc...");
      }

      if (files->write_co2) {
        printf("Open co2.out\n");
        files->fp_co2 = fopen("co2.out", "w"); 
        fprintf(files->fp_co2, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "CO2_ppm[0]", "CO2_ppm[1]", "CO2_ppm[2]",
                "CO2_ppm[3]", "CO2_ppm[4]");
        fprintf(files->fp_co2, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "CO2_ppm[5]", "CO2_ppm[6]", "CO2_ppm[7]", "CO2_ppm[8]",
                "CO2_ppm[9]", "etc...");
      }

      if (files->write_wflux) {
        printf("Open wflux.out\n");
        files->fp_wflux = fopen("wflux.out", "w"); 
        fprintf(files->fp_wflux, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "wflux[0]", "wflux[1]", "wflux[2]",
                "wflux[3]", "wflux[4]");
        fprintf(files->fp_wflux, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "wflux[5]", "wflux[6]", "wflux[7]", "wflux[8]", "wflux[9]",
                "etc...");
      }

      if (files->write_bio) {
        printf("Open bio.out\n");
        files->fp_bio = fopen("bio.out", "w"); 
        fprintf(files->fp_bio, "%8s  %8s  %10s  %10s  %10s  ",
                "time", "dayofyr", "aglivc", "bglivcj", "bglivcm");
        fprintf(files->fp_bio, "%10s  %10s  %10s  ",
                "aglivn", "bglivnj", "bglivnm");
        fprintf(files->fp_bio, "%10s  %10s  %10s  %10s  %10s  %10s",
                "rleavc", "frootcj", "frootcm", "fbrchc", "rlwodc", "crootc");
        fprintf(files->fp_bio, "%10s  %10s\n",
                "h2ogef(1)", "h2ogef(2)");
      }

      if (files->write_resp) {
        printf("Open resp.out\n");
        files->fp_resp = fopen("resp.out", "w");
        fprintf(files->fp_resp, "%4s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "oiresp", "oeresp", "slitrsp", "sminrlrsp",
                "hresp");
        fprintf(files->fp_resp, "%12s  %12s  %12s  %12s  %12s  ",
                "crtjresp", "crtmresp", "frtjresp", "frtmresp", "frtcresp");
        fprintf(files->fp_resp, "%12s  %12s  %12s  ",
                "sresp", "mresp", "gresp");
        fprintf(files->fp_resp, "%12s  %12s  %12s  %12s  %12s  ",
                "mrspflux(1)", "mrspflux(2)", "cmrspflux(1)", "cmrspflux(2)",
                "cmrspflux(3)");
        fprintf(files->fp_resp, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "fmrspflux(1)", "fmrspflux(2)", "fmrspflux(6)",
                "fmrspflux(3)", "fmrspflux(4)", "fmrspflux(5)");
        fprintf(files->fp_resp, "%12s  %12s  %12s  ",
                "mrspann(1)", "mrspann(2)", "tavedly");
        fprintf(files->fp_resp, "%12s  %12s  %12s  ",
                "mrspTempEffect(1,1)", "mrspTempEffect(1,2)",
                "mrspWaterEffect(1)");
        fprintf(files->fp_resp, "%12s  %12s  %12s  ",
                "mrspTempEffect(2,1)", "mrspTempEffect(2,2)",
                "mrspWaterEffect(2)");
        fprintf(files->fp_resp, "%12s  %12s  %12s  %12s  %12s  ",
                "grspflux(1)", "grspflux(2)", "cgrspflux(1)", "cgrspflux(2)",
                "cgrspflux(3)");
        fprintf(files->fp_resp, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "fgrspflux(1)", "fgrspflux(2)", "fgrspflux(6)",
                "fgrspflux(3)", "fgrspflux(4)", "fgrspflux(5)");
        fprintf(files->fp_resp, "%12s  %12s  ", "grspann(1)", "grspann(2)");
        fprintf(files->fp_resp, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "mcprd(1)", "mcprd(2)", "mcprd(3)", "mfprd(1)", "mfprd(2)",
                "mfprd(6)");
        fprintf(files->fp_resp, "%12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
                "mfprd(3)", "mfprd(4)", "mfprd(5)", "carbostg(1,1)",
                "carbostg(1,2)", "carbostg(2,1)", "carbostg(2,2)");
      }

      if (files->write_yearsum) {
        printf("Open year_summary.out\n");
        files->fp_yearsum = fopen("year_summary.out", "w"); 
        fprintf(files->fp_yearsum, "%8s  %12s  %12s  %12s  %12s  %12s",
                "time", "N2Oflux", "NOflux", "N2flux", "CH4_oxid", "NIT");
        fprintf(files->fp_yearsum, "%12s\n", "ANNPPT");
      }

      if (files->write_livec) {
        printf("Open livec.out\n");
        files->fp_livec = fopen("livec.out", "w"); 
        fprintf(files->fp_livec, "%8s  %8s  %10s  %10s  %10s  ", "time",
                "dayofyr", "aglivc", "bglivcj", "bglivcm");
        fprintf(files->fp_livec, "%10s  %10s  %10s  ",
                "rleavc", "frootcj", "frootcm");
        fprintf(files->fp_livec, "%10s  %10s  %10s\n",
                "fbrchc", "rlwodc", "crootc");
      }

      if (files->write_deadc) {
        printf("Open deadc.out\n");
        files->fp_deadc = fopen("deadc.out", "w"); 
        fprintf(files->fp_deadc, "%8s  %8s  %10s  %10s  %10s  %10s  ", "time",
                "dayofyr", "stdedc", "metabc(1)", "strucc(1)", "wood1c");
        fprintf(files->fp_deadc, "%10s  %10s\n",
                "wood2c", "wood3c");
      }

      if (files->write_soilc) {
        printf("Open soilc.out\n");
        files->fp_soilc = fopen("soilc.out", "w"); 
        fprintf(files->fp_soilc, "%8s  %8s  %10s  %10s  %10s  %10s  ", "time",
                "dayofyr", "metabc(2)", "strucc(2)", "som1c(1)", "som1c(2)");
        fprintf(files->fp_soilc, "%10s  %10s %10s\n",
                "som2c(1)", "som2c(2)", "som3c");
      }

      if (files->write_sysc) {
        printf("Open sysc.out\n");
        files->fp_sysc = fopen("sysc.out", "w"); 
        fprintf(files->fp_sysc, "%8s  %8s  %10s  %10s  %10s  %10s %10s\n",
                "time", "dayofyr", "livec", "deadc", "soilc", "sysc", "CO2resp");
      }

      if (files->write_tgmonth) {
        printf("Open tgmonth.out\n");
        files->fp_tgmonth = fopen("tgmonth.out", "w"); 
        fprintf(files->fp_tgmonth, "%8s  %12s  %12s  %12s  %12s  %12s",
                "time", "N2Oflux", "NOflux", "N2flux", "CH4_oxid", "NIT");
        fprintf(files->fp_tgmonth, "%12s\n", "PPT");
      }

      if (files->write_dN2lyr) {
        printf("Open dN2lyr.out\n");
        files->fp_dN2lyr = fopen("dN2lyr.out", "w"); 
        fprintf(files->fp_dN2lyr, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "dN2_g/m2[0]", "dN2_g/m2[1]", "dN2_g/m2[2]",
                "dN2_g/m2[3]", "dN2_g/m2[4]");
        fprintf(files->fp_dN2lyr, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "dN2_g/m2[5]", "dN2_g/m2[6]", "dN2_g/m2[7]", "dN2_g/m2[8]",
                "dN2_g/m2[9]", "etc...");
      }

      if (files->write_dN2Olyr) {
        printf("Open dN2Olyr.out\n");
        files->fp_dN2Olyr = fopen("dN2Olyr.out", "w"); 
        fprintf(files->fp_dN2Olyr, "%8s  %8s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "dN2O_g/m2[0]", "dN2O_g/m2[1]",
                "dN2O_g/m2[2]", "dN2O_g/m2[3]", "dN2O_g/m2[4]");
        fprintf(files->fp_dN2Olyr, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "dN2O_g/m2[5]", "dN2O_g/m2[6]", "dN2O_g/m2[7]",
                "dN2O_g/m2[8]", "dN2O_g/m2[9]", "etc...");
      }

      if (files->write_dels) {
        printf("Open dels.out\n");
        files->fp_dels = fopen("dels.out", "w"); 
        fprintf(files->fp_dels, "%8s  %8s  %12s  %12s  %12s  %12s  ",
                "time", "dayofyr", "deloi", "deloe", "dsrfclit", "dsmnrl");
        fprintf(files->fp_dels, "%12s  %12s  %12s  %12s  %12s  %12s  ",
                "dhetresp", "dsoilresp", "dcmresp", "dfmresp", "dcgresp",
                "dfgresp");
        fprintf(files->fp_dels, "%12s  %12s\n",
                "dccarbostg", "dfcarbostg");
      }

      if (files->write_dcsip) {
        printf("Open dc_sip.csv\n");
        files->fp_dcsip = fopen("dc_sip.csv", "w");
        fprintf(files->fp_dcsip, "%s,%s,", "time", "dayofyr");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,",
                "trandly", "evapdly", "intrcpt", "sublim", "drain", "runoff");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,%s,",
                "ppt", "accum", "melt", "snow", "snlq", "petdly", "stemp");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,",
                "wc_2cm", "wc_3cm", "wc_5cm", "wc_10cm", "wc_15cm");
        fprintf(files->fp_dcsip, "%s,%s,",
                "wc_30cm", "CO2resp");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,",
                "mcprd(1)", "mcprd(2)", "mcprd(3)", "mfprd(1)", "mfprd(2)",
                "mfprd(6)");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,",
                "mfprd(3)", "mfprd(4)", "mfprd(5)", "NPP", "NEE");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,%s,",
                "aglivc", "bglivcj", "bglivcm", "rleavc", "frootcj",
                "frootcm", "fbrchc");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,",
                "rlwodc", "crootc", "tlai", "stdedc", "wood1c", "wood2c");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,",
                "wood3c", "strucc(1)", "metabc(1)", "strucc(2)", "metabc(2)");
        fprintf(files->fp_dcsip, "%s,%s,%s,%s,%s,%s,\n",
                "som1c(1)", "som1c(2)", "som2c(1)", "som2c(2)", "som3c",
                "totsysc");
      }

      if (files->write_harvest) {
        printf("Open harvest.csv\n");
        files->fp_harv = fopen("harvest.csv", "w");
        fprintf(files->fp_harv, "%s,%s,", "time", "dayofyr");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,%s,",
                "crpval", "agcacc", "bgcjacc", "bgcmacc", "cgrain",
                "egrain(N)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,%s,%s,",
                "egrain(P)", "egrain(S)", "crmvst", "ermvst(N)", "ermvst(P)",
                "ermvst(S)", "cstraw");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,",
                "estraw(N)", "estraw(P)", "estraw(S)", "stdstraw",
                "estdstraw(N)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,%s,",
                "estdstraw(P)", "estdstraw(S)", "addsdc", "addsde(N)",
                "addsde(P)", "addsde(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,",
                "resid", "reside(N)", "reside(P)", "reside(S)", "irrapp");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,%s,%s,",
                "fertapp(N)", "fertapp(P)", "fertapp(S)", "omadapp",
                "omaeapp(N)", "omaeapp(P)", "omaeapp(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,%s,",
                "strmac(1)", "strmac(2)", "strmac(3)", "strmac(4)",
                "strmac(5)", "strmac(6)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,",
                "strmac(7)", "strmac(8)", "cgracc", "egracc(N)",
                "egracc(P)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,%s,",
                "egracc(S)", "accrst", "accrste(N)", "accrste(P)",
                "accrste(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "ctubesj", "etubesj(N)", "etubesj(P)", "etubesj(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "ctubesm", "etubesm(N)", "etubesm(P)", "etubesm(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "srfclittrj", "etubesm(N)", "esrfclittrj(P)",
                "esrfclittrj(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "soillittrj", "esoillittrj(N)", "esoillittrj(P)",
                "esoillittrj(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s,",
                "srfclittrm", "esrfclittrm(N)", "esrfclittrm(P)",
                "esrfclittrm(S)");
        fprintf(files->fp_harv, "%s,%s,%s,%s\n",
                "soillittrm", "esoillittrm(N)", "esoillittrm(P)",
                "esoillittrm(S)");
      }

      if (files->write_cflows) {
        files->fp_cflows = fopen("cflows.out", "w");
        fprintf(files->fp_cflows, "%8s  %8s  %10s  %10s  %10s  %10s  %10s  ",
                "time", "dayofyr", "som11tosom21", "som12tosom22",
                "som12tosom3", "som21tosom11", "som21tosom22");
        fprintf(files->fp_cflows, "%10s  %10s  %10s  %10s  %10s  ",
                "som22tosom12", "som22tosom3", "som3tosom12", "metc1tosom11",
                "metc2tosom12");
        fprintf(files->fp_cflows, "%10s  %10s  %10s  %10s  ",
                "struc1tosom11", "struc1tosom21", "struc2tosom12",
                "struc2tosom22");
        fprintf(files->fp_cflows, "%10s  %10s  %10s  %10s  ",
                "wood1tosom11", "wood1tosom21", "wood2tosom11",
                "wood2tosom21");
        fprintf(files->fp_cflows, "%10s  %10s\n",
                "wood3tosom12", "wood3tosom22");
      }

      if (files->write_yrcflows) {
        files->fp_yrcflows = fopen("year_cflows.out", "w"); 
        fprintf(files->fp_yrcflows, "%8s  %10s  %10s  %10s  %10s  %10s  ",
                "time", "asom11tosom21", "asom12tosom22", "asom12tosom3",
                "asom21tosom11", "asom21tosom22");
        fprintf(files->fp_yrcflows, "%10s  %10s  %10s  %10s  %10s  ",
                "asom22tosom12", "asom22tosom3", "asom3tosom12",
                "ametc1tosom11", "ametc2tosom12");
        fprintf(files->fp_yrcflows, "%10s  %10s  %10s  %10s  ",
                "astruc1tosom11", "astruc1tosom21", "astruc2tosom12",
                "astruc2tosom22");
        fprintf(files->fp_yrcflows, "%10s  %10s  %10s  %10s  ",
                "awood1tosom11", "awood1tosom21", "awood2tosom11",
                "awood2tosom21");
        fprintf(files->fp_yrcflows, "%10s  %10s\n",
                "awood3tosom12", "awood3tosom22");
      }

      if (files->write_daily) {
        files->fp_daily = fopen("daily.out", "w"); 
        fprintf(files->fp_daily, "%8s  %8s  %10s  %10s  %10s  %10s  ",
                "time", "dayofyr", "PET(cm)", "agdefac", "bgdefac",
                "stemp(C)");
        fprintf(files->fp_daily, "%10s  %10s  %10s  %10s  %10s  %10s  %10s\n",
                "snow", "snlq", "thermunits", "aglivc", "aggreenc",
                "hwstress", "scenfrac");
        fprintf(files->fp_daily, "%10s\n", "srad");
      }

      if (files->write_nflux) {
        files->fp_nflux = fopen("nflux.out", "w"); 
        fprintf(files->fp_nflux, "%8s  %8s  %10s  %10s  %10s  %10s  ",
                "time", "dayofyr", "nit_N2O-N", "dnit_N2O-N", "dnit_N2-N",
                "NO-N");
        fprintf(files->fp_nflux, "%10s  %10s\n",
                "CUM-N2O(gN/ha)", "CUM-NO(gN/ha)");
      }

      if (files->write_summary) {
        files->fp_summary = fopen("summary.out", "w"); 
        fprintf(files->fp_summary, "%8s  %8s  %10s  %10s  %10s  %10s  ",
                "time", "dayofyr", "tmax", "tmin", "ppt",
                "N2Oflux");
        fprintf(files->fp_summary, "%10s  %10s  %10s  %10s  %10s\n",
                "NOflux", "CH4_oxid", "NIT", "CO2resp", "ntCO2resp");
      }

      if (files->write_methane) {
        files->fp_methane = fopen("methane.out", "w"); 
        fprintf(files->fp_methane, "%4s  %8s  %10s  %10s  %10s  ",
                "time", "dayofyr", "aglivc", "bglivcj", "bglivcm");
        fprintf(files->fp_methane, "%11s  %11s  %11s  ",
                "prev_mcprd1", "prev_mcprd2", "prev_mcprd3");
        fprintf(files->fp_methane, "%10s  %10s  %10s  %10s  %10s  ",
                "COM", "ppt", "irri", "watr2sat", "avgst_10cm");
        fprintf(files->fp_methane, "%10s  %10s  %10s  %10s  %10s  ",
                "TI", "SI", "Cr", "Eh", "Feh");
        fprintf(files->fp_methane, "%10s  %10s  %10s  %10s\n",
                "CH4_prod", "CH4_Ep", "CH4_Ebl", "CH4_oxid");
      }

      if (files->write_psyn) {
        printf("Open psyn.out\n");
        files->fp_psyn = fopen("psyn.out", "w"); 
        fprintf(files->fp_psyn, "%8s  %8s  %10s  %10s  %10s  %10s   %10s  ",
                "time", "dayofyr", "tmindly", "tmaxdly", "prcann",
                "pptdly", "aetdly");
        fprintf(files->fp_psyn, "%10s  %10s  %10s  %10s  %10s  %10s  %10s  ",
                "petdly", "daylength", "srad", "avg_temp", "avg_vpd", "crpLAI",
                "crpdTemp");
        fprintf(files->fp_psyn, "%10s  %10s  %10s  %10s  %10s  ",
                "crpdVpd", "crpdWater", "crpLtEff", "crpPGrPsn", "crpGrPsn");
        fprintf(files->fp_psyn, "%10s  %10s  %10s  %10s  %10s  ",
                "forLAI", "fordTemp", "fordVpd", "fordWater", "forLtEff");
        fprintf(files->fp_psyn, "%10s  %10s\n", "forPGrPsn", "forGrPsn");
      }

/*!!
files->fp_snow = fopen("snow.out", "w"); 
fprintf(files->fp_snow, "%8s %8s %7s %7s %7s %7s %7s", "time", "dayofyr", "tave",
        "rain", "pet", "snow1", "snlq1");
fprintf(files->fp_snow,"%7s %7s %7s %7s %7s %7s", "snow2",
       "snlq2", "sublim", "snow3", "snlq3", "melt");
fprintf(files->fp_snow, "%7s %7s %7s %7s %7s %7s\n", "snow4",
       "snlq4", "snlq5", "pptrem", "petrem", "runoff");
!!*/

      fclose(files->fp_outf);

      if (flags->debug) {
        printf("Exitting initsw...\n");
      }

      return;
    }
