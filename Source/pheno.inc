
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      common/pheno/accumdd, basetemp(2), cgrwdys, clsgres, cpsndys,
     &             curgdys, dayhrs, ddbase, ddemerg, decidgrow, emerg,
     &             fgrwdys, flsgres, fpsndys, frstmth, furgdys, hrsinc,
     &             grnfill, grnfldys, mnddhrv, mxddhrv, plntkill,
     &             soiltavewk, tfstmth, thermunits,
     &             tmpkill, tmplff, tmplfs

c ... accumdd     - flag, set to true on FRST/PLTM event, is reset to
c ...               false when LAST occurs
c ... basetemp(1) - base temperature for crop growth, growing degree
c ...               days will accumulate only on days when the average
c ...               temperature is greater than the base temperature
c ...               for the crop (degrees C)
c ... basetemp(2) - ceiling on the maximum temperature used to
c ...               accumulate growing degree days (degrees C)
c ... cgrwdys     - number of days of grass/crop growth in the current
c ...               growing season
c ... clsgres     - late season crop/grass growth restriction factor
c ...               (0.0 - 1.0)
c ... cpsndys     - number of days photosynthesis has occurred during
c ...               the current growing season for crop/grass system
c ... curgdys     - number of days after the start of growing season
c ...               of unrestricted growth of grass/crop components
c ... dayhrs      - number of hours in a day based on day of year
c ... ddbase      - thermal units (growing degree days) required to
c ...               trigger a senescence event for a perennial plant
c ...               (FRTCINDX 3) or a non-grain filling annual
c ...               (FRTCINDX = 4) or anthesis for a grain filling
c ...               annual (FRTCINDX 5 or 6)
c ... ddemerg     - thermal units (growing degree days) required for
c ...               after planting for emergence to occur in annual
c ...               plants (FRTCINDX 4, 5, or 6)
c ... decidgrow   - true if in period between leaf out and leaf
c ...               drop in a deciduous forest
c ... emerg       - flag, indicates if emergence has occurred
c ... fgrwdys     - number of days that tree growth has occurred in the
c ...               current growing season
c ... flsgres     - late season forest growth restriction factor
c ...               (0.0 - 1.0)
c ... fpsndys     - number of days photosynthesis has occurred during
c ...               the current growing season for forest system
c ... frstmth     - month FRST/PLTM event occurs
c ... furgdys     - number of days after leaf out of unrestricted growth
c ...               of woody tree components for a decidious tree type
c ... grnfill     - true if in grain fill period for a grain filling
c ...               annual (FRTCINDX 5 or 6)
c ... grnfldys    - number of days that grain filling has been
c ...               occurring
c ... hrsinc      - indicates if the length of days in increasing
c ...               or decreasing
c ... mnddhrv     - minimum number of degree days from anthesis to
c ...               harvest for grain filling annuals (FRTCINDX 5 or 6)
c ... mxddhrv     - maximum number of degree days from anthesis to
c ...               harvest for grain filling annuals (FRTCINDX 5 or 6)
c ... plntkill    - true if a killing frost has occurred
c ... soiltavewk  - average soil temperature in second DayCent soils.in
c ...               layer over the previous 7 days, degrees C
c ... tfstmth     - month TFST event occurs
c ... thermunits  - accumulator of thermal units for growing degree
c ...               day implementation
c ... tmpkill     - air temperature at which plant growth stops,
c ...               degrees C
c ... tmplff      - temperature at which leaf drop will occur in a
c ...               deciduous tree type, degrees C
c ... tmplfs      - temperature at which leaf out will occur in a
c ...               deciduous tree type, degrees C

      integer cgrwdys, cpsndys, curgdys, fgrwdys, fpsndys, frstmth,
     &        furgdys, grnfldys, tfstmth
      real basetemp, clsgres, dayhrs, ddbase, ddemerg, flsgres,
     &     mnddhrv, mxddhrv, soiltavewk, thermunits,
     &     tmpkill, tmplff, tmplfs
      logical accumdd, decidgrow, emerg, grnfill, hrsinc, plntkill

      save /pheno/
