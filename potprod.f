
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine potprod(cancvr, tmaxdly, tmindly, tavedly, pptdly,
     &                   petdly, tfrac, tavemth, curday, scenfrac,
     &                   aetdly, srad, vpd, crpGrossPsn, forGrossPsn,
     &                   tminslope, tminintercept)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'isovar.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'pheno.inc'
      include 'photosyn.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Argument declarations
      real cancvr
      real tmaxdly, tmindly, tavedly
      real pptdly, petdly
      real tfrac, tavemth, tminslope, tminintercept
      real scenfrac
      real aetdly
      integer curday
      double precision crpGrossPsn, forGrossPsn, srad, vpd

c ... Function declarations
      real     line
      external line

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE calcphotosyn(srad, daytime, lai, potGrossPsn,
     &                          Psyn, dTemp, dVpd, dWater,
     &                          lightEff, grwdys, aMaxScalar1,
     &                          aMaxScalar2, aMaxScalar3, aMaxScalar4,
     &                          aMax, aMaxFrac, attenuation,
     &                          baseFolRespFrac, cFracLeaf, halfSatPar,
     &                          leafCSpWt, growthDays1, growthDays2,
     &                          growthDays3, growthDays4)
          !MS$ATTRIBUTES ALIAS:'_calcphotosyn' :: calcphotosyn
c !CMS_CHG add the bind(C, name = '') to handle case issue      
     &    bind(C, name = 'calcphotosyn') 
          DOUBLE PRECISION srad
          DOUBLE PRECISION daytime
          DOUBLE PRECISION lai
          DOUBLE PRECISION potGrossPsn
          DOUBLE PRECISION Psyn
          DOUBLE PRECISION dTemp
          DOUBLE PRECISION dVpd
          REAL             dWater
          DOUBLE PRECISION lightEff
          INTEGER          grwdys
          DOUBLE PRECISION aMaxScalar1
          DOUBLE PRECISION aMaxScalar2
          DOUBLE PRECISION aMaxScalar3
          DOUBLE PRECISION aMaxScalar4
          DOUBLE PRECISION aMax
          DOUBLE PRECISION aMaxFrac
          DOUBLE PRECISION attenuation
          DOUBLE PRECISION baseFolRespFrac
          DOUBLE PRECISION cFracLeaf
          DOUBLE PRECISION halfSatPar
          DOUBLE PRECISION leafCSpWt
          DOUBLE PRECISION growthDays1
          DOUBLE PRECISION growthDays2
          DOUBLE PRECISION growthDays3
          DOUBLE PRECISION growthDays4
        END SUBROUTINE calcphotosyn

        SUBROUTINE calcpsneffects(maxTemp, minTemp, etrans, potETrans,
     &                            average_temp, average_vpd, dTemp,
     &                            dVpd, dWater, daylength, vpd,
     &                            dVpdExp, dVpdSlope, psnTMin, psnTOpt,
     &                            tminslope, tminintercept, h2ogef3)
     &    bind(C, name = 'calcpsneffects') 
          !MS$ATTRIBUTES ALIAS:'_calcpsneffects' :: calcpsneffects
          DOUBLE PRECISION maxTemp
          DOUBLE PRECISION minTemp
          DOUBLE PRECISION etrans
          DOUBLE PRECISION potETrans
          DOUBLE PRECISION average_temp
          DOUBLE PRECISION average_vpd
          DOUBLE PRECISION dTemp
          DOUBLE PRECISION dVpd
          REAL             dWater
          REAL             daylength
          DOUBLE PRECISION vpd
          DOUBLE PRECISION dVpdExp
          DOUBLE PRECISION dVpdSlope
          DOUBLE PRECISION psnTMin
          DOUBLE PRECISION psnTOpt
          REAL             tminslope
          REAL             tminintercept
          REAL             h2ogef3
        END SUBROUTINE calcpsneffects

        SUBROUTINE wrtpsyn(time, curday, minTemp, maxTemp, annPrecip,
     &                     dailyPrecip, aetdly, petdly, daytime,
     &                     srad, average_temp, average_vpd, crpLAI,
     &                     crpdTemp, crpdVpd, crpdWater, crpLightEff,
     &                     crpPotGrossPsn, crpGrossPsn, forLAI,
     &                     fordTemp, fordVpd, fordWater, forLightEff,
     &                     forPotGrossPsn, forGrossPsn)
          !MS$ATTRIBUTES ALIAS:'_wrtpsyn' :: wrtpsyn
          REAL             time
          INTEGER          curday
          DOUBLE PRECISION minTemp
          DOUBLE PRECISION maxTemp
          DOUBLE PRECISION annPrecip
          DOUBLE PRECISION dailyPrecip
          REAL             aetdly
          REAL             petdly
          DOUBLE PRECISION daytime
          DOUBLE PRECISION srad
          DOUBLE PRECISION average_temp
          DOUBLE PRECISION average_vpd
          DOUBLE PRECISION crpLAI
          DOUBLE PRECISION crpdTemp
          DOUBLE PRECISION crpdVpd
          REAL             crpdWater
          DOUBLE PRECISION crpLightEff
          DOUBLE PRECISION crpPotGrossPsn
          DOUBLE PRECISION crpGrossPsn
          DOUBLE PRECISION forLAI
          DOUBLE PRECISION fordTemp
          DOUBLE PRECISION fordVpd
          REAL             fordWater
          DOUBLE PRECISION forLightEff
          DOUBLE PRECISION forPotGrossPsn
          DOUBLE PRECISION forGrossPsn
        END SUBROUTINE wrtpsyn

      END INTERFACE

c ... Local Variables
      real             accum(ISOS), grossPsn
      real             minTempEff, crpdWater, fordWater
      double precision average_temp, average_vpd, annPrecip, crpLAI
      double precision crpdTemp, crpdVpd, crpLightEff
      double precision crpPotGrossPsn, dailyPrecip, daytime 
      double precision etrans, forLAI, fordTemp, fordVpd
      double precision forLightEff, forPotGrossPsn, maxTemp, minTemp
      double precision potETrans

      accum(UNLABL) = 0.0
      accum(LABELD) = 0.0

c ... Since we are outputing the water stress term add initialization
c ... for the h2ogef array so that it will output zero during periods
c ... of no growth, cak - 04/28/2006
      h2ogef(1) = 0.0
      h2ogef(2) = 0.0

c ... Initialize varaibles used by the photosynthesis subroutine
      minTemp = tmindly
      maxTemp = tmaxdly
      annPrecip = prcann
      dailyPrecip = pptdly
      daytime = daylength(curday)/24.0
      etrans = aetdly
      potETrans = petdly
      crpLAI = 0.0
      crpdTemp = 0.0
      crpdVpd = 0.0
      crpdWater = 0.0
      crpLightEff = 0.0
      crpPotGrossPsn = 0.0
      crpGrossPsn = 0.0
      forLAI = 0.0
      fordTemp = 0.0
      fordVpd = 0.0
      fordWater = 0.0
      forLightEff = 0.0
      forPotGrossPsn = 0.0
      forGrossPsn = 0.0
      average_temp = 0.0
      average_vpd = 0.0

c ... For a Crop System...
      if (crpgrw .eq. 1) then
        if ((frtcindx .lt. 3) .or.
     &      ((frtcindx .ge. 3) .and. (.not. plntkill))) then
c ....... For crops and grasses a leaf area of 1 = 100 grams of biomass
          crplai = aglivc * 2.5 * 0.01
c ....... Call the photosynthesis submodel only on days when growth
c ....... occurs
c          if (cgrwdys .gt. prevcgrwdy) then
          if (cgrwdys .gt. 0) then
            cpsndys = cpsndys + 1
c ......... Calculate the photosynthesis effect on potential growth,
c ......... cak - 01/06/2009
            call calcpsneffects(maxTemp, minTemp, etrans, potETrans,
     &                          average_temp, average_vpd, crpdTemp,
     &                          crpdVpd, crpdWater, daylength(curday),
     &                          vpd, dVpdExp(CRPSYS),
     &                          dVpdSlope(CRPSYS), psnTMin(CRPSYS),
     &                          psnTOpt(CRPSYS), tminslope,
     &                          tminintercept, h2ogef(3))
            call calcphotosyn(srad, daytime, crpLAI, crpPotGrossPsn,
     &                        crpGrossPsn, crpdTemp, crpdVpd,
     &                        crpdWater, crpLightEff, cpsndys,
     &                        aMaxScalar1(CRPSYS), aMaxScalar2(CRPSYS),
     &                        aMaxScalar3(CRPSYS), aMaxScalar4(CRPSYS),
     &                        aMax(CRPSYS), aMaxFrac(CRPSYS),
     &                        attenuation(CRPSYS),
     &                        baseFolRespFrac(CRPSYS),
     &                        cFracLeaf(CRPSYS), halfSatPar(CRPSYS),
     &                        leafCSpWt(CRPSYS), growthDays1(CRPSYS),
     &                        growthDays2(CRPSYS), growthDays3(CRPSYS),
     &                        growthDays4(CRPSYS))
c ......... Modify the temperature effect on photosynthesis due to low
c ......... daily minimum temperature values.  When the plant gets cold
c ......... the stomatal conductance is reduced.  cak - 09/01/2009
            minTempEff = line(tmindly, -6.0, 0.4, 2.0, 1.0)
            minTempEff = min(minTempEff, 1.0)
            minTempEff = max(minTempEff, 0.4)
            crpGrossPsn = crpGrossPsn * minTempEff
            grossPsn = crpGrossPsn
c ......... Add the gross photosynthesis to the carbohydrate storage pool
c ......... cak - 08/12/2009
            call csched(grossPsn,cisofr,1.0,
     &                  csrsnk(UNLABL),carbostg(CRPSYS,UNLABL),
     &                  csrsnk(LABELD),carbostg(CRPSYS,LABELD),
     &                  1.0,accum)
          endif
          call potcrp(cancvr, tavedly, petdly, tfrac, scenfrac, srad,
     &                daylength(curday))
        endif
      endif

c ... For a Forest System...
      if (forgrw .eq. 1) then
        forlai = rleavc * 2.5 * btolai
c ..... Call the photosynthesis submodel only on days when growth
c ..... occurs
        if (fgrwdys .gt. 0) then
          fpsndys = fpsndys + 1
c ....... Calculate the photosynthesis effect on potential growth,
c ....... cak - 01/06/2009
          call calcpsneffects(maxTemp, minTemp, etrans, potETrans,
     &                        average_temp, average_vpd, fordTemp,
     &                        fordVpd, fordWater, daylength(curday),
     &                        vpd, dVpdExp(FORSYS), dVpdSlope(FORSYS),
     &                        psnTMin(FORSYS), psnTOpt(FORSYS),
     &                        tminslope, tminintercept, h2ogef(3))
          call calcphotosyn(srad, daytime, forLAI, forPotGrossPsn,
     &                      forGrossPsn, fordTemp, fordVpd, fordWater,
     &                      forLightEff, fpsndys,
     &                      aMaxScalar1(FORSYS), aMaxScalar2(FORSYS),
     &                      aMaxScalar3(FORSYS), aMaxScalar4(FORSYS),
     &                      aMax(FORSYS), aMaxFrac(FORSYS),
     &                      attenuation(FORSYS),
     &                      baseFolRespFrac(FORSYS), cFracLeaf(FORSYS),
     &                      halfSatPar(FORSYS), leafCSpWt(FORSYS),
     &                      growthDays1(FORSYS), growthDays2(FORSYS),
     &                      growthDays3(FORSYS), growthDays4(FORSYS))
c ....... Modify the temperature effect on photosynthesis due to low
c ....... daily minimum temperature values.  When the plant gets cold
c ....... the stomatal conductance is reduced.  cak - 09/01/2009
          minTempEff = line(tmindly, -6.0, 0.4, 2.0, 1.0)
          minTempEff = min(minTempEff, 1.0)
          minTempEff = max(minTempEff, 0.4)
          forGrossPsn = forGrossPsn * minTempEff
          grossPsn = forGrossPsn
c ....... Add the gross photosynthesis to the carbohydrate storage pool
c ....... cak - 08/12/2009
          call csched(grossPsn,cisotf,1.0,
     &                csrsnk(UNLABL),carbostg(FORSYS,UNLABL),
     &                csrsnk(LABELD),carbostg(FORSYS,LABELD),
     &                1.0,accum)
        endif
        call potfor(tavedly, petdly, tfrac, tavemth, srad,
     &              daylength(curday))
      endif

c ... Write the photosynthesis values to the output file
      if (time .ge. strplt) then
        call wrtpsyn(time, curday, minTemp, maxTemp, annPrecip,
     &               dailyPrecip, aetdly, petdly, daytime, srad,
     &               average_temp, average_vpd, crpLAI, crpdTemp,
     &               crpdVpd, crpdWater,  crpLightEff, crpPotGrossPsn,
     &               crpGrossPsn, forLAI, fordTemp, fordVpd, fordWater,
     &               forLightEff, forPotGrossPsn, forGrossPsn)
      endif

      return
      end
