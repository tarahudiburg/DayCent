
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine growth(tfrac, tavedly, month, crpGrossPsn)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'monprd.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'zztim.inc'

c ... Argument declarations
      real             tfrac, tavedly
      integer          month
      double precision crpGrossPsn

c ... Simulate production for the month.
c ...   tfrac   - fraction of month over which current production event
c ...             occurs (0-1)
c ...   tavedly - mean air temperature over production period (deg C)

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE cmpnfrac(lyr, ammonium, nitrate, minerl,
     &                      frac_nh4, frac_no3)
          !MS$ATTRIBUTES ALIAS:'_cmpnfrac' :: cmpnfrac
          INTEGER          lyr
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          REAL             minerl(*)
          DOUBLE PRECISION frac_nh4
          DOUBLE PRECISION frac_no3
        END SUBROUTINE cmpnfrac

        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &             ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER          clyr
          REAL             amt
          DOUBLE PRECISION frac_nh4
          DOUBLE PRECISION frac_no3
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          CHARACTER        subname*10
        END SUBROUTINE update_npool

      END INTERFACE

c ... Function declarations
      real      fsfunc, line, maxswpot, rtimp
      external  fsfunc, line, maxswpot, rtimp

c ... Local variables
      integer          iel, lyr, ipart
      real             accum(ISOS)
      real             agfrac, amt, availm(MAXIEL),
     &                 bgfracj, bgfracm, calcup, cfrac(CPARTS-1),
     &                 euf(CPARTS), fsol,
     &                 gnfrac, rimpct,
     &                 tm, uptake(4,MAXIEL), toler,
     &                 fraclblstg
      real             soilm(MAXIEL)
      real             cropNfix
      real             namt
      double precision frac_nh4, frac_no3, totmrspdyflux
      real             cprodcdy, eprodcdy(MAXIEL)

      real             agnpp, mrspReduce, rootadj
      real             cstor, mrspprev
      character        subname*10

      if (tfrac .lt. 0.0 .or. tfrac .gt. 1.0) then
        write(*,*) 'Error in growth, tfrac = ', tfrac
        STOP
      endif
      if (tavedly .le. -999.0) then
        write(*,*) 'Error in growth, tavedly = ', tavedly
        STOP
      endif

      subname = 'growth    '
      toler = 1.0E-30

      do  5 iel = 1, nelem
        uptake(ESTOR,iel) = 0.0
        uptake(ESOIL,iel) = 0.0
        uptake(ENFIX,iel) = 0.0
        uptake(EFERT,iel) = 0.0
5     continue

      accum(UNLABL) = 0.0
      accum(LABELD) = 0.0

c ... Use crop/grass specific favail(1) value
      favail(1) = sfavail(1)

c ... Compute fraction of labeled material in the crop/grass
c ... carbohydrate storage pool
      if (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD) .gt. 
     &        0.0) then
        fraclblstg = carbostg(CRPSYS,LABELD) /
     &              (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD))
      else
        fraclblstg = cisofr
      endif

c ... Temperature effect on maintenance respiration for aboveground
c ... components
      mrspTempEffect(CRPSYS,SRFC) = 0.1 * exp(0.07 * tavedly)
c ... Bound maintenance respiration temperature effect between 0.0 and 1.0,
c ... cak - 09/16/02
      mrspTempEffect(CRPSYS,SRFC) =
     &  min(1.0, mrspTempEffect(CRPSYS,SRFC))
      mrspTempEffect(CRPSYS,SRFC) =
     &  max(0.0, mrspTempEffect(CRPSYS,SRFC))

c ... Temperature effect on maintenance respiration for belowground
c ... components
      mrspTempEffect(CRPSYS,SOIL) = 0.1 * exp(0.07 * tavedly)
c ... Bound maintenance respiration temperature effect between 0.0 and 1.0,
c ... cak - 09/16/02
      mrspTempEffect(CRPSYS,SOIL) =
     &  min(1.0, mrspTempEffect(CRPSYS,SOIL))
      mrspTempEffect(CRPSYS,SOIL) =
     &  max(0.0, mrspTempEffect(CRPSYS,SOIL))

c ... Add a soil water term to the root maintenance respiration
c ... equation, cak - 06/27/2007
c ... Calculate the soil water potential of the wettest soil layer
c ... in the crop rooting zone
      mrspWaterEffect(CRPSYS) = maxswpot(claypg)
      if (mrspWaterEffect(CRPSYS) .le. 76.0) then
        mrspWaterEffect(CRPSYS) =
     &    (80.0 - mrspWaterEffect(CRPSYS)) / 80.0
      else if (mrspWaterEffect(CRPSYS) .gt. 76.0) then
        mrspWaterEffect(CRPSYS) = 0.05
      endif
      mrspWaterEffect(CRPSYS) = min(1.0, mrspWaterEffect(CRPSYS))
      mrspWaterEffect(CRPSYS) = max(0.0, mrspWaterEffect(CRPSYS))

c ... Reduce the amount of maintenance respiration based on
c ... predicted annual aboveground NPP and root biomass.
c ... Calculate a predicted NPP value (gmC/m2/yr) based on average
c ... annual precipitation, use 50% of this predicited value as an
c ... approximation of aboveground production
      agnpp = (-40.0 + 3.0 * prcann) * 0.5
c ... Adjust this value based on current root biomass
      if ((bglivcm+bglivcj) .gt. 150.0) then
        rootadj = 1.0
      elseif ((bglivcm+bglivcj) .lt. 0.0) then
        rootadj = 0.05
      else
        rootadj = line((bglivcm+bglivcj), 0.0, 0.05, 150.0, 1.0)
      endif
      agnpp = agnpp * rootadj
c ... Use two line functions to linearly decrease maintenance
c ... respiration as the amount of carbohydrate stored in the
c ... carbohydrate storage pool gets smaller based on the predicted
c ... annual aboveground NPP calculated above, cak - 01/08/2010
C     if (agnpp > 0.000000001) then
C       if (carbostg(CRPSYS,UNLABL)+carbostg(CRPSYS,LABELD) .lt.
C    &      cmrspnpp(3) * agnpp) then
C         mrspReduce =
C    &      line(carbostg(CRPSYS,UNLABL)+carbostg(CRPSYS,LABELD),
C    &           cmrspnpp(1) * agnpp, cmrspnpp(2),
C    &           cmrspnpp(3) * agnpp, cmrspnpp(4))
C       elseif (carbostg(CRPSYS,UNLABL)+carbostg(CRPSYS,LABELD) .gt.
C    &          cmrspnpp(5) * agnpp) then
C         mrspReduce = cmrspnpp(6)
C       else
C         mrspReduce =
C    &      line(carbostg(CRPSYS,UNLABL)+carbostg(CRPSYS,LABELD),
C    &           cmrspnpp(3) * agnpp, cmrspnpp(4),
C    &           cmrspnpp(5) * agnpp, cmrspnpp(6))
C       endif
C     else
C       mrspReduce = 0.0
C     endif
C
c BEGIN NEW CODE =====================================================
c ... This code is consistent with that in treegrowth. -mdh 1/7/2018
      cstor = carbostg(CRPSYS,UNLABL)+carbostg(CRPSYS,LABELD)
      if (agnpp > 0.000000001) then
        if (cstor .lt. cmrspnpp(1) * agnpp) then
c ....... Minimum value        
          mrspReduce = cmrspnpp(2)
        elseif (cstor .gt. cmrspnpp(5) * agnpp) then
c ....... Maximum value
          mrspReduce = cmrspnpp(6)
        elseif (cstor .lt. cmrspnpp(3) * agnpp) then
c ....... First segment
          mrspReduce =
     &      line(cstor, cmrspnpp(1) * agnpp, cmrspnpp(2),
     &                  cmrspnpp(3) * agnpp, cmrspnpp(4))
        else
c ....... Second segment
          mrspReduce =
     &      line(cstor, cmrspnpp(3) * agnpp, cmrspnpp(4),
     &                  cmrspnpp(5) * agnpp, cmrspnpp(6))
        endif
      else
        mrspReduce = 0.0
      endif
c END NEW CODE =======================================================

c ... Added maintenance respiration (mrspflux) calculation. -mdh 2/99
      cmrspdyflux(ABOVE) = ckmrspmx(ABOVE) *
     &                     mrspTempEffect(CRPSYS,SRFC) * aglivc *
     &                     tfrac * mrspReduce
      cmrspdyflux(BELOWJ) = ckmrspmx(BELOWJ) *
     &                      mrspTempEffect(CRPSYS,SOIL) *
     &                      mrspWaterEffect(CRPSYS) * bglivcj *
     &                      tfrac * mrspReduce
      cmrspdyflux(BELOWM) = ckmrspmx(BELOWM) *
     &                      mrspTempEffect(CRPSYS,SOIL) *
     &                      mrspWaterEffect(CRPSYS) * bglivcm *
     &                      tfrac * mrspReduce


c ... When photosynthesis occurs adjust the maintenance respiration
c ... based on photosynthate for the day and a weighted averaging of
c  .. the maintenance respiration values calculated above,
c ... cak - 12/22/2009
      if (crpGrossPsn .gt. 0.0) then
        totmrspdyflux = 0.0
        do 80 ipart = 1, CPARTS
          totmrspdyflux = totmrspdyflux + cmrspdyflux(ipart)
80      continue
        do 90 ipart = 1, CPARTS
          if (totmrspdyflux .gt. 0.0) then
            cmrspdyflux(ipart) =
     &        (cmrspdyflux(ipart) / totmrspdyflux) *
     &        (mrspReduce * ps2mrsp(CRPSYS) * crpGrossPsn)
          else
            cmrspdyflux(ipart) = 0.0
          endif
90      continue
      endif

c ... Maintenance respiration fluxes reduce carbohydrate storage
c ... pool, mdh - 9/4/01
      mrspdyflux(CRPSYS) = cmrspdyflux(ABOVE) +
     &                     cmrspdyflux(BELOWJ) +
     &                     cmrspdyflux(BELOWM)

c BEGIN NEW CODE =====================================================
c ... Don't allow total maintenance respiration to exceed total  
c ... carbohydrate storage (cstor) -mdh 10/7/2018
      mrspprev = mrspdyflux(CRPSYS)
      mrspdyflux(CRPSYS) = min(mrspdyflux(CRPSYS), cstor)
      if (mrspdyflux(CRPSYS) .gt. 0.0) then
        do 82 ipart = 1, CPARTS
          cmrspdyflux(ipart) = cmrspdyflux(ipart) * 
     &      mrspdyflux(CRPSYS)/mrspprev
82      continue
      else
        do 84 ipart = 1, FPARTS
          cmrspdyflux(ipart) = 0.0 
84      continue
       endif
c END NEW CODE ======================++===============================

      call csched(mrspdyflux(CRPSYS),fraclblstg,1.0,
     &            carbostg(CRPSYS,UNLABL),csrsnk(UNLABL),
     &            carbostg(CRPSYS,LABELD),csrsnk(LABELD),
     &            1.0,cautoresp)

c ... Determine actual production values, restricting the C/E ratios
      if (crpgrw .eq. 1 .and. pcropc .gt. 0.0 .and.
     &    .not. (senecnt .gt. 0)) then

c ..... Calculate impact of root biomass on available nutrients
        rimpct = rtimp(riint, rictrl, bglivcj+bglivcm)

c ..... Calculate carbon fraction in each part
        cfrac(ABOVE) = agp / tgprod
        cfrac(BELOW) = 1.0 - cfrac(ABOVE)

c ..... Determine nutrients available to plants for growth.
        do 10 iel = 1, nelem
          availm(iel) = 0.0
c ....... Nutrients available to grasses/crops are in the top claypg layers,
c ....... cak 01/29/03
          do 15 lyr = 1, claypg
            if (minerl(lyr,iel) .gt. toler) then
              availm(iel) = availm(iel) + minerl(lyr, iel)
            endif
15        continue
c ....... save the soil total minerals; savannas change the availability
c ....... KLK - 03/27/2009
          soilm(iel) = availm(iel)
10      continue

c ..... Calculate savanna available fractions
        if (cursys .eq. SAVSYS) then
          tm = MIN(availm(N), 1.5)
          gnfrac = exp(-1.664*exp(-.00102*tm*sitpot)*basfc2*trbasl)
c ....... Bound GNFRAC between 0 and 1
          gnfrac = MIN(gnfrac, 1.0)
          gnfrac = MAX(gnfrac, 0.0)
          do 70 iel = 1, nelem
            availm(iel) = availm(iel) * gnfrac
70        continue
        endif

c ..... Determine actual production values, restricting the C/E ratios
        call restrp(elimit, nelem, availm, cercrp, CPARTS-1, cfrac,
     &              pcropc, rimpct, crpstg, snfxmx(CRPSYS), cprodcdy,
     &              eprodcdy, uptake, crop_a2drat, cropNfix, relyld)

c ..... If growth occurs...
        if (cprodcdy .gt. 0.) then

c ....... If the carbohydrate storage pool falls below a critical value
c ....... add a minimal amount of carbon from the csrsnk to allow plant
c ....... growth.  This should only occur when the plants are small.
          if (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD) .lt.
     &        15.0) then
            write(*,*) 'Warning, carbostg pool below minimal in growth'
            write(*,*) 'time =', time, ' carbostg =', 
     &                  carbostg(CRPSYS,UNLABL)+carbostg(CRPSYS,LABELD)
            carbostg(CRPSYS,UNLABL) = carbostg(CRPSYS,UNLABL) +
     &                                (15.0 * (1.0 - fraclblstg))
            carbostg(CRPSYS,LABELD) = carbostg(CRPSYS,LABELD) +
     &                                (15.0 * fraclblstg)
            csrsnk(UNLABL) = csrsnk(UNLABL) -
     &                       (15.0 * (1.0 - fraclblstg))
            csrsnk(LABELD) = csrsnk(LABELD) - (15.0 * fraclblstg)
          endif

c ....... Increment the counter that is tracking the number of days to
c ....... growth for the current growing season, cak - 03/11/2010
          cgrwdys = cgrwdys + 1

c ....... Calculations for symbiotic N fixation accumulators moved
c ....... from nutrlm subroutine, cak - 10/17/02
c ....... Compute N fixation which actually occurs and add to the
c ....... N fixation accumulator.
          nfix = nfix + cropNfix
          snfxac(CRPSYS) = snfxac(CRPSYS) + cropNfix
c ....... Add computation for nfixac -mdh 1/16/02
          nfixac = nfixac + cropNfix

c ....... C/N ratio for production
          tcnpro = cprodcdy/eprodcdy(N)

c ....... Calculate production for each grass/crop part
          agfrac = agp/tgprod
          mcprd(ABOVE) = cprodcdy * agfrac
          bgfracj = (1.0 - agfrac) * (1.0 - mrtfrac)
          mcprd(BELOWJ) = cprodcdy * bgfracj
          bgfracm = 1.0 - (agfrac + bgfracj)
          mcprd(BELOWM) = cprodcdy * bgfracm

c ....... Crop/grass growth
c ....... All growth comes from the carbohydrate pool, cak - 08/12/2009
c ....... For a perennial, if we have reached the late growing season
c ....... use the late season growth restriction parameter value to
c ....... determine how much carbohydrate to flow out of the crop/grass
c .....,, carbohydrate storage pool, cak - 03/11/2010
          if ((cgrwdys .gt. curgdys) .and.
     &        ((frtcindx .le. 1) .or. (frtcindx .eq. 3))) then
c ......... Late season growth restriction of shoots
            call csched(mcprd(ABOVE)*(1.0-clsgres),fraclblstg,1.0,
     &                  carbostg(CRPSYS,UNLABL),aglcis(UNLABL),
     &                  carbostg(CRPSYS,LABELD),aglcis(LABELD),
     &                  1.0,agcisa)
          else
c ......... Unrestricted growth of shoots
            call csched(mcprd(ABOVE),fraclblstg,1.0,
     &                  carbostg(CRPSYS,UNLABL),aglcis(UNLABL),
     &                  carbostg(CRPSYS,LABELD),aglcis(LABELD),
     &                  1.0,agcisa)
          endif
c ....... Growth of juvenile roots
          call csched(mcprd(BELOWJ),fraclblstg,1.0,
     &                carbostg(CRPSYS,UNLABL),bglcisj(UNLABL),
     &                carbostg(CRPSYS,LABELD),bglcisj(LABELD),
     &                1.0,bgcisja)
c ....... Growth of mature roots
          call csched(mcprd(BELOWM),fraclblstg,1.0,
     &                carbostg(CRPSYS,UNLABL),bglcism(UNLABL),
     &                carbostg(CRPSYS,LABELD),bglcism(LABELD),
     &                1.0,bgcisma)

c ....... Growth respiration
          if ((cgrwdys .gt. curgdys) .and.
     &        ((frtcindx .le. 1) .or. (frtcindx .eq. 3))) then
            cgrspdyflux(ABOVE) = mcprd(ABOVE) * (1.0 - clsgres) *
     &                           cgresp(ABOVE)
          else
            cgrspdyflux(ABOVE) = mcprd(ABOVE) * cgresp(ABOVE)
          endif
          cgrspdyflux(BELOWJ) = mcprd(BELOWJ) * cgresp(BELOWJ)
          cgrspdyflux(BELOWM) = mcprd(BELOWM) * cgresp(BELOWM)
c ....... Growth respiration is subtracted from the carbohydrate
c ....... storage pool.
          grspdyflux(CRPSYS) = cgrspdyflux(ABOVE) +
     &                         cgrspdyflux(BELOWJ) +
     &                         cgrspdyflux(BELOWM)
          call csched(grspdyflux(CRPSYS),fraclblstg,1.0,
     &                carbostg(CRPSYS,UNLABL),csrsnk(UNLABL),
     &                carbostg(CRPSYS,LABELD),csrsnk(LABELD),
     &                1.0,cautoresp)

c ....... Actual uptake
          do 40 iel = 1, nelem
            euf(ABOVE) = eup(ABOVE,iel) / eprodcdy(iel)
            euf(BELOWJ) = (eup(BELOW,iel) * (1.0 - mrtfrac)) /
     &                    eprodcdy(iel)
            euf(BELOWM) = (eup(BELOW,iel) * mrtfrac) / eprodcdy(iel)
c ......... Reset eprodcdy(iel) to track the actual uptake which can be
c ......... restricted late in the growing season
            eprodcdy(iel) = 0.0

c ......... Take up nutrients from internal storage pool
c ......... Don't allow uptake from storage if crpstg is negative,
c ......... cak 07/21/03
c ......... If we have reached the late growing season use the late
c ......... season growth restriction parameter value to determine how
c ......... much nutrients to flow out of the crop/grass nutrient
c .......,, storage pool, cak - 03/11/2010
            if (crpstg(iel) .gt. 0.0) then
              if ((cgrwdys .gt. curgdys) .and.
     &            ((frtcindx .le. 1) .or. (frtcindx .eq. 3))) then
                amt = uptake(ESTOR,iel) * euf(ABOVE) * (1.0 - clsgres)
                call flow(crpstg(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
              else
                amt = uptake(ESTOR,iel) * euf(ABOVE)
                call flow(crpstg(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
              endif
              amt = uptake(ESTOR,iel) * euf(BELOWJ)
              call flow(crpstg(iel),bglivej(iel),time,amt)
              eupbga(iel) = eupbga(iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodcdy(iel) = eprodcdy(iel) + amt
              amt = uptake(ESTOR,iel) * euf(BELOWM)
              call flow(crpstg(iel),bglivem(iel),time,amt)
              eupbga(iel) = eupbga(iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodcdy(iel) = eprodcdy(iel) + amt
            endif

c ......... Take up nutrients from soil
c ......... Nutrients for uptake are available in the top claypg layers,
c ......... cak 01/29/03
c ......... If we have reached the late growing season flow nutrients to
c ......... the crop/grass storage pool rather than to the component
c ......... nutrient pools based on the crop/grass late season growth
c ......... restriction parameter value, cak - 03/11/2010
            do 30 lyr = 1, claypg
              if (minerl(lyr,iel) .gt. toler) then
                fsol = 1.0
                if (iel .eq. P) then
                  fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
                endif
                call cmpnfrac(lyr,ammonium,nitrate,minerl,
     &                        frac_nh4,frac_no3)
c ............. Changed availm(iel) to soilm(iel) to get the correct
c ............. normalization, KLK - 03/27/2009
c                calcup = uptake(ESOIL,iel) *
c     &                   minerl(lyr,iel) * fsol / availm(iel)
                calcup = uptake(ESOIL,iel) *
     &                   minerl(lyr,iel) * fsol / soilm(iel)
c ............. Aboveground live
                namt = 0.0
                if ((cgrwdys .gt. curgdys) .and.
     &              ((frtcindx .le. 1) .or. (frtcindx .eq. 3))) then
                  amt = calcup * euf(ABOVE) * clsgres
                  namt = namt + amt
                  call flow(minerl(lyr,iel),crpstg(iel),time,amt)
                  eupaga(iel) = eupaga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                  eprodcdy(iel) = eprodcdy(iel) + amt
                  amt = calcup * euf(ABOVE) * (1.0 - clsgres)
                  namt = namt + amt
                  call flow(minerl(lyr,iel),aglive(iel),time,amt)
                  eupaga(iel) = eupaga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                  eprodcdy(iel) = eprodcdy(iel) + amt
                else
                  amt = calcup * euf(ABOVE)
                  namt = namt + amt
                  call flow(minerl(lyr,iel),aglive(iel),time,amt)
                  eupaga(iel) = eupaga(iel) + amt
                  eupacc(iel) = eupacc(iel) + amt
                  eprodcdy(iel) = eprodcdy(iel) + amt
                endif
                if (iel .eq. N) then
                  namt = -1.0*namt
                  call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                              ammonium, nitrate, subname)
                endif
c ............. Juvenile fine roots
                namt = 0.0
                amt = calcup * euf(BELOWJ)
                namt = namt + amt
                call flow(minerl(lyr,iel),bglivej(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
                if (iel .eq. N) then
                  namt = -1.0*namt
                  call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                              ammonium, nitrate, subname)
                endif
c ............. Mature fine roots
                namt = 0.0
                amt = calcup * euf(BELOWM)
                namt = namt + amt
                call flow(minerl(lyr,iel),bglivem(iel),time,amt)
                eupbga(iel) = eupbga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
                if (iel .eq. N) then
                  namt = -1.0*namt
                  call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                              ammonium, nitrate, subname)
                endif
              endif
30          continue

c ......... Take up nutrients from nitrogen fixation
            if (iel .eq. N .and. cropNfix .gt. 0) then
c ........... If we have reached the late growing season flow nutrients to
c ........... the crop/grass storage pool rather than to the component
c ........... nutrient pools based on the crop/grass late season growth
c ........... restriction parameter value, cak - 03/11/2010
              if ((cgrwdys .gt. curgdys) .and.
     &            ((frtcindx .le. 1) .or. (frtcindx .eq. 3))) then
                amt = uptake(ENFIX,iel) * euf(ABOVE) * clsgres
                call flow(esrsnk(iel),crpstg(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
                amt = uptake(ENFIX,iel) * euf(ABOVE) * (1.0 - clsgres)
                call flow(esrsnk(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
              else
                amt = uptake(ENFIX,iel) * euf(ABOVE)
                call flow(esrsnk(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
              endif
              amt = uptake(ENFIX,iel) * euf(BELOWJ)
              call flow(esrsnk(iel),bglivej(iel),time,amt)
              eupbga(iel) = eupbga(iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodcdy(iel) = eprodcdy(iel) + amt
              amt = uptake(ENFIX,iel) * euf(BELOWM)
              call flow(esrsnk(iel),bglivem(iel),time,amt)
              eupbga(iel) = eupbga(iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodcdy(iel) = eprodcdy(iel) + amt
            endif

c ......... Take up nutrients from automatic fertilizer
            if (aufert .ne. 0 .and. uptake(EFERT,iel) .gt. 0.0) then
c ........... Automatic fertilizer added to plant pools
c ........... If we have reached the late growing season flow nutrients to
c ........... the crop/grass storage pool rather than to the component
c ........... nutrient pools based on the crop/grass late season growth
c ........... restriction parameter value, cak - 03/11/2010
              if ((cgrwdys .gt. curgdys)  .and.
     &            ((frtcindx .le. 1) .or. (frtcindx .eq. 3))) then
                amt = uptake(EFERT,iel) * euf(ABOVE) * clsgres
                call flow(esrsnk(iel),crpstg(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
                amt = uptake(EFERT,iel) * euf(ABOVE) * (1.0 - clsgres)
                call flow(esrsnk(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(ABOVE)
                call flow(esrsnk(iel),aglive(iel),time,amt)
                eupaga(iel) = eupaga(iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodcdy(iel) = eprodcdy(iel) + amt
              endif
              amt = uptake(EFERT,iel) * euf(BELOWJ)
              call flow(esrsnk(iel),bglivej(iel),time,amt)
              eupbga(iel) = eupbga(iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodcdy(iel) = eprodcdy(iel) + amt
              amt = uptake(EFERT,iel) * euf(BELOWM)
              call flow(esrsnk(iel),bglivem(iel),time,amt)
              eupbga(iel) = eupbga(iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodcdy(iel) = eprodcdy(iel) + amt
c ........... Automatic fertilizer added to mineral pool
              amt = uptake(EFERT,iel) * (1./favail(iel) - 1.)  
              fertot(iel) = fertot(iel) + uptake(EFERT,iel) + amt
              fertac(iel) = fertac(iel) + uptake(EFERT,iel) + amt
              fertmth(month, iel) = fertmth(month,iel) +
     &                              uptake(EFERT,iel) + amt
              if (iel .eq. N) then
                lyr = SRFC
                call update_npool(lyr, amt, frac_nh4_fert, 
     &                            frac_no3_fert, ammonium, nitrate,
     &                            subname)
              endif
              call flow(esrsnk(iel),minerl(SRFC,iel),time,amt)
            endif
40        continue
c ..... Not enough nutrients for production this time step
        else
          cprodcdy = 0.0
          do 55 iel = 1, nelem
            eprodcdy(iel) = 0.0
55        continue
        endif

c ... Else no production this time step
      else
        cprodcdy = 0.0
        do 50 iel = 1, nelem
          eprodcdy(iel) = 0.0
50      continue
      endif

c ... Accumulate monthly output
      mrspflux(CRPSYS) = mrspflux(CRPSYS) + mrspdyflux(CRPSYS)
      grspflux(CRPSYS) = grspflux(CRPSYS) + grspdyflux(CRPSYS)
      cmrspflux(ABOVE) = cmrspflux(ABOVE) + cmrspdyflux(ABOVE)
      cmrspflux(BELOWJ) = cmrspflux(BELOWJ) + cmrspdyflux(BELOWJ)
      cmrspflux(BELOWM) = cmrspflux(BELOWM) + cmrspdyflux(BELOWM)
      mrspmth(CRPSYS)  = mrspmth(CRPSYS) + cmrspdyflux(ABOVE) +
     &                   cmrspdyflux(BELOWJ) + cmrspdyflux(BELOWM)
      cgrspflux(ABOVE) = cgrspflux(ABOVE) + cgrspdyflux(ABOVE)
      cgrspflux(BELOWJ) = cgrspflux(BELOWJ) + cgrspdyflux(BELOWJ)
      cgrspflux(BELOWM) = cgrspflux(BELOWM) + cgrspdyflux(BELOWM)
      grspmth(CRPSYS)  = grspmth(CRPSYS) + cgrspdyflux(ABOVE) +
     &                   cgrspdyflux(BELOWJ) + cgrspdyflux(BELOWM)
      srspmth(CRPSYS)  = srspmth(CRPSYS) +
     &                   cmrspdyflux(BELOWJ) + cmrspdyflux(BELOWM) +
     &                   cgrspdyflux(BELOWJ) + cgrspdyflux(BELOWM)

c ... Accumulate annual output
      mrspann(CRPSYS) = mrspann(CRPSYS) + cmrspdyflux(ABOVE) +
     &                  cmrspdyflux(BELOWJ) + cmrspdyflux(BELOWM)
      grspann(CRPSYS) = grspann(CRPSYS) + cgrspdyflux(ABOVE) +
     &                  cgrspdyflux(BELOWJ) + cgrspdyflux(BELOWM)
      srspann(CRPSYS) = srspann(CRPSYS) +
     &                  cmrspdyflux(BELOWJ) + cmrspdyflux(BELOWJ) +
     &                  cgrspdyflux(BELOWM) + cgrspdyflux(BELOWM)

c ... If we are restricting production late in the growing season update the
c ... production output varaiables
      if ((cgrwdys .gt. curgdys) .and.
     &    ((frtcindx .le. 1) .or. (frtcindx .eq. 3))) then
        mcprd(ABOVE) = mcprd(ABOVE) * (1.0 - clsgres)
        cprodcdy = mcprd(ABOVE) + mcprd(BELOWJ) + mcprd(BELOWM)
      endif

c ... Sum the daily production variables for output to the monthly *.bin file
      cprodc = cprodc + cprodcdy
      do 60 iel = 1, nelem
        eprodc(iel) = eprodc(iel) + eprodcdy(iel)
60    continue

      return
      end
