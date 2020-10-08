
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine treegrow (tfrac, tavedly, forGrossPsn)

      implicit none
      include 'const.inc'
      include 'dynam.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'monprd.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'zztim.inc'

c ... Argument declarations
      real             tfrac, tavedly
      double precision forGrossPsn

c ... Simulate forest production associated with tree growth.
c ... This function was a part of TREES.F.
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
      real      leafa, line, maxswpot, rtimp
      external  leafa, line, maxswpot, rtimp

c ... Local variables
      integer   iel, ipart, iptr, lyr
      real      accum(ISOS), availm(MAXIEL), amt
      real      calcup, cprodfLeft
      real      euf(FPARTS), fsol
      real      namt
      real      remCfrac, rimpct
      real      sum_cfrac, toler, totCup
      real      uptake(4,MAXIEL)
      real      treeNfix
      double precision frac_nh4, frac_no3, totmrspdyflux
      real      cprodfdy, eprodfdy(MAXIEL)
      real      lai, rleavc_opt, mrspReduce, fraclblstg
      character subname*10

      if (tfrac .lt. 0.0 .or. tfrac .gt. 1.0) then
        write(*,*) 'Error in treegrow, tfrac = ', tfrac
        STOP
      endif
      if (tavedly .le. -999.0) then
        write(*,*) 'Error in treegrow, tavedly = ', tavedly
        STOP
      endif

      subname = 'treegrow  '
      toler = 1.0E-30

      accum(UNLABL) = 0.0
      accum(LABELD) = 0.0

      do 15 iel = 1, nelem
        uptake(ESTOR,iel) = 0.0
        uptake(ESOIL,iel) = 0.0
        uptake(ENFIX,iel) = 0.0
        uptake(EFERT,iel) = 0.0
15    continue

c ... Use tree specific favail(1) value
      favail(1) = sfavail(2)

c ... Compute fraction of labeled material in the forest carbohydrate
c ... storage pool
      if (carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD) .gt.
     &        0.0) then
        fraclblstg = carbostg(FORSYS,LABELD) /
     &              (carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD))
      else
        fraclblstg = cisotf
      endif

c ... Previous flowup, in trees, should have updated mineral pools.
c ... Determine nutrients available to plants for growth.
      do 20 iel = 1, nelem
        availm(iel) = 0.0
c ..... Nutrients available to trees are in the top tlaypg layers,
c ..... cak 01/29/03
c        do 30 lyr = 1, nlayer
        do 30 lyr = 1, tlaypg
          if (minerl(lyr,iel) .gt. toler) then
            availm(iel) = availm(iel) + minerl(lyr, iel)
          endif
30      continue
20    continue

c ... Determine old or new forest
c ... iptr points to new forest carbon allocation fractions (iptr = 1) or
c ... mature forest carbon allocation fractions (iptr = 2) for each of the
c ... tree parts; leaves, fine roots, fine branches, large wood, and
c ... coarse roots.  Switch from new forest allocation fractions to old
c ... forest allocation fractions at time = swold
      if (time .le. swold) then
c ..... Use juvenile forest C allocation fractions
        iptr = 1
      else
c ..... Use mature forest C allocation fractions
        iptr = 2
      endif

c ... Temperature effect on maintenance respiration for aboveground
c ... components
      mrspTempEffect(FORSYS,SRFC) = 0.1 * exp(0.07 * tavedly)
c ... Bound maintenance respiration temperature effect between 0.0 and 1.0,
c ... cak - 09/16/02
      mrspTempEffect(FORSYS,SRFC) =
     &  min(1.0, mrspTempEffect(FORSYS,SRFC))
      mrspTempEffect(FORSYS,SRFC) =
     &  max(0.0, mrspTempEffect(FORSYS,SRFC))

c ... Temperature effect on maintenance respiration for belowground
c ... components
      mrspTempEffect(FORSYS,SOIL) = 0.1 * exp(0.07 * tavedly)
c ... Bound maintenance respiration temperature effect between 0.0 and 1.0,
c ... cak - 09/16/02
      mrspTempEffect(FORSYS,SOIL) =
     &  min(1.0, mrspTempEffect(FORSYS,SOIL))
      mrspTempEffect(FORSYS,SOIL) =
     &  max(0.0, mrspTempEffect(FORSYS,SOIL))

c ... Add a soil water term to the fine root maintenance respiration
c ... equation, cak - 06/27/2007
c ... Calculate the soil water potential of the wettest soil layer
c ... in the tree rooting zone
      mrspWaterEffect(FORSYS) = maxswpot(tlaypg)
      if (mrspWaterEffect(FORSYS) .le. 76.0) then
        mrspWaterEffect(FORSYS) =
     &    (80.0 - mrspWaterEffect(FORSYS)) / 80.0
      else if (mrspWaterEffect(FORSYS) .gt. 76.0) then
        mrspWaterEffect(FORSYS) = 0.05
      endif
      mrspWaterEffect(FORSYS) = min(1.0, mrspWaterEffect(FORSYS))
      mrspWaterEffect(FORSYS) = max(0.0, mrspWaterEffect(FORSYS))

c ... Linearly decrease maintenance respiration as the amount of
c ... carbohydrate stored in the carbohydrate storage pool gets
c ... smaller based on the carbon level that is required for optimal
c ... leaf area, cak - 10/20/2006
c ... We will be using two line functions for this calculation
c ... depending on the optimal leaf carbon value, cak - 08/13/2009
c ... Calculate the amount of carbon in the leaves if the tree is at
c ... its optimal LAI, cak - 10/20/2006
      call lacalc(lai, fbrchc, rlwodc, maxlai, klai)
      rleavc_opt = lai / (2.5 * btolai)
      if (rleavc_opt > 0.000000001) then
        if (carbostg(FORSYS,UNLABL)+carbostg(FORSYS,LABELD) .lt.
     &      fmrsplai(3) * rleavc_opt) then
          mrspReduce =
     &      line(carbostg(FORSYS,UNLABL)+carbostg(FORSYS,LABELD),
     &           fmrsplai(1) * rleavc_opt, fmrsplai(2),
     &           fmrsplai(3) * rleavc_opt, fmrsplai(4))
        elseif (carbostg(FORSYS,UNLABL)+carbostg(FORSYS,LABELD) .gt.
     &          fmrsplai(5) * rleavc_opt) then
          mrspReduce = fmrsplai(6)
        else
          mrspReduce =
     &      line(carbostg(FORSYS,UNLABL)+carbostg(FORSYS,LABELD),
     &           fmrsplai(3) * rleavc_opt, fmrsplai(4),
     &           fmrsplai(5) * rleavc_opt, fmrsplai(6))
        endif
      else
        mrspReduce = 0.0
      endif

c ... Maintenance respiration flux calculation added, mdh - 9/4/01
      fmrspdyflux(LEAF) = fkmrspmx(LEAF) *
     &                    mrspTempEffect(FORSYS,SRFC) * rleavc *
     &                    tfrac * mrspReduce
      fmrspdyflux(FROOTJ) = fkmrspmx(FROOTJ) *
     &                      mrspTempEffect(FORSYS,SOIL) * 
     &                      mrspWaterEffect(FORSYS) * frootcj * tfrac *
     &                      mrspReduce
      fmrspdyflux(FROOTM) = fkmrspmx(FROOTM) *
     &                      mrspTempEffect(FORSYS,SOIL) *
     &                      mrspWaterEffect(FORSYS) * frootcm * tfrac *
     &                      mrspReduce
      fmrspdyflux(FBRCH) = fkmrspmx(FBRCH) *
     &                     mrspTempEffect(FORSYS,SRFC) * fbrchc *
     &                     tfrac * mrspReduce
      fmrspdyflux(LWOOD) = fkmrspmx(LWOOD) *
     &                     mrspTempEffect(FORSYS,SRFC) * rlwodc *
     &                     tfrac * mrspReduce
      fmrspdyflux(CROOT) = fkmrspmx(CROOT) *
     &                     mrspTempEffect(FORSYS,SOIL) * crootc *
     &                     tfrac * mrspReduce

c ... When photosynthesis occurs adjust the maintenance respiration
c ... based on photosynthate for the day and a weighted averaging of
c  .. the maintenance respiration values calculated above,
c ... cak - 12/22/2009
      if (forGrossPsn .gt. 0.0) then
        totmrspdyflux = 0.0
        do 40 ipart = 1, FPARTS
          totmrspdyflux = totmrspdyflux + fmrspdyflux(ipart)
40      continue
        do 50 ipart = 1, FPARTS
          if (totmrspdyflux .gt. 0.0) then
            fmrspdyflux(ipart) =
     &        (fmrspdyflux(ipart) / totmrspdyflux) *
     &        (mrspReduce * ps2mrsp(FORSYS) * forGrossPsn)
          else
            fmrspdyflux(ipart) = 0.0
          endif
50      continue
      endif

c ... Maintenance respiration is subtracted from the carbohydrate
c ... storage pool, cak - 08/12/2009
      mrspdyflux(FORSYS) = 0.0
      do 55 ipart = 1, FPARTS
        mrspdyflux(FORSYS) = mrspdyflux(FORSYS) + fmrspdyflux(ipart)
55    continue
      call csched(mrspdyflux(FORSYS),fraclblstg,1.0,
     &            carbostg(FORSYS,UNLABL),csrsnk(UNLABL),
     &            carbostg(FORSYS,LABELD),csrsnk(LABELD),
     &            1.0,fautoresp)

c ... If growth can occur
      if (forgrw .eq. 1 .and. pforc .gt. 0.0) then

c ..... Calculate actual production values
c ..... Calculate impact of root biomass on available nutrients
        rimpct = rtimp(riint, rictrl, frootcj+frootcm)
c ..... Determine actual production values, restricting the C/E ratios
c ..... When calling restrp we are only looking at allocation to fine roots
c ..... and leaves, cak - 07/02/02
        call restrp(elimit, nelem, availm, ccefor, 2, tree_cfrac,
     &              pforc, rimpct, forstg, snfxmx(FORSYS), cprodfdy,
     &              eprodfdy, uptake, tree_a2drat, treeNfix, relyld)
      else
        cprodfdy = 0.0
        do 58 iel = 1, MAXIEL
          eprodfdy(iel) = 0.0
58      continue
      endif

c ... If growth occurs...
      if (cprodfdy .gt. 0.) then

c ..... If the carbohydrate storage pool falls below a critical value
c ..... add a minimal amount of carbon from the csrsnk to allow plant
c ..... growth.  This should only occur when the tree is small.
        if (carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD) .lt. 15.0)
     &    then
          write(*,*) 'Warning, carbostg pool below minimal in treegrow'
          write(*,*) 'time =' , time, '  carbostg =', 
     &                carbostg(FORSYS,UNLABL)+carbostg(FORSYS,LABELD)
          carbostg(FORSYS,UNLABL) = carbostg(FORSYS,UNLABL) +
     &                              (15.0 * (1.0 - fraclblstg))
          carbostg(FORSYS,LABELD) = carbostg(FORSYS,LABELD) +
     &                              (15.0 * fraclblstg)
          csrsnk(UNLABL) = csrsnk(UNLABL) - (15.0 * (1.0 - fraclblstg))
          csrsnk(LABELD) = csrsnk(LABELD) - (15.0 * fraclblstg)
        endif

c ..... Increment the counter that is tracking the number of days to
c ..... allow woody growth, cak - 10/20/2006
        fgrwdys = fgrwdys + 1

c ..... Compute carbon allocation fractions for each tree part
c ..... Calculate how much of the carbon the roots use
        cprodfLeft = cprodfdy - (cprodfdy * tree_cfrac(FROOT))

c ..... Calculate how much of the carbon the leaves use, we allocate leaves
c ..... up to a optimal LAI
        if (rleavc .lt. 0.0) then
          write(*,*) 'Warning in treegrow, rleavc < 0.0'
          write(*,*) 'time=', time, '  rleavc=', rleavc
          rlvcis(UNLABL) = 0.0
          rlvcis(LABELD) = 0.0
          rleavc = rlvcis(UNLABL) + rlvcis(LABELD) 
        endif
        tree_cfrac(LEAF) = leafa(rleavc, fbrchc, rlwodc, cprodfLeft,
     &                           cprodfdy)
        remCfrac = 1.0 - tree_cfrac(FROOT) - tree_cfrac(LEAF)

c ..... If we have leftover carbon allocate it to the woody plant parts
c ..... using a weighted average
        if (remCfrac .lt. 1.0E-05) then
c ..... for FBRCH, LWOOD, and CROOT ...
          do 60 ipart = FBRCH, CROOT
            tree_cfrac(ipart) = 0.0
60        continue
        else
c ....... for FBRCH, LWOOD, and CROOT ...
          totCup = 0.0
          do 70 ipart = FBRCH, CROOT
            tree_cfrac(ipart) = fcfrac(ipart, iptr)
            totCup = totCup + tree_cfrac(ipart)
70        continue
          if (totCup .gt. 0.0) then
c ....... for FBRCH, LWOOD, and CROOT ...
            do 80 ipart = FBRCH, CROOT
              tree_cfrac(ipart) = tree_cfrac(ipart) / totCup * remCfrac
80          continue
          else
            write(*,*) 'Error in treegrow'
            write(*,*) 'fcfrac(FBRCH)+fcfrac(LWOOD)+fcfrac(CROOT) <= 0'
            write(*,*) 'should be > 0 for tree'
            STOP
          endif
        endif

c ... Error checking
        sum_cfrac = 0.0
        do 90 ipart = 1, FPARTS-1
          if (tree_cfrac(ipart) .lt. 0.0) then
            write(*,*) 'Error in treegrow, tree_cfrac(ipart) < 0'
            STOP
          else if (tree_cfrac(ipart) .gt. 1.0) then
            write(*,*) 'Error in treegrow, tree_cfrac(ipart) > 1'
            STOP
          else
            sum_cfrac = sum_cfrac + tree_cfrac(ipart)
          endif
90      continue
        if (abs(1.0 - sum_cfrac) .gt. 0.001) then
          write(*,*) "Error in tree carbon allocation fractions!"
          write(*,*) "sum_cfrac = ", sum_cfrac
c          STOP
        endif

c ..... Recalculate actual production values with updated C-allocation
c ..... fractions, restricting the C/E ratios -mdh 5/11/01
c ..... Calculate impact of root biomass on available nutrients
        rimpct = rtimp(riint, rictrl, frootcj+frootcm)
c ..... Determine actual production values, restricting the C/E ratios
        call restrp(elimit, nelem, availm, ccefor, FPARTS-1,
     &              tree_cfrac, pforc, rimpct, forstg, snfxmx(FORSYS),
     &              cprodfdy, eprodfdy, uptake, tree_a2drat, treeNfix,
     &              relyld)

c ..... Calculations for symbiotic N fixation accumulators moved from
c ..... nutrlm subroutine, cak - 10/17/02
c ..... Compute N fixation which actually occurs and add to the
c ..... N fixation accumulator.
        nfix = nfix + treeNfix
        snfxac(FORSYS) = snfxac(FORSYS) + treeNfix
c ..... Add computation for nfixac -mdh 1/16/02
        nfixac = nfixac + treeNfix

c ..... C/N ratio for production
C       if (eprodfdy(N) .eq. 0.0) then
        if (eprodfdy(N) .le. 0.0) then
          write(*,*) 'Error in treegrowth, eprodfdy(N) = 0.0)'
          STOP
        endif
        tcnpro = cprodfdy/eprodfdy(N)

c ..... Calculate production for each tree part
c ..... New variable MFPRD added for gridded output - 6/96 rm
        do 95 ipart = 1, FPARTS
          if (ipart .eq. LEAF .or. ipart .eq. FBRCH .or.
     &        ipart .eq. LWOOD .or. ipart .eq. CROOT) then
            mfprd(ipart) = tree_cfrac(ipart) * cprodfdy
          else if (ipart .eq. FROOTJ) then
            mfprd(ipart) = tree_cfrac(FROOTJ) * cprodfdy *
     &                     (1.0 - wmrtfrac)
          else if (ipart .eq. FROOTM) then
            mfprd(ipart) = tree_cfrac(FROOTJ) * cprodfdy * wmrtfrac
          else
            write(*,*) 'Error in treegrow, ipart out of bounds = ',
     &                  ipart
            STOP
          endif
95      continue

c ..... Forest Growth
c ..... All growth comes from the carbohydrate pool, cak - 08/12/2009
c ..... Growth of leaves split into labeled & unlabeled parts
        call csched(mfprd(LEAF),fraclblstg,1.0,
     &              carbostg(FORSYS,UNLABL),rlvcis(UNLABL),
     &              carbostg(FORSYS,LABELD),rlvcis(LABELD),
     &              1.0,alvcis)
c ..... Growth juvenile fine roots; split into labeled & unlabeled 
c ..... parts
        call csched(mfprd(FROOTJ),fraclblstg,1.0,
     &              carbostg(FORSYS,UNLABL),frtcisj(UNLABL),
     &              carbostg(FORSYS,LABELD),frtcisj(LABELD),
     &              1.0,afrcisj)
c ..... Growth mature fine roots; split into labeled & unlabeled parts
        call csched(mfprd(FROOTM),fraclblstg,1.0,
     &              carbostg(FORSYS,UNLABL),frtcism(UNLABL),
     &              carbostg(FORSYS,LABELD),frtcism(LABELD),
     &              1.0,afrcism)
c ..... If the tree is decidious and the time alloted to grow the
c ..... woody components of the trees has passed use the late season
c ..... growth restriction parameter value to determine how much
c ..... carbohydrate to flow out of the forest carbohydrate storage
c ..... pool for the woody components, cak - 03/11/2010
        if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys)) then
c ....... Growth of fine branches; split into labeled & unlabeled parts
          call csched(mfprd(FBRCH)*(1.0-flsgres),fraclblstg,1.0,
     &                carbostg(FORSYS,UNLABL),fbrcis(UNLABL),
     &                carbostg(FORSYS,LABELD),fbrcis(LABELD),
     &                1.0,afbcis)
c ....... Growth of large wood; split into labeled & unlabeled parts
          call csched(mfprd(LWOOD)*(1.0-flsgres),fraclblstg,1.0,
     &                carbostg(FORSYS,UNLABL),rlwcis(UNLABL),
     &                carbostg(FORSYS,LABELD),rlwcis(LABELD),
     &                1.0,alwcis)
c ....... Growth of coarse roots; split into labeled & unlabeled parts
          call csched(mfprd(CROOT)*(1.0-flsgres),fraclblstg,1.0,
     &                carbostg(FORSYS,UNLABL),crtcis(UNLABL),
     &                carbostg(FORSYS,LABELD),crtcis(LABELD),
     &                1.0,acrcis)
        else
c ....... Growth of fine branches; split into labeled & unlabeled parts
          call csched(mfprd(FBRCH),fraclblstg,1.0,
     &                carbostg(FORSYS,UNLABL),fbrcis(UNLABL),
     &                carbostg(FORSYS,LABELD),fbrcis(LABELD),
     &                1.0,afbcis)
c ....... Growth of large wood; split into labeled & unlabeled parts
          call csched(mfprd(LWOOD),fraclblstg,1.0,
     &                carbostg(FORSYS,UNLABL),rlwcis(UNLABL),
     &                carbostg(FORSYS,LABELD),rlwcis(LABELD),
     &                1.0,alwcis)
c ....... Growth of coarse roots; split into labeled & unlabeled parts
          call csched(mfprd(CROOT),fraclblstg,1.0,
     &                carbostg(FORSYS,UNLABL),crtcis(UNLABL),
     &                carbostg(FORSYS,LABELD),crtcis(LABELD),
     &                1.0,acrcis)
        endif

c ..... Growth respiration
        fgrspdyflux(LEAF) = mfprd(LEAF) * fgresp(LEAF)
        fgrspdyflux(FROOTJ) = mfprd(FROOTJ) * fgresp(FROOTJ)
        fgrspdyflux(FROOTM) = mfprd(FROOTM) * fgresp(FROOTM)
        if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys)) then
          fgrspdyflux(FBRCH) = mfprd(FBRCH) * (1.0 - flsgres) *
     &                         fgresp(FBRCH)
          fgrspdyflux(LWOOD) = mfprd(LWOOD) * (1.0 - flsgres) *
     &                         fgresp(LWOOD)
          fgrspdyflux(CROOT) = mfprd(CROOT) * (1.0 - flsgres) *
     &                         fgresp(CROOT)
        else
          fgrspdyflux(FBRCH) = mfprd(FBRCH) * fgresp(FBRCH)
          fgrspdyflux(LWOOD) = mfprd(LWOOD) * fgresp(LWOOD)
          fgrspdyflux(CROOT) = mfprd(CROOT) * fgresp(CROOT)
        endif
c ..... Growth respiration is subtracted from the carbohydrate
c ..... storage pool.
        grspdyflux(FORSYS) = 0.0
        do 105 ipart = 1, FPARTS
          grspdyflux(FORSYS) = grspdyflux(FORSYS) + fgrspdyflux(ipart)
105     continue
        call csched(grspdyflux(FORSYS),fraclblstg,1.0,
     &              carbostg(FORSYS,UNLABL),csrsnk(UNLABL),
     &              carbostg(FORSYS,LABELD),csrsnk(LABELD),
     &              1.0,fautoresp)

c ..... Actual Uptake
        do 110 iel = 1, nelem
C         if (eprodfdy(iel) .eq. 0.0) then
          if (eprodfdy(iel) .le. 0.0) then
            write(*,*) 'Divide by zero in treegrow, eprodfdy(iel) = 0'
            STOP
          endif
          euf(LEAF) = eup(LEAF,iel) / eprodfdy(iel)
          euf(FROOTJ) = (eup(FROOT,iel) * (1.0 - wmrtfrac)) /
     &                  eprodfdy(iel)
          euf(FROOTM) = (eup(FROOT,iel) * wmrtfrac) / eprodfdy(iel)
          euf(FBRCH) = eup(FBRCH,iel) / eprodfdy(iel)
          euf(LWOOD) = eup(LWOOD,iel) / eprodfdy(iel)
          euf(CROOT) = eup(CROOT,iel) / eprodfdy(iel)
c ....... Reset eprodcdy(iel) to track the actual uptake which can be
c ....... restricted late in the growing season
          eprodfdy(iel) = 0.0

c ....... Take up nutrients from internal storage pool
c ....... Don't allow uptake from storage if forstg is negative -mdh 8/8/00
c ....... If the tree is decidious and the time alloted to grow the
c ....... woody components of the trees has passed use the late season
c ....... growth restriction parameter value to determine how much
c ....... nutrients to flow out of the forest nutrient storage pool for
c ....... the woody components, cak - 03/11/2010
          if (forstg(iel) .gt. 0.0) then
            amt = uptake(ESTOR,iel) * euf(LEAF)
            call flow(forstg(iel),rleave(iel),time,amt)
            eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
            amt = uptake(ESTOR,iel) * euf(FROOTJ)
            call flow(forstg(iel),frootej(iel),time,amt)
            eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
            amt = uptake(ESTOR,iel) * euf(FROOTM)
            call flow(forstg(iel),frootem(iel),time,amt)
            eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
            if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys)) then
              amt = uptake(ESTOR,iel) * euf(FBRCH) * (1.0 - flsgres)
              call flow(forstg(iel),fbrche(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ESTOR,iel) * euf(LWOOD) * (1.0 - flsgres)
              call flow(forstg(iel),rlwode(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ESTOR,iel) * euf(CROOT) * (1.0 - flsgres)
              call flow(forstg(iel),croote(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            else
              amt = uptake(ESTOR,iel) * euf(FBRCH)
              call flow(forstg(iel),fbrche(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ESTOR,iel) * euf(LWOOD)
              call flow(forstg(iel),rlwode(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ESTOR,iel) * euf(CROOT)
              call flow(forstg(iel),croote(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            endif
          endif

c ....... Take up nutrients from soil
c ....... Nutrients for uptake are available in the top tlaypg layers,
c ....... cak 01/29/03
c ....... If the tree is decidious and the time alloted to grow the
c ....... woody components of the trees has passed flow nutrients to
c ....... the forest storage pool rather than to the component nutrient
c ....... pools based on the forest late season growth restriction
c ....... parameter value, cak - 03/11/2010
          do 100 lyr = 1, tlaypg
            if (minerl(lyr,iel) .gt. toler) then
              fsol = 1.0
c ........... The fsol calculation for P is not needed here to compute
c ........... the weighted average, cak - 04/05/02
c              if (iel .eq. P) then
c                fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
c              endif
              call cmpnfrac(lyr,ammonium,nitrate,minerl,
     &                      frac_nh4,frac_no3)
              calcup = uptake(ESOIL,iel)*minerl(lyr,iel)*
     &                 fsol/availm(iel)
c ........... Leaves
              amt = calcup * euf(LEAF)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),rleave(iel),time,amt)
              eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
c ........... Juvenile fine roots
              amt = calcup * euf(FROOTJ)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),frootej(iel),time,amt)
              eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
c ........... Mature fine roots
              amt = calcup * euf(FROOTM)
              if (iel .eq. N) then
                namt = -1.0*amt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
              call flow(minerl(lyr,iel),frootem(iel),time,amt)
              eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
c ........... Fine branch
              namt = 0.0
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = calcup * euf(FBRCH) * flsgres
                namt = namt + amt
                call flow(minerl(lyr,iel),forstg(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = calcup * euf(FBRCH) * (1.0 - flsgres)
                namt = namt + amt
                call flow(minerl(lyr,iel),fbrche(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = calcup * euf(FBRCH)
                namt = namt + amt
                call flow(minerl(lyr,iel),fbrche(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
              if (iel .eq. N) then
                namt = -1.0*namt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
c ........... Large wood
              namt = 0.0
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = calcup * euf(LWOOD) * flsgres
                namt = namt + amt
                call flow(minerl(lyr,iel),forstg(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = calcup * euf(LWOOD) * (1.0 - flsgres)
                namt = namt + amt
                call flow(minerl(lyr,iel),rlwode(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = calcup * euf(LWOOD)
                namt = namt + amt
                call flow(minerl(lyr,iel),rlwode(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
              if (iel .eq. N) then
                namt = -1.0*namt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
c ........... Coarse roots
              namt = 0.0
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = calcup * euf(CROOT) * flsgres
                namt = namt + amt
                call flow(minerl(lyr,iel),forstg(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = calcup * euf(CROOT) * (1.0 - flsgres)
                namt = namt + amt
                call flow(minerl(lyr,iel),croote(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = calcup * euf(CROOT)
                namt = namt + amt
                call flow(minerl(lyr,iel),croote(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
              if (iel .eq. N) then
                namt = -1.0*namt
                call update_npool(lyr, namt, frac_nh4, frac_no3, 
     &                            ammonium, nitrate, subname)
              endif
            endif
100       continue

c ....... Take up nutrients from nitrogen fixation
c ....... If the tree is decidious and the time alloted to grow the
c ....... woody components of the trees has passed flow nutrients to
c ....... the forest storage pool rather than to the component nutrient
c ....... pools based on the forest late season growth restriction
c ....... parameter value, cak - 03/11/2010
          if (iel .eq. N .and. treeNfix .gt. 0) then
c ......... Leaves
            amt = uptake(ENFIX,iel) * euf(LEAF)
            call flow(esrsnk(iel),rleave(iel),time,amt)
            eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
c ......... Juvenile fine roots
            amt = uptake(ENFIX,iel) * euf(FROOTJ)
            call flow(esrsnk(iel),frootej(iel),time,amt)
            eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
c ......... Mature fine roots
            amt = uptake(ENFIX,iel) * euf(FROOTM)
            call flow(esrsnk(iel),frootem(iel),time,amt)
            eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
c ......... Fine branch
            if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
              amt = uptake(ENFIX,iel) * euf(FBRCH) * flsgres
              call flow(esrsnk(iel),forstg(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ENFIX,iel) * euf(FBRCH) * (1.0 - flsgres)
              call flow(esrsnk(iel),fbrche(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            else
              amt = uptake(ENFIX,iel) * euf(FBRCH)
              call flow(esrsnk(iel),fbrche(iel),time,amt)
              eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            endif
c ......... Large wood
            if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
              amt = uptake(ENFIX,iel) * euf(LWOOD) * flsgres
              call flow(esrsnk(iel),forstg(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ENFIX,iel) * euf(LWOOD) * (1.0 - flsgres)
              call flow(esrsnk(iel),rlwode(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            else
              amt = uptake(ENFIX,iel) * euf(LWOOD)
              call flow(esrsnk(iel),rlwode(iel),time,amt)
              eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            endif
c ......... Coarse roots
            if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
              amt = uptake(ENFIX,iel) * euf(CROOT) * flsgres
              call flow(esrsnk(iel),forstg(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ENFIX,iel) * euf(CROOT) * (1.0 - flsgres)
              call flow(esrsnk(iel),croote(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            else
              amt = uptake(ENFIX,iel) * euf(CROOT)
              call flow(esrsnk(iel),croote(iel),time,amt)
              eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            endif
          endif

c ....... Take up nutrients from automatic fertilizer
c ....... If the tree is decidious and the time alloted to grow the
c ....... woody components of the trees has passed flow nutrients to
c ....... the forest storage pool rather than to the component nutrient
c ....... pools based on the forest late season growth restriction
c ....... parameter value, cak - 03/11/2010
          if (aufert .ne. 0) then
            if (uptake(EFERT,iel) .gt. 0.) then
c ........... Automatic fertilizer added to plant pools
c ........... Leaves
              amt = uptake(EFERT,iel) * euf(LEAF)
              call flow(esrsnk(iel),rleave(iel),time,amt)
              eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
c ........... Juvenile fine roots
              amt = uptake(EFERT,iel) * euf(FROOTJ)
              call flow(esrsnk(iel),frootej(iel),time,amt)
              eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
c ........... Mature fine roots
              amt = uptake(EFERT,iel) * euf(FROOTM)
              call flow(esrsnk(iel),frootem(iel),time,amt)
              eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
c ........... Fine branch
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = uptake(EFERT,iel) * euf(FBRCH) * flsgres
                call flow(esrsnk(iel),forstg(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = uptake(EFERT,iel) * euf(FBRCH) * (1.0 - flsgres)
                call flow(esrsnk(iel),fbrche(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(FBRCH)
                call flow(esrsnk(iel),fbrche(iel),time,amt)
                eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
c ........... Large wood
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = uptake(EFERT,iel) * euf(LWOOD) * flsgres
                call flow(esrsnk(iel),forstg(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = uptake(EFERT,iel) * euf(LWOOD) * (1.0 - flsgres)
                call flow(esrsnk(iel),rlwode(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(LWOOD)
                call flow(esrsnk(iel),rlwode(iel),time,amt)
                eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
c ........... Coarse roots
              if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys))then
                amt = uptake(EFERT,iel) * euf(CROOT) * flsgres
                call flow(esrsnk(iel),forstg(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = uptake(EFERT,iel) * euf(CROOT) * (1.0 - flsgres)
                call flow(esrsnk(iel),croote(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(CROOT)
                call flow(esrsnk(iel),croote(iel),time,amt)
                eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
c ........... Automatic fertilizer added to mineral pool
C             if (favail(iel) .eq. 0.0) then
              if (favail(iel) .le. 0.0) then
                write(*,*) 'Error in treegrow, favail(iel) = 0'
                STOP
              endif
              amt = uptake(EFERT,iel) * (1.0/favail(iel) - 1.0)
c ........... For now we are only accumulating fertilizer additions
c ........... for crop/grass, cak - 06/06/2008
c              fertot(iel) = fertot(iel) + uptake(EFERT,iel) + amt
c              fertac(iel) = fertac(iel) + uptake(EFERT,iel) + amt
              if (iel .eq. N) then
                lyr = SRFC
                call update_npool(lyr, amt, frac_nh4_fert, 
     &                            frac_no3_fert, ammonium, nitrate,
     &                            subname)
              endif
              call flow(esrsnk(iel),minerl(SRFC,iel),time,amt)
            endif
          endif
110     continue

c ... Else there is no production this time step due to nutrient
c ... limitation
      else  
        cprodfdy = 0.0
        do 140 iel = 1, MAXIEL
          eprodfdy(iel) = 0.0
          do 130 ipart = 1, FPARTS
            eup(ipart,iel) = 0.0
130       continue
140     continue
      endif

c ... Accumulate monthly output
      mrspflux(FORSYS) = mrspflux(FORSYS) + mrspdyflux(FORSYS)
      grspflux(FORSYS) = grspflux(FORSYS) + grspdyflux(FORSYS)
      do 150 ipart = 1, FPARTS
        sumrsp = sumrsp + fmrspdyflux(ipart)
        fmrspflux(ipart) = fmrspflux(ipart) + fmrspdyflux(ipart)
        mrspmth(FORSYS) = mrspmth(FORSYS) + fmrspdyflux(ipart)
        fgrspflux(ipart) = fgrspflux(ipart) + fgrspdyflux(ipart)
        grspmth(FORSYS) = grspmth(FORSYS) + fgrspdyflux(ipart)
c ..... Accumulate annual output
        mrspann(FORSYS) = mrspann(FORSYS) + fmrspdyflux(ipart)
        grspann(FORSYS) = grspann(FORSYS) + fgrspdyflux(ipart)
150   continue

c ... Accumulate monthly and annual output for soil respiration
      srspmth(FORSYS) = srspmth(FORSYS) + fmrspdyflux(FROOTJ) +
     &                  fmrspdyflux(FROOTM) + fmrspdyflux(CROOT) +
     &                  fgrspdyflux(FROOTJ) + fgrspdyflux(FROOTM) +
     &                  fgrspdyflux(CROOT)
      srspann(FORSYS) = srspann(FORSYS) + fmrspdyflux(FROOTJ) +
     &                  fmrspdyflux(FROOTM) + fmrspdyflux(CROOT) +
     &                  fgrspdyflux(FROOTJ) + fgrspdyflux(FROOTM) +
     &                  fgrspdyflux(CROOT)

c ... If the tree is decidious and the time alloted to grow the
c ... woody components of the trees has passed and the late season
c ... growth restriction parameter value is used to determine how much
c ... carbohydrate to flow out of the forest carbohydrate storage
c ... pool for the woody components reset the wood production
c ... components so the output accumulators are tracking correctly
      if ((decid .eq. 1) .and. (fgrwdys .gt. furgdys)) then
        mfprd(FBRCH) = mfprd(FBRCH) * (1.0 - flsgres)
        mfprd(LWOOD) = mfprd(LWOOD) * (1.0 - flsgres)
        mfprd(CROOT) = mfprd(CROOT) * (1.0 - flsgres)
        cprodfdy = mfprd(LEAF) + mfprd(FROOTJ) + mfprd(FROOTM) +
     &             mfprd(FBRCH) + mfprd(LWOOD) + mfprd(CROOT)
      endif

c ... Sum the daily production variables for output to the monthly *.bin file
      cprodf = cprodf + cprodfdy
      do 160 iel = 1, nelem
        eprodf(iel) = eprodf(iel) + eprodfdy(iel)
160   continue

      return
      end
