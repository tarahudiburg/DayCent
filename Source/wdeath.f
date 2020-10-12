
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... WDEATH

      subroutine wdeath (tavewk, bgwfunc, tfrac, avgstemp)

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Argument declarations
      real tavewk, bgwfunc, tfrac, avgstemp

c ... Death of leaves, fine branches, large wood, fine roots, and coarse roots.
c ... Modifications:
c ...
c ... Corrected a bug in the death of lableled fine branches, coarse wood, 
c ... and coarse roots. Use the actual labled fraction instead of the growth 
c ... ratio, cisotf. Prevents the removal of to much labeled material when 
c ... non equilibrium, under ratio, wood dies.  7/2007  K. Killian

c ... Fortran to C prototype
      INTERFACE
        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow
      END INTERFACE

c ... Function declarations
      real      carctanf, gpdf, maxswpot
      external  carctanf, gpdf, maxswpot

c ... Local variables
      integer iel, drpdys
      logical drpdlv
      real    accum(ISOS), ctodie, etodie, fr14, recres(MAXIEL),
     &        tostore, srfclittr, soillittr,
     &        rtdh, tempeff, watreff, cturn, eturn, temp
      data    drpdys /0/

c ... Saved variables
      save drpdlv
      save drpdys

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

      if ((time - strtyr) .le. 0.00001) drpdlv = .FALSE.

c ... Death of leaves
c ... NOTE:  WOODDR(1)   - the death rate in fall for deciduous forests
c ...        LEAFDR(MTH) - the monthly death rate for leaves in every
c ...                      case except for fall in deciduous forests.
      if (rleavc .gt. 0.0001) then
        if (decid .ge. 1) then

c ....... Deciduous forest
c ....... If the daylight hours are increasing - it must be spring
          if ((hrsinc) .and. (tavewk .gt. tmplff)) drpdlv = .FALSE.

c ....... If daylight hours are decreasing and the temperature is low
c ....... enough drop leaves for fall season.
c ....... Add check for number of daylight hours to conditional for
c ....... determining if leaf drop should occur, cak - 06/30/03
c ....... If leaf drop has not occurred by the time the winter solstice
c ....... is reached force leaf drop to occur, cak - 10/28/04
          if (decid .eq. 1) then
            if (((tavewk .lt. tmplff) .and. (.not. drpdlv) .and.
     &           (.not. hrsinc) .and. (dayhrs .lt. 12.0)) .or.
     &          ((.not.drpdlv) .and. (sitlat .ge. 0) .and.
     &           (month .eq. 12)) .or.
     &          ((.not.drpdlv) .and. (sitlat .lt. 0) .and.
     &           (month .eq. 6))) then
c ........... Allow leaf drop to occur over a 30 day period, dropping
c ........... all of the remaining leaves on the last day.
              if (drpdys .lt. 30) then
                ctodie = rleavc * wooddr(LEAF)*tfrac
                drpdys = drpdys + 1
                decidgrow = .FALSE.
               else
                ctodie = rleavc * wooddr(LEAF)
                drpdlv = .TRUE.
                drpdys = 0
              endif
            else
              ctodie = rleavc * leafdr(month)*tfrac
            endif
          elseif (decid .eq. 2) then
c ......... Drought deciduous forest
c ......... Compute death for drought deciduous forests
            ctodie = rleavc * (1. - bgwfunc) * wooddr(LEAF)*tfrac
          endif
        else
c ....... Continuous forest
c ....... Use leaf death rate multiplier from EACHYR
          ctodie = rleavc * leafdr(month)*tfrac * ldrmlt
        endif

c ..... Compute E/C ratios
        do 10 iel = 1, nelem
          recres(iel) = rleave(iel) / rleavc

c ....... Compute flow to retranslocation storage
          tostore = recres(iel) * ctodie * forrtf(iel)
          call flow(rleave(iel), forstg(iel), time, tostore)

c ....... Decrease E/C by the amount that is retranslocated
          recres(iel) = recres(iel) * (1 - forrtf(iel))
10      continue

c ..... If evntyp is greater than 1 the leaves go to the source/sink
c ..... rather than to the litter, cak - 02/07/2006
c ..... ATTENTION: What is evntyp >= 2? -mdh 9/19/2018
c ..... Updated so evergreen leaves that die go to attached
c ..... dead leaves.  Deciduous leaves still fall to the ground
c ..... and become litter. -mdh 9/19/2018

c       write(*,*) "wdeath:"
c       write(*,*) 'time=', time, 'evntyp=', evntyp
        fr14 = rlvcis(LABELD) / rleavc
        if (evntyp .lt. 2) then
          if (decid .ge. 1) then
            call partit(ctodie, recres, 1, rlvcis, rleave, wdlig(LEAF),
     &                  fr14)
          else
c ......... Evergreen. Live leaves become dead attached leaves.
c           write(*,*) 'wdeath: time=', time
c           write(*,*) 'wdeath: ctodie=', ctodie
c           write(*,*) 'wdeath: rleavc=', rleavc
c           write(*,*) 'wdeath: dleavc=', dleavc
            call csched(ctodie, fr14, 1.0,
     &                  rlvcis(UNLABL), dlvcis(UNLABL),
     &                  rlvcis(LABELD), dlvcis(LABELD),
     &                  1.0, accum)
            do 15 iel = 1, nelem
              etodie = ctodie * (rleave(iel) / rleavc)
              call flow(rleave(iel), dleave(iel), time,
     &                  etodie - tostore)
15        continue
          endif
        else
c ....... evntyp .ge. 2 
c         write(*,*) 'wdeath: evntyp =', evntyp
          call csched(ctodie, fr14, 1.0,
     &                rlvcis(UNLABL), csrsnk(UNLABL),
     &                rlvcis(LABELD), csrsnk(LABELD),
     &                1.0, accum)
          do 16 iel = 1, nelem
            etodie = ctodie * (rleave(iel) / rleavc)
            call flow(rleave(iel), esrsnk(iel), time,
     &                etodie - tostore)
16        continue
        endif
      endif

c ... Add code to age fine roots, juvenile fine roots age to the mature
c ... fine root pool.  Modify this subroutine so that the death rate of
c ... roots is a function of soil water potential and soil temperature,
c ... cak - 06/28/2007
c ... See:  A Model of Production and Turnover of Roots in Shortgrass Prairie
c ...       Parton, Singh, and Coleman, 1978
c ...       Journal of Applied Ecology
c ... Cap the temperature effect on fine roots at -2 and +28 degrees C
      if (avgstemp .gt. 28.0) then
        temp = 28.0
      else if (avgstemp .lt. -2.0) then
        temp = -2.0
      else
        temp = avgstemp
      endif

c ... Soil temperature effect on aging of juvenile roots
      if (frootcj .gt. 0.0) then
        tempeff = gpdf(temp, 37.0, 0.0, 3.0, 3.0)
        cturn = tmxturn * tempeff * frootcj * tfrac
        fr14 = frtcisj(LABELD) / frootcj
        call csched(cturn, fr14, 1.0,
     &              frtcisj(UNLABL), frtcism(UNLABL),
     &              frtcisj(LABELD), frtcism(LABELD),
     &              1.0, accum)
        do 60 iel = 1, nelem
          eturn = cturn * (frootej(iel) / frootcj)
          call flow(frootej(iel), frootem(iel), time, eturn)
60      continue
      endif

c ... Soil temperature effect on root death rate
      tempeff = (temp - 10.0)**2 / 4.0 * 0.00175 + 0.1
      tempeff = min(tempeff, 0.5)
c ... Soil water potential effect on root death rate
      watreff = maxswpot(tlaypg)
      watreff = carctanf(watreff, 35.0, 0.5, 1.0, 0.05)
      rtdh = max(tempeff, watreff)

c ... Death of juvenile fine roots
      if (frootcj .gt. 0.0) then
        ctodie = frootcj * wooddr(FROOTJ) * rtdh * tfrac
        do 20 iel = 1, nelem
          recres(iel) = frootej(iel) / frootcj
20      continue
        fr14 = frtcisj(LABELD) / frootcj
c ..... A fraction of the dead roots are transferred to the surface
c ..... litter layer, the remainder goes to the soil litter layer
c ..... cak - 05/14/2007
        srfclittr = ctodie * wrdsrfc
        soillittr = ctodie - srfclittr
        call partit(srfclittr, recres, SRFC, frtcisj, frootej,
     &              wdlig(FROOTJ), fr14)
        call partit(soillittr, recres, SOIL, frtcisj, frootej,
     &              wdlig(FROOTJ), fr14)
      endif

c ... Death of mature fine roots
      if (frootcm .gt. 0.0) then
        ctodie = frootcm * wooddr(FROOTM) * rtdh * tfrac
        do 25 iel = 1, nelem
          recres(iel) = frootem(iel) / frootcm
25      continue
        fr14 = frtcism(LABELD) / frootcm
c ..... A fraction of the dead roots are transferred to the surface
c ..... litter layer, the remainder goes to the soil litter layer
c ..... cak - 05/14/2007
        srfclittr = ctodie * wrdsrfc
        soillittr = ctodie - srfclittr
        call partit(srfclittr, recres, SRFC, frtcism, frootem,
     &              wdlig(FROOTM), fr14)
        call partit(soillittr, recres, SOIL, frtcism, frootem,
     &              wdlig(FROOTM), fr14)
      endif

c ... Fine Branches, Large Wood, and Coarse Roots go to the dead wood
c ... compartments: WOOD1, WOOD2, WOOD3

C ... Death of fine branches
C     if (fbrchc .gt. 0.0) then
C       ctodie = fbrchc * wooddr(FBRCH)*tfrac
C ..... remove labled fraction instead of the growth ratio, cisotf.  
C ..... KLK 7/2007
C        call csched(ctodie, cisotf, 1.0,
C       call csched(ctodie, fbrcis(LABELD), fbrchc,
C    &              fbrcis(UNLABL), wd1cis(UNLABL),
C    &              fbrcis(LABELD), wd1cis(LABELD),
C    &              1.0, accum)
C
C       do 30 iel = 1, nelem
C         etodie = ctodie * (fbrche(iel) / fbrchc)
C         call flow(fbrche(iel), wood1e(iel), time, etodie)
C30      continue
C     endif
C
C ... Death of large wood
C     if (rlwodc .gt. 0.0) then
C       ctodie = rlwodc * wooddr(LWOOD)*tfrac
C ..... remove labled fraction instead of the growth ratio, cisotf.  
C ..... KLK 7/2007
C        call csched(ctodie, cisotf, 1.0,
C       call csched(ctodie, rlwcis(LABELD), rlwodc,
C    &              rlwcis(UNLABL), wd2cis(UNLABL),
C    &              rlwcis(LABELD), wd2cis(LABELD),
C    &              1.0, accum)
C
C       do 40 iel = 1, nelem
C         etodie = ctodie * (rlwode(iel) / rlwodc)
C         call flow(rlwode(iel), wood2e(iel), time, etodie)
C40      continue
C     endif
 
c ------------------------------------------------------------------
c ... Death of fine branches and large wood flows into 
c ... standing wood pools instead of surface wood pools.
c ... -mdh 9/19/2018

c ... Death of fine branches
      if (fbrchc .gt. 0.0) then
        ctodie = fbrchc * wooddr(FBRCH)*tfrac
        call csched(ctodie, fbrcis(LABELD), fbrchc,
     &              fbrcis(UNLABL), dfbrcis(UNLABL),
     &              fbrcis(LABELD), dfbrcis(LABELD),
     &              1.0, accum)

        do 30 iel = 1, nelem
          etodie = ctodie * (fbrche(iel) / fbrchc)
          call flow(fbrche(iel), dfbrche(iel), time, etodie)
30      continue
      endif

c ... Death of large wood
      if (rlwodc .gt. 0.0) then
        ctodie = rlwodc * wooddr(LWOOD)*tfrac
        call csched(ctodie, rlwcis(LABELD), rlwodc,
     &              rlwcis(UNLABL), dlwcis(UNLABL),
     &              rlwcis(LABELD), dlwcis(LABELD),
     &              1.0, accum)

        do 40 iel = 1, nelem
          etodie = ctodie * (rlwode(iel) / rlwodc)
          call flow(rlwode(iel), dlwode(iel), time, etodie)
40      continue
      endif
c ------------------------------------------------------------------

c ... Death of coarse roots
      if (crootc .gt. 0.0) then
        ctodie = crootc * wooddr(CROOT)*tfrac
c ..... remove labled fraction instead of the growth ratio, cisotf.  
c ..... KLK 7/2007
c        call csched(ctodie, cisotf, 1.0,
        call csched(ctodie, crtcis(LABELD), crootc,
     &              crtcis(UNLABL), wd3cis(UNLABL),
     &              crtcis(LABELD), wd3cis(LABELD),
     &              1.0, accum)

        do 50 iel = 1, nelem
          etodie = ctodie * (croote(iel) / crootc)
          call flow(croote(iel), wood3e(iel), time, etodie)
50      continue
      endif

      return
      end
