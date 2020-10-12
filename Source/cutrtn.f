
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... CUTRTN

      subroutine cutrtn(accum)

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      accum(ISOS)

c ... Return of leaves as litter, and wood to surface dead wood 
c ....during a TREM cutting event.
c ... Added dead attached leaves, dead attached fine branches,
c ... and standing dead wood. -mdh 9/19/2018

c ... Called from:  frem

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

c ... Local Variables
      integer   iel
c     real      cgain, egain(MAXIEL)
      real      frc14, recres(MAXIEL)
      real      cret, eret(MAXIEL)

      cret = 0.0
      do 5 iel = 1, MAXIEL
        eret(iel) = 0.0
5     continue

c ... LIVE LEAVES are returned to LITTER
      if (rleavc .gt. 0.001) then
c       cgain = remf(1) * retf(1,1) * rleavc
        cret = remf(1) * (1.0 - lv2std(1)) * retf(1,1) * rleavc
        tcreta = tcreta + cret
        if (cret .gt. 0.0) then
          do 10 iel = 1, nelem
c           egain(iel) = remf(1) * retf(1,iel+1) * rleave(iel)
            eret(iel) = remf(1) * (1.0 - lv2std(1)) * retf(1,iel+1) 
     &                  * rleave(iel)
            tereta(iel) = tereta(iel) + eret(iel)
            recres(iel) = eret(iel) / cret
10        continue
          frc14 = rlvcis(LABELD) / rleavc
          call partit(cret,recres,1,csrsnk,esrsnk,wdlig(LEAF),frc14)
        endif
      endif

c ... DEAD ATTACHED LEAVES are returned to LITTER
      if (dleavc .gt. 0.001) then
        cret = remf(6) * retf(4,1) * dleavc
        tcreta = tcreta + cret
        if (cret .gt. 0.0) then
          do 15 iel = 1, nelem
            eret(iel) = remf(6) * retf(4,iel+1) * dleave(iel)
            tereta(iel) = tereta(iel) + eret(iel)
            recres(iel) = eret(iel) / cret
15        continue
          frc14 = dlvcis(LABELD) / dleavc
          call partit(cret,recres,1,csrsnk,esrsnk,wdlig(LEAF),frc14)
        endif
      endif

c ... LIVE FINE BRANCHES go to DEAD SURFACE FINE BRANCHES
      if (fbrchc .gt. 0.001) then
c       cgain = remf(2) * retf(2,1) * fbrchc
        cret = remf(2) * (1.0 - lv2std(2)) * retf(2,1) * fbrchc
        tcreta = tcreta + cret
        call csched(cret,fbrcis(LABELD),fbrchc,
     &              csrsnk(UNLABL),wd1cis(UNLABL),
     &              csrsnk(LABELD),wd1cis(LABELD),
     &              1.0,accum)
        do 20 iel = 1, nelem
c         egain(iel) = remf(2) * retf(2,iel+1) * fbrche(iel)
          eret(iel) = remf(2) * (1.0 - lv2std(2)) * retf(2,iel+1) 
     &                * fbrche(iel)
          tereta(iel) = tereta(iel) + eret(iel)
          call flow(esrsnk(iel),wood1e(iel),time,eret(iel))
20      continue
      endif

c ... DEAD ATTACHED FINE BRANCHES go to DEAD SURFACE FINE BRANCHES
      if (dfbrchc .gt. 0.001) then
        cret = remf(7) * retf(5,1) * dfbrchc
        tcreta = tcreta + cret
        call csched(cret,dfbrcis(LABELD),dfbrchc,
     &              csrsnk(UNLABL),wd1cis(UNLABL),
     &              csrsnk(LABELD),wd1cis(LABELD),
     &              1.0,accum)
        do 25 iel = 1, nelem
          eret(iel) = remf(5) * retf(5,iel+1) * dfbrche(iel)
          tereta(iel) = tereta(iel) + eret(iel)
          call flow(esrsnk(iel),wood1e(iel),time,eret(iel))
25      continue
      endif

c ... LIVE LARGE WOOD goes to DEAD SURFACE LARGE WOOD
      if (rlwodc .gt. 0.001) then
c       cgain = remf(3) * retf(3,1) * rlwodc
        cret = remf(3) * (1.0 - lv2std(3)) * retf(3,1) * rlwodc
        tcreta = tcreta + cret
        call csched(cret,rlwcis(LABELD),rlwodc,
     &              csrsnk(UNLABL),wd2cis(UNLABL),
     &              csrsnk(LABELD),wd2cis(LABELD),
     &              1.0,accum)
        do 30 iel = 1, nelem
c         egain(iel) = remf(3) * retf(3,iel+1) * rlwode(iel)
          eret(iel) = remf(3) * (1.0 - lv2std(3)) * retf(3,iel+1) 
     &                * rlwode(iel)
          tereta(iel) = tereta(iel) + eret(iel)
          call flow(esrsnk(iel),wood2e(iel),time,eret(iel))
30      continue
      endif

c ... DEAD STANDING LARGE WOOD goes to DEAD SURFACE LARGE WOOD
      if (dlwodc .gt. 0.001) then
        cret = remf(8) * retf(6,1) * dlwodc
        tcreta = tcreta + cret
        call csched(cret,dlwcis(LABELD),dlwodc,
     &              csrsnk(UNLABL),wd2cis(UNLABL),
     &              csrsnk(LABELD),wd2cis(LABELD),
     &              1.0,accum)
        do 35 iel = 1, nelem
          eret(iel) = remf(8) * retf(6,iel+1) * dlwode(iel)
          tereta(iel) = tereta(iel) + eret(iel)
          call flow(esrsnk(iel),wood2e(iel),time,eret(iel))
35      continue
      endif

c ... Add STORAGE back
      do 40 iel = 1, nelem
c       egain(iel) = remf(3) * retf(3,iel+1) * forstg(iel)
c       call flow(esrsnk(iel),metabe(SRFC,iel),time,egain(iel))
        eret(iel) = remf(3) * retf(3,iel+1) * forstg(iel)
        tereta(iel) = tereta(iel) + eret(iel)
        call flow(esrsnk(iel),metabe(SRFC,iel),time,eret(iel))
40    continue

      return
      end
