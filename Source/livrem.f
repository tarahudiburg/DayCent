
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... LIVREM

      subroutine livrem(accum)

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      accum(ISOS)

c ... Removal of live biomass due to cutting or fire in a forest.

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

c ... Local variables
      integer   iel
      real      closs, eloss(MAXIEL),
     &          liveCtotal, liveCremoved, liveCfrac, mrspstgLoss,
     &          mRespStorage

      liveCtotal = fbrchc + rlwodc
      liveCremoved = 0.0

c ... Remove live LEAVES

      if (rleavc .gt. 0.0) then
        closs = remf(1) * (1.0 - lv2std(1)) * rleavc
        tcrem = tcrem + closs
        call csched(closs,rlvcis(LABELD),rleavc,
     &              rlvcis(UNLABL),csrsnk(UNLABL),
     &              rlvcis(LABELD),csrsnk(LABELD),
     &              1.0,accum)

        do 10 iel = 1, nelem
          eloss(iel) = closs * (rleave(iel) / rleavc)
          terem(iel) = terem(iel) + eloss(iel)
          call flow(rleave(iel),esrsnk(iel),time,eloss(iel))
10      continue
      endif

c ... Remove live FINE BRANCHES

      if (fbrchc .gt. 0.0) then
        closs = remf(2) * (1.0 - lv2std(2)) * fbrchc
        tcrem = tcrem + closs
c ..... Add calculation of liveCremoved for maintenance respiration,
c ..... mdh - 05/01
        liveCremoved = liveCremoved + closs
        call csched(closs,fbrcis(LABELD),fbrchc,
     &              fbrcis(UNLABL),csrsnk(UNLABL),
     &              fbrcis(LABELD),csrsnk(LABELD),
     &              1.0,accum)

        do 20 iel = 1, nelem
          eloss(iel) = closs * (fbrche(iel) / fbrchc)
          terem(iel) = terem(iel) + eloss(iel)
          call flow(fbrche(iel),esrsnk(iel),time,eloss(iel))
20      continue
      endif

c ... Remove live LARGE WOOD

      if (rlwodc .gt. 0.0) then
        closs = remf(3) * (1.0 - lv2std(3)) * rlwodc
        tcrem = tcrem + closs
c ..... Add calculation of liveCremoved for maintenance respiration,
c ..... mdh - 05/01
        liveCremoved = liveCremoved + closs
        call csched(closs,rlwcis(LABELD),rlwodc,
     &              rlwcis(UNLABL),csrsnk(UNLABL),
     &              rlwcis(LABELD),csrsnk(LABELD),
     &              1.0,accum)

        do 30 iel = 1, nelem
          eloss(iel) = closs * (rlwode(iel) / rlwodc)
          terem(iel) = terem(iel) + eloss(iel)
          call flow(rlwode(iel),esrsnk(iel),time,eloss(iel))
30      continue

c ..... Remove C from carbohydrate storage pool based on fraction of
c ..... total live wood removed, mdh - 7/9/01
        liveCfrac = liveCremoved / liveCtotal
        if (carbostg(FORSYS,UNLABL) .lt. 0.0) then
          write(*,*) 'Warning in livrem, carbostg(FORSYS,UNLABL) < 0.0'
c         STOP
        endif
        if (carbostg(FORSYS,LABELD) .lt. 0.0) then
          write(*,*) 'Warning in livrem, carbostg(FORSYS,LABELD) < 0.0'
c         STOP
        endif
        mRespStorage = carbostg(FORSYS,UNLABL)+carbostg(FORSYS,LABELD)
        mrspstgLoss = max(liveCfrac * mRespStorage, 0.0)
        if (mRespStorage .gt. 0.001) then
          call csched(mrspstgLoss, carbostg(FORSYS,LABELD),mRespStorage,
     &                carbostg(FORSYS,UNLABL), csrsnk(UNLABL),
     &                carbostg(FORSYS,LABELD), csrsnk(LABELD),
     &                1.0, accum)
        endif

c ..... Remove from STORAGE pool based on fraction of large wood
c ..... based on fraction killed (remf(3)) rather than just removed.

        do 40 iel = 1, nelem
          eloss(iel) = MAX(remf(3) * forstg(iel), 0.0)
          call flow(forstg(iel),esrsnk(iel),time,eloss(iel))
40      continue
      endif

      return
      end
