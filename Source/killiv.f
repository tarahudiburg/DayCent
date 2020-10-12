
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... KILLIV

      subroutine killiv(accum)

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'zztim.inc'

c ... Argument declarations
      real     accum(ISOS)

c ... New subroutine to transfer above ground live tree parts to their 
c ... standing dead counterparts when trees die during a TREM fire or 
c ... cutting event. -mdh 9/19/2018

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
      real      closs, eloss(MAXIEL)

c ... Live leaves are transferred to dead attached leaves
      if (rleavc .gt. 0.001) then
        closs = remf(1) * lv2std(1) * rleavc
c ..... tcrem = tcrem + closs
        call csched(closs,rlvcis(LABELD),rleavc,
     &              rlvcis(UNLABL),dlvcis(UNLABL),
     &              rlvcis(LABELD),dlvcis(LABELD),
     &              1.0,accum)

        do 10 iel = 1, nelem
          eloss(iel) = closs * (rleave(iel) / rleavc)
c ....... terem(iel) = terem(iel) + eloss(iel)
          call flow(rleave(iel),dleave(iel),time,eloss(iel))
10      continue
      endif

c ... Live fine branches are transferred to dead attached fine branches
      if (fbrchc .gt. 0.001) then
        closs = remf(2) * lv2std(2) * fbrchc
c ..... tcrem = tcrem + closs
        call csched(closs,fbrcis(LABELD),fbrchc,
     &              fbrcis(UNLABL),dfbrcis(UNLABL),
     &              fbrcis(LABELD),dfbrcis(LABELD),
     &              1.0,accum)

        do 20 iel = 1, nelem
          eloss(iel) = closs * (fbrche(iel) / fbrchc)
c ....... terem(iel) = terem(iel) + eloss(iel)
          call flow(fbrche(iel),dfbrche(iel),time,eloss(iel))
20      continue
      endif


c ... Live large wood is transferred to standing dead large wood
      if (rlwodc .gt. 0.001) then
        closs = remf(3) * lv2std(3) * rlwodc
c ..... tcrem = tcrem + closs
        call csched(closs,rlwcis(LABELD),rlwodc,
     &              rlwcis(UNLABL),dlwcis(UNLABL),
     &              rlwcis(LABELD),dlwcis(LABELD),
     &              1.0,accum)

        do 30 iel = 1, nelem
          eloss(iel) = closs * (rlwode(iel) / rlwodc)
c ....... terem(iel) = terem(iel) + eloss(iel)
          call flow(rlwode(iel),dlwode(iel),time,eloss(iel))
30      continue
      endif

      return
      end
