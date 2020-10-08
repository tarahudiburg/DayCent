
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... KILLRT

      subroutine killrt(accum)

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'plot3.inc'
      include 'zztim.inc'

c ... Argument declarations
      real     accum(ISOS), srfclittr, soillittr

c ... Death of roots due to cutting or fire in a forest.

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
      integer  iel
      real     crd, dethe, frc14, frd, recres(MAXIEL)

c ... Death of juvenile FINE ROOTS

      if (frootcj .gt. 0.001) then
        frd = frootcj * fd(1)
        do 10 iel = 1, nelem
          recres(iel) = frootej(iel) / frootcj
10      continue

        frc14 = frtcisj(LABELD) / frootcj
        srfclittr = frd * wrdsrfc
        soillittr = frd - srfclittr
        call partit(srfclittr,recres,SRFC,frtcisj,frootej,wdlig(FROOTJ),
     &              frc14)
        call partit(soillittr,recres,SOIL,frtcisj,frootej,wdlig(FROOTJ),
     &              frc14)
      endif

c ... Death of mature FINE ROOTS

      if (frootcm .gt. 0.001) then
        frd = frootcm * fd(1)
        do 15 iel = 1, nelem
          recres(iel) = frootem(iel) / frootcm
15      continue

        frc14 = frtcism(LABELD) / frootcm
        srfclittr = frd * wrdsrfc
        soillittr = frd - srfclittr
        call partit(srfclittr,recres,SRFC,frtcism,frootem,wdlig(FROOTM),
     &              frc14)
        call partit(soillittr,recres,SOIL,frtcism,frootem,wdlig(FROOTM),
     &              frc14)
      endif

c ... Death of COARSE ROOTS

      if (crootc .gt. 0.001) then
        crd = crootc * fd(2)
        do 20 iel = 1, nelem
          dethe = crd * (croote(iel)/crootc)
          call flow(croote(iel),wood3e(iel),time,dethe)
20      continue

        call csched(crd,crtcis(LABELD),crootc,
     &              crtcis(UNLABL),wd3cis(UNLABL),
     &              crtcis(LABELD),wd3cis(LABELD),
     &              1.0,accum)
      endif

      return
      end
