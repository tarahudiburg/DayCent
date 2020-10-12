
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... STDEDREM

      subroutine stdedrem(accum)

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      accum(ISOS)

c ... New subroutine to simulate the removal of dead attached leaves and 
c ... standing dead wood due to a TREM cutting or fire in a forest. 
c ... -mdh 9/19/2018

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


c ... Remove standing dead LARGE WOOD

      if (dlwodc .gt. 0.001) then
        closs = remf(8) * dlwodc
        tcrem = tcrem + closs
        call csched(closs,dlwcis(LABELD),dlwodc,
     &              dlwcis(UNLABL),csrsnk(UNLABL),
     &              dlwcis(LABELD),csrsnk(LABELD),
     &              1.0,accum)

        do 10 iel = 1, nelem
          eloss(iel) = closs * (dlwode(iel) / dlwodc)
          terem(iel) = terem(iel) + eloss(iel)
          call flow(dlwode(iel),esrsnk(iel),time,eloss(iel))
10      continue
      endif

c ... Remove attached dead FINE BRANCHES

      if (dfbrchc .gt. 0.001) then
        closs = remf(7) * dfbrchc
        tcrem = tcrem + closs
        call csched(closs,dfbrcis(LABELD),dfbrchc,
     &              dfbrcis(UNLABL),csrsnk(UNLABL),
     &              dfbrcis(LABELD),csrsnk(LABELD),
     &              1.0,accum)

        do 20 iel = 1, nelem
          eloss(iel) = closs * (dfbrche(iel) / dfbrchc)
          terem(iel) = terem(iel) + eloss(iel)
          call flow(dfbrche(iel),esrsnk(iel),time,eloss(iel))
20      continue
      endif

c ... Remove attached dead LEAVES
c ... ATTENTION: check that there is standing live or dead biomass that  
c ... leaves can remain attached to? -mdh 9/19/2018
      if (dleavc .gt. 0.001) then
        closs = remf(6) * dleavc
        tcrem = tcrem + closs
        call csched(closs,dlvcis(LABELD),dleavc,
     &              dlvcis(UNLABL),csrsnk(UNLABL),
     &              dlvcis(LABELD),csrsnk(LABELD),
     &              1.0,accum)

        do 30 iel = 1, nelem
          eloss(iel) = closs * (dleave(iel) / dleavc)
          terem(iel) = terem(iel) + eloss(iel)
          call flow(dleave(iel),esrsnk(iel),time,eloss(iel))
30      continue
      endif

      return
      end
