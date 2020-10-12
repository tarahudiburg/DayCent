
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine wfalstd(tfrac)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'zztim.inc'

c ... Argument declarations
      real tfrac

c ... New subroutine to simulate the fall of attached dead leaves and 
c ... standing dead wood. -mdh 9/19/2018
c ... tfrac - fraction of month over which this event occurs (~1/30 for daily events)
c ... dlvfalrt - fall rate of dead attached leaves (fraction per month, in absence of disturbance, tree.100)
c ... dfbfalrt - fall rate of dead attached fine branches (fraction per month, in absence of disturbance, tree.100)
c ... dlwfalrt - fall rate of dead standing large wood (fraction per month, in absence of disturbance, tree.100)

c ... Local variables
      integer   iel
      real      fr14, fsdc, fsde, recres(MAXIEL), accum(ISOS)

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Fall of dead attached leaves
      if (dleavc .gt. 0) then
        fsdc = dleavc * dlvfalrt * tfrac
        do 10 iel = 1, nelem
          recres(iel) = dleave(iel)/dleavc
10      continue
        fr14 = dlvcis(LABELD)/dleavc
        call partit(fsdc,recres,1,dlvcis,dleave,wdlig(LEAF),fr14)
      endif

c ... Fall of dead attached fine branches
      if (dfbrchc .gt. 0) then
        fsdc = dfbrchc * dfbfalrt * tfrac
        call csched(fsdc, dfbrcis(LABELD), dfbrchc,
     &              dfbrcis(UNLABL), wd1cis(UNLABL),
     &              dfbrcis(LABELD), wd1cis(LABELD),
     &              1.0, accum)

        do 20 iel = 1, nelem
          fsde = fsdc * (dfbrche(iel) / dfbrchc)
          call flow(dfbrche(iel), wood1e(iel), time, fsde)
20      continue
      endif

c ... Fall of standing dead large wood
      if (dlwodc .gt. 0) then
        fsdc = dlwodc * dlwfalrt * tfrac
        call csched(fsdc, dlwcis(LABELD), dlwodc,
     &              dlwcis(UNLABL), wd2cis(UNLABL),
     &              dlwcis(LABELD), wd2cis(LABELD),
     &              1.0, accum)

        do 30 iel = 1, nelem
          fsde = fsdc * (dlwode(iel) / dlwodc)
          call flow(dlwode(iel), wood2e(iel), time, fsde)
30      continue
      endif

      return
      end
