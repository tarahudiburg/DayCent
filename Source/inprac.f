
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine inprac(system)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... Argument declarations
      integer system

c ... Initialize annual production accumulators.

c ... Local variables
      integer ii, jj

c ... In the CROP system or SAVANNA, if it is the last month
c ... of the growing season reset the accumulators for the grasses.
      if (dolast .and. crpgrw .eq. 0 .and. system .eq. CRPSYS) then
c ..... Aboveground carbon production
        agcisa(UNLABL) = 0.0
        agcisa(LABELD) = 0.0
        agcprd = agcacc
        agcacc = 0.0
        ptagc = 0.0
c ..... Belowground carbon production
        bgcisja(UNLABL) = 0.0
        bgcisja(LABELD) = 0.0
        bgcjprd = bgcjacc
        bgcjacc = 0.0
        bgcisma(UNLABL) = 0.0
        bgcisma(LABELD) = 0.0
        bgcmprd = bgcmacc
        bgcmacc = 0.0
        ptbgc = 0.0
        do 10 ii = 1, MAXIEL
c ....... N, P, and S uptake by plants
          eupaga(ii) = 0.0
          eupbga(ii) = 0.0
10      continue
c ..... Add code to reset the growing season accumulators for fertilizer
c ..... and organic matter addition and N2O flux, cak - 06/06/2008
c ..... Fertilizer addition
        do 20 ii = 1, MAXIEL
          fertprd(ii) = fertac(ii)
          fertac(ii) = 0.0
          do 30 jj = 1, MONTHS
            fertmth(jj,ii) = 0.0
30        continue
20      continue
c ..... Organic matter addition
        omadprd = omadac
        omadac = 0.0
        do 40 ii = 1, MAXIEL
          omadpre(ii) = omadae(ii)
          omadae(ii) = 0.0
          do 50 jj = 1, MONTHS
            omadmth(jj) = 0.0
            omadmte(jj,ii) = 0.0
50        continue
40      continue
c ..... N2O flux
        n2oprd = n2oacc
        n2oacc = 0.0
        do 60 ii = 1, MONTHS
          n2omth(ii) = 0.0
60      continue
      endif

c ... In the FOREST system or SAVANNA, if it is the last month
c ... of the growing season reset the accumulators for the trees.
      if (doflst .and. forgrw .eq. 0 .and. system .eq. FORSYS) then
c ..... Total forest carbon
        fcprd = fcacc
        fcacc = 0
c ..... Leaf carbon production
        alvcis(UNLABL) = 0.0
        alvcis(LABELD) = 0.0
        rlvprd = rlvacc
        rlvacc = 0.0
c ..... Fine root carbon production
        afrcisj(UNLABL) = 0.0
        afrcisj(LABELD) = 0.0
        frtjprd = frtjacc
        frtjacc = 0.0
        afrcism(UNLABL) = 0.0
        afrcism(LABELD) = 0.0
        frtmprd = frtmacc
        frtmacc = 0.0
c ..... Fine branch carbon production
        afbcis(UNLABL) = 0.0
        afbcis(LABELD) = 0.0
        fbrprd = fbracc
        fbracc = 0.0
c ..... Large wood carbon production
        alwcis(UNLABL) = 0.0
        alwcis(LABELD) = 0.0
        rlwprd = rlwacc
        rlwacc = 0.0
c ..... Coarse root carbon production
        acrcis(UNLABL) = 0.0
        acrcis(LABELD) = 0.0
        crtprd = crtacc
        crtacc = 0.0
c ..... N, P, and S uptake by plants
        do 70 ii = 1, MAXIEL
          do 80 jj = 1, FPARTS-1
            eupprt(jj,ii) = 0.0
80        continue
70      continue
      endif

c ... N, P, and S uptake by plants, reset if no growth is occurring
      if (crpgrw .eq. 0 .and. forgrw .eq. 0) then
        do 90 ii = 1, MAXIEL
          eupprd(ii) = eupacc(ii)
          eupacc(ii) = 0.0
90      continue
      endif

      return
      end
