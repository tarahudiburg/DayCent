
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine eachyr

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'isovar.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ... Perform tasks that only need to be done once a year.

c ... Function declarations
      real     fracis, line
      external fracis, line

c ... Local variables
      integer  iel, ipart, mon
      real     lfncmax, lfncmin, lfncon

c ... Added for savanna model (plot3.inc) BO

c ... Correct for the rounding problem with time. The value for time
c ... drifts downward during long runs since dt=1/12 cannot be represented
c ... precisely.  At this point in the run, time should be a whole number.

c ... Changed increment value from .001 to .5 to correct error in time calculation
c ... occuring after year 8192. (mse 3/95).  New code if from Kendrick Killian.

      time = sign(int(abs(time)+.5),int(time))

c ... Reset annual accumulators to zero
      call annacc

c ... Set climate scalars as indicated, cak - 10/18/05
      if (wthinput .gt. 0) then
        call climscale(time)
      endif

c ... Weather data
c ... This call removed for RAMS/Century linkage -mdh 12/96
c      call weathr(precip,prcstd,prcskw,mintmp,maxtmp)

c ... Wet-dry fixation of N

c ... Determine annual precipitation and annual PET
c ... For RAMS/Daily Century, prcann = average annual precip. 
c ... petann (output var) is computed in dailymoist. -mdh 12/96

c      do 10 mon = 1, MONTHS
c        prcann = prcann + prcurr(mon)
c        petann = petann + pevap(mon)
c10    continue

      do 10 mon = 1, MONTHS
        prcann = prcann + precip(mon) * precscalar(mon)
        agdefacm(mon) = -1.
        bgdefacm(mon) = -1.
10    continue

c ... Scale the OMAD inputs as indicated, cak - 04/06/04
      if (OMADinput .gt. 0) then
        call omadscale(time, OMADstart)
      endif

c ... Scale the N inputs as indicated, cak - 04/06/04
      if (Ninput .gt. 0) then
        call nscale(time, Nstart)
      endif

c ... N fixation in atmosphere
      baseNdep = epnfa(INTCPT)+epnfa(SLOPE)*MIN(prcann,80.0)
      if (baseNdep .lt. 0.) then
        baseNdep = 0.0
      endif

c      wdfxs = epnfs(INTCPT)+epnfs(SLOPE)*MIN(prcann,100.0)
c ... Now using annual ET in wdfxs calculation, cak - 02/21/02
c ... Use annual ET unless it is the first timestep
c ... No longer using the intercept in the calculation.
      if (annet .le. 0.0) then
        wdfxs = epnfs(SLOPE)*MIN(prcann,100.0)
      else 
        wdfxs = epnfs(2) * (annet - epnfs(1))
      endif
c ... Reset annual accumulator for evapotranspiration
      annet = 0
      if (wdfxs .lt. 0.)  then
        wdfxs = 0.0
      endif

c ... This output varible represents only the non-symbiotic soil
c ... N-fixation, the atmospheric N deposition is added to this output
c ... variable in the simsom subroutine, cak - 04/05/04
c      wdfx = wdfxa+wdfxs
      wdfx = wdfxs

c ... Atmospheric S deposition
      satmt = max(0.0, satmos(1) + satmos(2)*prcann)

c ... Determine what fraction of the carbon in new plant tissue is labeled
      if (labtyp .eq. 0) then
        cisofr = 0.0
        cisotf = 0.0
      elseif (labtyp .eq. 1) then
        cisofr = fracis(time,labyr)
        cisotf = cisofr
c      elseif (labtyp .eq. 2) then
c ..... cropin has set cisofr
c ..... treein has set cisotf
      endif

c ... Initialize co2 effects
      call co2eff(time)

c ... Implement pH shift as indicated, cak - 08/02/02
      if (phsys .gt. 0) then
        call phshift(time)
      endif

c ... Added effect of co2 for forest; done here because not calcualted
c ... dynamically based on biomass like grassland/crop
c ... Direct CO2 effects only C/E ratio of leaves.
      do 30 iel = 1, nelem
        ccefor(IMIN,LEAF,iel) = cerfor(IMIN,LEAF,iel) *
     &                          co2cce(FORSYS,IMIN,iel)
        ccefor(IMAX,LEAF,iel) = cerfor(IMAX,LEAF,iel) *
     &                          co2cce(FORSYS,IMAX,iel)
30    continue

      do 50 ipart = 2, FPARTS-1
        do 40 iel = 1, nelem 
          ccefor(IMIN,ipart,iel) = cerfor(IMIN,ipart,iel)
          ccefor(IMAX,ipart,iel) = cerfor(IMAX,ipart,iel)
40      continue 
50    continue

c ... Calculate leaf death rate multiplier for continuous forests 11/20/92
c ... Initialize LDRMLT to 1.0
      ldrmlt = 1.0

c ... Change leaf death rate multiplier if you have floating C/E ratios.
      if (ccefor(IMIN,LEAF,N) .ne. ccefor(IMAX,LEAF,N)) then
        if (rleavc .gt. 0) then
          lfncon = rleave(N) / rleavc
          lfncmin = 1 / ccefor(IMIN,LEAF,N)
          lfncmax = 1 / ccefor(IMAX,LEAF,N)
          ldrmlt = 1 + (maxldr - 1) *
     &             (lfncon - lfncmin) / (lfncmax - lfncmin)
        endif
      endif

      if (cursys .ne. FORSYS) then

c ..... Determine what fraction of plant residue added this year
c ..... will be lignin.
        call cmplig(cursys,fligni,wdlig,pltlig)
      endif

c ... a2drat - available E to plant demand for E.  -mdh 8/25/00
c ... Is this still necessary?
c ... Create separte a2drat arrays for crops and trees, mdh 5/11/01
      do 60 iel = 1, MAXIEL
        crop_a2drat(iel) = 1.0
        tree_a2drat(iel) = 1.0
60    continue

c ... Compute SITPOT as a function of the annual precipitation, cak - 05/02/03
c ... sitpot_m is the SITPOT parameter value as read from the tree.100 file
c ... for the current tree, cak - 11/21/01
c      if (prcann .lt. 30.0) then
c        sitpot = 1000.0
      if (prcann .lt. 20.0) then
        sitpot = 1500.0
c      else if (prcann .gt. 70.0) then
c        sitpot = 3000.0
c      else if (prcann .gt. 90.0) then
c        sitpot = 4000.0
c      else if (prcann .gt. 80.0) then
c        sitpot = 3500.0
      else if (prcann .gt. 90.0) then
        sitpot = 3250.0
      else
c        sitpot = line(prcann, 30.0, 1000.0, 70.0, 3000.0)
c        sitpot = line(prcann, 30.0, 1000.0, 90.0, 4000.0)
c        sitpot = line(prcann, 30.0, 1000.0, 80.0, 3500.0)
        sitpot = line(prcann, 20.0, 1500.0, 90.0, 3250.0)
      endif
      sitpot = sitpot * sitpot_m

      return
      end
