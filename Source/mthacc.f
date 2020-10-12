
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine mthacc(agdefacsum, bgdefacsum)

      implicit none
      include 'const.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'

c ... Argument declarations
      real      agdefacsum, bgdefacsum

c ... Reset monthly accumulators.

c ... Local variables
      integer iel, ii, ilyr, iso

c ... Initialize monthly accumulators
      do 20 ii = 1, 8
        stream(ii) = 0
20    continue
      pet = 0
      evap = 0
      tran = 0
      pttr = 0
      rain = 0
      agdefacsum = 0.0
      bgdefacsum = 0.0
      irract = 0.0
      runoff = 0.0
      do 102 ilyr = 1,nlayer
        amov(ilyr) = 0
102   continue

c ... Initialize monthly co2 accumlators (10/92)
      do 25 iso = 1, 2
        mt1c2(iso) = 0.0
        mt2c2(iso) = 0.0
        st1c2(iso) = 0.0
        st2c2(iso) = 0.0
        st1uvc2(iso) = 0.0
        s11c2(iso) = 0.0
        s12c2(iso) = 0.0
        s21c2(iso) = 0.0
        s22c2(iso) = 0.0
        s3c2(iso)  = 0.0
        stduvc2(iso) = 0.0
        wstduvc2(iso) = 0.0
        wd1c2(iso) = 0.0
        wd2c2(iso) = 0.0
        wd3c2(iso) = 0.0
        dwd1c2(iso) = 0.0
        dwd2c2(iso) = 0.0
25    continue

c ... Initialize monthly accumulator for volatilization of N during
c ... harvest, senescence, and return from grazing animal waste,
c ... cak 01/02
      volpl = 0.0

c ... Initialize monthly accumulator for symbiotic N fixation to track
c ... fixation for both grasses and trees as necessary, cak - 10/15/02
      nfix = 0.0

c ... Initialize monthly C production, cak - 11/20/03
      cprodc = 0.0
      cprodf = 0.0
      do 30 iel = 1, MAXIEL
        eprodc(iel) = 0.0
        eprodf(iel) = 0.0
30    continue

c ... Initialize monthly accumulator for soil surface temperature,
c ... cak - 11/20/03
      stempmth = 0.0

c ... Initialize monthly accumulators for maintenance respiration
c ... cak - 05/14/04
      mrspflux(CRPSYS) = 0.0
      mrspflux(FORSYS) = 0.0
      mrspmth(CRPSYS) = 0.0
      mrspmth(FORSYS) = 0.0
      do 40 ii = 1, CPARTS
        cmrspflux(ii) = 0.0
40    continue
      do 45 ii = 1, FPARTS
        fmrspflux(ii) = 0.0
45    continue
      sumrsp = 0.0

c ... Initialize monthly accumulators for growth respiration
c ... cak - 01/16/2007
      grspflux(CRPSYS) = 0.0
      grspflux(FORSYS) = 0.0
      grspmth(CRPSYS) = 0.0
      grspmth(FORSYS) = 0.0
      do 50 ii = 1, CPARTS
        cgrspflux(ii) = 0.0
50    continue
      do 55 ii = 1, FPARTS
        fgrspflux(ii) = 0.0
55    continue

c ... Initialize monthly respiration from decomposition output variables
c ... cak - 02/21/2007
      respmth(UNLABL) = 0.0
      respmth(LABELD) = 0.0

c ... Initialize monthly soil respiration output variables
c ... cak - 02/21/2007
      srspmth(CRPSYS) = 0.0
      srspmth(FORSYS) = 0.0

c ... Initialize monthly autotrophic respiration output variables
c ... cak - 03/27/2007
      arspmth(CRPSYS,UNLABL) = 0.0
      arspmth(CRPSYS,LABELD) = 0.0
      arspmth(FORSYS,UNLABL) = 0.0
      arspmth(FORSYS,LABELD) = 0.0

c ... Initialize monthly trace gas accumulator output variables,
c ... cak - 05/14/04
      N2O_month = 0.0
      NO_month = 0.0
      N2_month = 0.0
      ch4_oxid_month = 0.0
      nit_amt_month = 0.0
      pptmonth = 0.0

      return
      end
