
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine annacc

      implicit none
      include 'cflows.inc'
      include 'const.inc'
      include 'monprd.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'param.inc'

c ... Reset annual accumulators.
c ... NOTE: The annet annual accumulator is reset in eachyr as it is
c ...       before used in the calculation for non-symbiotic soil N
c ...       fixation being reset

c ... Local variables
      integer iel, ii

c ... Initialize annual removal accumulators
      prcann = 0.0
      petann = 0.0
      nfixac = 0.0
      cgracc = 0.0
      snfxac(CRPSYS) = 0.0
      snfxac(FORSYS) = 0.0
      accrst = 0.0
      shrema = 0.0
      shrmai(UNLABL) = 0.0
      shrmai(LABELD) = 0.0
      sdrema = 0.0
      sdrmai(UNLABL) = 0.0
      sdrmai(LABELD) = 0.0
      creta = 0.0
      resp(UNLABL) = 0.0
      resp(LABELD) = 0.0
      cautoresp(UNLABL) = 0.0
      cautoresp(LABELD) = 0.0
      fautoresp(UNLABL) = 0.0
      fautoresp(LABELD) = 0.0
      do 10 iel = 1, nelem
        accrste(iel) = 0.0
        ereta(iel) = 0.0
        shrmae(iel) = 0.0
        sdrmae(iel) = 0.0
        egracc(iel) = 0.0
c ..... Initialize mineralization accumulators
        tnetmn(iel) = 0.0
        sumnrs(iel) = 0.0
        soilnm(iel) = 0.0
10    continue

c ... Initialize annual C production
      cproda = 0.0

c ... Initialize cinputs
      cinput = 0.0

c ... Reset minimum total non-living C, an annual value
      totc = 1000000
      
c ... Initialize co2 accumulators (10/92)
      ast1c2 = 0.0
      ast2c2 = 0.0
      amt1c2 = 0.0
      amt2c2 = 0.0
      as11c2 = 0.0
      as12c2 = 0.0
      as21c2 = 0.0
      as22c2 = 0.0
      as3c2 = 0.0
      ast1uvc2 = 0.0
      astduvc2 = 0.0

c ... Initialize the annual maintenance respiration accumulators
      mrspann(CRPSYS) = 0.0
      mrspann(FORSYS) = 0.0

c ... Initialize the annual growth respiration accumulators
      grspann(CRPSYS) = 0.0
      grspann(FORSYS) = 0.0

c ... Initialize the annual soil respiration accumulators
      srspann(CRPSYS) = 0.0
      srspann(FORSYS) = 0.0

c ... Initialize stream accumulators (04/03)
      do 20 ii = 1, 8
        strmac(ii) = 0.0
20    continue

c ... Initialize annual accumulators for N volatilization (04/14/04)
      voleac = 0.0
      volgac = 0.0
      volpac = 0.0

c ... Initialize E return from grazing accumulators (04/14/04)
      do 40 ii = 1, 3
        tgzrte(ii) = 0.0
40    continue

c ... Initialize accumulators for N deposition and non-symbiotic soil N
c ... fixation
      wdfxas = 0.0
      wdfxaa = 0.0
      wdfxa = 0.0

c ... Reset accumulators for yearly trace gas output, cak - 09/23/02
      N2O_year = 0.0
      NO_year = 0.0
      N2_year = 0.0
      CH4_oxid_year = 0.0
      nit_amt_year = 0.0
      annppt = 0.0

c ... Modify the irrtot output variable so that it is an annual
c ... accumulator rather than a a simulation long accumulator,
c ... cak - 06/25/2008
      irrtot = 0.0

c ... Modify the fertot output variable so that it is an annual
c ... accumulator for mineral fertilizer additions rather than a
c ... simulation long accumulator, add new annual accumulators for
c ... organic matter addition carbon, omadtot, and minerals,
c ... omaetot(1..3), cak - 10/20/2008
      omadtot = 0.0
      do 50 iel = 1, nelem
        fertot(iel) = 0.0
        omaetot(iel) = 0.0
50    continue

c ... Initialize annual accumulator output variables that are tracking
c ... the carbon flows due to decomposition.
      ametc1tosom11  = 0.0
      ametc2tosom12  = 0.0
      astruc1tosom11 = 0.0
      astruc1tosom21 = 0.0
      astruc2tosom12 = 0.0
      astruc2tosom22 = 0.0
      asom11tosom21  = 0.0
      asom12tosom22  = 0.0
      asom12tosom3   = 0.0
      asom21tosom11  = 0.0
      asom21tosom22  = 0.0
      asom22tosom12  = 0.0
      asom22tosom3   = 0.0
      asom3tosom12   = 0.0
      awood1tosom11  = 0.0
      awood1tosom21  = 0.0
      awood2tosom11  = 0.0
      awood2tosom21  = 0.0
      awood3tosom12  = 0.0
      awood3tosom22  = 0.0

      return
      end
