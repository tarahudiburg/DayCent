
c               Copyright 1993 Colorado State University
c                       All Rights Reserved
c
c  File: plot4.inc
c...............................................................................
c... 18 new state variables, 4 aggregate pools, and 17 fluxes for standing dead 
c... tree pools. 
c... Melannie Hartman 
c... September 26, 2018
c... DEAD ATTACHED LEAVES
c...   dleavc = carbon in attached dead leaf biomass (g C m-2)
c...   dleave(MAXIEL) = N, P, S in attached dead leaf biomass (g E m-2)
c...   dlvcis(ISOS) = unlabeled and labeled carbon in attached dead leaf biomass (g C m-2)
c... DEAD ATTACHED FINE BRANCHES
c...   dfbrchc = carbon in attached dead fine branch biomass (g C m-2)
c...   dfbrche(MAXIEL) = N, P, S in attached dead fine branch biomass (g E m-2)
c...   dfbrcis(ISOS) = unlabeled and labeled carbon in attached dead fine branch biomass (g C m-2)
c... DEAD STANDING LARGE WOOD
c...   dlwodc = carbon in attached dead standing large wood biomass (g C m-2)
c...   dlwode(MAXIEL) = N, P, S in attached dead standing large wood biomass (g E m-2)
c...   dlwcis(ISOS) = unlabeled and labeled carbon in standing large wood biomass (g C m-2)
c... AGGREGATED POOLS
c...   dfrstc - total C in dead standing biomass (gC/m2)
c...   dfrste(3) - total N,P,S in dead standing biomass (gE/m2)
c.. FLUXES
c...   dw1mnr(MAXIEL) - mineralization from the decomposition of attached dead fine branches (gE/m2/mo)
c...   dw2mnr(MAXIEL) - mineralization from the decomposition of standing dead large wood (gE/m2/mo)
c...   dwd1c2(ISOS) - heterotrophic respiration from the decomposition of attached dead fine branches (gCO2-C/m2/mo)
c...   dwd2c2(ISOS) - heterotrophic respiration from the decomposition of standing dead large wood (gCO2-C/m2/mo)
c...   wstduvc2(ISOS) - CO2 released when dead attached leaves undergo photodegradation (gC/m2/mo)
c...   awstduvc2 - annual accumulator for CO2 released when dead attached leaves undergo photodegradation (gCO2-C/m2/yr)
c...   tcreta- annual accumulator of C returned to system as litter, dead surface wood, or charcoal during a 
c...       TREM (tree removal) event for a forest/savanna system (gC/m2/yr)
c...   tereta(MAXIEL)- annual accumulator of N, P, S returned to system as litter, dead surface wood, or 
c...       elemental return to mineral soil during a TREM (tree removal) event for a forest/savanna system (gE/m2/yr)


c...............................................................................
c... 111 new variables for grasstree here. Additional counterparts in other 
c ...plot*.inc files.
c... Melannie Hartman
c... March 1, 2014
c...............................................................................
c
c... LIVE LEAVES
c...   agtlvcis(ISOS) = annual accumulator for leaf production
c...   gtlvcis(ISOS)
c...   gtleavc = gtlvcis(1) + gtlvcis(2)
c...   gtleave(MAXIEL)
c...   gtlvacc = annual accumulator: agtlvcis(1) + agtlvcis(2)
c...   gtlvprd = gtlvacc in the last month of the growing season
c... DEAD LEAVES
c...   gtdlvcis(ISOS)
c...   gtdleavc = gtdlvcis(1) + gtdlvcis(2)
c...   gtdleave(MAXIEL)
c... LIVE STEMS
c...   agtstmcis(ISOS) = annual accumulator for stem production
c...   gtstmcis(ISOS) 
c...   gtstemc = gtstmcis(1) + gtstmcis(2)
c...   gtsteme(MAXIEL)
c...   gtstemacc = annual accumulator: agtstmcis(1) + agtstmcis(2)
c...   gtstemprd = gtstemacc in the last month of the growing season
c... DEAD STEMS
c...   gtdstmcis(ISOS)
c...   gtdstemc = gtdstmcis(1) + gtdstmcis(2)
c...   gtdsteme(MAXIEL)
c... LIVE JUVENILE AND MATURE FINE ROOTS
c...   agtfrcisj(ISOS) = annual accumulator for juvenile fine root production
c...   agtfrcism(ISOS) = annual accumulator for mature fine root production
c...   gtfrootcj = gtfrtcisj(1) + gtfrtcisj(2)
c...   gtfrootcm = gtfrtcism(1) + gtfrtcism(2)
c...   gtfrtcisj(ISOS)
c...   gtfrtcism(ISOS)
c...   gtfrootej(MAXIEL) 
c...   gtfrootem(MAXIEL)
c...   gtfrtjacc = annual accumulator: agtfrcisj(1) + agtfrcisj(2)
c...   gtfrtmacc = annual accumulator: agtfrcism(1) + agtfrcism(2)
c...   gtfrtjprd = gtfrtjacc in the last month of the growing season
c...   gtfrtmprd = gtfrtmacc in the last month of the growing season
c... LIVE COARSE ROOTS
c...   agtcrtcis(ISOS) = annual accumulator for coarse root production
c...   gtcrtcis(ISOS)
c...   gtcrootc = gtcrtcis(1) + gtcrtcis(2)
c...   gtcroote(MAXIEL)
c...   gtcrtacc = annual accumulator: agtcrtcis(1) + agtcrtcis(2)
c...   gtcrtprd = gtcrtacc in the last month of the growing season
c...
c...   eupgtprt(5,MAXIEL) = annual accumulator of E uptake (Note: GTLIVPARTS=5)
c...   gtstg(MAXIEL) = internal N,P,S storage
c...
c...   plot1.inc: gtgrspflux(5) = monthly growth respiration accumulator (GTLIVPARTS=5)
c...   plot1.inc: gtmrspflux(5) = monthly maintenance respiration accumulator (GTLIVPARTS=5)
c...
c...   gtpltc = gtleavc + gtfrootcj + gtfrootcm + gtstemc + gtcrootc + gtdleavc + gtdstemc
c...   gtplte(3) = gtleave(iel) + gtfrootej(iel) + gtfrootem(iel) + gtsteme(iel) + gtcroote(iel) 
c...             + gtdleave(iel) + gtdsteme(iel)
c...   gtsysc = gtpltc + litter and soil C
c...   gtsyse(MAXIEL) = gtplte(iel) + litter and soil E
c...   gtcmth(12) = monthly production for grasstree system
c...   gtcacc = gtlvacc + gtfrtjacc + gtfrtmacc + gtstemacc + gtcrtacc
c...   gtcprd = gtcacc in the last month of the growing season
c...
c...   gtcrem = C removed from harvest of GTLEAF, GTDLEAF, GTLSTEM, GTDSTEM
c...   gterem(MAXIEL) = E removed from harvest of GTLEAF, GTDLEAF, GTLSTEM, GTDSTEM
c...
c... POTENTIAL PRODUCTION
c...   potgtacc - potential growing season production for a grasstree (gC/m2/yr)

      common/plot4/agtcrtcis(2), agtfrcisj(2), agtfrcism(2), 
     & agtlvcis(2), agtstmcis(2), awstduvc2,
     & dfbrchc, dfbrche(3), dfbrcis(2),
     & dfrstc, dfrste(3), 
     & dleavc, dleave(3), dlvcis(2),
     & dlwodc, dlwode(3), dlwcis(2),
     & dw1mnr(3), dw2mnr(3), dwd1c2(2), dwd2c2(2), 
     & eupgtprt(5,3), 
     & gtcacc, gtcmth(12), gtcprd, gtcrem, gtcrootc, gtcroote(3),
     & gtcrtacc, gtcrtcis(2), gtcrtprd, 
     & gtdleavc, gtdleave(3), gtdlvcis(2),  
     & gtdstemc, gtdsteme(3), gtdstmcis(2), gterem(3),
     & gtfrootcj, gtfrootcm, gtfrootej(3), gtfrootem(3),
     & gtfrtcisj(2), gtfrtcism(2), 
     & gtfrtjacc, gtfrtjprd, gtfrtmacc, gtfrtmprd,
     & gtleavc, gtleave(3), gtlvacc, gtlvcis(2), gtlvprd, 
     & gtpltc, gtplte(3), gtstc, gtste(3), 
     & gtstemacc, gtstemc, gtsteme(3), gtstemprd, 
     & gtstg(3), gtstmcis(2), gtsysc, gtsyse(3), potgtacc,
     & tcreta, tereta(3), wstduvc2(2)

      real agtcrtcis, agtfrcisj, agtfrcism, 
     & agtlvcis, agtstmcis, awstduvc2, 
     & dfbrchc, dfbrche, dfbrcis,
     & dfrstc, dfrste, 
     & dleavc, dleave, dlvcis,
     & dlwodc, dlwode, dlwcis,
     & dw1mnr, dw2mnr, dwd1c2, dwd2c2, 
     & eupgtprt, 
     & gtcacc, gtcmth, gtcprd, gtcrem, gtcrootc, gtcroote,
     & gtcrtacc, gtcrtcis, gtcrtprd, 
     & gtdleavc, gtdleave, gtdlvcis,  
     & gtdstemc, gtdsteme, gtdstmcis, gterem,
     & gtfrootcj, gtfrootcm, gtfrootej, gtfrootem,
     & gtfrtcisj, gtfrtcism, 
     & gtfrtjacc, gtfrtjprd, gtfrtmacc, gtfrtmprd,
     & gtleavc, gtleave, gtlvacc, gtlvcis, gtlvprd, 
     & gtpltc, gtplte, gtstc, gtste, 
     & gtstemacc, gtstemc, gtsteme, gtstemprd, 
     & gtstg, gtstmcis, gtsysc, gtsyse, potgtacc,
     & tcreta, tereta, wstduvc2

      save /plot4/
