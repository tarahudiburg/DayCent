
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... The variables in /plot2/ are written to the plot file in the UNIX
c ... version of the model, but not the Microsoft version.

      common /plot2/ aagdefac, abgdefac, accrst, accrste(3), agcisa(2),
     &    aglcn, bgcisja(2), bgcisma(2), bglcnj, bglcnm,
     &    carbostg(2,2), cgracc, cisgra(2), cltfac(4), csrsnk(2),
     &    dautoresp(2), dbglivc, dbglivcj, dbglivcm, dblit,
     &    dcarbostg(2), deloe, deloi, dfrootc, dfrootcj, dfrootcm,
     &    dmetc(2), dhetresp, dslit, dsoilresp, dsom1c(2), dsom2c(2),
     &    dsom3c, dsomtc, dstruc(2),
     &    egracc(3), ereta(3),
     &    esrsnk(3),
     &    gromin(3),
     &    lhzcac, lhzeac(3),
     &    metmnr(2,3),
     &    n2oacc, n2omth(12), n2oprd,
     &    omadac, omadae(3), omadmte(12,3), omadmth(12),
     &    omadprd, omadpre(3), omadtot, omaetot(3), prcfal, rnpml1,
     &    sdrmae(3), sdrmai(2), shrmai(2), 
     &    shrmae(3), somsci(2), somtci(2), somte(3), 
     &    strmnr(2,3), s1mnr(2,3), s2mnr(2,3), s3mnr(3),
     &    tcerat(3), tcnpro, tgzrte(3), tomres(2), totalc, totale(3),
     &    totsysc, totsyse(3),
     &    voleac, volex, volgac, volgm, volpac, volpl,
     &    wdfx, wdfxa, wdfxma, wdfxms, wdfxs

      real aagdefac, abgdefac, accrst, accrste, agcisa,
     &     aglcn, bgcisja, bgcisma, bglcnj, bglcnm,
     &     carbostg, cgracc, cisgra, cltfac, csrsnk,
     &     dautoresp, dbglivc, dbglivcj, dbglivcm, dblit,
     &     dcarbostg, deloe, deloi, dfrootc, dfrootcj, dfrootcm,
     &     dmetc, dhetresp, dslit, dsoilresp, dsom1c, dsom2c,
     &     dsom3c, dsomtc, dstruc,
     &     egracc, ereta,
     &     esrsnk,
     &     gromin,
     &     lhzcac, lhzeac,
     &     metmnr,
     &     n2oacc, n2omth, n2oprd,
     &     omadac, omadae, omadmte, omadmth,
     &     omadprd, omadpre, omadtot, omaetot, prcfal, rnpml1,
     &     sdrmae, sdrmai, shrmai,
     &     shrmae, somsci, somtci, somte,
     &     strmnr, s1mnr, s2mnr, s3mnr,
     &     tcerat, tcnpro, tgzrte, tomres, totalc, totale,
     &     totsysc, totsyse,
     &     voleac, volex, volgac, volgm, volpac, volpl,
     &     wdfx, wdfxa, wdfxma, wdfxms, wdfxs

      save /plot2/
