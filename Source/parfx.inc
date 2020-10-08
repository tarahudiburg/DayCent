
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... sdco2sum - soil decomposition CO2 flow sum. Needed to calculate the
c                effective decomposition effect.
c ... ntdco2sm - no-till soil decomposition CO2 flow sum
c                KLK 14Aug2012

      common/parfx/psecoc2,
     &    adep(10), agppa, agppb, aneref(3), animpt, 
     &    awtl(10), bgppa, bgppb, co2ppm(2), co2rmp,
     &    damr(2,3), damrmn(3),
     &    dec1(2), dec2(2), dec3(2), dec4, dec5(2), deck5, dligdf,
     &    dresp, sdco2sum, ntdco2sm, edepth, elitst, enrich,
     &    favail(6), fleach(5), fwloss(4),
     &    fxmca, fxmcb, fxmxs, fxnpb, gremb, lhzf(3),
     &    minlch,
     &    omlech(3), 
     &    p1co2a(2), p1co2b(2), p2co2(2), p3co2, pabres,
     &    peftxa, peftxb, phesp(4),
     &    pligst(2), pmco2(2), pmnsec(3),
     &    pmntmp, pmxbio, pmxtmp, pparmn(3), pprpts(3),
     &    ps1co2(2), ps1s3(2), ps2s3(2),
     &    psecmn(3), psecoc1, rad1p(3,3), rcestr(3),
     &    rictrl, riint, rsplig,
     &    sfavail(2), spl(2), strmax(2),
     &    texepp(5), texesp(3), teff(4), tmelt(2),
     &    varat11(3,3), varat12(3,3), varat21(3,3), varat22(3,3),
     &    varat3(3,3), vlosse, vlossg, vlossg_m,
     &    idef, nsnfix, ntspm, seed

      double precision psecoc2

      real adep, agppa, agppb, aneref, animpt,
     &     awtl, bgppa, bgppb, co2ppm, co2rmp,
     &     damr, damrmn,
     &     dec1, dec2, dec3, dec4, dec5, deck5, dligdf,
     &     dresp, sdco2sum, ntdco2sm, edepth, elitst, enrich,
     &     favail, fleach, fwloss,
     &     fxmca, fxmcb, fxmxs, fxnpb, gremb, lhzf,
     &     minlch, omlech,
     &     p1co2a, p1co2b, p2co2, p3co2, pabres,
     &     peftxa, peftxb, phesp,
     &     pligst, pmco2, pmnsec,
     &     pmntmp, pmxbio, pmxtmp, pparmn, pprpts,
     &     ps1co2, ps1s3, ps2s3,
     &     psecmn, psecoc1, rad1p, rcestr,
     &     rictrl, riint, rsplig,
     &     sfavail, spl, strmax,
     &     texepp, texesp, teff, tmelt,
     &     varat11, varat12, varat21, varat22,
     &     varat3, vlosse, vlossg, vlossg_m

      integer idef, nsnfix, ntspm, seed

      save /parfx/
