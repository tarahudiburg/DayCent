
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine harvst(month, pltlig, curday)

      implicit none
      include 'const.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      integer   month, curday
      real      pltlig(3)

c ... Harvest the crop

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow

        SUBROUTINE flowup(time)
          !MS$ATTRIBUTES ALIAS:'_flowup' :: flowup
          REAL time
        END SUBROUTINE flowup

        SUBROUTINE flowup_double(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double' :: flowup_double
          REAL time
        END SUBROUTINE flowup_double

        SUBROUTINE flowup_double_in(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double_in' :: flowup_double_in
          REAL time
        END SUBROUTINE flowup_double_in

        SUBROUTINE flowup_double_out(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double_out' :: flowup_double_out
          REAL time
        END SUBROUTINE flowup_double_out

        SUBROUTINE wrtharvest(time, curday, crpval, agcacc, bgcjacc,
     &                        bgcmacc, cgrain, egrain1, egrain2,
     &                        egrain3, crmvst, ermvst1, ermvst2,
     &                        ermvst3, cstraw, estraw1, estraw2,
     &                        estraw3, stdstraw, estdstraw1,
     &                        estdstraw2, estdstraw3, addsdc, addsde1,
     &                        addsde2, addsde3, resid, reside1,
     &                        reside2, reside3, irrapp, fertapp1,
     &                        fertapp2, fertapp3, omadapp, omaeapp1,
     &                        omaeapp2, omaeapp3, strmac1, strmac2,
     &                        strmac3, strmac4, strmac5, strmac6,
     &                        strmac7, strmac8, cgracc, egracc1,
     &                        egracc2, egracc3, accrst, accrste1,
     &                        accrste2, accrste3, ctubesj,
     &                        etubesj1, etubesj2, etubesj3, ctubesm,
     &                        etubesm1, etubesm2, etubesm3, srfclittrj,
     &                        esrfclittrj1, esrfclittrj2, esrfclittrj3,
     &                        soillittrj, esoillittrj1, esoillittrj2,
     &                        esoillittrj3, srfclittrm, esrfclittrm1,
     &                        esrfclittrm2, esrfclittrm3, soillittrm,
     &                        esoillittrm1, esoillittrm2, esoillittrm3)
          !MS$ATTRIBUTES ALIAS:'_wrtharvest' :: wrtharvest
          REAL    time
          INTEGER curday
          REAL    crpval
          REAL    agcacc
          REAL    bgcjacc
          REAL    bgcmacc
          REAL    cgrain
          REAL    egrain1
          REAL    egrain2
          REAL    egrain3
          REAL    crmvst
          REAL    ermvst1
          REAL    ermvst2
          REAL    ermvst3
          REAL    cstraw
          REAL    estraw1
          REAL    estraw2
          REAL    estraw3
          REAL    stdstraw
          REAL    estdstraw1
          REAL    estdstraw2
          REAL    estdstraw3
          REAL    addsdc
          REAL    addsde1
          REAL    addsde2
          REAL    addsde3
          REAL    resid
          REAL    reside1
          REAL    reside2
          REAL    reside3
          REAL    irrapp
          REAL    fertapp1
          REAL    fertapp2
          REAL    fertapp3
          REAL    omadapp
          REAL    omaeapp1
          REAL    omaeapp2
          REAL    omaeapp3
          REAL    strmac1
          REAL    strmac2
          REAL    strmac3
          REAL    strmac4
          REAL    strmac5
          REAL    strmac6
          REAL    strmac7
          REAL    strmac8
          REAL    cgracc
          REAL    egracc1
          REAL    egracc2
          REAL    egracc3
          REAL    accrst
          REAL    accrste1
          REAL    accrste2
          REAL    accrste3
          REAL    ctubesj
          REAL    etubesj1
          REAL    etubesj2
          REAL    etubesj3
          REAL    ctubesm
          REAL    etubesm1
          REAL    etubesm2
          REAL    etubesm3
          REAL    srfclittrj
          REAL    esrfclittrj1
          REAL    esrfclittrj2
          REAL    esrfclittrj3
          REAL    soillittrj
          REAL    esoillittrj1
          REAL    esoillittrj2
          REAL    esoillittrj3
          REAL    srfclittrm
          REAL    esrfclittrm1
          REAL    esrfclittrm2
          REAL    esrfclittrm3
          REAL    soillittrm
          REAL    esoillittrm1
          REAL    esoillittrm2
          REAL    esoillittrm3
        END SUBROUTINE wrtharvest

      END INTERFACE

c ... Local variables
      integer   iel, mm, hmonth(2)
      real      accum(ISOS), addsdc, addsde(MAXIEL), bgd,
     &          cisbgdj(ISOS), cisbgdm(ISOS), cstraw, ctubesj, ctubesm,
     &          etubes, etubesj(MAXIEL), etubesm(MAXIEL), fr14,
     &          recres(MAXIEL), recresj(MAXIEL), recresm(MAXIEL),
     &          resid, reside(MAXIEL), sumpttr, sumtran, harv_volpl
      real      stdstraw, srfclittrj, srfclittrm, soillittrj,
     &          soillittrm, esrfclittrj(MAXIEL), esoillittrj(MAXIEL),
     &          esrfclittrm(MAXIEL), esoillittrm(MAXIEL)
      real      estraw(MAXIEL), estdstraw(MAXIEL)
      real      fertapp(MAXIEL), sfertot(MAXIEL)
      real      irrapp, sirrtot
      real      omadapp, somadtot, omaeapp(MAXIEL), somaetot(MAXIEL)

      data sfertot  /3*0.0/
      data sirrtot  /0.0/
      data somadtot /0.0/
      data somaetot /3*0.0/
      save sfertot, sirrtot, somadtot, somaetot

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Initialization
      cgrain = 0.0
      crmvst = 0.0
      cstraw = 0.0
      stdstraw = 0.0
      addsdc = 0.0
      resid = 0.0
      ctubesj = 0.0
      ctubesm = 0.0
      etubes = 0.0
      srfclittrj = 0.0
      srfclittrm = 0.0
      soillittrj = 0.0
      soillittrm = 0.0
      do 5 iel = 1, MAXIEL
        egrain(iel) = 0.0
        ermvst(iel) = 0.0
        estraw(iel) = 0.0
        estdstraw(iel) = 0.0
        addsde(iel) = 0.0
        reside(iel) = 0.0
        etubesj(iel) = 0.0
        etubesm(iel) = 0.0
        esrfclittrj(iel) = 0.0
        esrfclittrm(iel) = 0.0
        esoillittrj(iel) = 0.0
        esoillittrm(iel) = 0.0
        fertapp(iel) = 0.0
        omaeapp(iel) = 0.0
5     continue

c ... Check that there is material to harvest
      if ((aglivc .le. 0.001 .and. flghrv .eq. 1) .or.
     &    (aglivc .le. 0.001 .and. stdedc .le. 0.001 .and.
     &     flghrv .eq. 0)) then
        goto 999
      endif

c ... Carbon

c ... Grain
c ... (Alister's new way of calculating:)
      if (flghrv .eq. 1) then
        if (frtcindx .lt. 5) then
          sumtran = 0
          sumpttr = 0
          hmonth(1) = month - himon(1)
          hmonth(2) = month - himon(2)
          if (hmonth(1) .lt. 1) then
            hmonth(1) = hmonth(1) + MONTHS
          endif
          if (hmonth(2) .lt. 1)  then
            hmonth(2) = hmonth(2) + MONTHS
          endif
          if (hmonth(2) .ge. hmonth(1)) then
            do 10 mm = hmonth(1), hmonth(2)
              sumtran = sumtran + htran(mm)
              sumpttr = sumpttr + hpttr(mm)
10          continue
          else
            do 15 mm = hmonth(1), MONTHS
              sumtran = sumtran + htran(mm)
              sumpttr = sumpttr + hpttr(mm)
15          continue
            do 16 mm = 1, hmonth(2)
              sumtran = sumtran + htran(mm)
              sumpttr = sumpttr + hpttr(mm)
16          continue
          endif
C         if (sumpttr .eq. 0.0) then
          if (sumpttr .le. 0.0001) then
            sumpttr = 0.0001
          endif
          hi = himax * (1.0 - hiwsf * (1.0 - (sumtran / sumpttr)))
        else if (frtcindx .ge. 5) then
          if (grnfldys .gt. 0) then
            hi = himax * (1.0 - hiwsf * (1.0 - (gwstress/grnfldys)))
          else
            hi = 0.0
          endif
        endif
        cgrain = hi * aglivc * (1.0 - aglrem)
        call csched(cgrain,aglcis(LABELD),aglivc,
     &              aglcis(UNLABL),csrsnk(UNLABL),
     &              aglcis(LABELD),csrsnk(LABELD),
     &              1.0,cisgra)
c ..... Straw
        cstraw = aglivc * (1.0 - aglrem) - cgrain
c ..... Straw removal
        crmvst = rmvstr * cstraw
        accrst = accrst + crmvst
        call csched(crmvst,aglcis(LABELD),aglivc,
     &              aglcis(UNLABL),csrsnk(UNLABL),
     &              aglcis(LABELD),csrsnk(LABELD),
     &              1.0,accum)
c ..... Some straw will remain as standing dead
        addsdc = remwsd * (cstraw-crmvst)
        call csched(addsdc,aglcis(LABELD),aglivc,
     &              aglcis(UNLABL),stdcis(UNLABL),
     &              aglcis(LABELD),stdcis(LABELD),
     &              1.0,accum)
c ... Non-grain harvest
      else
        cgrain = 0.0
c ..... Straw can come from aboveground live carbon and/or standing
c ..... dead carbon
        if (aglivc .gt. 0.001) then
          cstraw = aglivc * (1.0 - aglrem)
          call csched(cstraw,aglcis(LABELD),aglivc,
     &                aglcis(UNLABL),csrsnk(UNLABL),
     &                aglcis(LABELD),csrsnk(LABELD),
     &                1.0,accum)
        else
          cstraw = 0.0
        endif
        if (stdedc .gt. 0.001) then
          stdstraw = stdedc * (1.0 - aglrem)
          call csched(stdstraw,stdcis(LABELD),stdedc,
     &                stdcis(UNLABL),csrsnk(UNLABL),
     &                stdcis(LABELD),csrsnk(LABELD),
     &                1.0,accum)
        else
          stdstraw = 0.0
        endif
        crmvst = cstraw + stdstraw
        accrst = accrst + crmvst
      endif

c ... Other elements
      do 20 iel = 1, nelem

c ..... Grain
        if (flghrv .eq. 1) then
          egrain(iel) = efrgrn(iel) * aglive(iel) * (1.0 - aglrem) *
     &                  sqrt(hi/himax)
          call flow(aglive(iel),esrsnk(iel),time,egrain(iel))
c ....... Volatilization of N from plants
c ....... Use the local variable harv_volpl so that volatilization that
c ....... occurs at harvest and senescence can both be tracked, see dshoot
c ....... cak - 01/02
          if (iel .eq. N) then
            harv_volpl = vlossp * aglive(iel)
            call flow(aglive(iel),esrsnk(iel),time,harv_volpl)
            volpl = volpl + harv_volpl
            volpla = volpla + harv_volpl
            volpac = volpac + harv_volpl
c ......... N/C ratio in straw
            recres(iel) = ((aglive(iel) - harv_volpl) * (1.0 - aglrem) -
     &                      egrain(iel)) / cstraw
          else
c ......... P/C, or S/C ratio in straw
            recres(iel) = (aglive(iel) * (1.0 - aglrem) -
     &                     egrain(iel)) / cstraw
          endif
c ....... Straw removal
          ermvst(iel) = crmvst * recres(iel)
          call flow(aglive(iel),esrsnk(iel),time,ermvst(iel))
c ....... Some straw remains as standing dead
          addsde(iel) = addsdc * recres(iel)
          call flow(aglive(iel),stdede(iel),time,addsde(iel))
        else
          egrain(iel) = 0.0
c ....... E/C ratio in aboveground live removed as straw/hay
          if (aglivc .gt. 0.001) then
            recres(iel) = aglive(iel) / aglivc
c ......... Straw/hay removal from aboveground live
            estraw(iel) = cstraw * recres(iel)
            call flow(aglive(iel),esrsnk(iel),time,estraw(iel))
          endif
c ....... E/C ratio in standing dead removed as straw/hay
          if (stdedc .gt. 0.001) then
            recres(iel) = stdede(iel) / stdedc
c ......... Straw/hay removal from standing dead
            estdstraw(iel) = stdstraw * recres(iel)
            call flow(stdede(iel),esrsnk(iel),time,estdstraw(iel))
          endif
          ermvst(iel) = estraw(iel) + estdstraw(iel)
          accrste(iel) = accrste(iel) + ermvst(iel)
        endif
20    continue

      if (flghrv .eq. 1) then
c ..... Partition c, n, p, and s in remaining straw into top layer
c ..... of structural and metabolic
        resid = cstraw - crmvst - addsdc
        if (aglivc .gt. 0.0) then
          fr14 = aglcis(LABELD)/aglivc
        else
          fr14 = cisofr
        endif
        call partit(resid,recres,1,aglcis,aglive,pltlig(ABOVE),fr14)
        do 25 iel = 1, nelem
          reside(iel) = resid * recres(iel)
25      continue
      endif

c ... Below ground removal (root harvest) -lh 8/91
c ... Harvest both juvenile and mature roots, cak - 05/18/2007

c ... Juvenile root harvest
      if (bglivcj .gt. 0.001) then
        ctubesj = hibg * bglivcj * (1.0 - bglrem)
        cgrain = cgrain + ctubesj
        call csched(ctubesj,bglcisj(LABELD),bglivcj,
     &              bglcisj(UNLABL), csrsnk(UNLABL),
     &              bglcisj(LABELD), csrsnk(LABELD),
     &              1.0,cisgra)
        if (bglivcj .gt. 0.0001) then
          cisbgdj(LABELD) = (bglivcj * (1.0 - bglrem) - ctubesj) *
     &                      (bglcisj(LABELD) / bglivcj)
        else
          cisbgdj(LABELD) = 0.0
        endif
        cisbgdj(UNLABL) = (bglivcj * (1.0 - bglrem) - ctubesj) -
     &                    cisbgdj(LABELD)
        do 30 iel = 1, nelem
          etubesj(iel) = hibg * (1.0 - bglrem) * bglivej(iel)
          call flow(bglivej(iel), esrsnk(iel), time, etubesj(iel))
          egrain(iel) = egrain(iel) + etubesj(iel)
          recresj(iel) = bglivej(iel) / bglivcj
30      continue
      else
        cisbgdj(LABELD) = 0.0
        cisbgdj(UNLABL) = 0.0
        recresj(N) = 0.0
        recresj(P) = 0.0
        recresj(S) = 0.0
      endif

c ... Mature root harvest
      if (bglivcm .gt. 0.001) then
        ctubesm = hibg * bglivcm * (1.0 - bglrem)
        cgrain = cgrain + ctubesm
        call csched(ctubesm,bglcism(LABELD),bglivcm,
     &              bglcism(UNLABL), csrsnk(UNLABL),
     &              bglcism(LABELD), csrsnk(LABELD),
     &              1.0,cisgra)
        if (bglivcm .gt. 0.0001) then
          cisbgdm(LABELD) = (bglivcm * (1.0 - bglrem) - ctubesm) *
     &                      (bglcism(LABELD) / bglivcm)
        else
          cisbgdm(LABELD) = 0.0
        endif
        cisbgdm(UNLABL) = (bglivcm * (1.0 - bglrem) - ctubesm) -
     &                    cisbgdm(LABELD)
        do 35 iel = 1, nelem
          etubesm(iel) = hibg * (1.0 - bglrem) * bglivem(iel)
          call flow(bglivem(iel), esrsnk(iel), time, etubesm(iel))
          egrain(iel) = egrain(iel) + etubesm(iel)
          recresm(iel) = bglivem(iel) / bglivcm
35      continue
      else
        cisbgdm(LABELD) = 0.0
        cisbgdm(UNLABL) = 0.0
        recresm(N) = 0.0
        recresm(P) = 0.0
        recresm(S) = 0.0
      endif

c ... When harvesting roots remove E from crop storage as well
      if (bglivcj+bglivcm .gt. 0.001) then
        do 37 iel = 1, nelem
          etubes = hibg * (1.0 - bglrem) * crpstg(iel)
          call flow(crpstg(iel), esrsnk(iel), time, etubes)
37      continue
      endif

c ... Calculation of accumulator for grain production
      cgracc = cgracc + cgrain
      do 40 iel = 1, nelem
        egracc(iel) = egracc(iel)+egrain(iel)
40    continue

c ... Partition c, n, p, and s in remaining roots into bottom layer of
c ... structural and metabolic
c ... Juvenile roots
      bgd = cisbgdj(LABELD) + cisbgdj(UNLABL)
      if (bglivcj .gt. 0.0) then
        fr14 = bglcisj(LABELD)/bglivcj
      else
        fr14 = cisofr
      endif
c ... A fraction of the dead roots are transferred to the surface
c ... litter layer, the remainder goes to the soil litter layer
c ... cak - 05/14/2007
      srfclittrj = bgd * rdsrfc
      soillittrj = bgd - srfclittrj
      call partit(srfclittrj, recresj, SRFC, bglcisj, bglivej,
     &           pltlig(BELOWJ), fr14)
      call partit(soillittrj, recresj, SOIL, bglcisj, bglivej,
     &           pltlig(BELOWJ), fr14)
      do 90 iel = 1, nelem
          esrfclittrj(iel) = srfclittrj * recresj(iel)
          esoillittrj(iel) = soillittrj * recresj(iel)
90      continue
c ... Mature roots
      bgd = cisbgdm(LABELD) + cisbgdm(UNLABL)
      if (bglivcm .gt. 0.0) then
        fr14 = bglcism(LABELD)/bglivcm
      else
        fr14 = cisofr
      endif
c ... A fraction of the dead roots are transferred to the surface
c ... litter layer, the remainder goes to the soil litter layer
c ... cak - 05/14/2007
      srfclittrm = bgd * rdsrfc
      soillittrm = bgd - srfclittrm
      call partit(srfclittrm, recresm, SRFC, bglcism, bglivem,
     &           pltlig(BELOWM), fr14)
      call partit(soillittrm, recresm, SOIL, bglcism, bglivem,
     &           pltlig(BELOWM), fr14)
      do 100 iel = 1, nelem
          esrfclittrm(iel) = srfclittrm * recresm(iel)
          esoillittrm(iel) = soillittrm * recresm(iel)
100   continue

c ... Write output to the harvest.csv file
      if (flghrv .eq. 1) then
        cstraw = crmvst
        estraw(N) = ermvst(N)
        estraw(P) = ermvst(P)
        estraw(S) = ermvst(S)
      endif
      do 105 iel = 1, nelem
        fertapp(iel) = gfertot(iel) - sfertot(iel)
        omaeapp(iel) = gomaetot(iel) - somaetot(iel)
        sfertot(iel) = gfertot(iel)
        somaetot(iel) = gomaetot(iel)
105   continue
      irrapp = girrtot - sirrtot
      omadapp = gomadtot - somadtot
      sirrtot = girrtot
      somadtot = gomadtot
      call wrtharvest(time, curday, crpval, agcacc, bgcjacc, bgcmacc,
     &                cgrain, egrain(N), egrain(P), egrain(S), crmvst,
     &                ermvst(N), ermvst(P), ermvst(S), cstraw,
     &                estraw(N), estraw(P), estraw(S), stdstraw,
     &                estdstraw(N), estdstraw(P), estdstraw(S), addsdc,
     &                addsde(N), addsde(P), addsde(S), resid,
     &                reside(N), reside(P), reside(S), irrapp,
     &                fertapp(N), fertapp(P), fertapp(S), omadapp,
     &                omaeapp(N), omaeapp(P), omaeapp(S), strmac(1),
     &                strmac(2), strmac(3), strmac(4), strmac(5),
     &                strmac(6), strmac(7), strmac(8), cgracc,
     &                egracc(N), egracc(P), egracc(S), accrst,
     &                accrste(N), accrste(P), accrste(S), ctubesj,
     &                etubesj(N), etubesj(P), etubesj(S), ctubesm,
     &                etubesm(N), etubesm(P), etubesm(S),
     &                srfclittrj, esrfclittrj(N), esrfclittrj(P),
     &                esrfclittrj(S), soillittrj, esoillittrj(N),
     &                esoillittrj(P), esoillittrj(S), srfclittrm,
     &                esrfclittrm(N), esrfclittrm(P), esrfclittrm(S),
     &                soillittrm, esoillittrm(N), esoillittrm(P),
     &                esoillittrm(S))

c ... Update state variables and accumulators.
      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

c ... Check status of values to make sure that everything
c ... has been reset correctly

      if (aglcis(UNLABL)+aglcis(LABELD) .lt. 1.e-05) then
        aglcis(UNLABL) = 0.0
        aglcis(LABELD) = 0.0
      endif
C     if (aglcis(UNLABL)+aglcis(LABELD) .eq. 0) then
      if (aglcis(UNLABL)+aglcis(LABELD) .le. 0.0) then
        do 50 iel = 1, MAXIEL
          aglive(iel) = 0.0
50      continue
      endif
      
      if (bglcisj(UNLABL)+bglcisj(LABELD) .lt. 1.e-05) then
        bglcisj(UNLABL) = 0.0
        bglcisj(LABELD) = 0.0
      endif  
C     if (bglcisj(UNLABL)+bglcisj(LABELD) .eq. 0) then
      if (bglcisj(UNLABL)+bglcisj(LABELD) .le. 0.0) then
        do 60 iel = 1, MAXIEL
          bglivej(iel) = 0.0
60      continue
      endif

      if (bglcism(UNLABL)+bglcism(LABELD) .lt. 1.e-05) then
        bglcism(UNLABL) = 0.0
        bglcism(LABELD) = 0.0
      endif  
C     if (bglcism(UNLABL)+bglcism(LABELD) .eq. 0) then
      if (bglcism(UNLABL)+bglcism(LABELD) .le. 0.0) then
        do 65 iel = 1, MAXIEL
          bglivem(iel) = 0.0
65      continue
      endif

      if (stdcis(UNLABL)+stdcis(LABELD) .lt. 1.e-05)then
        stdcis(UNLABL) = 0.0
        stdcis(LABELD) = 0.0
      endif  
C     if (stdcis(UNLABL)+stdcis(LABELD) .eq. 0) then
      if (stdcis(UNLABL)+stdcis(LABELD) .le. 0.0) then
        do 70 iel = 1, MAXIEL
          stdede(iel) = 0
70      continue
      endif
          
      do 80 iel = 1, MAXIEL
        if (aglive(iel) .lt. 1.e-05) then
          aglive(iel) = 0.0
        endif
        if (bglivej(iel) .lt. 1.e-05) then
          bglivej(iel) = 0.0
        endif
        if (bglivem(iel) .lt. 1.e-05) then
          bglivem(iel) = 0.0
        endif
        if (stdede(iel) .lt. 1.e-05) then
          stdede(iel) = 0.0
        endif
80    continue

999   continue

      return
      end
