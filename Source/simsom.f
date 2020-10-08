
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine simsom()

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'jday.inc'
      include 'ligvar.inc'
      include 'monprd.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

c ... Simulate flow of carbon, nitrogen, phosphorous, and sulfur.
c ... This routine is executed each time step.  It calls the decomposition
c ... submodel and a producer submodel.  It also includes a bunch of
c ... N fixation stuff that needs to be rewritten and put in its own routine.
c ...
c ... Added new local variable FSOL.  Added calls to new function FSFUNC
c ... to calculate the amount of mineral P that is in solution.  Added
c ... call to new subroutine PSCHEM, which calculates and schedules the
c ... Phosophorus and Sulfur flows during decomposition.  Previously
c ... this was calculated in the DECOMP routine.  -rm 6/91
c ...
c ... Removed decomposition submodel and put in call to dailymoist for
c ... the daily water budget version of Century.  -mdh 8/94

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

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &                          ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER          clyr
          REAL             amt
          DOUBLE PRECISION frac_nh4
          DOUBLE PRECISION frac_no3
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          CHARACTER        subname*10
        END SUBROUTINE update_npool

        SUBROUTINE wrtbio(time, curday, aglivc, bglivcj, bglivcm,
     &                    aglive, bglivej, bglivem, rleavc, frootcj,
     &                    frootcm, fbrchc, rlwodc, crootc, h2ogef1,
     &                    h2ogef2)
          !MS$ATTRIBUTES ALIAS:'_wrtbio' :: wrtbio
          REAL    time
          INTEGER curday
          REAL    aglivc
          REAL    bglivcj
          REAL    bglivcm
          REAL    aglive
          REAL    bglivej
          REAL    bglivem
          REAL    rleavc
          REAL    frootcj
          REAL    frootcm
          REAL    fbrchc
          REAL    rlwodc
          REAL    crootc
          REAL    h2ogef1
          REAL    h2ogef2
        END SUBROUTINE wrtbio

        SUBROUTINE wrtdcsip(time, curday, trandly, evapdly, intrcpt,
     &                      sublim, outflow, runoffdly, ppt, saccum,
     &                      melt, snow, snlq, petdly, stemp, wc_2cm,
     &                      wc_3cm, wc_5cm, wc_10cm, wc_15cm, wc_30cm,
     &                      CO2resp, mcprd1, mcprd2, mcprd3, mfprd1,
     &                      mfprd2, mfprd6, mfprd3, mfprd4, mfprd5,
     &                      npptot, nee, aglivc, bglivcj, bglivcm,
     &                      rleavc, frootcj, frootcm, fbrchc, rlwodc,
     &                      crootc, tlai, stdedc, wood1c, wood2c,
     &                      wood3c, strucc1, metabc1, strucc2, metabc2,
     &                      som1c1, som1c2, som2c1, som2c2, som3c,
     &                      systemc)
          !MS$ATTRIBUTES ALIAS:'_wrtdcsip' :: wrtdcsip
          REAL    time
          INTEGER curday
          REAL    trandly
          REAL    evapdly
          REAL    intrcpt
          REAL    sublim
          REAL    outflow
          REAL    runoffdly
          REAL    ppt
          REAL    saccum
          REAL    melt
          REAL    snow
          REAL    snlq
          REAL    petdly
          REAL    stemp
          REAL    wc_2cm
          REAL    wc_3cm
          REAL    wc_5cm
          REAL    wc_10cm
          REAL    wc_15cm
          REAL    wc_30cm
          REAL    CO2resp
          REAL    mcprd1
          REAL    mcprd2
          REAL    mcprd3
          REAL    mfprd1
          REAL    mfprd2
          REAL    mfprd6
          REAL    mfprd3
          REAL    mfprd4
          REAL    mfprd5
          REAL    npptot
          REAL    nee
          REAL    aglivc
          REAL    bglivcj
          REAL    bglivcm
          REAL    rleavc
          REAL    frootcj
          REAL    frootcm
          REAL    fbrchc
          REAL    rlwodc
          REAL    crootc
          REAL    tlai
          REAL    stdedc
          REAL    wood1c
          REAL    wood2c
          REAL    wood3c
          REAL    strucc1
          REAL    metabc1
          REAL    strucc2
          REAL    metabc2
          REAL    som1c1
          REAL    som1c2
          REAL    som2c1
          REAL    som2c2
          REAL    som3c
          REAL    systemc
        END SUBROUTINE wrtdcsip

        SUBROUTINE wrtdeadc(time, curday, stdedc, metabc1, strucc1,
     &                      wood1c, wood2c, wood3c)
          !MS$ATTRIBUTES ALIAS:'_wrtdeadc' :: wrtdeadc
          REAL    time
          INTEGER curday
          REAL    stdedc
          REAL    metabc1
          REAL    strucc1
          REAL    wood1c
          REAL    wood2c
          REAL    wood3c
        END SUBROUTINE wrtdeadc

        SUBROUTINE wrtdels(time, curday, ddeloi, ddeloe, ddsrfclit,
     &                     ddsmnrl, ddhetresp, ddsoilresp, ddcmresp,
     &                     ddfmresp, ddcgresp, ddfgresp, ddccarbostg,
     &                     ddfcarbostg)
          !MS$ATTRIBUTES ALIAS:'_wrtdels' :: wrtdels
          REAL    time
          INTEGER curday
          REAL    ddeloi
          REAL    ddeloe
          REAL    ddsrfclit
          REAL    ddsmnrl
          REAL    ddhetresp
          REAL    ddsoilresp
          REAL    ddcmresp
          REAL    ddfmresp
          REAL    ddcgresp
          REAL    ddfgresp
          REAL    ddccarbostg
          REAL    ddfcarbostg
        END SUBROUTINE wrtdels

        SUBROUTINE wrtlivec(time, curday, aglivc, bglivcj, bglivcm,
     &                      rleavc, frootcj, frootcm, fbrchc, rlwodc,
     &                      crootc)
          !MS$ATTRIBUTES ALIAS:'_wrtlivec' :: wrtlivec
          REAL    time
          INTEGER curday
          REAL    aglivc
          REAL    bglivcj
          REAL    bglivcm
          REAL    rleavc
          REAL    frootcj
          REAL    frootcm
          REAL    fbrchc
          REAL    rlwodc
          REAL    crootc
        END SUBROUTINE wrtlivec

        SUBROUTINE wrtmethane(time, curday, aglivc, bglivcj, bglivcm,
     &                        prev_mcprd1, prev_mcprd2, prev_mcprd3,
     &                        Com, ppt, irractwk, watr2sat, avgst_10cm,
     &                        TI, SI, Cr, Eh, Feh, CH4_prod, CH4_Ep,
     &                        CH4_Ebl, CH4_oxid)
          !MS$ATTRIBUTES ALIAS:'_wrtmethane' :: wrtmethane
          REAL             time
          INTEGER          curday
          REAL             aglivc
          REAL             bglivcj
          REAL             bglivcm
          REAL             prev_mcprd1
          REAL             prev_mcprd2
          REAL             prev_mcprd3
          REAL             Com
          REAL             ppt
          REAL             irractwk
          REAL             watr2sat
          REAL             avgst_10cm
          REAL             TI
          REAL             SI
          REAL             Cr
          REAL             Eh
          REAL             Feh
          REAL             CH4_prod
          REAL             CH4_Ep
          REAL             CH4_Ebl
          DOUBLE PRECISION CH4_oxid
        END SUBROUTINE wrtmethane

        SUBROUTINE wrtresp(time, curday, doiresp, doeresp, dslitrsp,
     &                     dmnrlrsp, dhresp, dcrtjresp, dcrtmresp,
     &                     dfrtjresp, dfrtmresp, dfrtcresp, dsresp,
     &                     dmresp, dgresp, mrspdyflux1, mrspdyflux2,
     &                     cmrspdyflux1, cmrspdyflux2, cmrspdyflux3,
     &                     fmrspdyflux1, fmrspdyflux2, fmrspdyflux6,
     &                     fmrspdyflux3, fmrspdyflux4, fmrspdyflux5,
     &                     mrspann1, mrspann2, tavedly,
     &                     mrspTempEffect11, mrspTempEffect12,
     &                     mrspWaterEffect1, mrspTempEffect21,
     &                     mrspTempEffect22, mrspWaterEffect2,
     &                     grspdyflux1, grspdyflux2,
     &                     cgrspdyflux1, cgrspdyflux2, cgrspdyflux3,
     &                     fgrspdyflux1, fgrspdyflux2, fgrspdyflux6,
     &                     fgrspdyflux3, fgrspdyflux4, fgrspdyflux5,
     &                     grspann1, grspann2, mcprd1, mcprd2, mcprd3,
     &                     mfprd1, mfprd2, mfprd6, mfprd3, mfprd4,
     &                     mfprd5, carbostg11, carbostg12, carbostg21,
     &                     carbostg22)
          !MS$ATTRIBUTES ALIAS:'_wrtresp' :: wrtresp
          REAL    time
          INTEGER curday
          REAL    doiresp
          REAL    doeresp
          REAL    dslitrsp
          REAL    dmnrlrsp
          REAL    dhresp
          REAL    dcrtjresp
          REAL    dcrtmresp
          REAL    dfrtjresp
          REAL    dfrtmresp
          REAL    dfrtcresp
          REAL    dsresp
          REAL    dmresp
          REAL    dgresp
          REAL    mrspdyflux1
          REAL    mrspdyflux2
          REAL    cmrspdyflux1
          REAL    cmrspdyflux2
          REAL    cmrspdyflux3
          REAL    fmrspdyflux1
          REAL    fmrspdyflux2
          REAL    fmrspdyflux6
          REAL    fmrspdyflux3
          REAL    fmrspdyflux4
          REAL    fmrspdyflux5
          REAL    mrspann1
          REAL    mrspann2
          REAL    tavedly
          REAL    mrspTempEffect11
          REAL    mrspTempEffect12
          REAL    mrspWaterEffect1
          REAL    mrspTempEffect21
          REAL    mrspTempEffect22
          REAL    mrspWaterEffect2
          REAL    grspdyflux1
          REAL    grspdyflux2
          REAL    cgrspdyflux1
          REAL    cgrspdyflux2
          REAL    cgrspdyflux3
          REAL    fgrspdyflux1
          REAL    fgrspdyflux2
          REAL    fgrspdyflux6
          REAL    fgrspdyflux3
          REAL    fgrspdyflux4
          REAL    fgrspdyflux5
          REAL    grspann1
          REAL    grspann2
          REAL    mcprd1
          REAL    mcprd2
          REAL    mcprd3
          REAL    mfprd1
          REAL    mfprd2
          REAL    mfprd6
          REAL    mfprd3
          REAL    mfprd4
          REAL    mfprd5
          REAL    carbostg11
          REAL    carbostg12
          REAL    carbostg21
          REAL    carbostg22
        END SUBROUTINE wrtresp

        SUBROUTINE wrtsoilc(time, curday, metabc2, strucc2, som1c1,
     &                      som1c2, som2c1, som2c2, som3c)
          !MS$ATTRIBUTES ALIAS:'_wrtsoilc' :: wrtsoilc
          REAL    time
          INTEGER curday
          REAL    metabc2
          REAL    strucc2
          REAL    som1c1
          REAL    som1c2
          REAL    som2c1
          REAL    som2c2
          REAL    som3c
        END SUBROUTINE wrtsoilc

        SUBROUTINE wrtsysc(time, curday, livec, deadc, soilc, sysc,
     &                     CO2resp)
          !MS$ATTRIBUTES ALIAS:'_wrtsysc' :: wrtsysc
          REAL    time
          INTEGER curday
          REAL    livec
          REAL    deadc
          REAL    soilc
          REAL    sysc
          REAL    CO2resp
        END SUBROUTINE wrtsysc

      END INTERFACE

c ... Function declarations
      real      fsfunc, del13out, del14out, irrigt, ramp
      external  fsfunc, del13out, del14out, irrigt, ramp

c ... Local variables
      integer   iel, lyr, imo, ii, clyr, ipart, iso
      real      basf, biof, cancvr, cmn
      real      frlech(MAXIEL), fsol, fwdfx, fxbiom
      real      satm, sirr, stot, tbiom, texeff
      real      wdbmas, wdfnp, wdfxm
      integer   curday
      real      irractwk, agdefacsum, bgdefacsum, tfrac
      real      bgwfunc
      real      prevstream(8)
      real      tavedly, tavemth, tavewk
      double precision frac_nh4, frac_no3
      real      co2val, CO2resp
      character subname*10
      integer   plntd
      real      respsum(2), arspsum(2,2)
      real      saccum, evapdly, intrcpt, melt, outflow, petdly,
     &          runoffdly, sublim, trandly, aetdly
      real      ddeloi, ddeloe, ddsrfclit, ddsmnrl, ddhetresp,
     &          ddsoilresp
      real      oiresp1(2), oiresp2(2), newoiresp(2)
      real      oeresp1(2), oeresp2(2), newoeresp(2)
      real      slitrsp1(2), slitrsp2(2), newslitrsp(2)
      real      mnrlrsp1(2), mnrlrsp2(2), newmnrlrsp(2)
      real      hetresp1(2), hetresp2(2), newhetresp(2)
      real      ddccarbostg, ddfcarbostg
      real      doiresp, doeresp, dslitrsp, dmnrlrsp, dhresp,
     &          dcrtresp, dcrtjresp, dcrtmresp,
     &          dfrtresp, dfrtjresp, dfrtmresp, dfrtcresp
      real      dsresp, dmresp, dgresp, crplblstg, forlblstg
      real      dcmresp, dfmresp, ddcmresp, ddfmresp
      real      dcgresp, dfgresp, ddcgresp, ddfgresp
      real      npptot, systemc, wc_2cm, wc_3cm, wc_5cm,
     &          wc_10cm, wc_15cm, wc_30cm, tlai, nee
      real      tmns_mlt, tmxs_mlt
      real      avgstemp, watr2sat
      real      prev_mcprd1, prev_mcprd2, prev_mcprd3
      real      accum(ISOS)
      real      scenfrac
      real      crpeff, trpeff, rpeff
c      real      dmt1c2, dmt2c2, dst1c2,
c     &          dst2c2, ds11c2, ds12c2,
c     &          ds21c2, ds22c2, ds3c2,
c     &          dwd1c2, dwd2c2, dwd3c2
c      real      mt1c2_1(ISOS), mt2c2_1(ISOS), st1c2_1(ISOS),
c     &          st2c2_1(ISOS), s11c2_1(ISOS), s12c2_1(ISOS),
c     &          s21c2_1(ISOS), s22c2_1(ISOS), s3c2_1(ISOS),
c     &          wd1c2_1(ISOS), wd2c2_1(ISOS), wd3c2_1(ISOS)
c      real      mt1c2_2(ISOS), mt2c2_2(ISOS), st1c2_2(ISOS),
c     &          st2c2_2(ISOS), s11c2_2(ISOS), s12c2_2(ISOS),
c     &          s21c2_2(ISOS), s22c2_2(ISOS), s3c2_2(ISOS),
c     &          wd1c2_2(ISOS), wd2c2_2(ISOS), wd3c2_2(ISOS)
      double precision crpGrossPsn, forGrossPsn
      double precision CH4_oxid
      real             avgst_10cm, Com
      real             TI, SI, Cr, Eh, Feh, CH4_prod, CH4_Ep, CH4_Ebl

      data plntd / 0 /
      data curday / 0 /
      data scenfrac /1.0/
      data respsum /2*0.0/
      data arspsum /4*0.0/
      data Eh / 300.0 /

c ... Saved variables
      save plntd, CO2resp, respsum, arspsum, scenfrac, Eh

c ... Initialize local variables
      cancvr = 0.0
      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0
      ddeloi = 0.0
      ddeloe = 0.0
      ddsrfclit = 0.0
      ddsmnrl = 0.0
      ddhetresp = 0.0
      ddsoilresp = 0.0
      ddccarbostg = 0.0
      ddfcarbostg = 0.0
      ddcmresp = 0.0
      ddfmresp = 0.0
      ddcgresp = 0.0
      ddfgresp = 0.0

c ... Added below for savanna model (rm)
      if (cursys .eq. SAVSYS ) then
        wdbmas = (fbrchc + rlwodc) * 2.0
c ..... Can get a divide by zero error when there is no wood biomass,
c ..... add a minimum wood biomass so that trees can grow from nothing,
c ..... cak - 10/08/02
        if (wdbmas .le. 0.0) then
          wdbmas = 50
        endif
c ..... Change the way that tree basal area is calculated,
c ..... cak 12/19/01
c        trbasl = wdbmas / basfct
        basf = (wdbmas/(0.88 * ((wdbmas * 0.01)**0.635)))
        if (basf .lt. 250.0) then
          basf = basf * basfct
        endif
        trbasl = wdbmas / basf
        cancvr = 1 - exp(-0.064 * trbasl)
        if (trbasl .le. 1.0E-6) then
c          trbasl = 0.1
          trbasl = 0.3
        endif
      endif

c ... Set aminrl for use in routines called from decomp
      do 10 iel = 1, nelem
        if (iel .eq. P) then
          fsol = fsfunc(minerl(1,P), pslsrb, sorpmx)
        else
          fsol = 1.0
        endif
        aminrl(iel) = minerl(1,iel) * fsol
10    continue

c ... Initialize accumulators
      call cycle()

c ... N Fixation
c ... Does not take into account the effect of irrigation
      if (nsnfix .eq. 1 .and. nelem .ge. P) then

c ..... Compute mineral N:P ratio for N-Fixation (use suface layer only)
c ..... rnpml1 is used to control soil N-fixation using a regression
c ..... equation based on Kansas data. This ratio is used only if nelem = 2.
c ..... rnpml1 is flagged if either minerl(1,1) or minerl(1,2) is zero.

        rnpml1 = minerl(1,N)/minerl(1,P)*
     &           fsfunc(minerl(1,P),pslsrb,sorpmx)

c ..... Wet-dry fixation of nitrogen -- monthly flow
c ..... Atmospheric fixation is split between monthly dry fall and
c ..... wdfnp is the N:P ratio control function for non-symbiotic
c ..... soil fixation.
c ..... Both wdfnp and fxbiom are relative functions
c ..... which vary from 0 to 1.
c ..... wdfnp computed as a negative natural log of rnpml1
c ..... symnfx is the symbiotic N fixation by legumes derived from Cole and
c ..... Heil (1981) using data from Walker et al. 1959.
C       if (rnpml1 .eq. 0) then
        if (rnpml1 .le. 0) then
          wdfnp = 1.
        else
          wdfnp = min(1., ((-alog(rnpml1))/fxnpb)-.2)
        endif

c ..... The definitions of fxmca and fxmcb originally refered to water,
c ..... not biomass. (Alister 11/91)
        tbiom = aglivc+stdedc+strucc(SRFC)
        biof  = fxmca + tbiom * fxmcb
        fxbiom = 1 - biof
        fxbiom = min(1.,fxbiom)
        if (wdfnp .lt. 0 .or. fxbiom .lt. 0 .or. stemp .lt. 7.5) then
          fwdfx = 0.0
        else
          fwdfx = wdfnp * fxbiom
        endif

c ..... Compute N-fixation for the month

c ..... Wet fall depending on the monthly precipitation (wdfxma)
c        wdfxma = wdfxa *  prcurr(month)/prcann
c ..... Add an optional multiplier on N deposition, cak - 04/05/04
        if (Ninput .eq. 2 .or. Ninput .eq. 3) then
          wdfxma = baseNdep *
     &             ((precip(month) * precscalar(month)) / prcann) *
     &             Nscalar(month)
        else
          wdfxma = baseNdep *
     &             ((precip(month) * precscalar(month)) / prcann)
        endif
        wdfxms = fxmxs * fwdfx
        wdfxm  = wdfxma + wdfxms

c ..... Compute annual N-fixation accumulators for atmosphere and soils
        frac_nh4 = 0.5
        frac_no3 = 0.5
        wdfxas = wdfxas + wdfxms
        wdfxaa = wdfxaa + wdfxma
        clyr = 1
        subname = 'simsom1a  '
        call update_npool(clyr, wdfxm, frac_nh4, frac_no3, ammonium,
     &                    nitrate, subname)
        call flow(esrsnk(N),minerl(1,N),time,wdfxm)
c ..... nfixac should accumulate SYMBIOTIC N-fixation -mdh 11/16/01
c        nfixac = nfixac+wdfxm

c ... Monthly N-fixation based on annual parameters
      else
c ..... USE PRCURR/PRCANN INSTEAD OF DT
c        wdfxms = wdfxs * prcurr(month)/prcann
c        wdfxma = wdfxa * prcurr(month)/prcann

        frac_nh4 = 0.5
        frac_no3 = 0.5
        wdfxms = wdfxs * ((precip(month) * precscalar(month)) / prcann)
c ..... Add an optional multiplier on N deposition, cak - 04/05/04
        if (Ninput .eq. 2 .or. Ninput .eq. 3) then
          wdfxma = baseNdep *
     &             ((precip(month) * precscalar(month)) / prcann) *
     &             Nscalar(month)
        else
          wdfxma = baseNdep *
     &             ((precip(month) * precscalar(month)) / prcann)
        endif
        wdfxas = wdfxas + wdfxms
        wdfxaa = wdfxaa + wdfxma
        wdfxm = wdfxma + wdfxms
        clyr = 1
        subname = 'simsom2a  '
        call update_npool(clyr, wdfxm, frac_nh4, frac_no3, ammonium,
     &                    nitrate, subname)
        call flow(esrsnk(N),minerl(1,N),time,wdfxm)
c ..... nfixac should accumulate SYMBIOTIC N-fixation -mdh 11/16/01
c        nfixac = nfixac + wdfxm
      endif

c ... Accumulate values for annual N deposition output variables,
c ... cak - 04/05/04
      wdfx = wdfx + wdfxma
      wdfxa = wdfxa + wdfxma

c ... Monthly atmospheric S deposition
c ... Note: irract is last month's irrigation at this point -mdh 12/9/96
      if (nelem .eq. S) then
c        satm = satmt * prcurr(month) / prcann
        satm = satmt * ((precip(month) * precscalar(month)) / prcann)
        satmac = satmac + satm
        if (doirri) then
          sirr = sirri * irract * 0.01
        else
          sirr = 0
        endif
        sirrac = sirrac + sirr
        stot = satm + sirr
        call flow(esrsnk(S),minerl(SRFC,S),time,stot)
      endif

c ... BEGIN MONTHLY INITIALIZATION

c ... Determine the number of days in each month.  The leap year exception
c ... will be handled in getwth.  -mdh 12/10/96
      if (curday .le. 1) then
        do 110 imo = 1, 12
          dysimo(imo) = idysimo(imo)
          lstdy(imo) = ilstdy(imo)
          frstdy(imo) = ifrstdy(imo)
110     continue
      endif

      do 20 ii = 1, 8
        prevstream(ii) = 0
20    continue

c      prevgromin = 0.0

c ... Reset the monthly accumulator
      call mthacc(agdefacsum, bgdefacsum)

      CO2resp = 0

c ... BEGIN DAILY LOOP...

      do 200 curday = frstdy(month), lstdy(month)

c ..... Calculate the root priming effect on the som2c(2) decomposition
c ..... rate, CAK - 01/28/2014
        crpeff = 1.0
        trpeff = 1.0
c ..... Crop/grass root priming effect
        if (crpindx .eq. 0) then
c ....... No root priming effect on som2c(2) decomposition
          crpeff = 1.0
        elseif (crpindx .eq. 1) then
c ....... Root priming effect on som2c(2) decomposition based on total
c ....... soil respiration
          crpeff = ramp(dsresp, crpcmn, crpmnmul, crpcmx, crpmxmul)
        elseif (crpindx .eq. 2) then
c ....... Root priming effect on som2c(2) decompostion based on soil
c ....... hetertrophic respiration only
          crpeff = ramp(dhresp, crpcmn, crpmnmul, crpcmx, crpmxmul)
        elseif (crpindx .eq. 3) then
c ....... Root priming effect on som2c(2) decomposition based on fine
c ....... root production 
          crpeff = ramp(mcprd(BELOWJ)+mcprd(BELOWM), crpcmn, crpmnmul,
     &                  crpcmx, crpmxmul)
        endif
c ..... Tree root priming effect
        if (trpindx .eq. 0) then
c ....... No root priming effect on som2c(2) decomposition
         trpeff = 1.0
        elseif (trpindx .eq. 1) then
c ....... Root priming effect on som2c(2) decomposition based on total
c .......  soil respiration
          trpeff = ramp(dsresp, trpcmn, trpmnmul, trpcmx, trpmxmul)
        elseif (trpindx .eq. 2) then
c ....... Root priming effect on som2c(2) decompostion based on soil
c ....... hetertrophic respiration only
          trpeff = ramp(dhresp, trpcmn, trpmnmul, trpcmx, trpmxmul)
        elseif (trpindx .eq. 3) then
c ....... Root priming effect on som2c(2) decomposition based on fine
c ....... root production 
          trpeff = ramp(mfprd(FROOTJ)+mfprd(FROOTM), trpcmn, trpmnmul,
     &                 trpcmx, trpmxmul)
        endif
        if (cursys .eq. CRPSYS) then
          rpeff = crpeff
        elseif (cursys .eq. FORSYS) then
          rpeff = trpeff
        elseif (cursys .eq. SAVSYS) then
c ....... If there is fine root growth weight the root priming effect
c ....... based on crop/grass vs. tree fine root growth
          if (mcprd(BELOWJ)+mcprd(BELOWM) .gt. 0.0 .or.
     &        mfprd(FROOTJ)+mfprd(FROOTM) .gt. 0.0) then
            rpeff = crpeff * ((mcprd(BELOWJ) + mcprd(BELOWM)) /
     &                        (mcprd(BELOWJ) + mcprd(BELOWM) +
     &                         mfprd(FROOTJ) + mfprd(FROOTM))) +
     &              trpeff * ((mfprd(FROOTJ) + mfprd(FROOTM)) /
     &                        (mcprd(BELOWJ) + mcprd(BELOWM) +
     &                         mfprd(FROOTJ) + mfprd(FROOTM)))
          else
c ......... If there is no fine root production weight the root priming
c ......... effect based on crop/grass vs. tree fine root biomass
            rpeff = crpeff * ((bglivcj + bglivcm) /
     &                        (bglivcj + bglivcm + frootcj + frootcm))+
     &              trpeff * ((frootcj + frootcm) /
     &                        (bglivcj + bglivcm + frootcj + frootcm))
          endif
        endif

c ..... Initialize NPP and growth and maintenance respiration varaibles
        prev_mcprd1 = mcprd(ABOVE)
        prev_mcprd2 = mcprd(BELOWJ)
        prev_mcprd3 = mcprd(BELOWM)
        mcprd(ABOVE) = 0.0
        mcprd(BELOWJ) = 0.0
        mcprd(BELOWM) = 0.0
        mrspdyflux(CRPSYS) = 0.0
        cmrspdyflux(ABOVE) = 0.0
        cmrspdyflux(BELOWJ) = 0.0
        cmrspdyflux(BELOWM) = 0.0
        grspdyflux(CRPSYS) = 0.0
        cgrspdyflux(ABOVE) = 0.0
        cgrspdyflux(BELOWJ) = 0.0
        cgrspdyflux(BELOWM) = 0.0
        mrspdyflux(FORSYS) = 0.0
        grspdyflux(FORSYS) = 0.0
        do 210 ipart=1, FPARTS
          fgrspdyflux(ipart) = 0.0
          fmrspdyflux(ipart) = 0.0 
          mfprd(ipart) = 0.0
210     continue

c ..... Call schedl to determine scheduling options for this day of the
c ..... year
        call schedl(curday)

        tfrac = 1.0/real(dysimo(month))

c ..... Get a day's worth of weather from the weather file
        call getwth(curday, month, tempmax, tempmin, avgtemp, ppt,
     &              solrad, rhumid, windsp, tavewk, petdly, fwloss,
     &              sitlat, snow, tmn2m, tmx2m, tavemth, wthinput,
     &              precscalar, tmaxscalar, tminscalar, srad, vpd)
c ..... Set the value of tave for output to the *.bin file
        tave = tavemth

c ..... Calculate the average daily temperature for use in the growth
c ..... equations.  This code captures a diurnal temperature variation
c ..... by using the day length to compute the values for the multipliers
c ..... on daily minimum and maximum air temperature
        if (daylength(curday) .lt. 12.0) then
          tmns_mlt = ((12.0 - daylength(curday)) * 3.0 + 12.0) / 24.0
        else
          tmns_mlt = ((12.0 - daylength(curday)) * 1.2 + 12.0) / 24.0
        endif
        tmns_mlt = min(0.95, tmns_mlt)
        tmns_mlt = max(0.05, tmns_mlt)
        tmxs_mlt = 1.0 - tmns_mlt
        tavedly = tmxs_mlt * tempmax(curday) +
     &            tmns_mlt * tempmin(curday)

        irractwk = 0
        if (doirri .and. (irriday .eq. curday) .or.
     &      (irricnt .gt. 0. .and. irricnt .lt. 31)) then
          if (irriday .eq. curday) then
            irricnt = 0
          endif
          if (mod(irricnt, 7) .eq. 0) then
            irractwk = irrigt(ppt(curday), petdly, 7.0/30.0)
          endif
          irricnt = irricnt + 1
        else
          irractwk = 0
          irricnt = 0
        endif

        irract = irract + irractwk
        rain = rain + ppt(curday)
        ppt(curday) = ppt(curday) + irractwk
        annppt = annppt + ppt(curday)

c ..... Set the flags for growth start before other events based on the
c ..... day of the year the event is scheduled to occur
        if (dofone .and. (foneday .eq. curday)) then
          forgrw = 1
          fgrwdys = 0
          fpsndys = 0
          call initprod(FORSYS, month)
        endif
        if (dofrst .and. (frstday .eq. curday)) then
          seedl = 0
          crpgrw = 1
          cgrwdys = 0
          cpsndys = 0
          thermunits = 0.0
          accumdd = .true.
          plntkill = .false.
          call initprod(CRPSYS, month)
          scenfrac = 1.0
        endif
        if (doplnt .and. (plntday .eq. curday)) then
          seedl = 1
          plntd = 1
          crpgrw = 1
          cgrwdys = 0
          cpsndys = 0
          falprc = 0
          prcfal = 0.0
          thermunits = 0.0
          emerg = .false.
          accumdd = .true.
          plntcnt = 0
          plntkill = .false.
          grnfill = .false.
          gwstress = 0.0
          grnfldys = 0
          call initprod(CRPSYS, month)
          scenfrac = 1.0
        endif

c ..... Check if days since planting (plntcnt) needs to be updated
        if (plntd .eq. 1 .and. stemp .ge. rtdtmp) then
          plntcnt = plntcnt + 1
        endif

c ..... Check if days since the senescence event (senecnt) needs to
c ..... be updated, cak - 04/17/03
        if (senecnt .ge. 1) then
          senecnt = senecnt + 1
          if (senecnt .gt. 30) then
            senecnt = 0
          endif
        endif

c ..... Effect of cultivation on decomposition (used in decomp routine)
c ..... This code has been moved to the simsom subroutine, cak - 04/17/03
c ..... Determine effect of cultivation on decomposition. vk 03-13-91
c ..... cltfac is this month's effect of cultivation on decomposition
c ..... of som1, som2, som3, and structural.  It is set to clteff
c ..... in months when cultivation occurs; otherwise it equals 1.
c ..... clteff is the effect of cultivation on decomposition read from
c ..... the cult.100 file
        if (docult .and. (cultday .eq. curday) .or.
     &      (cultcnt .gt. 0. .and. cultcnt .lt. 31)) then
          do 34 ii = 1, 4
            cltfac(ii) = clteff(ii)
34        continue
          if (cultday .eq. curday) then
            cultcnt = 0
          endif
          cultcnt = cultcnt + 1
        else
          do 35 ii = 1, 4
            cltfac(ii) = 1.0
35        continue
          cultcnt = 0
        endif

c ..... Reset the reduction factor on nitrification rates due to
c ..... nitrification inhibitors as necessary, cak - 12/01/03
        if ((fertcnt .gt. 0) .and. (fertcnt .lt. (ninhtm * 7.0))) then
          fertcnt = fertcnt + 1
        else
          fertcnt = 0
          nreduce = 1.0
        endif

c ..... Save previous stream values from inorganic and organic leaching.  -mdh 1/97
        do 85 iel = 1, nelem
          prevstream(iel+1) = stream(iel+1)
          prevstream(iel+5) = stream(iel+5)
85      continue
        strm5l = 0.0
        strm5u = 0.0

c ..... Track respiration over the daily timestep
        do 30 iso = 1, ISOS
          autoresp1(iso) = cautoresp(iso) + fautoresp(iso)
          hetresp1(iso) = mt1c2(iso) + mt2c2(iso) + 
     &                    st1c2(iso) + st2c2(iso) +
     &                    s11c2(iso) + s12c2(iso) +
     &                    s21c2(iso) + s22c2(iso) +
     &                    s3c2(iso)  + wd1c2(iso) +
     &                    wd2c2(iso) + wd3c2(iso)
          mnrlrsp1(iso) = mt2c2(iso) + st2c2(iso) +
     &                    s12c2(iso) + s22c2(iso) +
     &                    s3c2(iso)  + wd3c2(iso)
          oiresp1(iso) = mt1c2(iso) + st1c2(iso) +
     &                   s11c2(iso)
          oeresp1(iso) = s21c2(iso)
          slitrsp1(iso) = mt1c2(iso) + st1c2(iso) +
     &                    s11c2(iso) + s21c2(iso) +
     &                    st1uvc2(iso) + stduvc2(iso)
c          mt1c2_1(iso) = mt1c2(iso)
c          mt2c2_1(iso) = mt2c2(iso)
c          st1c2_1(iso) = st1c2(iso)
c          st2c2_1(iso) = st2c2(iso)
c          s11c2_1(iso) = s11c2(iso)
c          s12c2_1(iso) = s12c2(iso)
c          s21c2_1(iso) = s21c2(iso)
c          s22c2_1(iso) = s22c2(iso)
c          s3c2_1(iso)  = s3c2(iso)
c          wd1c2_1(iso) = wd1c2(iso)
c          wd2c2_1(iso) = wd2c2(iso)
c          wd3c2_1(iso) = wd3c2(iso)
30      continue

c ..... Set frlech to leaching fraction vek june90
c ..... Compute normal value for frlech.  Recompute in flood routine
c ..... if flooding occurs.

        texeff = fleach(1) + fleach(2) * sand
        do 50 iel = 1, nelem
          if (iel .eq. P) then
            fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
          else
            fsol = 1.0
          endif
          frlech(iel) = texeff * fleach(iel+2) * fsol

50      continue

c ..... Determine co2 effect on transpiration, pass to dailymoist
        if (cursys .eq. SAVSYS) then
C         if (aglivc + rleavc .eq. 0.0) then
          if (aglivc + rleavc .le. 0.0) then
            co2val = 1.0
          else
            co2val = (co2ctr(CRPSYS)*aglivc + co2ctr(FORSYS)*rleavc) /
     &               (aglivc + rleavc)
          endif
        else
          co2val = co2ctr(cursys)
        endif
 
        call dailymoist(curday, agdefacsum, bgdefacsum, bgwfunc,
     &                  frlech, co2val, saccum, evapdly, intrcpt, melt,
     &                  outflow, petdly, runoffdly, sublim, trandly,
     &                  avgstemp, scenfrac, rpeff, watr2sat,
     &                  prev_mcprd2 + prev_mcprd3, Com, avgst_10cm, TI,
     &                  SI, Cr, Eh, Feh, CH4_prod, CH4_Ep, CH4_Ebl,
     &                  CH4_oxid, aetdly)

c ..... Compute potential production 
        call potprod(cancvr, tempmax(curday), tempmin(curday), tavedly,
     &               ppt(curday), petdly, tfrac, tavemth, curday,
     &               scenfrac, aetdly, srad(curday), vpd(curday),
     &               crpGrossPsn, forGrossPsn, tminslope,
     &               tminintercept)

c ..... Volatilization loss of nitrogen as a function of
c ..... gross mineralization
c ..... Compute the grosmin which occured since the last time step
c        grominwk = gromin(1) - prevgromin
c        prevgromin = gromin(1)
        volgm = 0.0
c        volgm = vlossg*gromin(1)
c        volgm = vlossg*grominwk
c        minerl(SRFC,N) = minerl(SRFC,N) - volgm
c        write(*,*) 'volgma (volatilization) = ', volgma
c        esrsnk(N) = esrsnk(N) + volgm

c ..... Soil erosion
        if (doerod .and. (erodday .eq. curday) .or.
     &     (erodcnt .gt. 0. .and. erodcnt .lt. 31)) then
          call erosn(psloss*tfrac,bulkd,edepth,enrich,lhzci,lhze,nelem)
          if (erodday .eq. curday) then
            erodcnt = 0
          endif
          erodcnt = erodcnt + 1
        else
          scloss = 0.0
          erodcnt = 0
        endif

c ..... Fertilization option
c ..... This code has been moved to dailymoist, cak - 04/17/03

c ..... Available nutrients
c ..... tminrl is the total amount of each element available in
c ..... mineral form.

        do 80 iel = 1, nelem
          tminrl(iel) = 0.
          if (iel .eq. P) then
            fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
          else 
            fsol = 1.0
          endif

          do 70 lyr = 1, nlayer

c ......... Plants can only uptake from a layer with a positive
c ......... value, so only the positive layers are summed here.

            if (minerl(lyr,iel) .gt. 0.)  then
              tminrl(iel) = tminrl(iel) + minerl(lyr,iel) * fsol
            endif
70        continue
80      continue

c        write(*,*) 'SIMSOM: available nutrients = ', tminrl(1)

c ..... Compute the fraction of labile (non-sorbed) P in the surface
c ..... layer available to plants

        favail(2) = max(favail(4),
     &              min(favail(4) + minerl(SRFC,N)*
     &              (favail(5) - favail(4)) / favail(6),
     &              favail(5)))

c ..... Add to fallow rain
        if (falprc .eq. 1) then
c          prcfal = prcfal + prcurr(month)
          prcfal = prcfal + ppt(curday)
        endif

c ..... Call the producer submodels

c ..... Determine if daylength is increasing or decreasing
        if (daylength(curday) .lt. dayhrs) then
          hrsinc = .FALSE.
        else if (daylength(curday) .gt. dayhrs) then
          hrsinc = .TRUE.
        endif

c ..... Save the today's day length for use in tomorrow's calculation
        dayhrs = daylength(curday)

c ..... Crop and Forest removal options - moved here from CROP
c ..... and TREES so the mineral pools are not radically changed
c ..... between the growth routines. - rm 7/94

        if (cursys .eq. CRPSYS) then
          call crop(time, bgwfunc, tfrac, tavedly, curday, avgstemp,
     &              crpGrossPsn)
          if ((dofire(CRPSYS) .and. (fireday .eq. curday)) .or.
     &        (dograz .and. (grazday .eq. curday)) .or.
     &        (grazcnt .gt. 1)) then
            if (dofire(CRPSYS)) then
c              firecnt = 1
            endif
            if (dograz) then
              grazcnt = 1
            endif
c ......... Burning of aboveground live, standing dead, and litter
c ......... layer or grazing
            call grem(tfrac)
          endif

        else if (cursys .eq. FORSYS) then
          call trees(bgwfunc, tavewk, tfrac, curday, tavedly, avgstemp,
     &               forGrossPsn)
          if (dotrem .and. (tremday .eq. curday)) then
c ......... Burning of live wood or cutting events
            call frem()
          endif
          if (dofire(FORSYS) .and. (fireday .eq. curday)) then
c            firecnt = 1
c ......... Burning of dead wood and litter layer, cak - 08/23/02
            call grem(tfrac)
          endif

        else if (cursys .eq. SAVSYS) then
          call crop(time, bgwfunc, tfrac, tavedly, curday, avgstemp,
     &              crpGrossPsn)
          call trees(bgwfunc, tavewk, tfrac, curday, tavedly, avgstemp,
     &               forGrossPsn)
          if (dotrem .and. (tremday .eq. curday)) then
c ......... Burning of live wood or cutting events
            call frem()
          endif
          if ((dofire(SAVSYS) .and. (fireday .eq. curday)) .or.
     &        (dograz .and. (grazday .eq. curday)) .or.
     &        (grazcnt .gt. 1)) then
            if (dofire(SAVSYS)) then
c              firecnt = 1
            endif
            if (dograz) then
              grazcnt = 1
            endif
c ......... Burning of aboveground live, standing dead, litter
c ......... layer, and dead wood or grazing
            call grem(tfrac)
          endif
        endif

c ..... Update state variables and accumulators and sum carbon isotopes.
        call flowup(time)
        call flowup_double(time)
        call flowup_double_in(time)
        call flowup_double_out(time)
        call sumcar

c ..... Track respiration over the daily timestep
        do 32 iso = 1, ISOS
          autoresp2(iso) = cautoresp(iso) + fautoresp(iso)
          hetresp2(iso) = mt1c2(iso) + mt2c2(iso) + 
     &                    st1c2(iso) + st2c2(iso) +
     &                    s11c2(iso) + s12c2(iso) +
     &                    s21c2(iso) + s22c2(iso) +
     &                    s3c2(iso)  + wd1c2(iso) +
     &                    wd2c2(iso) + wd3c2(iso)
          mnrlrsp2(iso) = mt2c2(iso) + st2c2(iso) +
     &                    s12c2(iso) + s22c2(iso) +
     &                    s3c2(iso)  + wd3c2(iso)
          oiresp2(iso) = mt1c2(iso) + st1c2(iso) +
     &                   s11c2(iso)
          oeresp2(iso) = s21c2(iso)
          slitrsp2(iso) = mt1c2(iso) + st1c2(iso) +
     &                    s11c2(iso) + s21c2(iso) +
     &                    st1uvc2(iso) + stduvc2(iso)
c          mt1c2_2(iso) = mt1c2(iso)
c          mt2c2_2(iso) = mt2c2(iso)
c          st1c2_2(iso) = st1c2(iso)
c          st2c2_2(iso) = st2c2(iso)
c          s11c2_2(iso) = s11c2(iso)
c          s12c2_2(iso) = s12c2(iso)
c          s21c2_2(iso) = s21c2(iso)
c          s22c2_2(iso) = s22c2(iso)
c          s3c2_2(iso)  = s3c2(iso)
c          wd1c2_2(iso) = wd1c2(iso)
c          wd2c2_2(iso) = wd2c2(iso)
c          wd3c2_2(iso) = wd3c2(iso)
          newautoresp(iso) = autoresp2(iso) - autoresp1(iso)
          newhetresp(iso) = hetresp2(iso) - hetresp1(iso)
          newmnrlrsp(iso) = mnrlrsp2(iso) - mnrlrsp1(iso)
          newoiresp(iso) = oiresp2(iso) - oiresp1(iso)
          newoeresp(iso) = oeresp2(iso) - oeresp1(iso)
          newslitrsp(iso) = slitrsp2(iso) - slitrsp1(iso)
32      continue

c ..... Calculate daily output values for the dels.out file
        dcmresp = cmrspdyflux(ABOVE)  + cmrspdyflux(BELOWJ) +
     &            cmrspdyflux(BELOWM)
        dfmresp = fmrspdyflux(LEAF)   + fmrspdyflux(FROOTJ) +
     &            fmrspdyflux(FROOTM) + fmrspdyflux(FBRCH)  +
     &            fmrspdyflux(LWOOD)  + fmrspdyflux(CROOT)
        dcgresp = cgrspdyflux(ABOVE)  + cgrspdyflux(BELOWJ) +
     &            cgrspdyflux(BELOWM)
        dfgresp = fgrspdyflux(LEAF)   + fgrspdyflux(FROOTJ) +
     &            fgrspdyflux(FROOTM) + fgrspdyflux(FBRCH)  +
     &            fgrspdyflux(LWOOD)  + fgrspdyflux(CROOT)

c ..... Calculate daily output values for the resp.out file
        doiresp = (oiresp2(UNLABL) + oiresp2(LABELD)) -
     &            (oiresp1(UNLABL) + oiresp1(LABELD))
        doeresp = (oeresp2(UNLABL) + oeresp2(LABELD)) -
     &            (oeresp1(UNLABL) + oeresp1(LABELD))
        dslitrsp = (slitrsp2(UNLABL) + slitrsp2(LABELD)) -
     &             (slitrsp1(UNLABL) + slitrsp1(LABELD))
        dmnrlrsp = (mnrlrsp2(UNLABL) + mnrlrsp2(LABELD)) -
     &             (mnrlrsp1(UNLABL) + mnrlrsp1(LABELD))
        dhresp = (hetresp2(UNLABL) + hetresp2(LABELD)) -
     &           (hetresp1(UNLABL) + hetresp1(LABELD))
        dcrtjresp = cmrspdyflux(BELOWJ) + cgrspdyflux(BELOWJ)
        dcrtmresp = cmrspdyflux(BELOWM) + cgrspdyflux(BELOWM)
        dfrtjresp = fmrspdyflux(FROOTJ) + fgrspdyflux(FROOTJ)
        dfrtmresp = fmrspdyflux(FROOTM) + fgrspdyflux(FROOTM)
        dfrtcresp = fmrspdyflux(CROOT)  + fgrspdyflux(CROOT)
        dmresp = cmrspdyflux(ABOVE)  + cmrspdyflux(BELOWJ) +
     &           cmrspdyflux(BELOWM) + fmrspdyflux(LEAF)   +
     &           fmrspdyflux(FROOTJ) + fmrspdyflux(FROOTM) +
     &           fmrspdyflux(FBRCH)  + fmrspdyflux(LWOOD)  +
     &           fmrspdyflux(CROOT)
        dgresp = cgrspdyflux(ABOVE)  + cgrspdyflux(BELOWJ) +
     &           cgrspdyflux(BELOWM) + fgrspdyflux(LEAF)   +
     &           fgrspdyflux(FROOTJ) + fgrspdyflux(FROOTM) +
     &           fgrspdyflux(FBRCH)  + fgrspdyflux(LWOOD)  +
     &           fgrspdyflux(CROOT)
        dcrtresp = dcrtjresp + dcrtmresp
        dfrtresp = dfrtjresp + dfrtmresp + dfrtcresp
        dsresp = dcrtresp + dfrtresp + dhresp
c        dmt1c2 = (mt1c2_2(UNLABL) + mt1c2_2(LABELD)) -
c     &           (mt1c2_1(UNLABL) + mt1c2_1(LABELD))
c        dmt2c2 = (mt2c2_2(UNLABL) + mt2c2_2(LABELD)) -
c     &           (mt2c2_1(UNLABL) + mt2c2_1(LABELD))
c        dst1c2 = (st1c2_2(UNLABL) + st1c2_2(LABELD)) -
c     &           (st1c2_1(UNLABL) + st1c2_1(LABELD))
c        dst2c2 = (st2c2_2(UNLABL) + st2c2_2(LABELD)) -
c     &           (st2c2_1(UNLABL) + st2c2_1(LABELD))
c        ds11c2 = (s11c2_2(UNLABL) + s11c2_2(LABELD)) -
c     &           (s11c2_1(UNLABL) + s11c2_1(LABELD))
c        ds12c2 = (s12c2_2(UNLABL) + s12c2_2(LABELD)) -
c     &           (s12c2_1(UNLABL) + s12c2_1(LABELD))
c        ds21c2 = (s21c2_2(UNLABL) + s21c2_2(LABELD)) -
c     &           (s21c2_1(UNLABL) + s21c2_1(LABELD))
c        ds22c2 = (s22c2_2(UNLABL) + s22c2_2(LABELD)) -
c     &           (s22c2_1(UNLABL) + s22c2_1(LABELD))
c        ds3c2  = (s3c2_2(UNLABL)  + s3c2_2(LABELD))  -
c     &           (s3c2_1(UNLABL)  + s3c2_1(LABELD))
c        dwd1c2 = (wd1c2_2(UNLABL) + wd1c2_2(LABELD)) -
c     &           (wd1c2_1(UNLABL) + wd1c2_1(LABELD))
c        dwd2c2 = (wd2c2_2(UNLABL) + wd2c2_2(LABELD)) -
c     &           (wd2c2_1(UNLABL) + wd2c2_1(LABELD))
c        dwd3c2 = (wd3c2_2(UNLABL) + wd3c2_2(LABELD)) -
c     &           (wd3c2_1(UNLABL) + wd3c2_1(LABELD))

c ..... Harvest may be performed after updating flows.  Put here for
c ..... consistency with the Savanna model - moved calls to flowup, 
c ..... sumcar and harvst from CROP routine to here. -rm 7/94
        if (dohrvt .and. (hrvtday .eq. curday)) then
          call harvst(month, pltlig, curday)
          plntd = 0
          falprc = 1
          prcfal = 0.0
          harmth = 1
          grnfill = .false.
          gwstress = 0.0
          grnfldys = 0
        endif

c ..... Set the flags for growth end after other events based on the
c ..... day of the year the event is scheduled to occur
        if (doflst .and. (flstday .eq. curday)) then
          forgrw = 0
          fgrwdys = 0
          fpsndys = 0
          call inprac(FORSYS)
        endif
        if (dolast .and. (lastday .eq. curday)) then
          crpgrw = 0
          thermunits = 0.0
          accumdd = .false.
          cgrwdys = 0
          cpsndys = 0
          plntcnt = 0
          call inprac(CRPSYS)
c ....... For an annual plant reset the rooting depth back to 1 on a
c ....... LAST event
          if ((frtcindx .eq. 2) .or. (frtcindx .ge. 4)) then
            claypg = 1
c ......... If this is an annual plant flow the carbohydrate storage pool
c ......... to the csrsnk
            csrsnk(UNLABL) = csrsnk(UNLABL) + carbostg(CRPSYS,UNLABL)
            csrsnk(LABELD) = csrsnk(LABELD) + carbostg(CRPSYS,LABELD)
            carbostg(CRPSYS,UNLABL) = 0.0
            carbostg(CRPSYS,LABELD) = 0.0
c ......... For an annual plant simulated using the growing degree day
c ......... submodel set the photosynthetic active carbon to zero
            if (frtcindx .ge. 4) then
              scenfrac = 0.0
            endif
          endif
        endif

c ..... Update state variables and accumulators and sum carbon isotopes.
        call flowup(time)
        call flowup_double(time)
        call flowup_double_in(time)
        call flowup_double_out(time)
        call sumcar

c ..... Accumulate leached C,N,P,S
        csrsnk(UNLABL) = csrsnk(UNLABL) + strm5u
        csrsnk(LABELD) = csrsnk(LABELD) + strm5l
        stream(5) = stream(5) + strm5u + strm5l
        do 90 iel = 1, nelem
          esrsnk(iel) = esrsnk(iel)+(stream(iel+1)-prevstream(iel+1))+
     &                  (stream(iel+5)-prevstream(iel+5))
90      continue

c ..... Volatilization loss as a function of the mineral n which
c ..... remains after uptake by plants

        volex = 0.0
c        if (minerl(SRFC,N) .gt. 0.0) then
c          volex = vlosse*minerl(SRFC,N)*dt
c          volex = vlosse*minerl(SRFC,N)*dt*tfrac
c          minerl(SRFC,N) = minerl(SRFC,N) - volex
c          write(*,*) 'volex (volatilization) = ', volex
c          esrsnk(N) = esrsnk(N) + volex
c        endif

c ..... Volatilization
c        volgma = volgma+volgm
c        volexa = volexa+volex
c        volgac = volgac+volgm
c        voleac = voleac+volex

        if (time .ge. strplt) then
          call wrtbio(time, curday, aglivc, bglivcj, bglivcm, aglive(N),
     &                bglivej(N), bglivem(N), rleavc, frootcj, frootcm,
     &                fbrchc, rlwodc, crootc, h2ogef(1), h2ogef(2))
          call wrtresp(time, curday, doiresp, doeresp, dslitrsp,
     &                 dmnrlrsp, dhresp, dcrtjresp, dcrtmresp,
     &                 dfrtjresp, dfrtmresp, dfrtcresp, dsresp,
     &                 dmresp, dgresp, mrspdyflux(CRPSYS),
     &                 mrspdyflux(FORSYS), cmrspdyflux(ABOVE),
     &                 cmrspdyflux(BELOWJ), cmrspdyflux(BELOWM),
     &                 fmrspdyflux(LEAF), fmrspdyflux(FROOTJ),
     &                 fmrspdyflux(FROOTM), fmrspdyflux(FBRCH),
     &                 fmrspdyflux(LWOOD), fmrspdyflux(CROOT),
     &                 mrspann(CRPSYS), mrspann(FORSYS), tavedly,
     &                 mrspTempEffect(CRPSYS,SRFC),
     &                 mrspTempEffect(CRPSYS,SOIL),
     &                 mrspWaterEffect(CRPSYS),
     &                 mrspTempEffect(FORSYS,SRFC),
     &                 mrspTempEffect(FORSYS,SOIL),
     &                 mrspWaterEffect(FORSYS), grspdyflux(CRPSYS),
     &                 grspdyflux(FORSYS), cgrspdyflux(ABOVE),
     &                 cgrspdyflux(BELOWJ), cgrspdyflux(BELOWM),
     &                 fgrspdyflux(LEAF), fgrspdyflux(FROOTJ),
     &                 fgrspdyflux(FROOTM), fgrspdyflux(FBRCH),
     &                 fgrspdyflux(LWOOD), fgrspdyflux(CROOT),
     &                 grspann(CRPSYS), grspann(FORSYS),
     &                 mcprd(ABOVE), mcprd(BELOWJ), mcprd(BELOWM),
     &                 mfprd(LEAF), mfprd(FROOTJ), mfprd(FROOTM),
     &                 mfprd(FBRCH), mfprd(LWOOD), mfprd(CROOT),
     &                 carbostg(CRPSYS,UNLABL),
     &                 carbostg(CRPSYS,LABELD),
     &                 carbostg(FORSYS,UNLABL),
     &                 carbostg(FORSYS,LABELD))
          call wrtlivec(time, curday, aglivc, bglivcj, bglivcm, rleavc,
     &                  frootcj, frootcm, fbrchc, rlwodc, crootc)
          call wrtdeadc(time, curday, stdedc, metabc(1), strucc(1), 
     &                  wood1c, wood2c, wood3c)
          call wrtsoilc(time, curday, metabc(2), strucc(2), som1c(1), 
     &                  som1c(2), som2c(1), som2c(2), som3c)
          call wrtsysc(time, curday,
     &                 aglivc + bglivcj + bglivcm + rleavc + frootcj +
     &                 frootcm + fbrchc + rlwodc + crootc,
     &                 stdedc + metabc(1) + strucc(1) + wood1c +
     &                 wood2c + wood3c,
     &                 metabc(2) + strucc(2) + som1c(1) + som1c(2) +
     &                 som2c(1) + som2c(2) + som3c,
     &                 aglivc + bglivcj + bglivcm + rleavc + frootcj +
     &                 frootcm +  fbrchc + rlwodc + crootc + stdedc +
     &                 metabc(1) + strucc(1) + wood1c + wood2c +
     &                 wood3c + metabc(2) + strucc(2) + som1c(1) +
     &                 som1c(2) + som2c(1) + som2c(2) + som3c,
     &                 (st1c2(UNLABL) + st1c2(LABELD) + st2c2(UNLABL) +
     &                  st2c2(LABELD) + mt1c2(UNLABL) + mt1c2(LABELD) +
     &                  mt2c2(UNLABL) + mt2c2(LABELD) + s11c2(UNLABL) +
     &                  s11c2(LABELD) + s12c2(UNLABL) + s12c2(LABELD) +
     &                  s21c2(UNLABL) + s21c2(LABELD) + s22c2(UNLABL) +
     &                  s22c2(LABELD) + s3c2(UNLABL)  + s3c2(LABELD)  +
     &                  wd1c2(UNLABL) + wd1c2(LABELD) + wd2c2(UNLABL) +
     &                  wd2c2(LABELD) + wd3c2(UNLABL) +
     &                  wd3c2(LABELD)) - CO2resp)
          CO2resp = st1c2(UNLABL) + st1c2(LABELD) + st2c2(UNLABL) +
     &              st2c2(LABELD) + mt1c2(UNLABL) + mt1c2(LABELD) +
     &              mt2c2(UNLABL) + mt2c2(LABELD) + s11c2(UNLABL) +
     &              s11c2(LABELD) + s12c2(UNLABL) + s12c2(LABELD) +
     &              s21c2(UNLABL) + s21c2(LABELD) + s22c2(UNLABL) +
     &              s22c2(LABELD) + s3c2(UNLABL)  + s3c2(LABELD)  +
     &              wd1c2(UNLABL) + wd1c2(LABELD) + wd2c2(UNLABL) +
     &              wd2c2(LABELD) + wd3c2(UNLABL) + wd3c2(LABELD) +
     &              st1uvc2(UNLABL) + st1uvc2(LABELD) +
     &              stduvc2(UNLABL) + stduvc2(LABELD)

c ....... Calculate fraction of crop/forest carbon that is labeled
          if (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD) .gt.
     &        0.0) then
            crplblstg = carbostg(CRPSYS,LABELD) /
     &        (carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD))
          else
            crplblstg = cisofr
          endif
          if (carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD) .gt.
     &        0.0) then
            forlblstg = carbostg(FORSYS,LABELD) /
     &        (carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD))
          else
            forlblstg = cisotf
          endif

c ....... Calculate the daily delta 13C/14C values for ouput
          if (labtyp .eq. 2) then
c ......... Delta 13C output
            if ((newoiresp(LABELD) + newoiresp(UNLABL)) .gt. 0.0) then
              ddeloi = del13out(newoiresp(LABELD), newoiresp(UNLABL),
     &                          ddeloi)
            endif
            if ((newoeresp(LABELD) + newoeresp(UNLABL)) .gt. 0.0) then
              ddeloe = del13out(newoeresp(LABELD), newoeresp(UNLABL),
     &                          ddeloe)
            endif
            if ((newslitrsp(LABELD) + newslitrsp(UNLABL)) .gt. 0.0) then
              ddsrfclit = del13out(newslitrsp(LABELD),
     &                             newslitrsp(UNLABL), ddsrfclit)
            endif
            if ((newmnrlrsp(LABELD) + newmnrlrsp(UNLABL)) .gt. 0.0) then
              ddsmnrl = del13out(newmnrlrsp(LABELD),
     &                           newmnrlrsp(UNLABL), ddsmnrl)
            endif
            if ((newhetresp(LABELD) + newhetresp(UNLABL)) .gt. 0.0) then
              ddhetresp = del13out(newhetresp(LABELD),
     &                             newhetresp(UNLABL), ddhetresp)
            endif
            if ((newhetresp(LABELD) + newhetresp(UNLABL) + dcrtresp +
     &           dfrtresp) .gt. 0.0) then
              ddsoilresp = del13out((dcrtresp * crplblstg) +
     &                              (dfrtresp * forlblstg) +
     &                              newhetresp(LABELD),
     &                              (dcrtresp * (1.0 - crplblstg)) +
     &                              (dfrtresp * (1.0 - forlblstg)) +
     &                              newhetresp(UNLABL), ddsoilresp)
            endif
            if (dcmresp .gt. 0.0) then
              ddcmresp = del13out((dcmresp * crplblstg),
     &                            (dcmresp * (1.0 - crplblstg)),
     &                            ddcmresp)
            endif
            if (dfmresp .gt. 0.0) then
              ddfmresp = del13out((dfmresp * forlblstg),
     &                            (dfmresp * (1.0 - forlblstg)),
     &                            ddfmresp)
            endif
            if (dcgresp .gt. 0.0) then
              ddcgresp = del13out((dcgresp * crplblstg),
     &                            (dcgresp * (1.0 - crplblstg)),
     &                            ddcgresp)
            endif
            if (dfgresp .gt. 0.0) then
              ddfgresp = del13out((dfgresp * forlblstg),
     &                            (dfgresp * (1.0 - forlblstg)),
     &                            ddfgresp)
            endif
            if ((carbostg(CRPSYS,LABELD)+carbostg(CRPSYS,UNLABL)) .gt.
     &          0.0) then
              ddccarbostg = del13out(carbostg(CRPSYS,LABELD),
     &                               carbostg(CRPSYS,UNLABL),
     &                               ddccarbostg)
            endif
            if ((carbostg(FORSYS,LABELD)+carbostg(FORSYS,UNLABL)) .gt.
     &          0.0) then
              ddfcarbostg = del13out(carbostg(FORSYS,LABELD),
     &                               carbostg(FORSYS,UNLABL),
     &                               ddfcarbostg)
            endif
            call wrtdels(time, curday, ddeloi, ddeloe, ddsrfclit,
     &                   ddsmnrl, ddhetresp, ddsoilresp, ddcmresp,
     &                   ddfmresp, ddcgresp, ddfgresp, ddccarbostg,
     &                   ddfcarbostg)
          elseif (labtyp .eq. 1) then
c ......... Delta 14C output
            if ((newoiresp(LABELD) + newoiresp(UNLABL)) .gt. 0.0) then
              ddeloi = del14out(newoiresp(LABELD), newoiresp(UNLABL),
     &                          ddeloi)
            endif
            if ((newoeresp(LABELD) + newoeresp(UNLABL)) .gt. 0.0) then
              ddeloe = del14out(newoeresp(LABELD), newoeresp(UNLABL),
     &                          ddeloe)
            endif
            if ((newslitrsp(LABELD) + newslitrsp(UNLABL)) .gt. 0.0) then
              ddsrfclit = del14out(newslitrsp(LABELD),
     &                             newslitrsp(UNLABL), ddsrfclit)
            endif
            if ((newmnrlrsp(LABELD) + newmnrlrsp(UNLABL)) .gt. 0.0) then
              ddsmnrl = del14out(newmnrlrsp(LABELD),
     &                           newmnrlrsp(UNLABL), ddsmnrl)
            endif
            if ((newhetresp(LABELD) + newhetresp(UNLABL)) .gt. 0.0) then
              ddhetresp = del14out(newhetresp(LABELD),
     &                             newhetresp(UNLABL), ddhetresp)
            endif
            if ((newhetresp(LABELD) + newhetresp(UNLABL) + dcrtresp +
     &           dfrtresp) .gt. 0.0) then
              ddsoilresp = del14out((dcrtresp * crplblstg) +
     &                              (dfrtresp * forlblstg) +
     &                              newhetresp(LABELD),
     &                              (dcrtresp * (1.0 - crplblstg)) +
     &                              (dfrtresp * (1.0 - forlblstg)) +
     &                              newhetresp(UNLABL), ddsoilresp)
            endif
            if (dcmresp .gt. 0.0) then
              ddcmresp = del14out((dcmresp * crplblstg),
     &                            (dcmresp * (1.0 - crplblstg)),
     &                            ddcmresp)
            endif
            if (dfmresp .gt. 0.0) then
              ddfmresp = del14out((dfmresp * forlblstg),
     &                            (dfmresp * (1.0 - forlblstg)),
     &                            ddfmresp)
            endif
            if (dcgresp .gt. 0.0) then
              ddcgresp = del14out((dcgresp * crplblstg),
     &                            (dcgresp * (1.0 - crplblstg)),
     &                            ddcgresp)
            endif
            if (dfgresp .gt. 0.0) then
              ddfgresp = del14out((dfgresp * forlblstg),
     &                            (dfgresp * (1.0 - forlblstg)),
     &                            ddfgresp)
            endif
            if ((carbostg(CRPSYS,LABELD)+carbostg(CRPSYS,UNLABL)) .gt.
     &          0.0) then
              ddccarbostg = del14out(carbostg(CRPSYS,LABELD),
     &                               carbostg(CRPSYS,UNLABL),
     &                               ddccarbostg)
            endif
            if ((carbostg(FORSYS,LABELD)+carbostg(FORSYS,UNLABL)) .gt.
     &          0.0) then
              ddfcarbostg = del14out(carbostg(FORSYS,LABELD),
     &                               carbostg(FORSYS,UNLABL),
     &                               ddfcarbostg)
            endif
            call wrtdels(time, curday, ddeloi, ddeloe, ddsrfclit,
     &                   ddsmnrl, ddhetresp, ddsoilresp, ddcmresp,
     &                   ddfmresp, ddcgresp, ddfgresp, ddccarbostg,
     &                   ddfcarbostg)
          endif
        endif

c ..... Write to dc_sip.out
        if (time .ge. strplt) then
          npptot = mcprd(1) + mcprd(2) + mcprd(3) + mfprd(1) +
     &             mfprd(2) + mfprd(3) + mfprd(4) + mfprd(5) + mfprd(6)
          nee = npptot - dhresp
          tlai = (rleavc * 2.5) * btolai
          systemc = aglivc + bglivcj + bglivcm + rleavc + frootcj +
     &              frootcm + fbrchc + rlwodc + crootc + stdedc +
     &              wood1c + wood2c + wood3c + strucc(1) + metabc(1) +
     &              strucc(2) + metabc(2) + som1c(1) + som1c(2) +
     &              som2c(1) + som2c(2) + som3c + 
     &              carbostg(CRPSYS,UNLABL) + carbostg(CRPSYS,LABELD) +
     &              carbostg(FORSYS,UNLABL) + carbostg(FORSYS,LABELD)
c ....... Since each layer of the soil profile for this exercise has
c ....... the same texture use the field capacity and wilting point of
c ....... the top Century soil layer to compute the water holding
c ....... capacity for the various soil widths
          wc_2cm = (afiel(1) - awilt(1)) * 2
          wc_3cm = (afiel(1) - awilt(1)) * 3
          wc_5cm = (afiel(1) - awilt(1)) * 5
          wc_10cm = (afiel(1) - awilt(1)) * 10
          wc_15cm = (afiel(1) - awilt(1)) * 15
          wc_30cm = (afiel(1) - awilt(1)) * 30
          call wrtdcsip(time, curday, trandly, evapdly, intrcpt,
     &                  sublim, outflow, runoffdly, ppt(curday),
     &                  saccum, melt, snow, snlq, petdly, stemp,
     &                  wc_2cm, wc_3cm, wc_5cm, wc_10cm, wc_15cm,
     &                  wc_30cm, dhresp, mcprd(1), mcprd(2), mcprd(3),
     &                  mfprd(1), mfprd(2), mfprd(6), mfprd(3),
     &                  mfprd(4), mfprd(5), npptot, nee, aglivc,
     &                  bglivcj, bglivcm, rleavc, frootcj, frootcm,
     &                  fbrchc, rlwodc, crootc, tlai, stdedc, wood1c,
     &                  wood2c, wood3c, strucc(1), metabc(1), 
     &                  strucc(2), metabc(2), som1c(1), som1c(2),
     &                  som2c(1), som2c(2), som3c, systemc)

c ..... CH4 oxidation is converted from g/ha to g/m^2
        call wrtmethane(time, curday, aglivc, bglivcj, bglivcm,
     &                  prev_mcprd1, prev_mcprd2, prev_mcprd3,
     &                  Com, ppt(curday) - watr2sat - irractwk,
     &                  irractwk, watr2sat, avgst_10cm, TI, SI, Cr, Eh,
     &                  Feh, CH4_prod, CH4_Ep, CH4_Ebl,
     &                  CH4_oxid/10000.0)
        endif
200   continue
c ... END DAILY LOOP

c ... If we are at the end of the year reset the curday back to 1
      if (month .eq. 12) then
        curday = 1
      endif

c ... Annual co2 accumulators (10/92)
      ast1c2 = ast1c2 + st1c2(UNLABL) + st1c2(LABELD)
      ast2c2 = ast2c2 + st2c2(UNLABL) + st2c2(LABELD)
      amt1c2 = amt1c2 + mt1c2(UNLABL) + mt1c2(LABELD)
      amt2c2 = amt2c2 + mt2c2(UNLABL) + mt2c2(LABELD)
      as11c2 = as11c2 + s11c2(UNLABL) + s11c2(LABELD)
      as12c2 = as12c2 + s12c2(UNLABL) + s12c2(LABELD)
      as21c2 = as21c2 + s21c2(UNLABL) + s21c2(LABELD)
      as22c2 = as22c2 + s22c2(UNLABL) + s22c2(LABELD)
      as3c2  = as3c2  + s3c2(UNLABL)  + s3c2(LABELD)
      ast1uvc2 = ast1uvc2 + st1uvc2(UNLABL) + st1uvc2(LABELD)
      astduvc2 = astduvc2 + stduvc2(UNLABL) + stduvc2(LABELD)

c ... Monthly output averaged for *.bin file
      stemp = stempmth/real(dysimo(month))
      agdefacm(month) = agdefacsum/real(dysimo(month))
      bgdefacm(month) = bgdefacsum/real(dysimo(month))
      agdefac = agdefacm(month)
      bgdefac = bgdefacm(month)

      htran(month) = tran
      hpttr(month) = pttr

c ... Annual production accumulator
      cproda = cproda + cprodc + cprodf

c ... Net Mineralization
      do 100 iel = 1, nelem

c ..... Net mineralization for the mineralizing compartments
c ..... The structural component of litter and the wood compartments
c ..... are not mineralizers.  They should not be added into cmn or
c ..... sumnrs.
        cmn = metmnr(SRFC,iel) + metmnr(SOIL,iel) +
     &        s1mnr(SRFC,iel) + s1mnr(SOIL,iel) +
     &        s2mnr(SRFC,iel) + s2mnr(SOIL,iel) + s3mnr(iel)
        sumnrs(iel) = sumnrs(iel) + cmn

c ..... soilnm is net mineralization in the soil.
        soilnm(iel) = soilnm(iel) + s1mnr(SOIL,iel) +
     &                s2mnr(SOIL,iel) + s3mnr(iel) +
     &                metmnr(SOIL,iel) + strmnr(SOIL,iel) + w3mnr(iel)

c ..... Total net mineralization
        tnetmn(iel) = tnetmn(iel) + cmn + 
     &                strmnr(SRFC,iel) + strmnr(SOIL,iel) +
     &                w1mnr(iel) + w2mnr(iel) + w3mnr(iel)
100   continue

c ... Add calculation for annet which is used in the N deposition equations
c ... in eachyr, cak - 06/25/02
c ... Compute annual actual evapotranspiration
      annet = annet + evap + tran

c ... Stream flow accumulators, cak - 04/08/03
      do 105 ii = 1, 8
        strmac(ii) = strmac(ii) + stream(ii)
105   continue

c ... Calculate monthly respiration from decomposition for output
      if (month .eq. 1) then
        respmth(1) = resp(1)
        respmth(2) = resp(2)
        respsum(1) = respmth(1)
        respsum(2) = respmth(2)
      else
        respmth(1) = resp(1) - respsum(1)
        respmth(2) = resp(2) - respsum(2)
        respsum(1) = respsum(1) + respmth(1)
        respsum(2) = respsum(2) + respmth(2)
      endif          

c ... Calculate monthly autotrophic respiration for output
      if (month .eq. 1) then
        arspmth(1,1) = cautoresp(1)
        arspmth(1,2) = cautoresp(2)
        arspmth(2,1) = fautoresp(1)
        arspmth(2,2) = fautoresp(2)
        arspsum(1,1) = arspmth(1,1)
        arspsum(1,2) = arspmth(1,2)
        arspsum(2,1) = arspmth(2,1)
        arspsum(2,2) = arspmth(2,2)
      else
        arspmth(1,1) = cautoresp(1) - arspsum(1,1)
        arspmth(1,2) = cautoresp(2) - arspsum(1,2)
        arspmth(2,1) = fautoresp(1) - arspsum(2,1)
        arspmth(2,2) = fautoresp(2) - arspsum(2,2)
        arspsum(1,1) = arspsum(1,1) + arspmth(1,1)
        arspsum(1,2) = arspsum(1,2) + arspmth(1,2)
        arspsum(2,1) = arspsum(2,1) + arspmth(2,1)
        arspsum(2,2) = arspsum(2,2) + arspmth(2,2)
      endif          

c ... Compute output variables for printing or plotting.
      call savarp

      return
      end
