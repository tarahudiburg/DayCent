
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine dailymoist(curday, agdefacsum, bgdefacsum, bgwfunc,
     &                      frlech, co2val, accum, evapdly, intrcpt,
     &                      melt, outflow, petdly, runoffdly, sublim,
     &                      trandly, avgstemp, scenfrac, rpeff, 
     &                      watr2sat, prev_bgprd, Com, avgst_10cm, TI,
     &                      SI, Cr, Eh, Feh, CH4_prod, CH4_Ep, CH4_Ebl,
     &                      CH4_oxid, aetdly)

      implicit none
      include 'cflows.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'doubles.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'jday.inc'
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
      include 'timvar.inc'
      include 't0par.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

c ... FORMAL PARAMETERS
      integer curday
      real    bgwfunc, agdefacsum, bgdefacsum
      real    frlech(MAXIEL)
      real    co2val
      real    accum, evapdly, intrcpt, melt, outflow, petdly,
     &        runoffdly, sublim, trandly, aetdly
      real    avgstemp, scenfrac, rpeff
      real    watr2sat, prev_bgprd, Com, avgst_10cm, TI, SI, Cr, Eh,
     &        Feh, CH4_prod, CH4_Ep, CH4_Ebl

c ... This routine loops calls the water budget routine (h2oflux), the
c ... decomposition routine (decomp) and the trace gas routines at a
c ... daily timestep.
c ...
c ... wfluxout[] - total net flux thru the bottom of layer each day (cm H2O)
c ...              (positive is downward, negative is upward)
c ... nitrate[]  - layers of the nitrogen pool (gN/m2)
c ... ammonium   - the ammonium pool, no layer structure (gN/m2)
c ... newminrl   - mineralization that has occurred in the current day (gN/m2)
c ... co2val     - CO2 effect on transpiration.   Added 8/14/98 -mdh
c ... dstr       - Soil structural CO2 resp for the day (gC/m2. Added 11/9/2014 -mdh
c ... efscltef   - Effective decomposition effect. Added 11/9/2014 -mdh
c ...              Ratio of the soil decomposition CO2 flow sum (sdco2sum) 
c ...              to the no-till soil decomposition CO2 flow sum (ntdco2sm)

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE bal_npool(nlayer, minerl, ammonium, nitrate, 
     &                       inorglch)
          !MS$ATTRIBUTES ALIAS:'_bal_npool' :: bal_npool
          INTEGER          nlayer
          REAL             minerl(*)
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          DOUBLE PRECISION inorglch
        END SUBROUTINE bal_npool

        SUBROUTINE calcdefac(texture, tfunc, bgwfunc, agdefac, bgdefac,
     &                       avgwfps, teff, rprpet, idef, ppt, snow,
     &                       avgstemp)
          !MS$ATTRIBUTES ALIAS:'_calcdefac' :: calcdefac
          INTEGER texture
          REAL    tfunc
          REAL    bgwfunc
          REAL    agdefac
          REAL    bgdefac
          REAL    avgwfps
          REAL    teff(*)
          REAL    rprpet
          INTEGER idef
          REAL    ppt
          REAL    snow
          REAL    avgstemp
        END SUBROUTINE calcdefac

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

        SUBROUTINE trace_gas_model(newminrl, ammonium, nitrate,
     &                             texture, sand, silt, clay, afiel,
     &                             bulkd, maxt, ppt, snow, avgwfps,
     &                             stormf, basef, frlechd, stream,
     &                             inorglch, critflow, wfluxout,
     &                             newCO2, co2_conc, efscltef, time,
     &                             NOflux, Nn2oflux, Dn2oflux, Dn2flux,
     &                             CH4_oxid, isdecid, isagri, aglivc,
     &                             rleavc, btolai, crpstore, forstore,
     &                             nit_amt, nreduce, curday, pHscale,
     &                             dN2lyr, dN2Olyr, prev_bgprd, Com,
     &                             avgst_10cm, TI, SI, Cr, Eh, Feh,
     &                             CH4_prod, CH4_Ep, CH4_Ebl,
     &                             watertable, watrflag, bglivc,
     &                             tmxbio)
          !MS$ATTRIBUTES ALIAS:'_trace_gas_model' :: trace_gas_model
          DOUBLE PRECISION newminrl
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          INTEGER          texture
          REAL             sand
          REAL             silt
          REAL             clay
          REAL             afiel
          REAL             bulkd
          REAL             maxt
          REAL             ppt
          REAL             snow
          REAL             avgwfps
          REAL             stormf
          REAL             basef
          REAL             frlechd(*)
          REAL             stream(*)
          DOUBLE PRECISION inorglch
          REAL             critflow
          REAL             wfluxout(*)
          REAL             newCO2
          REAL             co2_conc(*)
          REAL             efscltef
          REAL             time
          DOUBLE PRECISION NOflux
          DOUBLE PRECISION Nn2oflux
          DOUBLE PRECISION Dn2oflux
          DOUBLE PRECISION Dn2flux
          DOUBLE PRECISION CH4_oxid
          INTEGER          isdecid
          INTEGER          isagri
          REAL             aglivc
          REAL             rleavc
          REAL             btolai
          REAL             crpstore
          REAL             forstore
          DOUBLE PRECISION nit_amt
          REAL             nreduce
          INTEGER          curday
          REAL             pHscale
          DOUBLE PRECISION dN2lyr(*)
          DOUBLE PRECISION dN2Olyr(*)
          REAL             prev_bgprd
          REAL             Com
          REAL             avgst_10cm
          REAL             TI
          REAL             SI
          REAL             Cr
          REAL             Eh
          REAL             Feh
          REAL             CH4_prod
          REAL             CH4_Ep
          REAL             CH4_Ebl
          INTEGER          watertable
          INTEGER          watrflag
          REAL             bglivc
          REAL             tmxbio
        END SUBROUTINE trace_gas_model

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

        SUBROUTINE watrbal(curday, time, ppt, accum, melt, wbswc1, 
     &                     wbswc2, evapdly, trandly, sublim, intrcpt,
     &                     outflow, snlq1, snlq2, snow, runoffdly)
          !MS$ATTRIBUTES ALIAS:'_watrbal' :: watrbal
          INTEGER curday
          REAL    time
          REAL    ppt
          REAL    accum
          REAL    melt
          REAL    wbswc1
          REAL    wbswc2
          REAL    evapdly
          REAL    trandly
          REAL    sublim
          REAL    intrcpt
          REAL    outflow
          REAL    snlq1
          REAL    snlq2
          REAL    snow
          REAL    runoffdly
        END SUBROUTINE watrbal

        SUBROUTINE watrflow(curday, month, nlayer, nlaypg, watertable,
     &                      watrflag, avgtemp, tempmin, tempmax,
     &                      solrad, rhumid, windsp, ppt, aglivb,
     &                      sfclit, stdead, rwcf, avh2o, asmos, snow,
     &                      snlq, amovdly, petdly, evapdly, trandly,
     &                      stream1, basef, pottransp, baseflow, accum,
     &                      melt, intrcpt, outflow, tmelt, sublim,
     &                      wfluxout, time, strplt, co2val, tmns, tmxs,
     &                      runoffdly, trandep, soiltavewk, daylength,
     &                      woodb, elitst, pmxtmp, pmntmp, pmxbio,
     &                      stemp, stsys, ststart, stamt, litrcrbn,
     &                      aggreenc, watr2sat, aetdly, h2ogef3,
     &                      stcrlai, srad)
          !MS$ATTRIBUTES ALIAS:'_watrflow' :: watrflow
          INTEGER          curday
          INTEGER          month
          INTEGER          nlayer
          INTEGER          nlaypg
          INTEGER          watertable
          INTEGER          watrflag
          REAL             avgtemp
          REAL             tempmin
          REAL             tempmax
          REAL             solrad
          REAL             rhumid
          REAL             windsp
          REAL             ppt
          REAL             aglivb
          REAL             sfclit
          REAL             stdead
          REAL             rwcf(*)
          REAL             avh2o(*)
          REAL             asmos(*)
          REAL             snow
          REAL             snlq
          REAL             amovdly(*)
          REAL             petdly
          REAL             evapdly
          REAL             trandly
          REAL             stream1
          REAL             basef
          REAL             pottransp
          REAL             baseflow
          REAL             accum
          REAL             melt
          REAL             intrcpt
          REAL             outflow
          REAL             tmelt(*)
          REAL             sublim
          REAL             wfluxout(*)
          REAL             time
          REAL             strplt
          REAL             co2val
          REAL             tmns
          REAL             tmxs
          REAL             runoffdly
          REAL             trandep
          REAL             soiltavewk
          REAL             daylength
          REAL             woodb
          REAL             elitst
          REAL             pmxtmp
          REAL             pmntmp
          REAL             pmxbio
          REAL             stemp
          REAL             stsys
          REAL             ststart
          REAL             stamt
          REAL             litrcrbn
          REAL             aggreenc
          REAL             watr2sat
          REAL             aetdly
          REAL             h2ogef3
          REAL             stcrlai
          DOUBLE PRECISION srad
        END SUBROUTINE watrflow

        SUBROUTINE wrtcflows(time, jday, som11tosom21, som12tosom22,
     &                       som12tosom3, som21tosom11, som21tosom22,
     &                       som22tosom12, som22tosom3, som3tosom12,
     &                       metc1tosom11, metc2tosom12, struc1tosom11,
     &                       struc1tosom21, struc2tosom12,
     &                       struc2tosom22, wood1tosom11, wood1tosom21,
     &                       wood2tosom11, wood2tosom21, wood3tosom12,
     &                       wood3tosom22)
          !MS$ATTRIBUTES ALIAS:'_wrtcflows' :: wrtcflows
          REAL    time
          INTEGER jday
          REAL    som11tosom21
          REAL    som12tosom22
          REAL    som12tosom3
          REAL    som21tosom11
          REAL    som21tosom22
          REAL    som22tosom12
          REAL    som22tosom3
          REAL    som3tosom12
          REAL    metc1tosom11
          REAL    metc2tosom12
          REAL    struc1tosom11
          REAL    struc1tosom21
          REAL    struc2tosom12
          REAL    struc2tosom22
          REAL    wood1tosom11
          REAL    wood1tosom21
          REAL    wood2tosom11
          REAL    wood2tosom21
          REAL    wood3tosom12
          REAL    wood3tosom22
        END SUBROUTINE wrtcflows

        SUBROUTINE wrtco2(time, curday, co2_conc)
          !MS$ATTRIBUTES ALIAS:'_wrtco2' :: wrtco2
          REAL    time
          INTEGER curday
          REAL    co2_conc(*)
        END SUBROUTINE wrtco2

        SUBROUTINE wrtdaily(time, curday, petdly, agdefac, bgdefac,
     &                      stemp, snow, snlq, thermunits, aglivc,
     &                      aggreenc, hwstress, scenfrac, srad)
          !MS$ATTRIBUTES ALIAS:'_wrtdaily' :: wrtdaily
          REAL             time
          INTEGER          curday
          REAL             petdly
          REAL             agdefac
          REAL             bgdefac
          REAL             stemp
          REAL             snow
          REAL             snlq
          REAL             thermunits
          REAL             aglivc
          REAL             aggreenc
          REAL             hwstress
          REAL             scenfrac
          DOUBLE PRECISION srad
        END SUBROUTINE wrtdaily

        SUBROUTINE wrtdn2lyr(time, jday, dN2lyr)
          !MS$ATTRIBUTES ALIAS:'_wrtdn2lyr' :: wrtdn2lyr
          REAL             time
          INTEGER          jday
          DOUBLE PRECISION dn2lyr(*)
        END SUBROUTINE wrtdN2lyr

        SUBROUTINE wrtdn2olyr(time, jday, dN2Olyr)
          !MS$ATTRIBUTES ALIAS:'_wrtdn2olyr' :: wrtdn2olyr
          REAL             time
          INTEGER          jday
          DOUBLE PRECISION dN2Olyr(*)
        END SUBROUTINE wrtdn2olyr

        SUBROUTINE wrtnflux(time, curday, Nn2oflux, Dn2oflux, Dn2flux,
     &                      NOflux, nflux_sum, nflux_sum2)
          !MS$ATTRIBUTES ALIAS:'_wrtnflux' :: wrtnflux
          REAL             time
          INTEGER          curday
          DOUBLE PRECISION Nn2oflux
          DOUBLE PRECISION Dn2oflux
          DOUBLE PRECISION Dn2flux
          DOUBLE PRECISION NOflux
          DOUBLE PRECISION nflux_sum
          DOUBLE PRECISION nflux_sum2
        END SUBROUTINE wrtnflux

        SUBROUTINE wrtsoiln(time, curday, ammonium, nitrate)
          !MS$ATTRIBUTES ALIAS:'_wrtsoiln' :: wrtsoiln
          REAL             time
          INTEGER          curday
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
        END SUBROUTINE wrtsoiln

        SUBROUTINE wrtsummary(time, curday, tempmax, tempmin, ppt,
     &                        N2Oflux, NOflux, CH4_oxid, nit_amt,
     &                        CO2resp, ntdco2sm)
          !MS$ATTRIBUTES ALIAS:'_wrtsummary' :: wrtsummary
          REAL             time
          INTEGER          curday
          REAL             tempmax
          REAL             tempmin
          REAL             ppt
          DOUBLE PRECISION N2Oflux
          DOUBLE PRECISION NOflux
          DOUBLE PRECISION CH4_oxid
          DOUBLE PRECISION nit_amt
          DOUBLE PRECISION CO2resp
          REAL             ntdco2sm
        END SUBROUTINE wrtsummary

        SUBROUTINE wrtwflux(time, curday, wfluxout)
          !MS$ATTRIBUTES ALIAS:'_wrtwflux' :: wrtwflux
          REAL    time
          INTEGER curday
          REAL    wfluxout(*)
        END SUBROUTINE wrtwflux

      END INTERFACE

c ... LOCAL VARIABLES
      integer          iel, ilyr, kts, ii, clyr
      integer          isdecid, isagri
      real             fsol 
      real             stream1
      real             amovdly(10)
      real             rprpet
      real             pottransp
      real             baseflow
      real             wbswc1, wbswc2
      real             snlq1, snlq2
      real             tfunc
      real             wfluxout(SWMAXLYR)
      real             co2_conc(SWMAXLYR)
      real             newCO2
      real             soilresp1, soilresp2
      double precision newminrl, inorglch
      double precision nflux_sum, SCALE, nflux_sum2
      double precision Nn2oflux, Dn2oflux, Dn2flux, NOflux
      real             critflow
      real             frlechd(MAXIEL)
      real             avgwfps
      real             tmns, tmxs
      double precision CH4_oxid, nit_amt
      character        subname*10
      real             croplai, treelai, totlai
      real             trandep
      real             CO2resp
      real             maxrwcf
      real             hwstress
      real             dayunits
      real             thermtemp, tmns_mlt, tmxs_mlt
      real             litrcrbn, aggreenc
      real             sradKJ, soilsrad
      double precision dN2lyr(SWMAXLYR), dN2Olyr(SWMAXLYR)
      real             dstr, efscltef 
      real             Com1, Com2

      data isagri /0/

      save isagri
      save nflux_sum, nflux_sum2

c ... FUNCTION DECLARATIONS
      real      anerob, fsfunc, line, ramp
      external  anerob, fsfunc, line, ramp

      subname = 'dailymoist'

c      write(*,*) 'ENTERING FUNCTION DAILYMOIST...'
c ... decodt is equivalent to 4 times a day (vs. 4 times a month)
c ... -mdh 8/31/94  (Formerly, decodt=dt/ntspm, monthly version)

      decodt = dt/real(dysimo(month)*ntspm)

      dstr = 0 
      efscltef = 1.0

c ... ***********************************************************************

c ... Soil moisture

c ... This code was pulled from subroutine cycle. -mdh  8/31/94.

      cwstress = 0.0
      twstress = 0.0
      hwstress = 0.0

      if (curday .eq. 1) then
        nflux_sum = 0.0
        nflux_sum2 = 0.0
      endif
      snlq1 = snlq
c ... Amount of water in the soil at the beginning of the day
      wbswc1 = 0.0
      do 106 ilyr = 1,nlayer+1
        wbswc1 = wbswc1 + asmos(ilyr)
106   continue
      newminrl = 0.0

c ... Initialize output variables that are tracking the daily carbon
c ... flows due to decomposition.
      metc1tosom11  = 0.0
      metc2tosom12  = 0.0
      struc1tosom11 = 0.0
      struc1tosom21 = 0.0
      struc2tosom12 = 0.0
      struc2tosom22 = 0.0
      som11tosom21  = 0.0
      som12tosom22  = 0.0
      som12tosom3   = 0.0
      som21tosom11  = 0.0
      som21tosom22  = 0.0
      som22tosom12  = 0.0
      som22tosom3   = 0.0
      som3tosom12   = 0.0
      wood1tosom11  = 0.0
      wood1tosom21  = 0.0
      wood2tosom11  = 0.0
      wood2tosom21  = 0.0
      wood3tosom12  = 0.0
      wood3tosom22  = 0.0

c ... Fertilization option
c ... This code has been relocated from simsom so that fertilization
c ... can occur on a specific day, cak - 04/17/03
c ... Add an optional multiplier on feramt for N, cak - 04/05/04
      if (dofert .and. curday .eq. fertday) then
        do 60 iel = 1, nelem
          if (iel .eq. N) then
            clyr = SRFC
            if (Ninput .eq. 1 .or. Ninput .eq. 3) then
              esrsnk(iel) = esrsnk(iel) - feramt(iel)*Nscalar(month)
              call update_npool(clyr, feramt(iel)*Nscalar(month),
     &                          frac_nh4_fert, frac_no3_fert,
     &                          ammonium, nitrate, subname)
              minerl(SRFC,iel) = minerl(SRFC,iel) +
     &                           feramt(iel)*Nscalar(month)
              fertot(iel) = fertot(iel) +
     &                      feramt(iel)*Nscalar(month)
              fertac(iel) = fertac(iel) +
     &                      feramt(iel)*Nscalar(month)
              fertmth(month,iel) = fertmth(month,iel) +
     &                             feramt(iel)*Nscalar(month)
              gfertot(iel) = gfertot(iel) + feramt(iel)*Nscalar(month)
            else
              esrsnk(iel) = esrsnk(iel) - feramt(iel)
              call update_npool(clyr, feramt(iel),
     &                          frac_nh4_fert, frac_no3_fert,
     &                          ammonium, nitrate, subname)
              minerl(SRFC,iel) = minerl(SRFC,iel) + feramt(iel)
              fertot(iel) = fertot(iel) + feramt(iel)
              fertac(iel) = fertac(iel) + feramt(iel)
              fertmth(month,iel) = fertmth(month,iel) + feramt(iel)
              gfertot(iel) = gfertot(iel) + feramt(iel)
            endif
          else
            esrsnk(iel) = esrsnk(iel) - feramt(iel)
            minerl(SRFC,iel) = minerl(SRFC,iel) + feramt(iel)
            fertot(iel) = fertot(iel) + feramt(iel)
            fertac(iel) = fertac(iel) + feramt(iel)
            fertmth(month,iel) = fertmth(month,iel) + feramt(iel)
            gfertot(iel) = gfertot(iel) + feramt(iel)
          endif
60      continue
        fertcnt = 1
        nreduce = ninhib
      endif

c ... Accumulate thermal units for the growing degree day
c ... implementation, cak - 04/17/03
      if (accumdd) then
c ..... Use the day length to calculate the thermal units,
c ..... cak - 08/29/05
        if (daylength(curday) .lt. 12.0) then
          tmns_mlt = ((12.0 - daylength(curday)) * 3.0 + 12.0) / 24.0
        else
          tmns_mlt = ((12.0 - daylength(curday)) * 1.2 + 12.0) / 24.0
        endif
        tmns_mlt = min(0.95, tmns_mlt)
        tmns_mlt = max(0.05, tmns_mlt)
        tmxs_mlt = 1.0 - tmns_mlt
c ..... Set an upper limit on the calculation for accumulating growing
c ..... degree days such that when the maximum temperature for the day
c ..... is capped at basetemp(2) degrees C, cak - 05/21/2008
        thermtemp = tmxs_mlt * min(tempmax(curday), basetemp(2)) +
     &              tmns_mlt * max(tempmin(curday), basetemp(1))
        dayunits = max(0.0, thermtemp - basetemp(1))
        thermunits = thermunits + dayunits
c ..... Set the emergence flag to true when enough growing degree days
c ..... have occurred after planting (DDEMERG) for a grain filling
c ..... annual
        if ((thermunits .ge. ddemerg) .and. (frtcindx .ge. 4)) then
          emerg = .true.
        endif
c ..... For a grain producing crop reaching ddbase starts the grain
c ..... filling period, cak 06/02/05
c ..... Add code to simulate the senescence of annual plants, beginning
c ..... at anthesis (DDBASE), the amount of photosynthetic active carbon
c ..... (aggreenc) decreases with age, based on water stress, affecting
c ..... potential growth and transpiration, 06/10/2014
c ..... scenfrac is a 0.0 to 1.0 multiplier used to indicate the
c ..... fraction of the above live carbon that is photosynthetic active
c ..... carbon
c .....   1.0 = no senescence has occurred, 100% photosynthetic active
c .....         carbon
c .....   0.0 = full senescence, 0% photosynthetic active carbon
        if ((thermunits .ge. ddbase) .and. (frtcindx .ge. 4)) then
          if (.not. grnfill) then
            grnfill = .true.
            gwstress = 0.0
            grnfldys = 1
          endif
c ....... Check to see if the plant has reached maturity
          if (thermunits .ge. (ddbase + mxddhrv)) then
c ......... Stop plant growth
            crpgrw = 0
            cgrwdys = 0
            scenfrac = 0.0
          else
c ......... Use the grain water stress term to determine if we
c ......... have reached maturity, full senescence
            hwstress = ramp(gwstress/grnfldys, 0.0, mnddhrv,
     &                      1.0, mxddhrv)
            scenfrac = ramp(thermunits-ddbase, 0.0, 1.0, hwstress,
     &                      0.0)
          endif
          grnfldys = grnfldys + 1
        endif
      endif
      aggreenc = aglivc * scenfrac

      petann = petann + petdly

      do 110 ii = 1,10
        amovdly(ii) = 0.0
110   continue

c ... If the grass/crop plant functional type is annual calculate a
c ... dynamic value for claypg from 1 to claypg_const based on the
c ... number of days since planting and FRTC(3), claypg_const is the
c ... value of claypg as read for the current CROP option
      if ((cursys .eq. CRPSYS) .or. (cursys .eq. SAVSYS)) then
        if ((frtcindx .eq. 2) .or. (frtcindx .ge. 4)) then
          claypg = nint(ramp(real(plntcnt), 0.0, 1.0, frtc(3),
     &                       real(claypg_const)))
        endif
      endif

c ... Calculate a dynamic value for nlaypg based on the crop and/or tree
c ... option used, cak - 01/29/03
      if (cursys .eq. CRPSYS) then
        nlaypg = claypg
      else if (cursys .eq. SAVSYS) then
c ..... For crops and grasses a leaf area of 1 = 100 grams of biomass
        croplai = aglivc * 2.5 * 0.01
        treelai = rleavc * 2.5 * btolai
        totlai = croplai + treelai
        if (totlai .gt. 0.0) then
          nlaypg = nint(line(treelai/totlai, 0.0, real(claypg), 1.0,
     &                       real(tlaypg)))
        else
          nlaypg = min(claypg, tlaypg)
        endif
        if (nlaypg .lt. min(claypg, tlaypg)) then
          nlaypg = min(claypg, tlaypg)
        endif
        if (nlaypg .gt. max(claypg, tlaypg)) then
          nlaypg = max(claypg, tlaypg)
        endif
      else
c ..... This is a tree system
        nlaypg = tlaypg
      endif

c ... Calculate depth for transpiration, cak - 01/29/03
      trandep = 0
      do 120 ii = 1,nlaypg
        trandep = trandep + adep(ii)
120   continue

c ... Calculate biomass values to be used in soil surface temperature
c ... calculations, this code has been moved from the potprod
c ... subroutine, cak - 01/28/2010
      if (cursys .eq. FORSYS) then
c ..... Live biomass
        aglivb = rleavc * 2.5
c ..... Surface litter biomass
c ..... Second mod to remove effect of woodc -rm 1/91
        sfclit = (strucc(SRFC) + metabc(SRFC)) * 2.0
c ..... Standing dead biomass
        stdead = 0.0
c ..... Wood biomass
        woodb = (rlwodc + fbrchc) * 2.0
      elseif (cursys .eq. SAVSYS) then
c ..... Live biomass
        aglivb = (rleavc + aglivc) * 2.5
c ..... Surface litter biomass
        sfclit = (strucc(SRFC) + metabc(SRFC)) * 2.0
c ..... Standing dead biomass
        stdead = stdedc * 2.5
c ..... Wood biomass
        woodb = (rlwodc + fbrchc) * 2.0
      else
c ..... Live biomass
        aglivb = aglivc * 2.5
c ..... Surface litter biomass
        sfclit = (strucc(SRFC) + metabc(SRFC)) * 2.5
c ..... Standing dead biomass
        stdead = stdedc * 2.5
c ..... Wood biomass
        woodb = 0.0
      endif

      litrcrbn = som1c(SRFC) + som2c(SRFC) + metabc(SRFC) +
     &           strucc(SRFC)
      call watrflow(curday, month, nlayer, nlaypg, watertable,
     &              watrflag, avgtemp(curday), tempmin(curday),
     &              tempmax(curday), solrad(curday), rhumid(curday),
     &              windsp(curday), ppt(curday), aglivb, sfclit,
     &              stdead, rwcf, avh2o, asmos, snow, snlq, amovdly,
     &              petdly, evapdly, trandly, stream1, basef,
     &              pottransp, baseflow, accum, melt, intrcpt, outflow,
     &              tmelt, sublim, wfluxout, time, strplt, co2val,
     &              tmns, tmxs, runoffdly, trandep, soiltavewk, 
     &              daylength(curday), woodb, elitst, pmxtmp, pmntmp,
     &              pmxbio, stemp, stsys, ststart, stamt, litrcrbn,
     &              aggreenc, watr2sat, aetdly, h2ogef(3), stcrlai,
     &              srad(curday))
c ... Sum water that is added to keep the soil at saturation into the
c ... water input output variables
      ppt(curday) = ppt(curday) + watr2sat
      annppt = annppt + watr2sat

c ... Accumulate the value of stemp, returned via wateflow, for the
c ... month
      stempmth = stempmth + stemp

c ... Find the relative water content in the wettest soil layer to use in
c ... the calculation of water stress on potential growth, cak - 12/06/04
c ... Grass/crop water stress
      maxrwcf = -9999
      do 130 ii = 1, claypg
        if (rwcf(ii) .gt. maxrwcf) then
          maxrwcf = rwcf(ii)
        endif
130   continue
      cwstress = cwstress + min(1.0, maxrwcf)
c ... If this is a grain fill crop and we are within the grain
c ... filling period calculate the water stress term to be used in
c ... the harvest water stress calculation
      if (grnfill) then
        gwstress = gwstress + 1.0/
     &             (1.0 + exp(wscoeff(1,2)*
     &                        (wscoeff(1,1)-min(1.0,maxrwcf))))
      endif
c ... Tree water stress
      maxrwcf = -9999
      do 140 ii = 1, tlaypg
        if (rwcf(ii) .gt. maxrwcf) then
          maxrwcf = rwcf(ii)
        endif
140   continue
      twstress = twstress + min(1.0, maxrwcf)

c ... If there is snow melting into the ground use the melt value returned
c ... from the watrflow subroutine to determine how much snow is melting
c ... into the ground, cak - 10/21/02
c ... ppt(curday) includes any irrigation
      if (melt .gt. 0.0) then
c ..... melt represents the amount of water draining into the soil when
c ..... there is snow on the ground, both precipitation and irrigation
c ..... amounts have been taken into account in the snow calculations,
c ..... cak - 12/13/02
        rprpet = melt / petdly
      else
        rprpet = (ppt(curday) + avh2o(3))/ petdly
      endif

c ... Accumulate daily stream flow, drainage from each layer, pet, 
c ... evaporation, and transpiration by month
      stream(1) = stream(1) + stream1
      pet = pet + petdly
      evap = evap + evapdly + intrcpt + sublim
      tran = tran + trandly
      pttr = pttr + pottransp
      do 104 ilyr = 1,nlayer
        amov(ilyr) = amov(ilyr) + amovdly(ilyr)
104   continue
c ... Accumulate runoff for month
      runoff = runoff + runoffdly
c ... Accumulate precipitation + irrigation for the month
      pptmonth = pptmonth + ppt(curday)

c ... Calculate the effect impact of anerobic conditions on decomposition
c ... Last parameter = 1 if microcosm

      anerb = anerob(aneref,drain,rprpet,petdly,0)

c ... Combined effects of temperature and moisture on decomposition
      call calcdefac(texture, tfunc, bgwfunc, agdefac, bgdefac,
     &               avgwfps, teff, rprpet, idef, ppt(curday), snow,
     &               avgstemp)

c ... calculate defacm(month) in subroutine simsom. -mdh 10/94
      agdefacsum = agdefacsum + agdefac
      bgdefacsum = bgdefacsum + bgdefac

c ... *********************************************************************

c ... Decomposition Submodel

c ... Call decomp routines ntspm times per month.
c ... Removed the P and S chemistry from decomp and created the
c ... subroutine pschem.  -rm  6/91

c ... Change CO2 respiration value passed to the trace gas submodel
c ... so that we are passing only soil respiration, cak - 10/22/2007
      soilresp1 = mt2c2(UNLABL) + mt2c2(LABELD) +
     &            st2c2(UNLABL) + st2c2(LABELD) +
     &            s12c2(UNLABL) + s12c2(LABELD) +
     &            s22c2(UNLABL) + s22c2(LABELD) +
     &            s3c2(UNLABL)  + s3c2(LABELD)
c ... The methane production submodel requires CO2 respiration from all
c ... of the crop/grass litter and SOM pools
      Com1 = soilresp1 +
     &       mt1c2(UNLABL) + mt1c2(LABELD) +
     &       st1c2(UNLABL) + st1c2(LABELD) +
     &       s11c2(UNLABL) + s11c2(LABELD) +
     &       s21c2(UNLABL) + s21c2(LABELD)

c ... Scale pH values if necessary
      if (phsys .gt. 0) then
        ph = phstart * pHscalar(month)
      endif

c ... dstr = pre decomposition structural
      dstr = st2c2(UNLABL) + st2c2(LABELD) 
c ... sdco2sum = sum of soil CO2 flow with cultivation
      sdco2sum = 0.0    
c ... ntdco2sm = sum of soil CO2 flow without cultivation
      ntdco2sm = 0.0    

c ... Solar radiation in KJ
      sradKJ = srad(curday)*W2KJ

c ... Solar radiation at the bottom of the plant canopy
      soilsrad = sradKJ * exp(-0.5 * stcrlai)

c ... Photodecomposition, cak - 12/10/2009
      call photodecomp(sradKJ, soilsrad)

      do 40 kts = 1, ntspm
        call decomp(decodt,decsys,amovdly,newminrl,rpeff,soilsrad)
        if (nelem .ge. P) then
          call pschem(decodt)
        endif

c ..... Update decomposition and nitrogen fixation flows.
        call flowup(time)
        call flowup_double(time)
        call flowup_double_in(time)
        call flowup_double_out(time)
        call sumcar

c ..... Update the occlud and secndy single precision variables using
c ..... the values from the double precision variables, cak - 03/20/02
        occlud = real(occlud_double)
        secndy(1) = real(secndy_double(1))
        secndy(2) = real(secndy_double(2))
        secndy(3) = real(secndy_double(3))
     
c ..... aminrl contains the average amount of N, P, and S 
c ..... available in the top layer for the time period covered by 
c ..... dt/ntspm.  minerl contains the current value of mineral N,
c ..... P, and S by layer.

        do 30 iel = 1, nelem
          if (iel .eq. P) then
            fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
          else
            fsol = 1.0
          endif
          aminrl(iel) = (aminrl(iel) + minerl(SRFC,iel)*fsol)/2.0
30      continue
40    continue

c ... Change CO2 respiration value passed to the trace gas submodel
c ... so that we are passing only soil respiration, cak - 10/22/2007
      soilresp2 = mt2c2(UNLABL) + mt2c2(LABELD) +
     &            st2c2(UNLABL) + st2c2(LABELD) +
     &            s12c2(UNLABL) + s12c2(LABELD) +
     &            s22c2(UNLABL) + s22c2(LABELD) +
     &            s3c2(UNLABL)  + s3c2(LABELD)
c ... The methane production submodel requires CO2 respiration from all
c ... of the crop/grass litter and SOM pools
      Com2 = soilresp2 +
     &       mt1c2(UNLABL) + mt1c2(LABELD) +
     &       st1c2(UNLABL) + st1c2(LABELD) +
     &       s11c2(UNLABL) + s11c2(LABELD) +
     &       s21c2(UNLABL) + s21c2(LABELD)
      Com = Com2 - Com1
      newCO2 = soilresp2 - soilresp1
      if (newCO2 .le. 0.0) then
        newCO2 = 0.000001
      endif
      CO2resp = newCO2
      if (avgwfps .gt. 0.60) then
        newCO2 = newCO2 / bgwfunc
      endif

c ... Difference the structural accumulator to avoid having modify declig
c ... and tell the difference between soil litter and other pools.
      dstr = st2c2(UNLABL) + st2c2(LABELD) - dstr
c ... sdco2sum - soil decomposition CO2 flow sum. Needed to calculate the
c                effective decomposition effect.
c ... ntdco2sm - no-till soil decomposition CO2 flow sum
      sdco2sum = sdco2sum + dstr
      if(cltfac(4) .gt. 0.0) then
        ntdco2sm = ntdco2sm + dstr/cltfac(4)
      endif
c ... Bill Parton requested to restrict the CO2 respiration effect 
c ... on dentrification by capping ntdco2sm at 2.5 gC/m2/d. -mdh 11/12/2014
      ntdco2sm = min(2.5,ntdco2sm)

c ... ratio the total and no-till CO2 flows to get the weighted clteff
      if(efscltef .gt. 0) then
        efscltef = sdco2sum / ntdco2sm
      else
        efscltef = 1.0
      endif

c      critflow = minlch/real(dysimo(month))
      critflow = minlch

      do 45 iel = 1, nelem
c        frlechd(iel) = frlech(iel)/real(dysimo(month))
        frlechd(iel) = frlech(iel)
45    continue 

c ... *********************************************************************

c ... Trace Gas Model
   
c      write(*,*) 'Time = ', time, ' curday = ', curday

c ... Are we running a decidious forest?
      if ((cursys .eq. FORSYS) .and. (decid .eq. 1)) then
        isdecid = 1
      else
        isdecid = 0
      endif
c ... Once cultivation, fertilization, or harvesting occurs in a system the 
c ... agricultural effect needs to be "turned on".  Once this effect on
c ... methane oxidation has been invoked it stays in place.
      if (isagri .eq. 0) then
c        if (dofert) isagri = 1
        if (docult) isagri = 1
c        if (dohrvt) isagri = 1
      endif

      call trace_gas_model(newminrl, ammonium, nitrate, texture, 
     &                     sand, silt, clay, afiel(1), bulkd, maxt,
     &                     ppt(curday), snow, avgwfps, stormf, basef,
     &                     frlechd, stream, inorglch, critflow,
     &                     wfluxout, newCO2, co2_conc, efscltef, time,
     &                     NOflux, Nn2oflux, Dn2oflux, Dn2flux,
     &                     CH4_oxid, isdecid, isagri, aglivc, rleavc,
     &                     btolai, crpstg(N), forstg(N), nit_amt,
     &                     nreduce, curday, pHscalar(month), dN2lyr,
     &                     dN2Olyr, prev_bgprd, Com, avgst_10cm, TI,
     &                     SI, Cr, Eh, Feh, CH4_prod, CH4_Ep, CH4_Ebl,
     &                     watertable, watrflag, bglivcj+bglivcm,
     &                     tmxbio)

      esrsnk(N) = esrsnk(N) + NOflux + Nn2oflux +Dn2oflux + Dn2flux

c ... *********************************************************************

c ... Write to output files

c ... SCALE: convert gN/m^2 to gN/ha
      SCALE = 10000.0

      nflux_sum = nflux_sum+(Nn2oflux+Dn2oflux)*SCALE
      nflux_sum2 = nflux_sum2 + NOflux*SCALE

c ... Accumulate yearly trace gas output, cak - 09/23/02
      N2O_year = N2O_year + (Nn2oflux+Dn2oflux)*SCALE
      NO_year = NO_year + NOflux*SCALE
      N2_year = N2_year + Dn2flux*SCALE
      CH4_oxid_year = CH4_oxid_year + CH4_oxid
      nit_amt_year = nit_amt_year + nit_amt*SCALE

c ... Accumulate output variables that are tracking the carbon flows
c ... due to decomposition, cak - 11/09/2010
      ametc1tosom11  = ametc1tosom11  + metc1tosom11
      ametc2tosom12  = ametc2tosom12  + metc2tosom12
      astruc1tosom11 = astruc1tosom11 + struc1tosom11
      astruc1tosom21 = astruc1tosom21 + struc1tosom21
      astruc2tosom12 = astruc2tosom12 + struc2tosom12
      astruc2tosom22 = astruc2tosom22 + struc2tosom22
      asom11tosom21  = asom11tosom21  + som11tosom21
      asom12tosom22  = asom12tosom22  + som12tosom22
      asom12tosom3   = asom12tosom3   + som12tosom3
      asom21tosom11  = asom21tosom11  + som21tosom11
      asom21tosom22  = asom21tosom22  + som21tosom22
      asom22tosom12  = asom22tosom12  + som22tosom12
      asom22tosom3   = asom22tosom3   + som22tosom3
      asom3tosom12   = asom3tosom12   + som3tosom12
      awood1tosom11  = awood1tosom11  + wood1tosom11
      awood1tosom21  = awood1tosom21  + wood1tosom21
      awood2tosom11  = awood2tosom11  + wood2tosom11
      awood2tosom21  = awood2tosom21  + wood2tosom21
      awood3tosom12  = awood3tosom12  + wood3tosom12
      awood3tosom22  = awood3tosom22  + wood3tosom22

c ... Accumulate monthly trace gas output, cak - 05/14/04
      N2O_month = N2O_month + (Nn2oflux+Dn2oflux)*SCALE
      NO_month = NO_month + NOflux*SCALE
      N2_month = N2_month + Dn2flux*SCALE
      CH4_oxid_month = CH4_oxid_month + CH4_oxid
      nit_amt_month = nit_amt_month + nit_amt*SCALE
c ... Growing season accumulator for N2O flux, cak - 06/06/2008
      n2oacc = n2oacc + (Nn2oflux+Dn2oflux)
      n2omth(month) = n2omth(month) + (Nn2oflux+Dn2oflux)

      if (time .ge. strplt) then

        call wrtdaily(time, curday, petdly, agdefac, bgdefac, stemp,
     &                snow, snlq, thermunits, aglivc, aggreenc,
     &                hwstress, scenfrac, srad(curday))

        call wrtnflux(time, curday, Nn2oflux*SCALE, Dn2oflux*SCALE,
     &                Dn2flux*SCALE, NOflux*SCALE, nflux_sum,
     &                nflux_sum2)

        call wrtsoiln(time, curday, ammonium, nitrate)

        call wrtco2(time, curday, co2_conc)

        call wrtdn2lyr(time, curday, dN2lyr)

        call wrtdn2olyr(time, curday, dN2Olyr)

        call wrtwflux(time, curday, wfluxout)

        call wrtsummary(time, curday, tempmax(curday), tempmin(curday),
     &                  ppt(curday),(Nn2oflux+Dn2oflux)*SCALE,
     &                  NOflux*SCALE, CH4_oxid, nit_amt*SCALE,
     &                  CO2resp*SCALE, ntdco2sm*10000.0)

        call wrtcflows(time, curday, som11tosom21, som12tosom22,
     &                 som12tosom3, som21tosom11, som21tosom22,
     &                 som22tosom12, som22tosom3, som3tosom12,
     &                 metc1tosom11, metc2tosom12, struc1tosom11,
     &                 struc1tosom21, struc2tosom12, struc2tosom22,
     &                 wood1tosom11, wood1tosom21, wood2tosom11,
     &                 wood2tosom21, wood3tosom12, wood3tosom22)

      endif

c ... *********************************************************************

c ... Update state variables and accumulators and sum carbon isotopes
      call flowup(time)
      call sumcar

c ... Now check for N balance and rebalance nh4 and no3 pools with minerl N

      minerl(1,N) = minerl(1,N) - Nn2oflux - Dn2oflux - Dn2flux -
     &              NOflux

      subname = 'dailymst2 '
c      call showminrl(nlayer,minerl,ammonium,nitrate,subname)
      call bal_npool(nlayer, minerl, ammonium, nitrate, inorglch)

c ... *********************************************************************
     
c ... Report the water balnce at the end of the day.
      snlq2 = snlq
      wbswc2 = 0.0
      do 108 ilyr = 1,nlayer+1
        wbswc2 = wbswc2 + asmos(ilyr)
108   continue

      if (time .ge. strplt) then
        call watrbal(curday, time, ppt(curday), accum, melt, wbswc1,
     &               wbswc2, evapdly, trandly, sublim, intrcpt,
     &               outflow, snlq1, snlq2, snow, runoffdly)
      endif

c ... *********************************************************************

c ... If the minimum temperature for the day is lower than the tmpkill
c ... parameter for the crop and the crop has accumulated at least 1/2 of
c ... the base thermal units a killing frost has occurred
      if ((tempmin(curday) .le. tmpkill) .and.
     &    (thermunits .ge. (ddbase/2.0) .and. (frtcindx .ge. 3))) then
        plntkill = .true.
      endif

      return
      end
