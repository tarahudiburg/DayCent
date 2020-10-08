
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine crop(time, bgwfunc, tfrac, tavedly, curday, avgstemp,
     &                crpGrossPsn)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'timvar.inc'

c ... Argument declarations
      real             time
      real             bgwfunc
      real             tfrac
      real             tavedly
      real             avgstemp
      double precision crpGrossPsn
      integer          curday

c ... Driver for calling all of crop code.

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

      END INTERFACE

c ... Local variables
      real    accum(ISOS)
      real    fraclabl, k2
      integer iel, ii

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0
      fraclabl = 0.0

c ... Organic matter addition
      if (doomad .and. (omadday .eq. curday)) then
c ..... For C13/C14 labeling simulations astlbl is the concentration of
c ..... C13/C14 in the labeled C
c ..... Calculate fraction of labeled C as 14C if 14C labeling
        if (labtyp .eq. 1) then
c          fraclabl = ((astlbl / 1000.0) + 1.0) / 1000.0
          k2 = FRAC_C14 * (1.0 + (astlbl / 1000.0))
          fraclabl = k2 / (1.0 + k2)
        endif
c ..... Calculate fraction of labeled C as 13C if 13C labeling
        if (labtyp .eq. 2) then
          fraclabl = astlbl * PEEDEE * 1.0e-03 + PEEDEE
          fraclabl = 1 / (1/fraclabl + 1)
        endif
c ..... Allow for two types of organic matter addition.  Litter, for
c ..... example wheat straw, is added to the structural and metabolic
c ..... pools by the partit subroutine.  Partially decomposed organic
c ..... matter, for example compost, is added directly into the surface
c ..... slow pool (som2c(1)).
        if (omadtyp .eq. 1) then
          call partit(astgc*OMADscalar(month),astrec,1,csrsnk,esrsnk,
     &                astlig,fraclabl)
        elseif (omadtyp .eq. 2) then
          call csched(astgc*OMADscalar(month), fraclabl, 1.0,
     &                csrsnk(UNLABL), som2ci(SRFC,UNLABL),
     &                csrsnk(LABELD), som2ci(SRFC,LABELD),
     &                1.0, accum)
          do 10 iel = 1, nelem
            call flow(esrsnk(iel), som2e(SRFC,iel), time,
     &                astgc*astrec(iel))
10        continue
c ..... For omadtyp values of 3 and 4 astlbl is the fraction of material that
c ..... is labeled rather than the concentration
        elseif (omadtyp .eq. 3) then
          call partit(astgc*OMADscalar(month),astrec,1,csrsnk,esrsnk,
     &                astlig,astlbl)
        elseif (omadtyp .eq. 4) then
          call csched(astgc*OMADscalar(month), astlbl, 1.0,
     &                csrsnk(UNLABL), som2ci(SRFC,UNLABL),
     &                csrsnk(LABELD), som2ci(SRFC,LABELD),
     &                1.0, accum)
          do 15 iel = 1, nelem
            call flow(esrsnk(iel), som2e(SRFC,iel), time,
     &                astgc*astrec(iel))
15        continue
        else
          write(*,*) 'omadtyp out of range omad.100.'
          write(*,*) 'Must be equal to 1, 2, 3, or 4'
          write(*,*) 'omadtyp currently set to: ', omadtyp
          STOP
        endif
c ..... Update OMAD accumulator output variables, cak - 07/13/2006
        omadac = omadac + astgc*OMADscalar(month)
        omadmth(month) = omadmth(month)  + astgc*OMADscalar(month)
        omadtot = omadtot + astgc*OMADscalar(month)
        gomadtot = gomadtot + astgc*OMADscalar(month)
        do 20 ii = 1, nelem
          omadae(ii) = omadae(ii) +
     &                 (astgc*OMADscalar(month) * astrec(ii))
          omadmte(month, ii) = omadmte(month, ii) +
     &                         (astgc*OMADscalar(month) * astrec(ii))
          omaetot(ii) = omaetot(ii) +
     &                  (astgc*OMADscalar(month) * astrec(ii))
          gomaetot(ii) = gomaetot(ii) +
     &                   (astgc*OMADscalar(month) * astrec(ii))
20      continue
c ..... don't let organic matter be added twice in savana
        doomad = .FALSE.
      endif

c ... If microcosm selected, skip the rest of the crop code
      if (micosm .eq. 1) then
        goto 999
      endif

c ... Update flows so direct absorption will be accounted for
c ... before plant uptake.
      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

c ... Grow (growth checks crpgrw and exactly what should be done)
      call growth(tfrac, tavedly, month, crpGrossPsn)

c ... Fall of standing dead
      call falstd(pltlig, tfrac)

c ... Death of roots
      call droot(pltlig, tfrac, avgstemp)

c ... Death of shoots
      call dshoot(bgwfunc, tfrac, curday)

c ... Cultivation
      if (docult .and. (cultday .eq. curday)) then
        call cultiv(pltlig)
      endif

c ... Update state variables and accumulators and sum carbon isotopes.
      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

999   continue

      return
      end
