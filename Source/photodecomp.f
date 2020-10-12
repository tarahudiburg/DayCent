
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... PHOTODECOMP.F

      subroutine photodecomp(sradKJ, soilsrad)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'parcp.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot4.inc'
      include 'zztim.inc'

c ... Argument declarations
      real soilsrad, sradKJ

c ... Photo Decomposition
c ... Decompose standing dead and structrual litter based on incoming
c ... solar radiation
c ... written by cak 12/2009
c ... Added attached dead leaves. -mdh 9/19/2018

c ... Fortran to C prototype
      INTERFACE
        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow
      END INTERFACE

c ... Function declarations
      real     agdrat, line
      logical  candec
      external agdrat, candec, line

c ... Local variables
      integer  iel
      real     accum(ISOS), amt, biocnv, co2left, co2los, litabs,
     &         rceto1(MAXIEL), tcflow

c ... Factor to convert C to biomass is 2.5 for everything but wood.
      parameter (biocnv = 2.5)

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Standing dead has CO2 loss without nutrient loss
      if (stdedc .gt. 1.e-07) then

c ..... C/E ratios for standing dead material
        do 10 iel = 1, nelem
          rceto1(iel) = stdedc / stdede(iel)
10      continue

c ..... Check to see if decomposition can occur
        if (candec(nelem,aminrl,stdedc,stdede,1,1,rceto1)) then

c ....... Litter absorption coefficient
          litabs = line(stdedc * 2.5, 0.0, 0.0, bioabsorp, 1.0)

c ....... Compute total C flow out of standing dead pool
          tcflow = litabs * sradKJ * maxphoto
c ....... Make sure stdedc does not go negative.
          if (tcflow .gt. stdedc) then
            tcflow = stdedc
          endif

c ....... CO2 loss
          co2los = tcflow * 0.48

c ....... C flow from standing dead to CO2
          if (labtyp .eq. 2) then
            call csched(co2los,stdcis(LABELD),stdedc,
     &                  stdcis(UNLABL),stduvc2(UNLABL),
     &                  stdcis(LABELD),stduvc2(LABELD),
     &                  dresp,resp)
          else
            call csched(co2los,stdcis(LABELD),stdedc,
     &                  stdcis(UNLABL),stduvc2(UNLABL),
     &                  stdcis(LABELD),stduvc2(LABELD),
     &                  1.0,resp)
          endif

c ....... Net carbon flow from standing dead to C source/sink
          co2left = tcflow - co2los
c ....... Remaining C flows to the C source/sink with no nutrient loss
          call csched(co2left,stdcis(LABELD),stdedc,
     &                stdcis(UNLABL),csrsnk(UNLABL),
     &                stdcis(LABELD),csrsnk(LABELD),
     &                1.0,accum)
        endif
      endif

c ... Surface structural material decomposes to metabolic with CO2 loss
      if(strucc(SRFC) .gt. 1.e-07) then

c ..... Store the C/E ratios for surface structural residue in an array of
c ..... the correct size to be passed to CANDEC
        do 20 iel = 1, nelem
          rceto1(iel) = rcestr(iel)
20      continue

c ..... Check to see if decomposition can occur
        if (candec(nelem,aminrl,strucc(SRFC),struce,1,SRFC,rceto1)) then

c ....... Litter absorption coefficient
          litabs = line(strucc(SRFC) * 2.5, 0.0, 0.0, bioabsorp, 1.0)

c ....... Compute total C flow out of surface structural
          tcflow = litabs * soilsrad * maxphoto

c ....... Make sure surface structural does not go negative.
          if (tcflow .gt. strucc(SRFC)) then
            tcflow = strucc(SRFC)
          endif

c ....... CO2 loss
          co2los = tcflow * 0.48
c ....... C flow from surface structural to CO2
          if (labtyp .eq. 2) then
            call csched(co2los,strcis(SRFC,LABELD),strucc(SRFC),
     &                  strcis(SRFC,UNLABL),st1uvc2(UNLABL),
     &                  strcis(SRFC,LABELD),st1uvc2(LABELD),
     &                  dresp,resp)
          else
            call csched(co2los,strcis(SRFC,LABELD),strucc(SRFC),
     &                  strcis(SRFC,UNLABL),st1uvc2(UNLABL),
     &                  strcis(SRFC,LABELD),st1uvc2(LABELD),
     &                  1.0,resp)
          endif

c ....... Net carbon flow from surface structural to surface metabolic
          co2left = tcflow - co2los
c ....... Carbon flow from surface structural to surface metabolic
          call csched(co2left,strcis(SRFC,LABELD),strucc(SRFC),
     &                strcis(SRFC,UNLABL),metcis(SRFC,UNLABL),
     &                strcis(SRFC,LABELD),metcis(SRFC,LABELD),
     &                1.0,accum)

c ....... Nutrient flow from surface structural to surface metabolic
          do 30 iel = 1, nelem
c ......... N, P, or S flowing out of the surface structural pool is
c ......... proportional to the total carbon flow.
            amt = tcflow / rceto1(iel)
            call flow(struce(SRFC,iel),metabe(SRFC,iel),time,amt)
30        continue
        endif
      endif

c ----------------------------------------------------------------
c ... Attached dead leaves have CO2 loss without nutrient loss.
c ... -mdh 10/24/2018

      if (dleavc .gt. 1.e-07) then

c ..... C/E ratios for attached dead material
        do 40 iel = 1, nelem
          rceto1(iel) = dleavc / dleave(iel)
40      continue

c ..... Check to see if decomposition can occur
        if (candec(nelem,aminrl,dleavc,dleave,1,1,rceto1)) then

c ....... Litter absorption coefficient
          litabs = line(dleavc * 2.5, 0.0, 0.0, bioabsorp, 1.0)

c ....... Compute total C flow out of dead attached leaf pool
          tcflow = litabs * sradKJ * maxphoto
c ....... Make sure stdedc does not go negative.
          if (tcflow .gt. dleavc) then
            tcflow = dleavc
          endif

c ....... CO2 loss
          co2los = tcflow * 0.48

c ....... C flow from dead attached leaves to CO2
          if (labtyp .eq. 2) then
            call csched(co2los,dlvcis(LABELD),dleavc,
     &                  dlvcis(UNLABL),wstduvc2(UNLABL),
     &                  dlvcis(LABELD),wstduvc2(LABELD),
     &                  dresp,resp)
          else
            call csched(co2los,dlvcis(LABELD),dleavc,
     &                  dlvcis(UNLABL),wstduvc2(UNLABL),
     &                  dlvcis(LABELD),wstduvc2(LABELD),
     &                  1.0,resp)
          endif

c ....... Net carbon flow from dead attached leaves to C source/sink
          co2left = tcflow - co2los
c ....... Remaining C flows to the C source/sink with no nutrient loss
          call csched(co2left,dlvcis(LABELD),dleavc,
     &                dlvcis(UNLABL),csrsnk(UNLABL),
     &                dlvcis(LABELD),csrsnk(LABELD),
     &                1.0,accum)

c ....... Nutrient flow from dead attached leaves to E source/sink
          do 50 iel = 1, nelem
c ......... N, P, or S flowing out of the dead attached leaf pool is
c ......... proportional to the total carbon flow.
            amt = tcflow / rceto1(iel)
            call flow(dleave(iel),esrsnk(iel),time,amt)
50        continue

        endif
      endif
c ----------------------------------------------------------------

      return
      end
