
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine calciv

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'ligvar.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ... Calculate initial values for temperature, water and live root
c ... carbon variables.
c ... Called from detiv.
c ... Note that variables which are functions of changeable parameters
c ... (i.e. read from 'site'.par) should be computed in prelim instead
c ... of calciv.

c ... Fortran to C prototype
      INTERFACE

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

        SUBROUTINE showminrl(nlayer, minerl, ammonium, nitrate,
     &                       subname)
          !MS$ATTRIBUTES ALIAS:'_showminrl' :: showminrl
          INTEGER          nlayer
          REAL             minerl(*)
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          CHARACTER        subname*10
        END SUBROUTINE showminrl

      END INTERFACE

c ... Local variables
      integer   iel, ilayer, iso, mm
      real      avtemp, arain, dumye(MAXIEL), dumyc(ISOS), fraclabl,
     &          frc14, k2, reclt1(MAXIEL), reclt2(MAXIEL), tcg,
c     &          storFrac,
     &          tcl, som21frac, som22frac, som3frac
      character string*80, char1*1
      character subname*10

c ... Function declarations
      real      ramp
      external  ramp

c ... Initialize soil C pools using Burke's equations.
c ...   ivauto = 0  the user has supplied the initial values
c ...   ivauto = 1  initialize using the grassland soil parameters
c ...   ivauto = 2  initialize using the crop soil parameters
c ...   ivauto = 3  initialize soil C pools for a forest soil

      subname = 'calciv    '

c ... Initialize dumyc and dumye variables.
      dumyc(LABELD) = 1000.
      dumyc(UNLABL) = 1000.
      dumye(N) = 100.
      dumye(P) = 100.
      dumye(S) = 100.

c ... Initialize the variable for computing the fraction of labeled
c ... material based on the type of labeling being done
      if (labtyp .eq. 1) then
c ..... C14 labeling
c        fraclabl = 0.0011
c ..... Calculate a fraction that is equivalent to a zero delta 14C
c ..... value, cak - 05/03/2007
        k2 = FRAC_C14 * (1.0 + (1.0 / 1000.0))
        fraclabl = k2 / (1.0 + k2)
      elseif (labtyp .eq. 2) then
c ..... C13 labeling
        fraclabl = 0.011
      else
c ..... No labeling
        fraclabl = 0.0
      endif

c ... Initialize irrtot, accumulator in irrigt.  -mdh 12/9/96
      irrtot = 0.0

c ... Compute mean annual temperature (avtemp) and mean annual 
c ... precipitation (arain)
      avtemp = 0.
      arain = 0.
      do 10 mm = 1, MONTHS
        avtemp = avtemp + (tmn2m(mm) + tmx2m(mm))/2.
        arain = arain + precip(mm)
10    continue
      avtemp = avtemp/12.
      if (avtemp .gt. 23.) then
        avtemp = 23.
      endif
      if (arain .gt. 120.) then
        arain = 120.
      endif

c ... Initialize soil C pools for a grassland soil
      if (ivauto .eq. 1) then

c ..... tcg = total soil carbon in grams (som1c + som2c + som3c)
        tcg = (-8.27E-01 * avtemp + 2.24E-02 * avtemp * avtemp +
     &         arain * 1.27E-01 - 9.38E-04 * arain * arain +
     &         arain * silt * 8.99E-02 +
     &         arain * clay * 6.00E-02 + 4.09) *1000.

c ..... Do not allow initial soil carbon values to fall below 500 g/m^2,
c ..... cak - 03/22/02
        if (tcg .lt. 500.0) then
          tcg = 500.0
        endif

c ..... Use the line function to determine the ratio of the soil carbon
c ..... in the som2 and som3 pools
        som3frac = ramp(clay, 0.08, 0.47, 0.50, 0.57)
        som21frac = (0.98 - som3frac) * 0.0
        som22frac = (0.98 - som3frac) * 1.0

c ..... Assign initial values to the labeled pools as well as unlabeled
c ..... pools using a percentage of the value from the unlabeled pool,
c ..... cak - 03/22/02
c ..... som2 pool has been split into surface and soil pools, initialize
c ..... both pools, cak - 06/14/05
c ..... Assign a fixed value to surface som1.   vek  08-91
        som1ci(SRFC,UNLABL) = 10.
        som1ci(SRFC,LABELD) = som1ci(SRFC,UNLABL) * fraclabl
        som1ci(SRFC,UNLABL) = som1ci(SRFC,UNLABL) - som1ci(SRFC,LABELD)
        som2ci(SRFC,UNLABL) = som21frac
        som2ci(SRFC,LABELD) = som2ci(SRFC,UNLABL) * fraclabl
        som2ci(SRFC,UNLABL) = som2ci(SRFC,UNLABL) - som2ci(SRFC,LABELD)

c ..... Burke's equations only apply to soil compartments.
        som1ci(SOIL,UNLABL) = tcg * .02
        som1ci(SOIL,LABELD) = som1ci(SOIL,UNLABL) * fraclabl
        som1ci(SOIL,UNLABL) = som1ci(SOIL,UNLABL) - som1ci(SOIL,LABELD)
        som2ci(SOIL,UNLABL) = tcg * som22frac
        som2ci(SOIL,LABELD) = som2ci(SOIL,UNLABL) * fraclabl
        som2ci(SOIL,UNLABL) = som2ci(SOIL,UNLABL) - som2ci(SOIL,LABELD)
        som3ci(UNLABL) = tcg * som3frac
        som3ci(LABELD) = som3ci(UNLABL) * fraclabl
        som3ci(UNLABL) = som3ci(UNLABL) - som3ci(LABELD)
        stdcis(UNLABL) = 80.
        stdcis(LABELD) = stdcis(UNLABL) * fraclabl
        stdcis(UNLABL) = stdcis(UNLABL) - stdcis(LABELD)
        stdede(N) = 1.6
        stdede(P) = .3
        stdede(S) = .3
        bglcisj(UNLABL) = 100.
        bglcisj(LABELD) = bglcisj(UNLABL) * fraclabl
        bglcisj(UNLABL) = bglcisj(UNLABL) - bglcisj(LABELD)
        bglivej(N) = 1.5
        bglivej(P) = .25
        bglivej(S) = .25
        bglcism(UNLABL) = 100.
        bglcism(LABELD) = bglcism(UNLABL) * fraclabl
        bglcism(UNLABL) = bglcism(UNLABL) - bglcism(LABELD)
        bglivem(N) = 1.5
        bglivem(P) = .25
        bglivem(S) = .25
        clittr(SRFC,UNLABL) = 100.
        clittr(SRFC,LABELD) = clittr(SRFC,UNLABL) * fraclabl
        clittr(SRFC,UNLABL) = clittr(SRFC,UNLABL) - clittr(SRFC,LABELD)
        clittr(SOIL,UNLABL) = 100.
        clittr(SOIL,LABELD) = clittr(SOIL,UNLABL) * fraclabl
        clittr(SOIL,UNLABL) = clittr(SOIL,UNLABL) - clittr(SOIL,LABELD)
      endif

c ... Initialize soil C pools for cultivated soils
      if (ivauto .eq. 2) then

c ..... tcg = total soil carbon in grams (som1c + som2c + som3c)
        tcg = (-7.50E-01 * avtemp + 2.10E-02 * avtemp * avtemp +
     &         5.81E-02 * arain -4.58E-04 * arain * arain +
     &         arain * silt * 4.94E-02 +
     &         arain * 5.82E-02 * clay + 5.15) * 1000.

c ..... Do not allow initial soil carbon values to fall below 500 g/m^2,
c ..... cak - 03/22/02
        if (tcg .lt. 500.0) then
          tcg = 500.0
        endif

c ..... Assign a fixed value to surface som1.   vek  08-91
        som1ci(SRFC,UNLABL) = 10.
        som1ci(SRFC,LABELD) = som1ci(SRFC,UNLABL) * fraclabl
        som1ci(SRFC,UNLABL) = som1ci(SRFC,UNLABL) - som1ci(SRFC,LABELD)
        som2ci(SRFC,UNLABL) = 0.0
        som2ci(SRFC,LABELD) = som2ci(SRFC,UNLABL) * fraclabl
        som2ci(SRFC,UNLABL) = som2ci(SRFC,UNLABL) - som2ci(SRFC,LABELD)

c ..... Burke's equations only apply to soil compartments. vek  08-91
        som1ci(SOIL,UNLABL) = tcg * .02
        som1ci(SOIL,LABELD) = som1ci(SOIL,UNLABL) * fraclabl
        som1ci(SOIL,UNLABL) = som1ci(SOIL,UNLABL) - som1ci(SOIL,LABELD)
        som2ci(SOIL,UNLABL) = tcg * .54
        som2ci(SOIL,LABELD) = som2ci(SOIL,UNLABL) * fraclabl
        som2ci(SOIL,UNLABL) = som2ci(SOIL,UNLABL) - som2ci(SOIL,LABELD)
        som3ci(UNLABL) = tcg * .44
        som3ci(LABELD) = som3ci(UNLABL) * fraclabl
        som3ci(UNLABL) = som3ci(UNLABL) - som3ci(LABELD)
        stdcis(UNLABL) = 20.
        stdcis(LABELD) = stdcis(UNLABL) * fraclabl
        stdcis(UNLABL) = stdcis(UNLABL) - stdcis(LABELD)
        stdede(N) = .40
        stdede(P) = .075
        stdede(S) = .075
        clittr(SRFC,UNLABL) = 10.
        clittr(SRFC,LABELD) = clittr(SRFC,UNLABL) * fraclabl
        clittr(SRFC,UNLABL) = clittr(SRFC,UNLABL) - clittr(SRFC,LABELD)
        clittr(SOIL,UNLABL) = 10.
        clittr(SOIL,LABELD) = clittr(SOIL,UNLABL) * fraclabl
        clittr(SOIL,UNLABL) = clittr(SOIL,UNLABL) - clittr(SOIL,LABELD)
      endif

c ... Initialize soil C pools for a forest soil
      if (ivauto .eq. 3) then

c ..... tcg = total soil carbon in grams (som1c + som2c + som3c)
        tcg = (-8.27E-01 * avtemp + 2.24E-02 * avtemp * avtemp +
     &         arain * 1.27E-01 - 9.38E-04 * arain * arain +
     &         arain * silt * 8.99E-02 +
     &         arain * clay * 6.00E-02 + 4.09) *1000.

c ..... Do not allow initial soil carbon values to fall below 500 g/m^2,
c ..... cak - 03/22/02
        if (tcg .lt. 500.0) then
          tcg = 500.0
        endif

c ..... Use the line function to determine the ratio of the soil carbon
c ..... in the som2 and som3 pools
        som3frac = ramp(clay, 0.08, 0.47, 0.50, 0.57)
        som21frac = (0.98 - som3frac) * 0.2
        som22frac = (0.98 - som3frac) * 0.8

c ..... Assign initial values to the labeled pools as well as unlabeled
c ..... pools using a percentage of the value from the unlabeled pool,
c ..... cak - 03/22/02
c ..... som2 pool has been split into surface and soil pools, initialize
c ..... both pools, cak - 06/14/05
c ..... Assign a fixed value to surface som1.   vek  08-91
        som1ci(SRFC,UNLABL) = 10.
        som1ci(SRFC,LABELD) = som1ci(SRFC,UNLABL) * fraclabl
        som1ci(SRFC,UNLABL) = som1ci(SRFC,UNLABL) - som1ci(SRFC,LABELD)
        som2ci(SRFC,UNLABL) = tcg * som21frac
        som2ci(SRFC,LABELD) = som2ci(SRFC,UNLABL) * fraclabl
        som2ci(SRFC,UNLABL) = som2ci(SRFC,UNLABL) - som2ci(SRFC,LABELD)

c ..... Burke's equations only apply to soil compartments.
        som1ci(SOIL,UNLABL) = tcg * .02
        som1ci(SOIL,LABELD) = som1ci(SOIL,UNLABL) * fraclabl
        som1ci(SOIL,UNLABL) = som1ci(SOIL,UNLABL) - som1ci(SOIL,LABELD)
        som2ci(SOIL,UNLABL) = tcg * som22frac
        som2ci(SOIL,LABELD) = som2ci(SOIL,UNLABL) * fraclabl
        som2ci(SOIL,UNLABL) = som2ci(SOIL,UNLABL) - som2ci(SOIL,LABELD)
        som3ci(UNLABL) = tcg * som3frac
        som3ci(LABELD) = som3ci(UNLABL) * fraclabl
        som3ci(UNLABL) = som3ci(UNLABL) - som3ci(LABELD)
        stdcis(UNLABL) = 80.
        stdcis(LABELD) = stdcis(UNLABL) * fraclabl
        stdcis(UNLABL) = stdcis(UNLABL) - stdcis(LABELD)
        stdede(N) = 1.6
        stdede(P) = .3
        stdede(S) = .3
        bglcisj(UNLABL) = 100.
        bglcisj(LABELD) = bglcisj(UNLABL) * fraclabl
        bglcisj(UNLABL) = bglcisj(UNLABL) - bglcisj(LABELD)
        bglivej(N) = 1.5
        bglivej(P) = .25
        bglivej(S) = .25
        bglcism(UNLABL) = 100.
        bglcism(LABELD) = bglcism(UNLABL) * fraclabl
        bglcism(UNLABL) = bglcism(UNLABL) - bglcism(LABELD)
        bglivem(N) = 1.5
        bglivem(P) = .25
        bglivem(S) = .25
        clittr(SRFC,UNLABL) = 100.
        clittr(SRFC,LABELD) = clittr(SRFC,UNLABL) * fraclabl
        clittr(SRFC,UNLABL) = clittr(SRFC,UNLABL) - clittr(SRFC,LABELD)
        clittr(SOIL,UNLABL) = 100.
        clittr(SOIL,LABELD) = clittr(SOIL,UNLABL) * fraclabl
        clittr(SOIL,UNLABL) = clittr(SOIL,UNLABL) - clittr(SOIL,LABELD)
      endif

c ... End of soil C pool initialization

c ... Starting values for nitrogen, phosphorous, and sulfur depend on
c ... carbon values and the ratios of carbon to each other element.
c ... Initialize structural and metabolic pools C, N, P, and S.
c ... First set them to zero and calculate N/C, P/C, & S/C ratios.

      do 40 ilayer = SRFC, SOIL
        do 20 iso = 1, ISOS
          strcis(ilayer,iso) = 0.
          metcis(ilayer,iso) = 0.
20      continue
        do 30 iel = 1, MAXIEL
          struce(ilayer,iel) = 0.
          metabe(ilayer,iel) = 0.
30      continue
40    continue

c ... Compute N/C, P/C, and S/C ratios from C/N, C/P, and C/S.
c ... This is for use in partit.
c ... Added the conditional set to zero if rcelit <= 0 -rm 7/98
      do 50 iel = 1, MAXIEL
        if (rcelit(SRFC, iel) .gt. 0.) then
          reclt1(iel) = 1. / rcelit(SRFC, iel)
        else
          reclt1(iel) = 0.0
        endif
50    continue
      do 55 iel = 1, MAXIEL
        if (rcelit(SOIL, iel) .gt. 0.) then
          reclt2(iel) = 1. / rcelit(SOIL, iel)
        else
          reclt2(iel) = 0.0
        endif
55    continue

c ... Sum carbon isotopes for use in partit.
      call sumcar

c ... Split litter C content into structural/metabolic based upon 
c ... litter C and litter lignin content and compute structural and
c ... metabolic N, P, & S based upon amount of C and the ratios
c ... computed above.
      if (initcp .ne. ' ' .and. initre .ne. ' ') then
        pltlig(ABOVE) = (wdlig(LEAF)+fligni(INTCPT,ABOVE) +
     &                  fligni(SLOPE,ABOVE) * arain) / 2.0
        pltlig(BELOWJ) = (wdlig(FROOTJ)+fligni(INTCPT,BELOWJ) +
     &                   fligni(SLOPE,BELOWJ) * arain) / 2.0
        pltlig(BELOWM) = (wdlig(FROOTM)+fligni(INTCPT,BELOWM) +
     &                   fligni(SLOPE,BELOWM) * arain) / 2.0
      else if (initcp .ne. ' ') then
        pltlig(ABOVE) = fligni(INTCPT,ABOVE)+fligni(SLOPE,ABOVE) *
     &                  arain
        pltlig(BELOWJ) = fligni(INTCPT,BELOWJ)+fligni(SLOPE,BELOWJ) *
     &                   arain
        pltlig(BELOWM) = fligni(INTCPT,BELOWM)+fligni(SLOPE,BELOWM) *
     &                   arain
      else if (initre .ne. ' ') then
        pltlig(ABOVE) = wdlig(LEAF)
        pltlig(BELOWJ) = wdlig(FROOTJ)
        pltlig(BELOWM) = wdlig(FROOTM)
      endif

c ... Total C in litter
      tcl = clittr(SRFC,UNLABL)+clittr(SRFC,LABELD)
      frc14 = clittr(SRFC,LABELD)/tcl
      call partit(tcl,reclt1,1,dumyc,dumye, 
     &            pltlig(SRFC),frc14)
      tcl = clittr(SOIL,UNLABL)+clittr(SOIL,LABELD)
      frc14 = clittr(SOIL,LABELD)/tcl
      call partit(tcl,reclt2,2,dumyc,dumye, 
     &            pltlig(SOIL),frc14)

      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

      call showminrl(nlayer,minerl,ammonium,nitrate,subname)

c ... If the C/E ratio as read from the site file is <= 0.0 throw an
c ... error message, cak - 06/14/05
      do 70 iel=1,MAXIEL
c ..... Compute N, P, and S for surface and soil som1, as well as for
c ..... som2 and som3.   vek  08-91
        if (rces1(SRFC,iel) .gt. 0.) then
          som1e(SRFC,iel)=som1c(SRFC)/rces1(SRFC,iel)
        else
          call message('   There is an error in your <site>.100 file.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES1(1,' // char1 // ') is <= 0.0'
          call message(string)
        endif
        if (rces1(SOIL,iel) .gt. 0.) then
          som1e(SOIL,iel)=som1c(SOIL)/rces1(SOIL,iel)
        else
          call message('   There is an error in your <site>.100 file.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES1(2,' // char1 // ') is <= 0.0'
          call message(string)
        endif
c ..... som2 pool has been split into surface and soil pools, compute N,
c ..... P, and S for both pools, cak - 06/14/05
        if (rces2(SRFC,iel) .gt. 0.) then
          som2e(SRFC,iel)=som2c(SRFC)/rces2(SRFC,iel)
        else
          call message('   There is an error in your <site>.100 file.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES2(1,' // char1 // ') is <= 0.0'
          call message(string)
        endif
        if (rces2(SOIL,iel) .gt. 0.) then
          som2e(SOIL,iel)=som2c(SOIL)/rces2(SOIL,iel)
        else
          call message('   There is an error in your <site>.100 file.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES2(2,' // char1 // ') is <= 0.0'
          call message(string)
        endif
        if (rces3(iel) .gt. 0.) then
          som3e(iel)=som3c/rces3(iel)
        else
          call message('   There is an error in your <site>.100 file.')
          char1 = char(ichar(char(iel)) + ichar('0'))
          string = '   RCES3(' // char1 // ') is <= 0.0'
          call message(string)
        endif
70    continue

      if (initre .ne. ' ') then
        do 80 iel = 1, MAXIEL
          if (cerfor(IVAL,FBRCH,iel) .gt. 0.) then
            wood1e(iel)=wood1c/cerfor(IVAL,FBRCH,iel)
          else
            call message('   There is an error in your tree.100 file.')
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   cerfor(1,3,' // char1 // ') is <= 0.0'
            call message(string)
          endif
          if (cerfor(IVAL,LWOOD,iel) .gt. 0.) then
            wood2e(iel)=wood2c/cerfor(IVAL,LWOOD,iel)
          else
            call message('   There is an error in your tree.100 file.')
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   cerfor(1,4,' // char1 // ') is <= 0.0'
            call message(string)
          endif
          if (cerfor(IVAL,CROOT,iel) .gt. 0.) then
            wood3e(iel)=wood3c/cerfor(IVAL,CROOT,iel)
          else
            call message('   There is an error in your tree.100 file.')
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   cerfor(1,5,' // char1 // ') is <= 0.0'
            call message(string)
          endif
80      continue
      endif

c ... Surface temperature and soil temperature
      tave = (tmn2m(1) + tmx2m(1)) / 2.0
      stemp = tave

c ... Make sure there is N, P, and S for roots
      if (bglcisj(UNLABL)+bglcisj(LABELD) .gt. 0.0) then
        do 90 iel = 1, nelem
          if (bglivej(iel) .le. 0.) then
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   Value for bglivej(' // char1
     &                // ') must be greater than 0.'
            call message(string)
            STOP
          endif
90      continue
      endif
      if (bglcism(UNLABL)+bglcism(LABELD) .gt. 0.0) then
        do 95 iel = 1, nelem
          if (bglivem(iel) .le. 0.) then
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   Value for bglivem(' // char1
     &                // ') must be greater than 0.'
            call message(string)
            STOP
          endif
95      continue
      endif

c ... Initialize grain pools
      cgrain = 0.0
      do 100 iel = 1, MAXIEL
        egrain(iel) = 0.0
100   continue

c ... Initialize crop and tree carbohydrate storage to 20-30% of live C
c ... (Bill Parton - 11/29/01)
c      storFrac = 0.25
c      carbostg(CRPSYS,UNLABL) = storFrac *
c     &                          (aglcis(UNLABL) + bglcis(UNLABL))
c      carbostg(CRPSYS,LABELD) = storFrac *
c     &                          (aglcis(LABELD) + bglcis(LABELD))
c      carbostg(FORSYS,UNLABL) = storFrac *
c     &                          (rlvcis(UNLABL) + frtcis(UNLABL) +
c     &                           crtcis(UNLABL) + rlwcis(UNLABL) +
c     &                           fbrcis(UNLABL))
c      carbostg(FORSYS,LABELD) = storFrac *
c     &                          (rlvcis(LABELD) + frtcis(LABELD) +
c     &                           crtcis(LABELD) + rlwcis(LABELD) +
c     &                           fbrcis(LABELD))

c ... To grow vegetation up from zero set an initial values for the
c ... crop and tree carbohydrate storage pools, cak - 02/23/2007
      if (cursys .eq. CRPSYS) then
        carbostg(CRPSYS,UNLABL) = 250.0
        carbostg(CRPSYS,LABELD) = carbostg(CRPSYS,UNLABL) * fraclabl
        carbostg(CRPSYS,UNLABL) = carbostg(CRPSYS,UNLABL) -
     &                            carbostg(CRPSYS,LABELD)
      else if (cursys .eq. FORSYS) then
        carbostg(FORSYS,UNLABL) = 250.0
        carbostg(FORSYS,LABELD) = carbostg(FORSYS,UNLABL) * fraclabl
        carbostg(FORSYS,UNLABL) = carbostg(FORSYS,UNLABL) -
     &                            carbostg(FORSYS,LABELD)
      else if (cursys .eq. SAVSYS) then
        carbostg(CRPSYS,UNLABL) = 250.0
        carbostg(CRPSYS,LABELD) = carbostg(CRPSYS,UNLABL) * fraclabl
        carbostg(CRPSYS,UNLABL) = carbostg(CRPSYS,UNLABL) -
     &                            carbostg(CRPSYS,LABELD)
        carbostg(FORSYS,UNLABL) = 250.0
        carbostg(FORSYS,LABELD) = carbostg(FORSYS,UNLABL) * fraclabl
        carbostg(FORSYS,UNLABL) = carbostg(FORSYS,UNLABL) -
     &                            carbostg(FORSYS,LABELD)
      endif

c ... a2drat - available E to plant demand for E.  -mdh 8/25/00
c ... Create separte a2drat arrays for crops and trees, mdh 5/11/01
      do 110 iel = 1, MAXIEL
        crop_a2drat(iel) = 1.0
        tree_a2drat(iel) = 1.0
110   continue

      return
      end
