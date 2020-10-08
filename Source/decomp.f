
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... DECOMP.F

      subroutine decomp(dtm,decsys,amovdly,newminrl,rpeff,soilsrad)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'

c ... Argument declarations
      real      dtm
      integer   decsys
      real      amovdly(MAXLYR), rpeff, soilsrad
      double precision newminrl

c ... Decomposition Submodel (rewritten by vek 04/91)

c ... Function declarations
      real     agdrat, bgdrat
      external agdrat, bgdrat

c ... Local variables
      integer   iel

c ... ratnew1(iel,1) - the C/E ratio for new material created when a
c ...                  lignin component decomposes to SRFC som1.
c ... ratnew1(iel,2) - the C/E ratio for new material created when a
c ...                  lignin component decomposes to SRFC som2.
c ... ratnew2(iel,1) - the C/E ratio for new material created when a
c ...                  lignin component decomposes to SOIL som1.
c ... ratnew2(iel,2) - the C/E ratio for new material created when a
c ...                  lignin component decomposes to SOIL som2.

c ... Determine C/E ratios for flows from structural material to
c ... surface som1 and surface som2
c ... The second index of ratnew(iel,*) should be 1 or 2
c ... (not SRFC or SOIL) for som1c or som2c. -MDH 1/23/2012
c ... Create ratnew1 and ratnew2 for surface and soil. -MDH 9/24/2012
      do 30 iel=1,nelem
c ..... ratnew1: SRFC som1 and som2
        ratnew1(iel,1) = agdrat(aminrl,varat11,iel)
        ratnew1(iel,2) = agdrat(aminrl,varat21,iel)
c ..... ratnew2: SOIL som1 and som2
        ratnew2(iel,1) = bgdrat(aminrl,varat12,iel)
        ratnew2(iel,2) = bgdrat(aminrl,varat22,iel)
30    continue

c                   ********** LITTER **********
c ... Decompose structural and metabolic components for surface and soil.
      call litdec(dtm, newminrl, soilsrad)

c                   *********** WOOD ***********
c ... If the system is a forest or savanna...
c ... Decompose dead fine branches, large wood, and coarse roots.
c ... Dead fine roots are in the soil structural compartment.
      if (decsys .eq. FORSYS) then
        call woodec(dtm, newminrl)
      endif

c                 ***** SOIL ORGANIC MATTER *****
c ... Decompose som1 and som2 (surface and soil) and som3.
c ... Added amovdly parameter for daily version. -mdh 10/10/94
      call somdec(amovdly, dtm, newminrl, rpeff, soilsrad)

      return
      end
