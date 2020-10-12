
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function agdrat (aminrl,varat,iel)

      implicit none

c ... Argument declarations
      integer    iel
      real       aminrl(3), varat(3,3)

c ... AboveGround Decomposition RATio computation.
c ... Determine C/E of new material entering 'Box B'.
c ... 'E' represents N, P, or S.

c ... The C/E ratios for structural and wood can be computed once; 
c ... they then remain fixed throughout the run.  The ratios for 
c ... metabolic and som may vary and must be recomputed each time
c ... step.

c ... Input:
c ...   aminrl = mineral E in top soil layer
c ...   varat : varat(1,iel) = maximum C/E ratio for material
c ...                          entering 'Box B'
c ...           varat(2,iel) = minimum C/E ratio for material
c ...           varat(3,iel) = amount of E present when minimum
c ...                          ratio applies
c ...   iel  = index to element for which ratio is being computed;
c ...          1 for N, 2 for P, 3 for S

c ... Output:
c ...   agdrat = C/E ratio of new material where E is N, P, or S
c ...            depending on the value of iel

c ... Ratio depends on available E 

      if (aminrl(iel) .le. 0.0) then

c ..... Set ratio to maximum allowed
        agdrat = varat(1,iel)

      elseif (aminrl(iel) .gt. varat(3,iel)) then

c ..... Set ratio to minimum allowed
        agdrat = varat(2,iel)

      else

c ..... aminrl(iel) > 0 and <= varat(3,iel)
        agdrat = (1.-aminrl(iel)/varat(3,iel))*
     &           (varat(1,iel)-varat(2,iel))+varat(2,iel)
c ..... where: varat(1,iel) = maximum C/E ratio for material
c .....                       entering 'Box B'
c .....        varat(2,iel) = minimum C/E ratio for material
c .....        varat(3,iel) = amount of E present when minimum
c .....                       ratio applies

      endif

      return
      end
