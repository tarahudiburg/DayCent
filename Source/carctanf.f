
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function carctanf(x,a,b,c,d)

      implicit none
c     include 'pi.inc'

c ... Argument declarations
      real x, a, b, c, d
      real PI
      parameter (PI=3.141592653589793)

c ... ******************** flowlib ********************* 
c ...
c ... (run-time sub-set of modaid, exclusive of modctl)
c ...
c ... release 1.0  (first formal release of modaid)
c ...
c ... james m. vevea
c ... natural resource ecology lab
c ... colorado state university
c ... fort collins, colorado  80523
c ...
c ... this routine is functionally equivalent to the routine of the
c ... same name, described in the publication:
c ...
c ... some graphs and their functional forms
c ... technical report no. 153
c ... william parton and georgo innis  (1972)
c ... natural resource ecology lab.
c ... colorado state university
c ... fort collins, colorado  80523

c     write(*,*) a, b, c, d, x
      carctanf = b + (c / PI) * atan(PI * d * (x - a))

      return
      end function carctanf
