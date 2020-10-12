
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function maxswpot(numlayers)

      implicit none
      include 'parfx.inc'
      include 'plot1.inc'
      include 'site.inc'

c ... Argument declarations
      integer numlayers

c ... This function calculates the soil water potential of the wettest
c ... soil layer in the plant rooting zone.

c ... Declaration explanations:
c ...   numlayers - number of soil layers in the plant rooting zone

c ... Local variable explanations:
c ...   b      - slope of retention curve
c ...   BAR2CM - conversion factor for bars to centimeters H2O
c ...   base   - base value for power function
c ...   expon  - exponent value for power function
c ...   lyr    - current soil layer
c ...   psis   - "saturation" matric potential of "ilyr" (cm H2O ?)
c ...   swptnl - soil water potential of the current layer (bars)
c ...   theta  - volumetric soil water content * 100
c ...   thetas - volumetric soil water content at saturation for layer
c ...            (% volume)

c ... Local variables
      integer          BAR2CM, lyr
      real             b, psis, swptnl, theta, thetas
      double precision base, expon 

c ... 
      BAR2CM = 1024

      maxswpot = 9999

      do 130 lyr = 1, numlayers
c ..... Calculate the soil water potential for the current soil layer
        if (asmos(lyr) .gt. 0) then
          theta = (asmos(lyr) / adep(lyr)) * 100
          thetas = (-14.2 * sand) - (3.7 * clay) + 50.5
          base =  theta / thetas
          b = (-0.3 * sand) + (15.7 * clay) + 3.10
          expon = b
          psis = 10.0**((-1.58 * sand) - (0.63 * clay) + 2.17)
          swptnl = (psis / (base**expon)) / BAR2CM
        else
          swptnl = 80.0
        endif

        if (swptnl .lt. maxswpot) then
          maxswpot = swptnl
        endif
130   continue

c ... Place a limit the maximum water potential of the wettest layer
c ... because roots are not able to extract water below this level,
c ... cak - 03/20/2008
      if (maxswpot .gt. 30.0) then
        maxswpot = 30.0
      endif

      return
      end
