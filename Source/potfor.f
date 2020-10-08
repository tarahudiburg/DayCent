
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine potfor(tavedly, petdly, tfrac, tavemth, srad,
     &                  daylength)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dynam.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'site.inc'

c ... Argument declarations
      real    daylength, tavedly, petdly, tfrac
      real    tavemth
      double precision srad

c ... Compute monthly potential production for forest
c ...
c ... Outputs
c ...   pforc gC/m^2/time
c ...   h2ogef(2)

c ... RESPPT is an array of respiration values for the forest system
c ... production components.

c ... Added savanna model, pcropc now pforc
c ...                      tgprod now tfprod (BO)

c ... Function declarations
      real     gpdf, frespr, laprod, pprdwc
      external gpdf, frespr, laprod, pprdwc

c ... Local variables
      real     fcmax, frlive, lai, langleys, potprd

c ... Estimate potential production based on temp & h2o

      if (tavedly .gt. 0.0) then

c ..... Calculate temperature effect on growth.
        potprd = gpdf(tavedly, ppdf(1,2), ppdf(2,2), ppdf(3,2), 
     &                ppdf(4,2))

c ..... Added to match version 3.0 -lh 4/93
        potprd = potprd * .8

c ..... Calculate moisture effect on growth
        if (petdly .ge. .01) then
c ....... Calculate potential growth based on the relative water content
c ....... of the wettest soil layer, cak - 12/06/04
          h2ogef(2) = 1.0/(1.0+exp(wscoeff(2,2)*
     &                             (wscoeff(2,1)-twstress)))
        else
          h2ogef(2) = 0.01
        endif

c ..... For large wood, calculate the percentage which is live (sapwood)
        frlive = sapk / (sapk + rlwodc)
 
c ..... Calculate LAI and then use in function for LAI reducer on
c ..... production.  -rm 5/91
c ..... Calculate theoretical maximum for LAI based on large wood biomass,
c ..... cak - 07/24/02
c ..... Include fine branch carbon as part of the woody component in the
c ..... LAI calculation, cak - 10/20/2006
        call lacalc(lai, fbrchc, rlwodc, maxlai, klai)

c ..... Use the solar radiation value as calculated by Peter Thornton's
c ..... subroutines in the calculation of potential production,
c ..... cak - 06/18/2009
c ..... Convert the solar radiation value from W/m^2 to langleys/day
        langleys = srad * ((daylength * 60 * 60) / (10000 * 4.184))
        fcmax = langleys * prdx(2) * tfrac * potprd * h2ogef(2) *
     &          laprod(lai,laitop) * co2cpr(FORSYS)
        pforc = fcmax

c ..... Compute carbon allocation fractions for each tree part
c ..... mdh 5/11/01
c ..... Call moved from treegrow subroutine, cak - 07/01/02
        if (pforc .gt. 0.01) then
          call treeDynC(pforc, tavemth, tree_cfrac)
        else
          pforc = 0.0
        endif

      else
        pforc = 0.
      endif

      return
      end
