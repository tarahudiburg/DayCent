
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine potcrp (cancvr, tavedly, petdly, tfrac, scenfrac,
     &                   srad, daylength)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'

c ... Argument declarations
      real    cancvr
      real    daylength, tavedly, petdly, tfrac, scenfrac
      double precision srad

c ... Compute monthly production potential based upon montly precip
c ... and restrict potential production based upon method specified
c ... by grzeff.

c ... Function declarations
      real     gpdf, pprdwc
      external gpdf, pprdwc

c ... Local variables
      real     agprod, aisc, bgp, bgprod, bioc, biof,
     &         bioprd, fracrc, langleys, potprd, ratlc, rtsh,
     &         sdlng, shdmod, subcan, temp1, temp2, temp3

c ... Compute shading modifier for savanna
      if (cursys .eq. SAVSYS) then
        if (cancvr .le. 0.001) then
          aisc = 0.
        else
          aisc = 5 * exp(-.0035 * (rleavc*2.5)/cancvr)
        endif
        subcan = aisc/(aisc + 1.)
        shdmod = (1.0-cancvr) + (cancvr*subcan)
      else
        shdmod = 1.0
      endif

c ... Calculate moisture effect on growth
      if (petdly .ge. .01) then
c ..... Calculate potential growth based on the relative water content
c ..... of the wettest soil layer, cak - 12/06/04
        h2ogef(1) = 1.0/(1.0+exp(wscoeff(1,2)*(wscoeff(1,1)-cwstress)))
      else
        h2ogef(1) = 0.01
      endif

c ... Estimate plant production:
      if (tavedly .gt. 0.0) then

c ..... Calculate temperature effect on growth
        potprd = gpdf(tavedly, ppdf(1,1), ppdf(2,1), ppdf(3,1), 
     &                ppdf(4,1))

c ..... Calculate biof
        if (bioflg .eq. 1) then

c ....... Calculate maximum potential effect of standing dead on plant growth
c ....... (the effect of physical obstruction of litter and standing dead)
          bioc = stdedc + .1*strucc(SRFC)
          if (bioc .le. 0.) then
            bioc = .01
          endif

          if (bioc .gt. pmxbio) then
            bioc = pmxbio
          endif
          bioprd = 1. - (bioc/(biok5+bioc))

c ....... Calculate the effect of the ratio of live biomass to dead biomass
c ....... on the reduction of potential growth rate.  The intercept of this 
c ....... equation ( highest negative effect of dead plant biomass ) is equal
c ....... to bioprd when the ratio is zero.
          temp1 = (1. - bioprd)
          temp2 = temp1*0.75
          temp3 = temp1*0.25
          ratlc = aglivc/bioc 
          if (ratlc .le. 1.0) then
            biof = bioprd+(temp2*ratlc)
          endif
          if (ratlc .gt. 1.0 .and. ratlc .le. 2.0) then
            biof = (bioprd+temp2)+temp3*(ratlc-1.)
          endif
          if (ratlc .gt. 2.0) then
            biof = 1.0
          endif
        else
          biof = 1.0
        endif

c ..... Restriction on seedling growth
c ..... sdlng is the fraction that prdx is reduced
        if (aglivc .gt. fulcan) then
          seedl = 0
        endif

        if (seedl .eq. 1) then
          sdlng = min(1.0, pltmrf + aglivc*(1-pltmrf) /fulcan)
        else
          sdlng = 1.0
        endif

c ..... Compute potential total production
c ..... Use the solar radiation value as calculated by Peter Thornton's
c ..... subroutines in the calculation of potential production,
c ..... cak - 06/18/2009
c ..... For an annual plant use the scenfrac value to slow the growth
c ..... rate in the period from anthesis to maturity
c ..... Convert the solar radiation value from W/m^2 to langleys/day
        langleys = srad * ((daylength * 60 * 60) / (10000 * 4.184))
        tgprod = langleys * prdx(1) * potprd * h2ogef(1) * biof *
     &           shdmod * sdlng * co2cpr(CRPSYS) * tfrac * scenfrac

c ..... Dynamic carbon allocation routines for crop/grsss, used to compute
c ..... root/shoot ratio, cak - 07/01/02
c ..... Do not call the dynamic carbon allocation routine when there is no
c ..... production, cak - 09/09/02
        if (tgprod .gt. 0.0) then
          call cropDynC(rtsh, fracrc)
        else
          tgprod = 0.0
          agp = 0.0
          pcropc = 0.0
          goto 40
        endif

cc ..... Change root/shoot ratio if burning occurs
c        if (firecnt .ge. 1) then
c          rtsh = rtsh + frtsh
c          firecnt = firecnt + 1
c          if (firecnt .gt. 5) then
c            firecnt = 0
c          endif
c        endif

c ..... When simulating an annual crop if emergence has not occurred
c ..... no production occurs, 06/05/2014
        if ((frtcindx .ge. 4) .and. (.not. emerg)) then
          tgprod = 0.0
          agp = 0.0
          pcropc = 0.0
          goto 40
        else
c ....... Use the fraction of carbon allocated to the roots rather than
c ....... root to shoot ratio to determine amount of aboveground and
c ....... belowground production, cak - 08/22/03
c          bgprod = agprod * rtsh
          bgprod = tgprod * fracrc
          agprod = tgprod - bgprod
        endif
        agp = agprod
        bgp = bgprod
        tgprod = agp + bgp
      else
c ..... No production this time step
        tgprod = 0.0
        agp = 0.0
        pcropc = 0.0
        goto 40
      endif

c ... Determine if grazing occurs
      if (grazcnt .ge. 1) then
        call grazrst(agp, bgp, flgrem, gremb, grzeff, rtsh, tgprod)
        grazcnt = grazcnt + 1
        if (grazcnt .gt. 31) then
          grazcnt = 0
        endif
      endif

c ... Update accumulators & compute potential C production
      ptagc = ptagc + agp/2.5
      ptbgc = ptbgc + bgp/2.5
      pcropc = tgprod / 2.5

40    continue

      return
      end
