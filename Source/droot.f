
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine droot(pltlig, tfrac, avgstemp)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      pltlig(3)
      real      tfrac
      real      avgstemp

c ... Simulate death of roots for the month.

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
      real      carctanf, gpdf, maxswpot
      external  carctanf, gpdf, maxswpot

c ... Local variables
      integer   iel
      real      fr14, recres(MAXIEL), rdeath, rtdh,
     &          srfclittr, soillittr, tempeff, watreff
      real      accum(ISOS), cturn, eturn, temp
      real      liveCtotal, liveCremoved, liveCfrac
      real      carboLoss, carboStorage

c ... Death of roots

c ... Add code to age fine roots, juvenile fine roots age to the mature
c ... fine root pool.  Modify this subroutine so that the death rate of
c ... roots is a function of soil water potential and soil temperature,
c ... cak - 06/28/2007
c ... See:  A Model of Production and Turnover of Roots in Shortgrass Prairie
c ...       Parton, Singh, and Coleman, 1978
c ...       Journal of Applied Ecology

c ... Cap the temperature effect on fine roots at -2 and +28 degrees C
      if (avgstemp .gt. 28.0) then
        temp = 28.0
      else if (avgstemp .lt. -2.0) then
        temp = -2.0
      else
        temp = avgstemp
      endif

c ... Soil temperature effect on aging of juvenile roots
      if (bglivcj .gt. 0.0) then
        tempeff = gpdf(temp, 37.0, 0.0, 3.0, 3.0)
        cturn = cmxturn * tempeff * bglivcj * tfrac
        fr14 = bglcisj(LABELD) / bglivcj
        call csched(cturn, fr14, 1.0,
     &              bglcisj(UNLABL), bglcism(UNLABL),
     &              bglcisj(LABELD), bglcism(LABELD),
     &              1.0, accum)
        do 60 iel = 1, nelem
          eturn = cturn * (bglivej(iel) / bglivcj)
          call flow(bglivej(iel), bglivem(iel), time, eturn)
60      continue
      endif

c ... Soil temperature effect on root death rate
      tempeff = (temp - 10.0)**2 / 4.0 * 0.00175 + 0.1
      tempeff = min(tempeff, 0.5)
c ... Soil water potential effect on root death rate
      watreff = maxswpot(claypg)
      watreff = carctanf(watreff, 35.0, 0.5, 1.0, 0.05)
c ... Root death is driven by the maximum of the soil temperature
c ... effect and the soil water potential effect on root death rate,
c ... cak - 06/28/2007
      rtdh = max(tempeff, watreff)

c ... Upon root death for perennial plants transfer carbon from the
c ... carbohydrate storage pool to the C source/sink, CAK - 05/14/2014
      liveCtotal = bglivcj + bglivcm
      liveCremoved = 0.0
      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

      if (bglivcj .gt. 0.0) then
c ..... Death of juvenile fine roots
        rdeath = rdrj * tfrac * rtdh
        if (rdeath .gt. 0.95) then
          rdeath = 0.95
        endif
        rdeath = rdeath * bglivcj
        do 10 iel = 1, nelem
          recres(iel) = bglivej(iel)/bglivcj
10      continue
        fr14 = bglcisj(LABELD)/bglivcj
c ..... A fraction of the dead roots are transferred to the surface
c ..... litter layer, the remainder goes to the soil litter layer
c ..... cak - 05/14/2007
        srfclittr = rdeath * rdsrfc
        soillittr = rdeath - srfclittr
        if (frtcindx .eq. 1 .or. frtcindx .eq. 3) then
          liveCremoved = liveCremoved + srfclittr + soillittr
        endif
        call partit(srfclittr,recres,SRFC,bglcisj,bglivej,
     &              pltlig(BELOWJ),fr14)
        call partit(soillittr,recres,SOIL,bglcisj,bglivej,
     &              pltlig(BELOWJ),fr14)
      endif

      if (bglivcm .gt. 0.0) then
c ..... Death of mature fine roots
        rdeath = rdrm * tfrac * rtdh
        if (rdeath .gt. 0.95) then
          rdeath = 0.95
        endif
        rdeath = rdeath * bglivcm
        do 15 iel = 1, nelem
          recres(iel) = bglivem(iel)/bglivcm
15      continue
        fr14 = bglcism(LABELD)/bglivcm
c ..... A fraction of the dead roots are transferred to the surface
c ..... litter layer, the remainder goes to the soil litter layer
c ..... cak - 05/14/2007
        srfclittr = rdeath * rdsrfc
        soillittr = rdeath - srfclittr
        if (frtcindx .eq. 1 .or. frtcindx .eq. 3) then
          liveCremoved = liveCremoved + srfclittr + soillittr
        endif
        call partit(srfclittr,recres,SRFC,bglcism,bglivem,
     &              pltlig(BELOWM),fr14)
        call partit(soillittr,recres,SOIL,bglcism,bglivem,
     &              pltlig(BELOWM),fr14)
      endif

c ... For a perennial plant, remove C from carbohydrate storage pool
c ... based on fraction of total live roots removed, transfer carbon
c ... from the carbohydrate storage pool to the C source sink
      if (frtcindx .eq. 1 .or. frtcindx .eq. 3) then
        if (liveCtotal > 0.001) then
          liveCfrac = liveCremoved / liveCtotal
          carboStorage = carbostg(CRPSYS,UNLABL)+carbostg(CRPSYS,LABELD)
          carboLoss = max(0.0, liveCfrac * carboStorage)
          call csched(carboLoss,carbostg(CRPSYS,LABELD),carboStorage,
     &                carbostg(CRPSYS,UNLABL), csrsnk(UNLABL),
     &                carbostg(CRPSYS,LABELD), csrsnk(LABELD),
     &                1.0, accum)
        endif
      endif

      return
      end
