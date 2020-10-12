
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine cultiv(pltlig)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      pltlig(3)

c ... Implement cultivation option

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
      integer   iel
      real      accum(ISOS), fr14, recres(MAXIEL), tagsfc, tagsoi,
     &          tbgsoi, trans, tsdsfc, tsdsoi,
     &          srfclittr, soillittr

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Update flows and sum carbon isotopes
      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

c ... Some standing dead goes into surface litter
      if (stdedc .gt. 0.0) then
        tsdsfc = stdedc * cultra(4)
        do 10 iel = 1, nelem
          recres(iel) = stdede(iel)/stdedc
10      continue
        fr14 = stdcis(LABELD)/stdedc
        call partit(tsdsfc,recres,1,stdcis,stdede,pltlig(ABOVE),fr14)
      endif

c ... Some surface litter goes into the top soil layer.

c ... Structural
      trans = strucc(SRFC) * cultra(6)

      if (trans .gt. 0.) then
        call csched(trans,strcis(SRFC,LABELD),strucc(SRFC),
     &              strcis(SRFC,UNLABL),strcis(SOIL,UNLABL),
     &              strcis(SRFC,LABELD),strcis(SOIL,LABELD),
     &              1.0,accum)

c ..... Recompute lignin fraction in structural soil C
        call adjlig(strucc(SOIL),strlig(SRFC),trans,strlig(SOIL))

        do 20 iel = 1, nelem
          trans = struce(SRFC,iel) * cultra(6)
          call flow(struce(SRFC,iel),struce(SOIL,iel),time,trans)
20      continue
      endif

c ... Metabolic
      trans = metabc(SRFC) * cultra(6)
      if (trans .gt. 0.) then
        call csched(trans,metcis(SRFC,LABELD),metabc(SRFC),
     &              metcis(SRFC,UNLABL),metcis(SOIL,UNLABL),
     &              metcis(SRFC,LABELD),metcis(SOIL,LABELD),
     &              1.0,accum)
        do 30 iel = 1, nelem
          trans = metabe(SRFC,iel) * cultra(6)
          call flow(metabe(SRFC,iel),metabe(SOIL,iel),time,trans)
30      continue
      endif

c ... Surface SOM1
      trans = som1c(SRFC) * cultra(6)
      if (trans .gt. 0.) then
        call csched(trans,som1ci(SRFC,LABELD),som1c(SRFC),
     &              som1ci(SRFC,UNLABL),som1ci(SOIL,UNLABL),
     &              som1ci(SRFC,LABELD),som1ci(SOIL,LABELD),
     &              1.0,accum)
        do 40 iel = 1, nelem
          trans = som1e(SRFC,iel) * cultra(6)
          call flow(som1e(SRFC,iel),som1e(SOIL,iel),time,trans)
40      continue
      endif

c ... Surface SOM2
      trans = som2c(SRFC) * cultra(6)
      if (trans .gt. 0.) then
        call csched(trans,som2ci(SRFC,LABELD),som2c(SRFC),
     &              som2ci(SRFC,UNLABL),som2ci(SOIL,UNLABL),
     &              som2ci(SRFC,LABELD),som2ci(SOIL,LABELD),
     &              1.0,accum)
        do 45 iel = 1, nelem
          trans = som2e(SRFC,iel) * cultra(6)
          call flow(som2e(SRFC,iel),som2e(SOIL,iel),time,trans)
45      continue
      endif

c ... Some standing dead goes to the top soil layer.
      if (stdedc .gt. 0.0) then
        tsdsoi = stdedc * cultra(5)
        fr14 = stdcis(LABELD)/stdedc
        call partit(tsdsoi,recres,2,stdcis,stdede,pltlig(ABOVE),fr14)
      endif

c ... Some above ground live goes into surface litter
      if (aglivc .gt. 0.0) then
        tagsfc = aglivc * cultra(2)
        do 50 iel = 1, nelem
          recres(iel) = aglive(iel)/aglivc
50      continue
        fr14 = aglcis(LABELD)/aglivc
        call partit(tagsfc,recres,1,aglcis,aglive,pltlig(ABOVE),fr14)
      endif

c ... Some storage pool goes to metabolic surface pool
      if (aglivc .gt. 0.0) then
        do 60 iel = 1, nelem
          trans = crpstg(iel) * cultra(2)
          call flow(crpstg(iel),metabe(SRFC,iel),time,trans)
60      continue
      endif

c ... Some above ground live goes to the top soil layer.
      if (aglivc .gt. 0.0) then
        tagsoi = aglivc * cultra(3)
        fr14 = aglcis(LABELD)/aglivc
        call partit(tagsoi,recres,2,aglcis,aglive,pltlig(ABOVE),fr14)
      endif

c ... Some storage pool goes to metabolic soil pool
      if (aglivc .gt. 0.0) then
        do 70 iel = 1, nelem
          trans = crpstg(iel) * cultra(3)
          call flow(crpstg(iel),metabe(SOIL,iel),time,trans)
70      continue
      endif

c ... Live roots go to the top soil layer.
c ... Juvenile fine roots
      if (bglivcj .gt. 0.0) then
        tbgsoi = bglivcj * cultra(7)
        do 80 iel = 1, nelem
          recres(iel) = bglivej(iel)/bglivcj
80      continue
        fr14 = bglcisj(LABELD)/bglivcj
c ..... A fraction of the live roots are transferred to the surface
c ..... litter layer, the remainder goes to the soil litter layer,
c ..... cak - 05/14/2007
        srfclittr = tbgsoi * rdsrfc
        soillittr = tbgsoi - srfclittr
        call partit(srfclittr,recres,SRFC,bglcisj,bglivej,
     &              pltlig(BELOWJ),fr14)
        call partit(soillittr,recres,SOIL,bglcisj,bglivej,
     &              pltlig(BELOWJ),fr14)
      endif
c ... Mature fine roots
      if (bglivcm .gt. 0.0) then
        tbgsoi = bglivcm * cultra(7)
        do 85 iel = 1, nelem
          recres(iel) = bglivem(iel)/bglivcm
85      continue
        fr14 = bglcism(LABELD)/bglivcm
c ..... A fraction of the live roots are transferred to the surface
c ..... litter layer, the remainder goes to the soil litter layer
c ..... cak - 05/14/2007
        srfclittr = tbgsoi * rdsrfc
        soillittr = tbgsoi - srfclittr
        call partit(srfclittr,recres,SRFC,bglcism,bglivem,
     &              pltlig(BELOWM),fr14)
        call partit(soillittr,recres,SOIL,bglcism,bglivem,
     &              pltlig(BELOWM),fr14)
      endif

c ... Some above ground live goes to standing dead
      trans = aglivc * cultra(1)
      if (trans .gt. 0.0) then
        call csched(trans,aglcis(LABELD),aglivc,
     &              aglcis(UNLABL),stdcis(UNLABL),
     &              aglcis(LABELD),stdcis(LABELD),
     &              1.0,accum)
        do 90 iel = 1, nelem
          trans = aglive(iel) * cultra(1)
          call flow(aglive(iel),stdede(iel),time,trans)
90      continue
      endif

c ... State variables and accumulators
      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar
 
c ... Check status of carbon values to make sure that everything
c ... has been reset correctly; if carbon = 0, reset elements to 0 as well
      if (aglcis(UNLABL)+aglcis(LABELD) .lt. 1.e-05) then
        aglcis(UNLABL) = 0.0
        aglcis(LABELD) = 0.0
      endif

      if (aglcis(UNLABL)+aglcis(LABELD) .le. 0.0) then
        do 100 iel = 1, MAXIEL
          aglive(iel) = 0.0
100     continue
      endif

      if (stdcis(UNLABL)+stdcis(LABELD) .lt. 1.e-05)then
        stdcis(UNLABL) = 0.0
        stdcis(LABELD) = 0.0
      endif

      if (stdcis(UNLABL)+stdcis(LABELD) .le. 0.0) then
        do 110 iel = 1, MAXIEL
          stdede(iel) = 0.0
110     continue
      endif

      do 120 iel = 1, MAXIEL
        if (aglive(iel) .lt. 1.e-05) then
          aglive(iel) = 0.0
        endif
        if (bglivej(iel) .lt. 1.e-05) then
          bglivej(iel) = 0.0
        endif
        if (bglivem(iel) .lt. 1.e-05) then
          bglivem(iel) = 0.0
        endif
        if (stdede(iel) .lt. 1.e-05) then
          stdede(iel) = 0.0
        endif
120   continue
       
      return
      end
