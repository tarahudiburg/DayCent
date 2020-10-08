
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine grem(tfrac)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'zztim.inc'

c ... Argument declarations
      real tfrac

c ... Simulate removal of crop/grass by fire or grazing for the month.
c ... Fire events in forest and savanna systems will burn the litter layer
c ... as well as the dead fine branches and dead large wood, cak - 08/23/02

c ... Function declarations
      real      line
      external  line

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &                          ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER          clyr
          REAL             amt
          DOUBLE PRECISION frac_nh4
          DOUBLE PRECISION frac_no3
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          CHARACTER        subname*10
        END SUBROUTINE update_npool

      END INTERFACE

c ... Local variables
      integer   ii, iel, lyr, clyr
      real      flrem, fdrem(2), fcret, shremc, shreme(MAXIEL),
     &          litrme(MAXIEL), sdremc, sdreme(MAXIEL),
     &          cret, eret(MAXIEL), feces, urine,
     &          recres(MAXIEL), ciso, friso
      real      accum(ISOS), cgain, closs, egain(MAXIEL),
     &          eloss(MAXIEL), cpass, epass(MAXIEL)
      real      esum
      double precision frac_nh4, frac_no3
      character        subname*10

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... NOTES:
c ...   ciso tells what fraction of the C returned is labeled (C14).
c ...   Initialize flrem, fdrem, and fcret based on fire or grazing.
c ...   Mod. fire routines, created 'litburn() subroutine  -mse 4-94.
c ...   If dofire(cursys)=1 -> CRPSYS burns standing dead and litter.

      subname = 'grem      '
      cret = 0.0

      if (dofire(cursys)) then
        flrem = flfrem
        fdrem(1) = fdfrem(1)
        fdrem(2) = fdfrem(2)
      else
        flrem = flgrem * tfrac
        fdrem(1) = fdgrem * tfrac
        fcret = gfcret
      endif

c ... Added for local initialization of variables which may not
c ... get initialized during a run. 8-31-90 -rm

      ciso = 0.0
      cpass = 0.0
      shremc = 0.0
      sdremc = 0.0
      do 10 ii = 1, MAXIEL
        sdreme(ii) = 0.0
        shreme(ii) = 0.0
        litrme(ii) = 0.0
        egain(ii) = 0.0
        epass(ii) = 0.0
10    continue

c ... Shoots removed
      if (aglivc .gt. 0.) then
c ..... carbon
        shremc = flrem * aglivc
        shrema = shrema + shremc
        call csched(shremc,aglcis(LABELD),aglivc,
     &              aglcis(UNLABL),csrsnk(UNLABL),
     &              aglcis(LABELD),csrsnk(LABELD),
     &              1.0,shrmai)
        ciso = ciso + (shremc*aglcis(LABELD)/aglivc)
c ..... elements
        do 20 iel = 1, nelem
          shreme(iel) = shremc*aglive(iel)/aglivc
          shrmae(iel) = shrmae(iel) + shreme(iel)
          call flow(aglive(iel),esrsnk(iel),time,shreme(iel))
20      continue
      endif

c ... Standing dead removed
      if (stdedc .gt. 0.) then
c ..... carbon
        sdremc = fdrem(1) * stdedc
        sdrema = sdrema + sdremc
        call csched(sdremc,stdcis(LABELD),stdedc,
     &              stdcis(UNLABL),csrsnk(UNLABL),
     &              stdcis(LABELD),csrsnk(LABELD),
     &              1.0,sdrmai)
        ciso = ciso + (sdremc*stdcis(LABELD)/stdedc)
c ..... elements
        do 30 iel = 1, nelem
          sdreme(iel) = sdremc*stdede(iel)/stdedc
          sdrmae(iel) = sdrmae(iel) + sdreme(iel)
          call flow(stdede(iel),esrsnk(iel),time,sdreme(iel))
30      continue
      endif

c ... FIRE
      if (dofire(cursys)) then

c ..... Residue (surface litter) removed by fire       vek 5/26/90
        call litburn(litrme)

        ciso = ciso + (fdrem(2)*strcis(SRFC,LABELD))
        ciso = ciso + (fdrem(2)*metcis(SRFC,LABELD))

c ..... Dead fine branches removed by fire, cak - 01/02
        if (wood1c .gt. 0.) then
c ....... carbon
          closs = fdfrem(3) * wood1c
          tcrem = tcrem + closs
          call csched(closs,wd1cis(LABELD),wood1c,
     &                wd1cis(UNLABL),csrsnk(UNLABL),
     &                wd1cis(LABELD),csrsnk(LABELD),
     &                1.0,accum)
c ....... elements
          do 40 iel = 1, nelem
            eloss(iel) = closs * (wood1e(iel) / wood1c)
            terem(iel) = terem(iel) + eloss(iel)
            call flow(wood1e(iel),esrsnk(iel),time,eloss(iel))
40        continue
        endif

c ..... Dead large wood removed by fire, cak - 01/02
        if (wood2c .gt. 0.) then
c ....... carbon
          closs = fdfrem(4) * wood2c
          tcrem = tcrem + closs
          call csched(closs,wd2cis(LABELD),wood2c,
     &                wd2cis(UNLABL),csrsnk(UNLABL),
     &                wd2cis(LABELD),csrsnk(LABELD),
     &                1.0,accum)
c ....... elements
          do 50 iel = 1, nelem
            eloss(iel) = closs * (wood2e(iel) / wood2c)
            terem(iel) = terem(iel) + eloss(iel)
            call flow(wood2e(iel),esrsnk(iel),time,eloss(iel))
50        continue
        endif

c ..... Carbon and nutrient return following removal by fire
c .....   fret()    - fraction of element returned by fire
c ..... The following variables have units g/m**2/month and are:
c .....   eret(iel) - total elemental return for aboveground removal

c ..... Return carbon from burning live shoots by the fire as charcoal
c ..... to the passive SOM pool
        if (aglivc .gt. 0.) then
          cgain = flrem * fret(1,1) * aglivc
          cpass = cpass + cgain
          call csched(cgain,aglcis(LABELD),aglivc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon removed from standing dead by the fire as charcoal
c ..... to the passive SOM pool
        if (stdedc .gt. 0.) then
          cgain = fdrem(1) * fret(1,1) * stdedc
          cpass = cpass + cgain
          call csched(cgain,stdcis(LABELD),stdedc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the structural component of surface
c ..... litter by the fire as charcoal to the passive SOM pool
        if (strucc(SRFC) .gt. 0.) then
          cgain = fdrem(2) * fret(1,1) * strucc(SRFC)
          cpass = cpass + cgain
          call csched(cgain,strcis(SRFC,LABELD),strucc(SRFC),
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the metabolic component of surface
c ..... litter by the fire as charcoal to the passive SOM pool
        if (metabc(SRFC) .gt. 0.) then
          cgain = fdrem(2) * fret(1,1) * metabc(SRFC)
          cpass = cpass + cgain
          call csched(cgain,metcis(SRFC,LABELD),metabc(SRFC),
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the surface component of active soil
c ..... organic matter by the fire as charcoal to the passive SOM pool
        if (som1c(SRFC) .gt. 0.) then
          cgain = fdrem(2) * fret(1,1) * som1c(SRFC)
          cpass = cpass + cgain
          call csched(cgain,som1ci(SRFC,LABELD),som1c(SRFC),
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the surface component of
c ..... intermediate soil organic matter by the fire as charcoal to the
c ..... passive SOM pool
        if (som2c(SRFC) .gt. 0.) then
          cgain = fdrem(2) * fret(1,1) * som2c(SRFC)
          cpass = cpass + cgain
          call csched(cgain,som2ci(SRFC,LABELD),som2c(SRFC),
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the dead fine branches by the fire
c ..... as charcoal to the passive SOM pool
        if (wood1c .gt. 0.) then
          cgain = fdfrem(3) * fret(2,1) * wood1c
          cpass = cpass + cgain
          call csched(cgain,wd1cis(LABELD),wood1c,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the dead large wood by the fire
c ..... as charcoal to the passive SOM pool
        if (wood2c .gt. 0.) then
          cgain = fdfrem(4) * fret(3,1) * wood2c
          cpass = cpass + cgain
          call csched(cgain,wd2cis(LABELD),wood2c,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Associated elements to passive pool based on max C/E ratio
        do 100 iel = 1, nelem
          epass(iel) = cpass / varat3(1,iel)
c ....... Add SOM3 return to erata accumulator
          ereta(iel) = ereta(iel) + epass(iel)
          call flow(esrsnk(iel),som3e(iel),time,epass(iel))
100     continue
        frac_nh4 = 0.5
        frac_no3 = 0.5
c ..... Return nutrients
        do 60 iel = 1, nelem
c ....... Burnt live shoots, standing dead, and surface litter
          sdreme(iel) = sdreme(iel) + litrme(iel)
          eret(iel) = fret(1,iel+1) * (shreme(iel) + sdreme(iel))
c ....... Burnt dead fine branches and dead large wood
          egain(iel) = egain(iel) +
     &                 fdfrem(3) * fret(2,iel+1) * wood1e(iel) +
     &                 fdfrem(4) * fret(3,iel+1) * wood2e(iel)
          esum = eret(iel) + egain(iel) - epass(iel)
          if(esum .gt. 0.) then
c ......... Add the mineral pool return to the accumulator
            ereta(iel) = ereta(iel) + esum
            if (iel .eq. N) then
              clyr = 1
              subname = 'grem1     '
              call update_npool(clyr, esum, frac_nh4, frac_no3,
     &                          ammonium, nitrate, subname)
            endif
            call flow(esrsnk(iel),minerl(1,iel),time,esum)
          endif
60      continue

c ..... END FIRE

c ... GRAZE
      else

c ..... NOTES:
c .....   Carbon and nutrient return following removal by grazing.
c .....   Grazing return with feces and urine explicitly separated.
c .....   All carbon returned by grazing is in the form of feces.
c .....     fcret     - fraction of carbon returned
c .....     gret(iel) - fraction of element returned by grazing
c .....   The following variables have units g/m**2/month and are:
c .....     cret      - amount of carbon returned to system
c .....     eret(iel) - total elemental return for aboveground removal
c .....     urine     - amount of urine returned
c .....     feces     - amount of fecal material returned (N, P, S)
c .....   To adjust for the changing lignin content of added material
c .....   strucc(1) and strlig are recomputed.
        cret = fcret * (shremc + sdremc)
        if (cret .le. 0.0) then
          cret = 0.0
          do 70 iel = 1, nelem
            eret(iel) = 0.0
70        continue
        else
          frac_nh4 = 1.0
          frac_no3 = 0.0
          do 80 iel = 1, nelem
c ......... Fraction of N that is returned is a function of clay
c ......... content, cak - 03/03/02
            if (iel .eq. N) then
              if (clay .lt. 0.0) then
                gret(iel) = 0.7
              else if (clay .gt. 0.30) then
                gret(iel) = 0.85
              else
                gret(iel) = line(clay, 0.0, 0.7, 0.30, 0.85)
              endif
            endif
            eret(iel) = gret(iel) * (shreme(iel) + sdreme(iel))
            tgzrte(iel) = tgzrte(iel) + eret(iel)
c ......... Add the mineral pool return to the accumulator
            ereta(iel) = ereta(iel) + eret(iel)
            urine= (1.-fecf(iel)) * eret(iel)
            feces= fecf(iel) * eret(iel)
            recres(iel) = feces/cret
            if (iel .eq. N) then
              clyr = 1
              subname = 'grem2     '
              call update_npool(clyr, urine, frac_nh4, frac_no3, 
     &                          ammonium, nitrate, subname)
            endif
            call flow(esrsnk(iel),minerl(1,iel),time,urine)

c ......... Add the amount of N that is volatilized from excreted
c ......... animal waste to the VOLPL and VOLPLA output variables,
c ......... cak - 03/31/04
            if (iel .eq. N) then
              volpl = volpl + ((1.0 - gret(N)) * shreme(N)) +
     &                        ((1.0 - gret(N)) * sdreme(N))
              volpla = volpla + ((1.0 - gret(N)) * shreme(N)) +
     &                          ((1.0 - gret(N)) * sdreme(N))
              volpac = volpac + ((1.0 - gret(N)) * shreme(N)) +
     &                          ((1.0 - gret(N)) * sdreme(N))
            endif
80        continue

c ....... Mod. to add structural & metabolic C into labeled (numerator)
c ....... and total (denominator) C removed.  (vek  05/26/90)
c ....... friso tells what fraction of the C returned is labeled
          lyr = 1
          friso = ciso / (shremc + sdremc)
          call partit(cret,recres,lyr,csrsnk,esrsnk,feclig,friso)
        endif

c ... END GRAZE
      endif

c ... Accumulate carbon amount returned
      creta= creta + cret + cpass

      return
      end
