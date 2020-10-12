
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... FIRRTN

      subroutine firrtn()

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'npool.inc'
      include 'parfx.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'zztim.inc'

c ... Elemental return from a fire event.

c ... Called from:  frem

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
          CHARACTER subname*10
        END SUBROUTINE update_npool

      END INTERFACE

c ... Local Variables
      integer   iel, clyr
      real      accum(ISOS)
c     real      cgain, egain(MAXIEL)
      real      cpass, epass(MAXIEL)
      real      cret, eret(MAXIEL)
      character subname*10
      double precision frac_nh4, frac_no3

      subname = 'firrtn    '

      cpass = 0.0
      cret = 0.0
      do 10 iel = 1, nelem
        epass(iel) = 0.0
        eret(iel) = 0.0
10    continue

c ... Litter is being burned for the forest systems in grem.f as well,
c ... cak - 08/23/02
c      if (dofire(FORSYS)) then
c        call litburn(egain)
c      endif

c ... Return from TREE compartments
c ... Carbon return is usually 0. It is ignored since it would
c ... be returned as charcoal. (No longer true, see update below).
c ....N, P, and S returns go to the top layer of minerl.  EGAIN
c ... will contain the total returns for N, P, and S across pools.
c ... No longer returning elements from the dead fine branch
c ... and dead large wood forest components since these
c ... components are no longer burned during a TREM event,
c ... but as a FIRE event instead (see subroutine grem).
c ... cak - 01/02

c ---------------------------------------------------------------------
c ... When live shoots, standing dead shoots, surface litter, 
c ... surface SOM, and dead wood are burned, some C can be returned 
c ... as charcoal (see subroutine grem).  However, a similar C
c ... return was not implemented for burned live tree biomass.  
c ... I added it here for consitency. -mdh 9/14/2018
c ... Added dead attached leaves, fine branches, and standing dead wood. 
c ... -mdh 9/19/2018

c ..... Return carbon from burnt live leaves as charcoal
c ..... to the passive SOM pool
        if (rleavc .gt. 0.) then
          cret = remf(1) * (1.0 - lv2std(1)) * retf(1,1) * rleavc
          tcreta = tcreta + cret
          cpass = cpass + cret
          call csched(cret,rlvcis(LABELD),rleavc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif

c ..... Return carbon from burnt live fine branches as charcoal
c ..... to the passive SOM pool
        if (fbrchc .gt. 0.) then
          cret = remf(2) * (1.0 - lv2std(2)) * retf(2,1) * fbrchc
          tcreta = tcreta + cret
          cpass = cpass + cret
          tcreta = tcreta + cret
          call csched(cret,fbrcis(LABELD),fbrchc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif

c ..... Return carbon from burnt live large wood as charcoal
c ..... to the passive SOM pool
        if (rlwodc .gt. 0.) then
          cret = remf(3) * (1.0 - lv2std(3)) * retf(3,1) * rlwodc
          tcreta = tcreta + cret
          cpass = cpass + cret
          call csched(cret,rlwcis(LABELD),rlwodc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif

c ..... Return carbon from burnt dead leaves as charcoal
c ..... to the passive SOM pool
        if (dleavc .gt. 0.) then
          cret = remf(6) * retf(4,1) * dleavc
          tcreta = tcreta + cret
          cpass = cpass + cret
          call csched(cret,dlvcis(LABELD),dleavc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif

c ..... Return carbon from burnt dead fine branches as charcoal
c ..... to the passive SOM pool
        if (dfbrchc .gt. 0.) then
          cret = remf(7) * retf(5,1) * fbrchc
          tcreta = tcreta + cret
          cpass = cpass + cret
          call csched(cret,dfbrcis(LABELD),dfbrchc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif

c ..... Return carbon from burnt dead large wood as charcoal
c ..... to the passive SOM pool
        if (dlwodc .gt. 0.) then
          cret = remf(8) * retf(6,1) * dlwodc
          tcreta = tcreta + cret
          cpass = cpass + cret
          call csched(cret,dlwcis(LABELD),dlwodc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif

c ..... Associated elements to passive pool based on max C/E ratio
        do 100 iel = 1, nelem
          epass(iel) = cpass / varat3(1,iel)
c ....... Add SOM3 return to erata accumulator
c ....... Add to tereta accumulator insted. -mdh 9/26/2018
c         ereta(iel) = ereta(iel) + epass(iel)
          tereta(iel) = tereta(iel) + epass(iel)
          call flow(esrsnk(iel),som3e(iel),time,epass(iel))
100     continue

c ---------------------------------------------------------------------

      do 20 iel = 1, nelem
        eret(iel) = eret(iel) +
     &    remf(1) * retf(1,iel+1) * (1.0 - lv2std(1)) * rleave(iel) +
     &    remf(2) * retf(2,iel+1) * (1.0 - lv2std(2)) * fbrche(iel) +
     &    remf(3) * retf(3,iel+1) * (1.0 - lv2std(3)) * rlwode(iel) +
     &    remf(3) * retf(3,iel+1) * forstg(iel) +
     &    remf(6) * retf(4,iel+1) * dleave(iel) +
     &    remf(7) * retf(5,iel+1) * dfbrche(iel) +
     &    remf(8) * retf(6,iel+1) * dlwode(iel) 

c... Dead surface wood returns are calculated in subroutine grem.
c     &              remf(4) * retf(2,iel+1) * wood1e(iel) +
c     &              remf(5) * retf(3,iel+1) * wood2e(iel) +

        frac_nh4 = 0.5
        frac_no3 = 0.5
        if (iel .eq. N) then 
          clyr = 1
          call update_npool(clyr, eret(iel), frac_nh4, frac_no3, 
     &                      ammonium, nitrate, subname)
        endif
        call flow(esrsnk(iel),minerl(1,iel),time,eret(iel))
20    continue

      return
      end
