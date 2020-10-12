
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine woodec (dtm, newminrl)

      implicit none
      include 'cflows.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'

c ... Argument declarations
      real             dtm
      double precision newminrl

c ... Wood decomposition       written by vek 04/91
c ... Added dead attached fine branches and dead standing wood. -mdh 9/19/2018

c ... defac  = decomposition factor based on water and temperature
c ...          (computed in prelim and in cycle)
c ... pligst = fixed parameter that represents the effect of
c ...          of lignin-to-structural-ratio on structural
c ...          decomposition

c ... Function declarations
      real      carctanf
      external  carctanf

c ... Local variables
      real      tcflow, pheff

c ... DEAD SURFACE FINE BRANCHES
c ...   wood1c       = C in dead surface fine branch component of forest system (gC/m2)
c ...   decw1        = intrinsic rate of decomposition of dead surface fine branches (yr-1)
c ...   wdlig(FBRCH) = lignin fraction for fine branches

      if (wood1c .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of dead surface fine branches
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = wood1c * defac * decw1 * exp(-pligst(SRFC) * 
c     &           wdlig(FBRCH)) * dtm
        tcflow = wood1c * agdefac * decw1 * exp(-pligst(SRFC) * 
     &           wdlig(FBRCH)) * dtm * pheff

c ..... Decompose fine branches into som1 and som2 with CO2 loss.
c ..... Use surface decomposition ratios, MDH Nov 2012
        call declig(aminrl,wdlig(FBRCH),SRFC,nelem,1,ps1co2,ratnew1,
     &              rsplig,tcflow,wood1c,wd1c2,wd1cis,wood1e,
     &              gromin,minerl,w1mnr,resp,som1ci,som1e,som2ci,som2e,
     &              wood1tosom11,wood1tosom21,newminrl)
      endif

c ... DEAD SURFACE LARGE WOOD
c ...   wood2c       = C in dead surface large wood component of forest system (gC/m2)
c ...   decw2        = intrinsic rate of decomposition of dead surface large wood (yr-1)
c ...   wdlig(LWOOD) = lignin fraction for large wood

      if (wood2c .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of large wood
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = wood2c * defac * decw2 * exp(-pligst(SRFC) * 
c     &           wdlig(LWOOD)) * dtm
        tcflow = wood2c * agdefac * decw2 * exp(-pligst(SRFC) * 
     &           wdlig(LWOOD)) * dtm * pheff

c ..... Decompose large wood into som1 and som2 with CO2 loss.
c ..... Use surface decomposition ratios, MDH Nov 2012
        call declig(aminrl,wdlig(LWOOD),SRFC,nelem,1,ps1co2,ratnew1,
     &              rsplig,tcflow,wood2c,wd2c2,wd2cis,wood2e,
     &              gromin,minerl,w2mnr,resp,som1ci,som1e,som2ci,som2e,
     &              wood2tosom11,wood2tosom21,newminrl)
      endif

c ... DEAD COARSE ROOTS
c ...   wood3c       = C in dead coarse root component of forest system (gC/m2)
c ...   decw3        = intrinsic rate of decomposition of dead coarse roots (yr-1)
c ...   wdlig(CROOT) = lignin fraction for coarse roots

      if (wood3c .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = carctanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of coarse roots.
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = wood3c * defac * decw3 * exp(-pligst(SOIL) *
c     &           wdlig(CROOT)) *  anerb * dtm
        tcflow = wood3c * bgdefac * decw3 * exp(-pligst(SOIL) *
     &           wdlig(CROOT)) *  anerb * dtm * pheff

c ..... Decompose coarse roots into som1 and som2 with CO2 loss.
c ..... Use soil decomposition ratios, MDH Nov 2012
        call declig(aminrl,wdlig(CROOT),SOIL,nelem,1,ps1co2,ratnew2,
     &              rsplig,tcflow,wood3c,wd3c2,wd3cis,wood3e,
     &              gromin,minerl,w3mnr,resp,som1ci,som1e,som2ci,som2e,
     &              wood3tosom12,wood3tosom22,newminrl)

      endif

c ... DEAD ATTACHED FINE BRANCHES. -mdh 9/19/2018
c ...   dfbrchc      = C in dead attached fine branch component of forest system (gC/m2)
c ...   decw4        = intrinsic rate of decomposition of dead attached fine branches (yr-1)
c ...   wdlig(FBRCH) = lignin fraction for fine branches

      if (dfbrchc .gt. 1.e-07) then

c ..... Compute total C flow out of attached dead fine branches
c ..... ATTENTION. Need substitute for agdefac? -mdh 9/19/2018
        tcflow = dfbrchc * agdefac * decw4 * exp(-pligst(SRFC) * 
     &           wdlig(FBRCH)) * dtm 

c ..... Decompose dead attached fine branches into som1 and som2 with CO2 loss.
        call declig(aminrl,wdlig(FBRCH),SRFC,nelem,1,ps1co2,ratnew1,
     &              rsplig,tcflow,dfbrchc,dwd1c2,dfbrcis,dfbrche,
     &              gromin,minerl,dw1mnr,resp,som1ci,som1e,som2ci,som2e,
     &              dwood1tosom11,dwood1tosom21,newminrl)
      endif

c ... STANDIND DEAD LARGE WOOD. -mdh 9/19/2018
c ...   dlwodc       = C in standing dead large wood component of forest system (g/m2)
c ...   decw5        = intrinsic rate of decomposition of standing dead large wood (yr-1)
c ...   wdlig(LWOOD) = lignin fraction for large wood

      if (dlwodc .gt. 1.e-07) then

c ..... Compute total C flow out of standing dead large wood
c ..... ATTENTION. Need substitute for agdefac? -mdh 9/19/2018
        tcflow = dlwodc * agdefac * decw5 * exp(-pligst(SRFC) * 
     &           wdlig(LWOOD)) * dtm 

c ..... Decompose standing dead large wood into som1 and som2 with CO2 loss.
        call declig(aminrl,wdlig(LWOOD),SRFC,nelem,1,ps1co2,ratnew1,
     &              rsplig,tcflow,dlwodc,dwd2c2,dlwcis,dlwode,
     &              gromin,minerl,dw2mnr,resp,som1ci,som1e,som2ci,som2e,
     &              dwood2tosom11,dwood2tosom21,newminrl)
      endif

      return
      end