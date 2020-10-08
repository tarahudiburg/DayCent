
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine cmplig(cursys,fligni,wdlig,pltlig)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'wth.inc'

c ... Argument declarations
      integer  cursys
      real     fligni(2,CPARTS), wdlig(FPARTS), pltlig(CPARTS)

c ... Compute plant lignin; returns the fraction of residue which will
c ... lignin.

c ... Local variables
      integer  mm
      real     arain

c ... Cursys tells whether a crop, forest or savanna system is being simulated.
      if (cursys .eq. CRPSYS .or. cursys .eq. SAVSYS) then
c ..... crop or savanna system: lignin contend depends on annual rainfall
        arain = 0.

c        do 10 mm = 1, MONTHS
c          arain = arain + prcurr(mm)
c10      continue
c        if (arain .eq. 0.) then
c          do 20 mm = 1, MONTHS
c            arain = arain + precip(mm)
c20        continue
c        endif

c ..... For RAMS/Daily Century, arain = average annual rainfall
        do 30 mm = 1, MONTHS
          arain = arain + precip(mm) * precscalar(mm)
30      continue
      endif

      if (cursys .eq. CRPSYS) then
c ..... Crop/grass
        pltlig(ABOVE)=fligni(INTCPT,ABOVE)+fligni(SLOPE,ABOVE)*arain
        pltlig(BELOWJ)=fligni(INTCPT,BELOWJ)+fligni(SLOPE,BELOWJ)*arain
        pltlig(BELOWM)=fligni(INTCPT,BELOWM)+fligni(SLOPE,BELOWM)*arain

      else if (cursys .eq. FORSYS) then
c ..... Forest system (leaves, fine roots = above ground and below ground)
        pltlig(ABOVE)=wdlig(LEAF)
        pltlig(BELOWJ)=wdlig(FROOTJ)
        pltlig(BELOWM)=wdlig(FROOTM)

      else if (cursys .eq. SAVSYS) then
c ..... Savanna
        pltlig(ABOVE) = (wdlig(LEAF)+fligni(INTCPT,ABOVE) +
     &                  fligni(SLOPE,ABOVE) * arain) / 2.0
        pltlig(BELOWJ) = (wdlig(FROOTJ)+fligni(INTCPT,BELOWJ) +
     &                  fligni(SLOPE,BELOWJ) * arain) / 2.0
        pltlig(BELOWM) = (wdlig(FROOTM)+fligni(INTCPT,BELOWM) +
     &                  fligni(SLOPE,BELOWM) * arain) / 2.0
      endif

c ... Check range for pltlig
      pltlig(ABOVE) = max(0.03, pltlig(ABOVE))
      pltlig(ABOVE) = min(0.4, pltlig(ABOVE))
      pltlig(BELOWJ) = max(0.05, pltlig(BELOWJ))
      pltlig(BELOWJ) = min(0.4, pltlig(BELOWJ))
      pltlig(BELOWM) = max(0.05, pltlig(BELOWM))
      pltlig(BELOWM) = min(0.4, pltlig(BELOWM))

      return
      end
