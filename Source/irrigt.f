
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... IRRIGT.F

      real function irrigt(pptdly, petdly, tfrac)

      implicit none
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'

c ... Argument declarations
      real     pptdly, petdly, tfrac

c ... Simulate irrigation

c ... Local variables
c ... irract (a plot1 var) was changed to local variable lirract.  -mdh 12/96
c ... Also, since this is a function, irrtot should not be modified 
c ... in this function (suggestion).
      integer layr
      real    toth2o
      real    lirract
      real    adepsum

c ... Check tavedly so that irrigation water 
c ... is not added as snow
c ... Allow irrigation to occur even if the temperature is below freezing,
c ... cak - 09/16/02
c      if (tavedly .le. 0.0) then
c        lirract = 0.0

c ... Add amount given by user
c      else
        if (auirri .eq. 0) then
          lirract = irramt * tfrac

c ..... Add amount automatically to field capacity
        else if (auirri .eq. 1) then
          if (avh2o(1)/awhc .le. fawhc) then
            toth2o = pptdly
            adepsum = 0.0
            do 10 layr = 1, nlayer
              adepsum = adepsum + adep(layr)
              if (adepsum .le. 30.0) then
                toth2o = toth2o + asmos(layr)
              endif
10          continue
            lirract = max(twhc - toth2o,0.0)
          else
            lirract = 0.0
          endif

c ..... Add amount automatically to nominated amount
        else if (auirri .eq. 2) then
          if (avh2o(1)/awhc .le. fawhc) then
            lirract = irraut
          else
            lirract = 0.0
          endif

c ..... Add amount automatically to field capacity plus PET
        else if (auirri .eq. 3) then
          if (avh2o(1)/awhc .le. fawhc) then
            toth2o = pptdly
            adepsum = 0.0
            do 20 layr = 1, nlayer
              adepsum = adepsum + adep(layr)
              if (adepsum .le. 30.0) then
                toth2o = toth2o + asmos(layr)
              endif
20          continue
            lirract = max(twhc + petdly - toth2o,0.0)
          else
            lirract = 0.0
          endif
        endif
c      endif

      irrtot = irrtot + lirract
      girrtot = girrtot + lirract
      irrigt = lirract

      return
      end
