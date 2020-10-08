
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c***********************************************************************
c**
c**  FILE:     getwth.f
c**
c**  PURPOSE:  Retrieve a day's worth of weather for the weather file.
c**            Compute weekly average temperature, weekly min and max
c**            temperature, weekly pet, weekly pecip, and weekly soil
c**            surface temperature using circular arrays over a 7 day
c**            period.  Compute monthly average temperature using a
c**            circular array over a 30 day period.
c**
c**  This routine was developed for the RAMS / Daily Century linkage
c**
c**  Melannie D. Hartman
c**  12/5/96
c**
c**  Add more robust checking for valid weather data values.
c**  CAK - 04/09/01
c**
c**  INPUTS:
c**     curday       - the day of the year to read weather from the file
c**     fwloss       - as read from the fix.100 file
c**     month        - current month of the year (1..12)
c**     precscalar   - monthly precipitation scalar values
c**     sitlat       - site latitude in decimal degrees as read from
c**                    <site>.100 file
c**     snow         - current snowpack (equiv. cm H2O)
c**     tmaxscalar   - monthly maximum temperature scalar values
c**     tminscalar   - monthly minimum temperature scalar values
c**     tmn2m        - average minimum air temperature for the month
c**                    as read from <site>.100 file (deg C - 2m)
c**     tmx2m        - average maximum air temperature for the month
c**                    as read from <site>.100 file (deg C - 2m)
c**     wthinput     - flag to indicate weather scalars, if any, to be
c**                    applied to data read from weather data file:
c**                      0 - No scalars used
c**                      1 - Use scalars for minimum temperature only
c**                      2 - Use scalars for maximum temperature only
c**                      3 - Use scalars for both minimum and maximum
c**                          temperatures
c**                      4 - Use scalars for precipitation only
c**                      5 - Use scalars for minimum and maximum
c**                          temperatures and precipitation
c**
c**  OUTPUTS:
c**     (From weather file):
c**     tempmax - maximum air temperature for the day (deg C)
c**     tempmin - minimum air temperature for the day (deg C)
c**     avgtemp - average air temperature for the day (deg C)
c**     ppt     - precipitation for the day (cm)
c**     solrad  - total incoming shortwave radiation (langleys/day)
c**     srad    - total incoming shortwave radiation (W/m^2/day)
c**     rhumid  - average relative humidity for the day (% 1..100)
c**     vpd     - vapor pressure deficit (kPa/day)
c**     windsp  - average daily windspeed at 2 meters (mph)
c** 
c**     tmaxwk  - average of tempmax for the last 7 days (deg C)
c**     tminwk  - average of tempmin for the last 7 days (deg C)
c**     tavewk  - average of avgtemp for the last 7 days (deg C)
c**     pptwk   - total precip for the last 7 days (cm)
c**     petwk   - total PET for the last 7 days (cm H2O)
c**     tavemth - average of avgtemp for the last 30 days (deg C)
c**     petdly  - potential evapotranspiration rate for day (cm H2O)
c**  
c**  Called by:  simsom.f
c**
c**  Calls:  none
c**
c***********************************************************************

      subroutine getwth(curday, month, tempmax, tempmin, avgtemp, ppt,
     &                  solrad, rhumid, windsp, tavewk, petdly, fwloss,
     &                  sitlat, snow, tmn2m, tmx2m, tavemth, wthinput,
     &                  precscalar, tmaxscalar, tminscalar, srad, vpd)

      implicit none
      include 'dconst.inc'
      include 'jday.inc'

c ... Formal parameters

      integer curday, month, wthinput
      real    tempmax(NDAY+1),tempmin(NDAY+1),avgtemp(NDAY+1)
      real    ppt(NDAY+1),solrad(NDAY+1),rhumid(NDAY+1),windsp(NDAY+1)
      real    petdly, tavemth, tavewk
      real    fwloss(4)
      real    sitlat
      real    snow
      real    tmn2m(NMONTH), tmx2m(NMONTH)
      real    precscalar(NMONTH), tmaxscalar(NMONTH),
     &        tminscalar(NMONTH)
      double precision srad(NDAY+1),vpd(NDAY+1)

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE calc_srad(maxTemp, minTemp, dailyPrecip, curday,
     &                       tdew, srad)
          !MS$ATTRIBUTES ALIAS:'_calc_srad' :: calc_srad
          DOUBLE PRECISION maxTemp
          DOUBLE PRECISION minTemp
          DOUBLE PRECISION dailyPrecip
          INTEGER          curday
          DOUBLE PRECISION tdew
          DOUBLE PRECISION srad
        END SUBROUTINE calc_srad

        SUBROUTINE calcpet(jday, month, tempmin, tempmax, avgtemp,
     &                     solrad, rhumid, windsp, snow, usexdrvrs,
     &                     fwloss, sitlat, tmn2m, tmx2m, petdly)
          !MS$ATTRIBUTES ALIAS:'_calcpet' :: calcpet
          INTEGER jday
          INTEGER month
          REAL    tempmin
          REAL    tempmax
          REAL    avgtemp
          REAL    solrad
          REAL    rhumid
          REAL    windsp
          REAL    snow
          INTEGER usexdrvrs
          REAL    fwloss(*)
          REAL    sitlat
          REAL    tmn2m(*)
          REAL    tmx2m(*)
          REAL    petdly
        END SUBROUTINE calcpet

      END INTERFACE

c ... Local Variables

      integer ndy, nyr, njday, nmth, imo
      integer arrayindx, ii, wkarrayindx
      real    tmaxwk, tminwk, pptwk, petwk
      real    temparray(30), tmaxarray(7), tminarray(7), tavearray(7),
     &        pptarray(7), petarray(7)
      double precision dailyPrecip, maxTemp, minTemp
      double precision tdew
      logical startofrun, debug

      data startofrun /.true./
      data debug /.false./
      save startofrun
      save arrayindx, wkarrayindx
      save temparray, tmaxarray, tminarray, tavearray, pptarray,
     &     petarray

c ... If it is the start of the run initialize all of the values in the
c ... temparray and the arrays being used to calculate the weekly
c ... averages.
      if (startofrun) then
        do 5 arrayindx = 1, 30
          temparray(arrayindx) = (tmn2m(month) + tmx2m(month)) / 2.0
5       continue
        arrayindx = mod(arrayindx, 30)
        do 7 wkarrayindx = 1, 7
          tmaxarray(wkarrayindx) = tmx2m(month)
          tminarray(wkarrayindx) = tmn2m(month)
          tavearray(wkarrayindx) = (tmn2m(month) + tmx2m(month)) / 2.0
          pptarray(wkarrayindx)  = 0.0
          petarray(wkarrayindx)  = 0.0
7       continue
        wkarrayindx = mod(wkarrayindx, 7)
        startofrun = .false.
      endif 

c ... Return here if we have reached the end of the weather data file
10    continue
c ... Use extra weather drivers for PET calculations
      if (usexdrvrs .eq. 1) then
        read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(curday),
     &                   tempmin(curday),ppt(curday),solrad(curday),
     &                   rhumid(curday), windsp(curday)
        if ((solrad(curday) .le. -99.0) .or. 
     &      (rhumid(curday) .le. -99.0) .or. 
     &      (windsp(curday) .le. -99.0)) then
          write(*,*) 'Invalid weather data, day ', curday
          STOP
        endif
        srad(curday) = -999.0
        vpd(curday) = -999.0
c ... Use extra weather drivers for photosynthesis calculations
      elseif (usexdrvrs .eq. 2) then
        read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(curday),
     &                   tempmin(curday),ppt(curday), srad(curday),
     &                   vpd(curday)
        solrad(curday) = -999.0
        rhumid(curday) = -999.0
        windsp(curday) = -999.0
c ... Use extra weather drivers for both PET and photosynthesis
c ... calculations
      elseif (usexdrvrs .eq. 3) then
        read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(curday),
     &                   tempmin(curday),ppt(curday),solrad(curday),
     &                   rhumid(curday), windsp(curday), srad(curday),
     &                   vpd(curday)
        if ((solrad(curday) .le. -99.0) .or. 
     &      (rhumid(curday) .le. -99.0) .or. 
     &      (windsp(curday) .le. -99.0)) then
          write(*,*) 'Invalid weather data, day ', curday
          STOP
        endif
c ... No extra weather drivers used
      else
        read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(curday),
     &                   tempmin(curday),ppt(curday)
        solrad(curday) = -999.0
        rhumid(curday) = -999.0
        windsp(curday) = -999.0
        srad(curday) = -999.0
        vpd(curday) = -999.0
      endif

c ... Checks for valid weather data
      if (tempmax(curday) .le. -99.0) then
        if (debug) then
          write(*,*) 'Warning: missing maximum temperature data, ',
     &               'day ', curday, ' year ', nyr
        endif
        tempmax(curday) = tmx2m(month)
      endif
      if (tempmin(curday) .le. -99.0) then
        if (debug) then
          write(*,*) 'Warning: missing minimum temperature data, ',
     &               'day ', curday, ' year ', nyr
        endif
        tempmin(curday) = tmn2m(month)
      endif
      if (ppt(curday) .le. -99.0) then
        if (debug) then
          write(*,*) 'Warning:  missing precipitation data, day ',
     &               curday, ' year ', nyr
        endif
        ppt(curday) = 0
      endif
      if (tempmax(curday) .lt. tempmin(curday)) then
        write(*,*) 'Warning:  invalid weather data, day ', curday,
     &             ' year ', nyr, ', tmax < tmin'
        tempmax(curday) = tmx2m(month)
        tempmin(curday) = tmn2m(month)
      endif

      goto 30

c ... If necessary, start reading the weather file again from the beginning      
20    rewind(9)
      goto 10

30    continue

      if (njday .ne. curday) then
        write(*,*) 'Expect day ', curday, ' got day ', njday
        write(*,*) 'in weather file.'
        STOP
      endif

      if (nmth .ne. month) then
        write(*,*) 'Expect month ', month, ' got month ', nmth
        write(*,*) 'in weather file.'
        STOP
      endif

c ... Check for leap year
      if (curday .eq. 1) then 
        if ((mod(nyr,400) .eq. 0)) then
          leapyr = .TRUE.
        else if ((mod(nyr,100) .eq. 0)) then
          leapyr = .FALSE.
        else if ((mod(nyr,4) .eq. 0)) then
          leapyr = .TRUE.
        else
          leapyr = .FALSE.
        endif
        if (leapyr) then
          dysimo(2) = idysimo(2)+1
          do 40 imo = 2, 12
            lstdy(imo) = ilstdy(imo)+1
            if (imo .gt. 2) frstdy(imo) = ifrstdy(imo)+1
40        continue
        endif
      endif

c ... Apply weather scalars to the climate values as indicated,
c ... cak - 10/18/05
      if (wthinput .eq. 4 .or. wthinput .eq. 5) then
        ppt(curday) = ppt(curday) * precscalar(month)
      endif
      if (wthinput .eq. 2 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        tempmax(curday) = tempmax(curday) + tmaxscalar(month)
      endif
      if (wthinput .eq. 1 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        tempmin(curday) = tempmin(curday) + tminscalar(month)
      endif
c ... Check that application of scalars to temperature values does
c ... not result in invalid weather values
      if (tempmin(curday) .gt. tempmax(curday)) then
        write(*,*)
        write(*,*) 'Warning: scaled minimum temperature: ',
     &             tempmin(curday)
        write(*,*) '  greater than scaled maximum temperature: ',
     &             tempmax(curday)
        write(*,*) '  on day: ', curday, ' year: ', nyr
        write(*,*) 'Setting maximum temperature to scaled minimum'
        write(*,*) 'and minimum temperature 1 degree lower.'
        write(*,*)
        tempmax(curday) = tempmin(curday)
        tempmin(curday) = tempmax(curday) - 1.0
      endif

      avgtemp(curday) = (tempmax(curday) + tempmin(curday)) / 2.0

c ... If the solar radiation value has not been read from the weather data
c ... file calculate the solar radiation for the day using Peter Thornton's
c ... code, cak - 06/18/2009
      if (srad(curday) .le. -99.0) then
        minTemp = tempmin(curday)
        maxTemp = tempmax(curday)
        dailyPrecip = ppt(curday)
        tdew = 0.0
        call calc_srad(maxTemp, minTemp, dailyPrecip, curday, tdew,
     &                 srad(curday))
c ..... Adjust the total daily radiation value returned from Peter
c ..... Thornton's code to reflect the the transmission coefficient
c ..... and cloud cover at the site.
        srad(curday) = srad(curday) * sradadj(month)
      endif

      call calcpet(curday, month, tempmin(curday), tempmax(curday),
     &             avgtemp(curday), solrad(curday), rhumid(curday),
     &             windsp(curday), snow, usexdrvrs, fwloss, sitlat,
     &             tmn2m, tmx2m, petdly)

c ... Compute weeky average values over the last 7 days, the last 7
c ... days worth of values are stored in the circular arrays,
c ... cak - 11/07/2208
      tmaxarray(wkarrayindx) = tempmax(curday)
      tminarray(wkarrayindx) = tempmin(curday)
      tavearray(wkarrayindx) = avgtemp(curday)
      pptarray(wkarrayindx)  = ppt(curday)
      petarray(wkarrayindx)  = petdly
      wkarrayindx = wkarrayindx + 1
      if (wkarrayindx .gt. 7) then
        wkarrayindx = 1
      endif

c ... Code added to compute average temperature over the last 30 days,
c ... store the last 30 days worth of temperatures in the circular
c ... array, cak - 06/10/02
      temparray(arrayindx) = avgtemp(curday)
      arrayindx = arrayindx + 1
      if (arrayindx .gt. 30) then
        arrayindx = 1
      endif

c      if (usexdrvrs .eq. 1) then
c        write(*,201) ndy,nmth,nyr,njday,tempmax(curday),tempmin(curday),
c     &               ppt(curday),petdly,solrad(curday),rhumid(curday),
c     &               windsp(curday)
c      else
c        write(*,201) ndy,nmth,nyr,njday,tempmax(curday),tempmin(curday),
c     &               ppt(curday),petdly
c      endif
c201   format(i2,2x,i2,2x,i4,2x,i3,7(f9.3,2x))

c ... Compute weekly values over the last 7 days, cak - 11/07/2007
      tmaxwk = 0.0
      tminwk = 0.0
      tavewk = 0.0
      pptwk  = 0.0
      petwk  = 0.0
      do 42 ii = 1, 7
        tmaxwk = tmaxwk + tmaxarray(ii)
        tminwk = tminwk + tminarray(ii)
        tavewk = tavewk + tavearray(ii)
        pptwk  = pptwk  + pptarray(ii)
        petwk  = petwk  + petarray(ii)
42    continue
      tmaxwk = tmaxwk / 7.0
      tminwk = tminwk / 7.0
      tavewk = tavewk / 7.0

c ... Code added to compute average temperature over the last 30 days,
c ... this value is used in the grochk subroutine for determining when
c ... leaf out should start, cak - 06/10/02
      tavemth = 0.0
      do 45 ii = 1, 30
        tavemth = tavemth + temparray(ii)
45    continue
      tavemth = tavemth / 30.0

c     write(*,202) 'tmaxwk', 'tminwk', 'tavewk', 'pptwk', 'petwk'
c202  format(5(a9,2x))
c     write(*,203) tmaxwk, tminwk, tavewk, pptwk, petwk
c203  format(5(f9.3x))
 
      return
      end
