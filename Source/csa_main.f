
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c                           DISCLAIMER
c
c        Neither the Great Plains System Research Unit - USDA (GPSR) nor
c     Colorado State University (CSU) nor any of their employees make
c     any warranty or assumes any legal liability or responsibility for
c     the accuracy, completeness, or usefulness of any information,
c     apparatus, product, or process disclosed, or represents that its
c     use would not infringe privately owned rights.  Reference to any
c     special commercial products, process, or service by tradename,
c     trademark, manufacturer, or otherwise, does not necessarily
c     constitute or imply endorsement, recommendation, or favoring by  
c     the GPSR or CSU.  The views and opinions of the authors do not
c     necessarily state or reflect those of GPSR or CSU and shall not 
c     be used for advertising or product endorsement. 

      program main

c ... Century Soil Organic Matter Model
c ... Simulation of carbon, nitrogen, phosphorous, and sulfur cycling
c ... As of Dec. 1991, uses a 1 month time step
c ... Project - Soil Fertility in the Great Plains
c ... Modeler - Bill Parton
c ... Programmers - Vicki Kirchner, Becky McKeown, Laura Harding,
c ...               Melannie Hartman

c ... State variables and flows are grams/m2.

      implicit none
      include 'cflows.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'jday.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ...               (unit 1) = plot/print file used by modaid (unformatted)
c ... <site>.100    (unit 7) = parameter values and initial values for 
c ...                          state variables; see subroutine sitein.
c ... fix.100       (unit 8) = fixed parameter values values for 
c ...                          state variables; see subroutine fixin.
c ...               (unit 9) = a file of weather data read in subroutines 
c ...                          wthini, weathr
c ... c14data      (unit 10) = a data file specifying years in which labeled
c ...                          carbon is added via plant growth and what 
c ...                          fraction of the growth is labeled.
c ... nscale.dat   (unit 20) = a data file specifying years in which N input
c ...                          scalars are used and the scalar values
c ... omadscale.dat(unit 30) = a data file specifying years in which organic
c ...                          matter input scalars are used and the scalar
c ...                          values
c ... phscale.dat  (unit 40) = a data file specifying years in which pH
c ...                          scalars are used and the scalar values
c ... precscale.dat(unit 50) = a data file specifying years in which
c ...                          precipitation scalars are used and the scalar
c ...                          values, precipitation scalar are multipliers
c ... tmaxscale.dat(unit 55) = a data file specifying years in which
c ...                          maximum temperature scalars are used and the
c ...                          scalar values, maximum temperature scalar are
c ...                          addends
c ... tminscale.dat(unit 60) = a data file specifying years in which
c ...                          minimum temperature scalars are used and the
c ...                          scalar values, minimum temperature scalar are
c ...                          addends
c ... nflux.out    (unit 70) = N2/N2O fluxes computed by Trace Gas Model
c ... daily.out    (unit 80) = pet, defac, stemp, and snowpack water content
c ...                          computed by Trace Gas Model
c ... summary.out  (unit 90) = tmax, tmin, prec, N2O flux, NO flux, CH4
c ...                          oxidation, and gross nitrification computed
c ...                          by Trace Gas Model

c ... If you're getting floating point errors mentioned after you exit
c ... Century, uncomment the following lines, recompile, run Century
c ... in dbx with the 'catch FPE' option to find the offending code.
c ... You can also run Century outside of dbx, in which case you will
c ... get messages on your screen giving you an fpe error code (see
c ... the Floating Point Programmer's Guide, p20) and a not-very-
c ... useful-to-3rd-or-4th-generation-language-programmers location. 
c ... The error handler 'mysigal' is a Fortran callable C routine 
c ... written by Martin Fowler; it can be replaced by any user written
c ... handler or any of several library handlers, the most useful 
c ... probably being SIGFPE_ABORT.  The calls to ieee_handler won't 
c ... compile using poa's binaries.

c      external mysigal
c      ieeer=ieee_handler('set','invalid',SIGFPE_ABORT)
c      ieeer=ieee_handler('set','division',mysigal)
c      ieeer=ieee_handler('set','overflow',mysigal)
c      ieeer=ieee_handler('set','underflow',SIGFPE_ABORT)

c ... You probably won't want to uncomment the following line; inexact
c ... floating point signals occur all over the place.

c      ieeer=ieee_handler('set','inexact',mysigal)

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE closefiles()
          !MS$ATTRIBUTES ALIAS:'_closefiles' :: closefiles
        END SUBROUTINE closefiles

        SUBROUTINE wrttgmonth(time, N2O_month, NO_month, N2_month,
     &                        CH4_oxid_month, nit_amt_month, pptmonth)
          !MS$ATTRIBUTES ALIAS:'_wrttgmonth' :: wrttgmonth
          REAL             time
          REAL             N2O_month
          REAL             NO_month
          REAL             N2_month
          REAL             CH4_oxid_month
          DOUBLE PRECISION nit_amt_month
          REAL             pptmonth
        END SUBROUTINE wrttgmonth

        SUBROUTINE wrtyearsum(time, N2O_year, NO_year, N2_year,
     &                        CH4_oxid_year, nit_amt_year, annppt)
          !MS$ATTRIBUTES ALIAS:'_wrtyearsum' :: wrtyearsum
          REAL             time
          REAL             N2O_year
          REAL             NO_year
          REAL             N2_year
          REAL             CH4_oxid_year
          DOUBLE PRECISION nit_amt_year
          REAL             annppt
        END SUBROUTINE wrtyearsum

        SUBROUTINE wrtyrcflows(time, asom11tosom21, asom12tosom22,
     &                         asom12tosom3, asom21tosom11,
     &                         asom21tosom22, asom22tosom12,
     &                         asom22tosom3, asom3tosom12,
     &                         ametc1tosom11, ametc2tosom12,
     &                         astruc1tosom11, astruc1tosom21,
     &                         astruc2tosom12, astruc2tosom22,
     &                         awood1tosom11, awood1tosom21,
     &                         awood2tosom11, awood2tosom21,
     &                         awood3tosom12, awood3tosom22)
          !MS$ATTRIBUTES ALIAS:'_wrtyrcflows' :: wrtyrcflows
          REAL time
          REAL asom11tosom21
          REAL asom12tosom22
          REAL asom12tosom3
          REAL asom21tosom11
          REAL asom21tosom22
          REAL asom22tosom12
          REAL asom22tosom3
          REAL asom3tosom12
          REAL ametc1tosom11
          REAL ametc2tosom12
          REAL astruc1tosom11
          REAL astruc1tosom21
          REAL astruc2tosom12
          REAL astruc2tosom22
          REAL awood1tosom11
          REAL awood1tosom21
          REAL awood2tosom11
          REAL awood2tosom21
          REAL awood3tosom12
          REAL awood3tosom22
        END SUBROUTINE wrtyrcflows

      END INTERFACE

c ... Local variables
      integer mon
      logical useprompts, wrtsite
      real month_vals(12), neg_month_vals(12)

      data month_vals /0.08, 0.17, 0.25, 0.33, 0.42, 0.50, 0.58, 0.67,
     &                 0.75, 0.83, 0.92, 1.0/
      data neg_month_vals /-0.92, -0.83, -0.75, -0.67, -0.58, -0.50,
     &                     -0.42, -0.34, -0.25, -0.17, -0.08, 0.0/

      wrtsite = .false.

      data (idysimo(k), k=1,12) /31,28,31,30,31,30,31,31,30,31,30,31/
      data (ilstdy(k), k=1,12) 
     &                 /31,59,90,120,151,181,212,243,273,304,334,365/
      data (ifrstdy(k), k=1,12) 
     &                 /1,32,60,91,121,152,182,213,244,274,305,335/

c ... Obtain startup information from user, do initializations based on
c ... answers to Modaid questions
      call detiv(useprompts, wrtsite)

c ... Adjust parameters from crop.100 and fix.100 for daily production
c ... if necessary. -mdh 1/95
      call adjustpar

c ... Write out starting values
      call wrtbin(time)

c ... Update month
20    continue
      month = mod(month,12) + 1

c ... If time is greater than the ending time for the current block,
c ... read the next block
      if ((abs(time - blktnd) .lt. (0.5 * dt)) .and.
     &    (abs(time - tend)   .gt. (0.5 * dt))) then
        call readblk()
      endif

c ... Perform annual tasks
      if (month .eq. 1) then
        call eachyr
      endif

c ... The main driver for the model; call decomposition, growth, etc.
      call simsom()

c ... Update monthly production output variables for grass/crop system,
c ... cak - 10/01/03
      agcmth(month) = agcacc
      bgcjmth(month) = bgcjacc
      bgcmmth(month) = bgcmacc

C ... Prevent index out of bounds for monhtly accumulators. -mdh 11/7/2017
C     do 50 mon = month-1, frstmth, -1
C       agcmth(month) = agcmth(month) - agcmth(mon)
C       bgcjmth(month) = bgcjmth(month) - bgcjmth(mon)
C       bgcmmth(month) = bgcmmth(month) - bgcmmth(mon)
C50    continue
C     if (month .lt. frstmth) then
C       do 60 mon = MONTHS, frstmth, -1
C         agcmth(month) = agcmth(month) - agcmth(mon)
C         bgcjmth(month) = bgcjmth(month) - bgcjmth(mon)
C         bgcmmth(month) = bgcmmth(month) - bgcmmth(mon)
C60      continue
C       do 70 mon = month-1, 1, -1
C         agcmth(month) = agcmth(month) - agcmth(mon)
C         bgcjmth(month) = bgcjmth(month) - bgcjmth(mon)
C         bgcmmth(month) = bgcmmth(month) - bgcmmth(mon)
C70      continue
C     endif

      do 50 mon = month-1, frstmth, -1
        if (mon .ge. 1) then
          agcmth(month) = agcmth(month) - agcmth(mon)
          bgcjmth(month) = bgcjmth(month) - bgcjmth(mon)
          bgcmmth(month) = bgcmmth(month) - bgcmmth(mon)
        endif
50    continue
      if (month .lt. frstmth) then
        do 60 mon = MONTHS, frstmth, -1
          agcmth(month) = agcmth(month) - agcmth(mon)
          bgcjmth(month) = bgcjmth(month) - bgcjmth(mon)
          bgcmmth(month) = bgcmmth(month) - bgcmmth(mon)
60      continue
        do 70 mon = month-1, 1, -1
          if (mon .ge. 1) then
            agcmth(month) = agcmth(month) - agcmth(mon)
            bgcjmth(month) = bgcjmth(month) - bgcjmth(mon)
            bgcmmth(month) = bgcmmth(month) - bgcmmth(mon)
          endif
70      continue
      endif

      agcmth(month) = max(0.0, agcmth(month))
      bgcjmth(month) = max(0.0, bgcjmth(month))
      bgcmmth(month) = max(0.0, bgcmmth(month))

c ... Update monthly production output variables for forest system,
c ... cak - 10/01/03
      fcmth(month)  = fcacc
      do 100 mon = month-1, tfstmth, -1
        fcmth(month)  = fcmth(month)  - fcmth(mon)
100   continue
      if (month .lt. tfstmth) then
        do 110 mon = MONTHS, tfstmth, -1
          fcmth(month)  = fcmth(month)  - fcmth(mon)
110     continue
        do 120 mon = month-1, 1, -1
          fcmth(month)  = fcmth(month)  - fcmth(mon)
120     continue
      endif
      fcmth(month)  = max(0.0, fcmth(month))

c ... Write yearly output
c ... Add output for the N2 flux for the year and convert fluxes to
c ... g/m^2, cak - 01/16/03
      if ((time .ge. strplt) .and. (month .eq. 12)) then
        call wrtyearsum(time, N2O_year/10000, NO_year/10000,
     &                  N2_year/10000, CH4_oxid_year/10000,
     &                  nit_amt_year/10000, annppt)
        call wrtyrcflows(time, asom11tosom21, asom12tosom22,
     &                   asom12tosom3, asom21tosom11, asom21tosom22,
     &                   asom22tosom12, asom22tosom3, asom3tosom12,
     &                   ametc1tosom11, ametc2tosom12, astruc1tosom11,
     &                   astruc1tosom21, astruc2tosom12,
     &                   astruc2tosom22, awood1tosom11, awood1tosom21,
     &                   awood2tosom11, awood2tosom21, awood3tosom12,
     &                   awood3tosom22)
      endif

c ... Write monthly trace gas output, cak - 05/14/42
      if (time .ge. strplt) then
        call wrttgmonth(time, N2O_month/10000, NO_month/10000,
     &                  N2_month/10000, CH4_oxid_month/10000,
     &                  nit_amt_month/10000, pptmonth)
      endif

c ... Update time
c ... Add calculation to handle time drift caused by inexact floating
c ... point addition, cak - 08/23/01
c      time = time + dt
c ... Add calculation to handle incrementing the month for negative years,
c ... cak - 03/29/02
      if (time .ge. 0.0) then
        time = int(time) + month_vals(month)
      else
        time = int(time) + neg_month_vals(month)
        if (month .eq. 1) then
          time = time + 1.0
        endif
      endif
      if (time .ge. -1.0e-07 .and. time .le. 1.0e-07) then
        time = 0.0
      endif

c ... Write out values
      if ((tplt - time) .lt. (dt * 0.5)) then
        call wrtbin(time)
        tplt = time + dtpl
      endif

c ... Run for tend years
      if ((tend-time) .gt. (dt*.5)) then
        goto 20
      endif

c ... Write out final values
      call wrtbin(time)

c ... Write out extended site file, if desired
      if (wrtsite) then
        call wrtextsite()
      endif

c ... Close data files

c ... Close the weather file
      close(unit=9)
c ... Close the c14data file if necessary
      if (labtyp .gt. 0) then
        close(unit=10)
      endif
c ... Close the nscale.dat file if necessary
      if (Ninput .gt. 0) then
        close(unit=20)
      endif
c ... Close the omadscale.dat file if necessary
      if (OMADinput .gt. 0) then
        close(unit=30)
      endif
c ... Close the phscale.dat file if necessary
      if (phsys .gt. 0) then
        close(unit=40)
      endif
c ... Close the precscale.dat file if necessary
      if (wthinput .eq. 4 .or. wthinput .eq. 5) then
        close(unit=50)
      endif
c ... Close the tmaxscale.dat file if necessary
      if (wthinput .eq. 2 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        close(unit=55)
      endif
c ... Close the tminscale.dat file if necessary
      if (wthinput .eq. 1 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        close(unit=60)
      endif
c ... Close the schedule file
      close(unit=15)
c ... Close N2/N2O flux file, nflux.out (unit=70)
      close(unit=70)
c ... Close the daily.out file (unit=80)
      close(unit=80)
c ... Close the summary.out file (unit=90)
      close(unit=90)
c ... Close *.out and *.csv files opened by the initsw subroutine
      call closefiles()

c ... Mark end of file
      endfile(unit=1)

c ... Close binary file
      close(unit=1)

      if (useprompts) then
        call message('')
        call message(' Execution success.')
        call message('')
        call message(' Press <Enter> to continue:')
        read(*,*)
      else
        STOP 'Execution success.'
      endif

      end
