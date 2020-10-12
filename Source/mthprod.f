      subroutine mthprod(month)
      implicit none
      include 'const.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot3.inc'

      integer month
      integer mon

c ... This code was moved from csa_main to this new subroutine. -mdh 9/13/2018

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

c ... Added check for month > 1 and frstmth > 0 so there will 
c ... be no "index 0" error. -mdh 3/28/2014
      if (month .gt. 1 .and. frstmth .gt. 0) then
        do 50 mon = month-1, frstmth, -1
          agcmth(month) = agcmth(month) - agcmth(mon)
          bgcjmth(month) = bgcjmth(month) - bgcjmth(mon)
          bgcmmth(month) = bgcmmth(month) - bgcmmth(mon)
50      continue
      endif
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
c ... Added check for month > 1 and tfstmth > 0 so there will 
c ... be no "index 0" error. -mdh 3/28/2014
      if (month .gt. 1 .and. tfstmth .gt. 0) then
        do 100 mon = month-1, tfstmth, -1
          fcmth(month)  = fcmth(month)  - fcmth(mon)
100     continue
      endif
      if (month .lt. tfstmth) then
        do 110 mon = MONTHS, tfstmth, -1
          if (mon .ge. 1) then
            fcmth(month)  = fcmth(month)  - fcmth(mon)
          endif
110     continue
        do 120 mon = month-1, 1, -1
          if (mon .ge. 1) then
            fcmth(month)  = fcmth(month)  - fcmth(mon)
          endif
120     continue
      endif
      fcmth(month)  = max(0.0, fcmth(month))

      return
      end
