
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine initprod(system, month)

      implicit none
      include 'const.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot3.inc'

c ... Argument declarations
      integer system, month

c ... Initialize growing season production output variables.

c ... Local variables
      integer mon
      logical grassprod, treeprod

c ... Update production output variables for grass/crop system
      if (system .eq. CRPSYS) then
c ..... If no production has occurred over the past year the growing
c ..... season production output variables get reset to zero,
c ..... initialize the monthly production values for the year,
c ..... cak - 10/02/03
        grassprod = .false.
        do 30 mon = 1, MONTHS
          if ((agcmth(mon).gt.0.001) .or.
     &        ((bgcjmth(mon)+bgcmmth(mon)).gt.0.001)) then
            grassprod = .true.
          endif
30      continue
        if (.not. grassprod) then
          agcprd = 0.0
          bgcjprd = 0.0
          bgcmprd = 0.0
        endif
c ..... Initialize monthly production values
        do 40 mon = 1, MONTHS
          agcmth(mon) = 0.0
          bgcjmth(mon) = 0.0
          bgcmmth(mon) = 0.0
40      continue
        frstmth = month
      endif

c ... Update production output variables for forest system
      if (system .eq. FORSYS) then
c ..... If no production has occurred over the past year the growing
c ..... season production output variables get reset to zero,
c ..... initialize the monthly production values for the year,
c ..... cak - 10/02/03
        treeprod = .false.
        do 80 mon = 1, MONTHS
          if (fcmth(mon) .gt. 0.001) then
            treeprod = .true.
          endif
80      continue
        if (.not. treeprod) then
          rlvprd = 0.0
          frtjprd = 0.0
          frtmprd = 0.0
          fbrprd = 0.0
          rlwprd = 0.0
          crtprd = 0.0
          fcprd = 0.0
        endif
c ..... Initialize monthly production values
        do 90 mon = 1, MONTHS
          fcmth(mon) = 0
90      continue
        tfstmth = month
      endif

      return
      end
