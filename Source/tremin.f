
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... TREMIN.F

      subroutine tremin(tomatch,curtrm)

      implicit none
      include 'const.inc'
      include 'forrem.inc'

c ... Argument declarations
      character*8 tomatch
      character*8 curtrm

c ... Read in the new tree removal type

c ... Local variables
      integer   ii, TREMLNS
      real      temp
      character fromdat*8, name*20, string*80

c ... Number of lines to read for each tree removal type
c     parameter (TREMLNS = 20)
      parameter (TREMLNS = 38)

      open(unit=11, file='trem.100',status='OLD')
      rewind(11)
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 ii = 1, TREMLNS 
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) temp, name
        evntyp = int(temp)
        call ckdata('tremin','evntyp',name)
        read(11, *) remf(LEAF), name 
        call ckdata('tremin','remf',name)
        read(11, *) remf(FROOT), name 
        call ckdata('tremin','remf',name)
        read(11, *) remf(FBRCH), name 
        call ckdata('tremin','remf',name)
        read(11, *) remf(LWOOD), name 
        call ckdata('tremin','remf',name)
        read(11, *) remf(CROOT), name 
        call ckdata('tremin','remf',name)

c ..... Three new parameters for removal of standing dead 
c ..... tree biomass. -mdh 9/21/2018
        read(11, *) remf(6), name 
        call ckdata('tremin','remf',name)
        read(11, *) remf(7), name 
        call ckdata('tremin','remf',name)
        read(11, *) remf(8), name 
        call ckdata('tremin','remf',name)

        read(11, *) fd(1), name 
        call ckdata('tremin','fd',name)
        read(11, *) fd(2), name 
        call ckdata('tremin','fd',name)
        read(11, *) retf(1,1), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(1,2), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(1,3), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(1,4), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(2,1), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(2,2), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(2,3), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(2,4), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(3,1), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(3,2), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(3,3), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(3,4), name 
        call ckdata('tremin','retf',name)

c ..... Twelve new parameters for C,N,P,S return of litter 
c ..... and ash for standing dead tree biomass. -mdh 9/21/2018
        read(11, *) retf(4,1), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(4,2), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(4,3), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(4,4), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(5,1), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(5,2), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(5,3), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(5,4), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(6,1), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(6,2), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(6,3), name 
        call ckdata('tremin','retf',name)
        read(11, *) retf(6,4), name 
        call ckdata('tremin','retf',name)

c ..... Three new parameters for transfer of live tree biomass
c ..... two standing dead counterparts. -mdh 9/21/2018
        read(11, *) lv2std(1), name 
        call ckdata('tremin','lv2std',name)
        read(11, *) lv2std(2), name 
        call ckdata('tremin','lv2std',name)
        read(11, *) lv2std(3), name 
        call ckdata('tremin','lv2std',name)

        close(11)
        curtrm = tomatch 
      endif

      return

100   format(a5)

200   continue
      call message('  Error reading in values from the trem.100 file.')
      string = '   Looking for type: ' // tomatch
      call message(string)
      STOP

      end
