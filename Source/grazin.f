
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... GRAZIN.F

      subroutine grazin(tomatch,curgraz)

      implicit none
      include 'const.inc'
      include 'parcp.inc'

c ... Argument declarations
      character*8 tomatch
      character*8 curgraz

c ... Read in the new graze type

c ... Local variables
      integer   ii, GRAZLNS
      real      temp
      character fromdat*5, name*20, string*80

c ... Number of lines to read for each graze type
      parameter (GRAZLNS = 11)

      open(unit=11, file='graz.100',status='OLD')
      rewind(11)
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 ii = 1, GRAZLNS 
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) flgrem, name 
        call ckdata('grazin','flgrem',name)
        read(11, *) fdgrem, name 
        call ckdata('grazin','fdgrem',name)
        read(11, *) gfcret, name 
        call ckdata('grazin','gfcret',name)
        read(11, *) gret(N), name 
        call ckdata('grazin','gret',name)
        read(11, *) gret(P), name 
        call ckdata('grazin','gret',name)
        read(11, *) gret(S), name 
        call ckdata('grazin','gret',name)
        read(11, *) temp, name
        grzeff = int(temp)
        call ckdata('grazin','grzeff',name)
        read(11, *) fecf(N), name 
        call ckdata('grazin','fecf',name)
        read(11, *) fecf(P), name 
        call ckdata('grazin','fecf',name)
        read(11, *) fecf(S), name 
        call ckdata('grazin','fecf',name)
        read(11, *) feclig, name 
        call ckdata('grazin','feclig',name)
        close(11)
        curgraz = tomatch 
      endif

      return

100   format(a5)

200   continue
      call message('  Error reading in values from the graz.100 file.')
      string = '   Looking for type: ' // tomatch
      call message(string)
      STOP

      end
