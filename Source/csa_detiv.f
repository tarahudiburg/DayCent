
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine detiv(useprompts, wrtsite)

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'doubles.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'jday.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ... Argument declarations
      logical useprompts, wrtsite

c ... Determine name of schedule file, which contains the
c ... name of the site file, values of timing variables,
c ... and order of events

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE initsw(sitlat, swcinit, usexdrvrs, numlyrs, texture,
     &                    daylength, sradadj, tminslope, tminintercept,
     &                    maxphoto, bioabsorp, ext_flag, swcextend)
          !MS$ATTRIBUTES ALIAS:'_initsw' :: initsw
          REAL             sitlat
          REAL             swcinit(*)
          INTEGER          usexdrvrs
          INTEGER          numlyrs
          INTEGER          texture
          REAL             daylength(*)
          REAL             sradadj(*)
          REAL             tminslope
          REAL             tminintercept
          REAL             maxphoto
          REAL             bioabsorp
          INTEGER          ext_flag
          REAL             swcextend(*)
        END SUBROUTINE initsw

        SUBROUTINE setasmos(asmos, nlayer, swcinit, numlyrs, avh2o, 
     &                      rwcf)
          !MS$ATTRIBUTES ALIAS:'_setasmos' :: setasmos
          REAL    asmos(*)
          INTEGER nlayer
          REAL    swcinit(*)
          INTEGER numlyrs
          REAL    avh2o(*)
          REAL    rwcf(*)
        END SUBROUTINE setasmos

        SUBROUTINE setlyrs(adep,nlayer,numlyrs, sand, silt, clay, 
     &                     bulkd, ph, awilt, afiel, swflag)
          !MS$ATTRIBUTES ALIAS:'_setlyrs' :: setlyrs
          REAL    adep(*)
          INTEGER nlayer
          INTEGER numlyrs
          REAL    sand
          REAL    silt
          REAL    clay
          REAL    bulkd
          REAL    ph
          REAL    awilt(*)
          REAL    afiel(*)
          INTEGER swflag
        END SUBROUTINE setlyrs

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &             ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER          clyr
          REAL             amt
          DOUBLE PRECISION frac_nh4
          DOUBLE PRECISION frac_no3
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          CHARACTER        subname*10
        END SUBROUTINE update_npool

      END INTERFACE

c ... Function declarations
      integer iargc

      common /libpath/filpath
      character*100 filpath
      character*100 sitnam
      character*10 subname

c ... Local variables
      integer          ii, nargs, ext_flag
      character*100    newbin, oldbin
      logical          ext, ext_grid, ext_site
      character*80     string
      integer          iel, numlyrs
      real             swcinit(SWMAXLYR)
      real             swcextend(SWMAXLYR)

c ... Add new routine to do a "brute force" initialization of all common block
c ... variables, cak - 06/04/02
      call default

      subname = 'detiv     '

c ... Initialize potential command line arguments
      ext = .false.
      ext_grid = .false.
      ext_site = .false.
      schnam = ' '
      newbin = ' '
      oldbin = ' '
      filpath = ' '
      useprompts = .false.

c ... Determine command line arguments or interactive
      nargs = iargc()
      if (nargs .eq. 0) then
        useprompts = .true.
        call interactive(schnam, newbin, ext, oldbin, wrtsite,
     &                   extsiteout, ext_site, extsitein)
      else
        useprompts = .false.
        call getcmdline(nargs, schnam, newbin, ext, oldbin, ext_grid,
     &                  wrtsite, extsiteout, ext_site, extsitein,
     &                  filpath)
      endif

c ... Check that minimal information was entered
      if (schnam .eq. ' ') then
        call message(' ')
        call message('   No schedule file name was given.')
        call message(' ')
        STOP 'Execution error.'
      endif

      if (newbin .eq. ' ') then
        if (ext) then
          newbin = oldbin
        else
          call message(' ')
          call message('   No binary output file name was given.')
          call message(' ')
          STOP 'Execution error.'
        endif
      endif

c ... Open binary file to write to
      if (ext .and. newbin .eq. oldbin) then
        open(unit=1,file=newbin,form='UNFORMATTED',status='OLD')
      else
        open(unit=1,file=newbin,form='UNFORMATTED',status='NEW')
      endif

c ... Open the schedule file and read the header lines
      open(unit=15,file=schnam,status='OLD')

c ... Allow comments at the top of a schedule file for documentation
c ... purposes.  All of the comment lines have a # character at the
c ... start and are stored at the top of the schedule file and are 
c ... ignored here.  There is no blank line permitted between the last
c ... comment line and the line of the schedule file containing the
c ... start year information.  cak - 12/30/02
      read(15,100) string
100   format(a80)
      do while (string(1:1) .eq. '#') 
        read(15,100) string
      end do

c      read(15,*) strtyr
      read(string,*) strtyr

      read(15,*) tend
      tend = tend + 1

      read(15,*) sitnam
      if (ext_site) then
        sitnam = extsitein
      endif

      read(15,*) labtyp
      read(15,*) labyr

      read(15,*) mctemp
      micosm = 0
      if (mctemp .ge. 0) then
        call message(' ')
        call message('Microcosms are not implemented in this version')
        call message('of Daily Century')
        micosm = 1
        STOP
      endif

      read(15,*) co2sys
      if (co2sys .gt. 0) then
        read(15,*) co2tm(1), co2tm(2)
      endif

c ... New header line in schedule file to handle pH shift, cak - 08/02/02
c ... Change pH shift implementation to use scalar values, cak - 10/17/05
      read(15,*) phsys
      if (phsys .gt. 0) then
        read(15,*) phtm
      endif

c ... New header lines in schedule file to handle soil temperature warming
c ... experiments, cak - 07/02/03
      read(15,*) stsys
      if (stsys .gt. 0) then
        read(15,*) ststart
        read(15,*) stamt
      endif

c ... New header lines in schedule file to handle the N input scalars,
c ... cak - 04/06/04
      read(15,*) Ninput
      if (Ninput .gt. 0) then
        read(15,*) Nstart
      endif

c ... New header lines in schedule file to handle the OMAD input scalars,
c ... cak - 04/06/04
      read(15,*) OMADinput
      if (OMADinput .gt. 0) then
        read(15,*) OMADstart
      endif

c ... New header lines in schedule file to handle the weather input scalars,
c ... cak - 10/18/05
      read(15,*) wthinput
      if (wthinput .gt. 0) then
        read(15,*) wthstart
      endif

      read(15,*) decsys
      if (decsys .eq. SAVSYS) then
        decsys = FORSYS
      endif
      read(15,40) initcp
40    format(a5)
      if (initcp .eq. 'Initi') then
        initcp = ' '
      endif
      read(15,40) initre
      if (initre .eq. 'Initi') then
        initre = ' '
      endif

      read(15,*)
      read(15,*)

c ... Read starting values from fixed parameter file
      call fixin

c ... Read starting values from site-specific file
      open(unit=7,file=sitnam,status='OLD',err=1000)
      if (ext_grid) then
        call sitein_grid(ext)
      elseif (ext_site) then
        call sitein_ext(swcextend)
      else
        call sitein(ext)
      endif

c ... Moved the read calls for the initial tree and crop to inside the extend
c ... if statement.  This is done to prevent a rather subtle bug that occurs
c ... when the initial crop/tree do not match the final values in the 
c ... original schedule.  In that case, the derived output values, (crpval ...)
c ... do not match the current crop values.
c ... The crop/tree reads must occur before the calciv call on a normal run.
c ... 7/20/95  K. Killian

c ... Determine initial values
      if (ext) then
        if (oldbin .ne. newbin) then
          open(unit=3,file=oldbin,form='UNFORMATTED',status='OLD')
          call extend(3,.TRUE.)
          close(unit=3)
        else
          call extend(1,.FALSE.)
        endif
c ..... Save the single precision secndy and occlud variables read into their
c ..... double precision counterparts, cak - 09/29/2005
        secndy_double(1) = secndy(1)
        secndy_double(2) = secndy(2)
        secndy_double(3) = secndy(3)
        occlud_double = occlud
      endif

c ... Add new routine to initialize common block variables to other than
c ... default values as necessary
      call initialize(ext)

c ... Added call to initsw for Daily water budget version of Century
c ... -mdh 9/94
      if (ext_site) then
        ext_flag = 1
      else
        ext_flag = 0
      endif
      call initsw(sitlat, swcinit, usexdrvrs, numlyrs, texture,
     &            daylength, sradadj, tminslope, tminintercept,
     &            maxphoto, bioabsorp, ext_flag, swcextend)

c ... Initialize soil properties based on structure of Daily Soil Water
c ... Model -mdh 10/96
      call setlyrs(adep,nlayer,numlyrs, sand, silt, clay, bulkd, ph, 
     &             awilt, afiel, swflag)

c ... Set the initial pH value based on the value returned from the setlyrs
c ... subroutine, cak - 08/02/02
      phstart = ph

      if (.not. ext_site) then
        ammonium = 0.0
        do 45 ii = 1,SWMAXLYR
          nitrate(ii) = 0.0
45      continue
c ..... Initialize the layer beyond the last one the used for safety
c ..... Zero out ALL layers below nlayer. -mdh 7/27/01
        do 19 iel = 1, 3
          do 18 ii = nlayer+1, MAXLYR
            minerl(ii, iel) = 0.0
18        continue
19      continue
 
c ..... When initializing the simulation all of the initial mineral N
c ..... is assumed to be nitrate
        frac_nh4_fert = 0.0
        frac_no3_fert = 1.0
        do 110 ii=1,nlayer
          if (minerl(ii,N) .lt. 0.05) then
            minerl(ii,N) = 0.1
          endif
          call update_npool(ii, minerl(ii,N), frac_nh4_fert,
     &                      frac_no3_fert, ammonium, nitrate, subname)
110     continue
        frac_nh4_fert = 0.5
        frac_no3_fert = 0.5
      endif

c ... Obtain the initial values for the crop or forest system
c ... Initialize the fine root pools based on the information read from
c ... the site file and the crop/tree parameterization, cak - 05/24/2007
      cursys = 0
      if (initcp .ne. ' ') then
        call cropin(initcp)
        cursys = CRPSYS
        if (.not. ext .and. .not. ext_site) then
          bglcism(UNLABL) = bglcisj(UNLABL) * mrtfrac
          bglcism(LABELD) = bglcisj(LABELD) * mrtfrac
          bglcisj(UNLABL) = bglcisj(UNLABL) * (1.0 - mrtfrac)
          bglcisj(LABELD) = bglcisj(LABELD) * (1.0 - mrtfrac)
          do 120 ii = 1, MAXIEL
            bglivem(ii) = bglivej(ii) * mrtfrac
            bglivej(ii) = bglivej(ii) * (1.0 - mrtfrac)
120       continue
        endif
      endif
      if (initre .ne. ' ') then
        call treein(initre)
        cursys = cursys + FORSYS
        if (.not. ext .and. .not. ext_site) then
          frtcism(UNLABL) = frtcisj(UNLABL) * wmrtfrac
          frtcism(LABELD) = frtcisj(LABELD) * wmrtfrac
          frtcisj(UNLABL) = frtcisj(UNLABL) * (1.0 - wmrtfrac)
          frtcisj(LABELD) = frtcisj(LABELD) * (1.0 - wmrtfrac)
          do 130 ii = 1, MAXIEL
            frootem(ii) = frootej(ii) * wmrtfrac
            frootej(ii) = frootej(ii) * (1.0 - wmrtfrac)
130       continue
        endif
      endif
      if (.not. ext .and. .not. ext_site) then
        call calciv
      endif

c ... Sum up isotopes
      call sumcar

c ... Do preliminary initializations and calculations
      call prelim

c ... Initialize initial soil moisture (asmos) based on structure of Daily
c ... Soil Water Model -mdh 10/96
      call setasmos(asmos, nlayer, swcinit, numlyrs, avh2o, rwcf)

c ... Read the first block of events
      call readblk

      call message(' ')
      call message('   Model is running...')

      return
1000    call message(' Fatal error: unknown site file :'//sitnam)
        stop ' Abnormal Termination'
      end


      integer function getlen(name)

      implicit none
      character*(*) name
      integer jj

C -----------------------------------------------------------------------------
C     this subroutine left justifies the file name and determines the length
C
C Variables
C      Input
C   name    character (*)  the input and processed file name
C
C  Modified by K. Killian 8/11/94
C              included the left justification on a subroutine coded by Laura
C
C -----------------------------------------------------------------------------
 
15    getlen = index(name,' ')-1

      if (getlen .eq. -1) then
        getlen = len(name)
      else if (getlen .eq. 0) then
        do 20 jj= 1,len(name)
          if (name(jj:jj) .ne. ' ') then
            name = name(jj:)
            goto 15
          endif
20      continue
        getlen = 0
      endif

      return
      end


c ... SUBROUTINE INTERACTIVE
      subroutine interactive(schnam, newbin, ext, oldbin, wrtsite,
     &                       extsiteout, ext_site, extsitein)

      implicit none
      character*100 newbin, oldbin, schnam, extsiteout, extsitein
      logical       ext, ext_site, wrtsite

c ... Function declarations
      integer getlen
 
c ... Local variables
      integer clen
      logical goahead

c ... Print title
      call message(' ')
      call message('              DAYCENT SOIL ORGANIC MATTER')
      call message('                 AND TRACEGAS MODEL')
      call message('                 with Photosynthesis,')
      call message('                 UV litter degradation,')
      call message('                 and methanogenesis')
      call message('                  DailyDayCent_muvp')
      call message('                     09/15/2017')

c ... Obtain name of the schedule file
      call message(' ')
      call message('   Enter name of schedule file (no .sch):')
      read(*,*) schnam
      clen = getlen(schnam)
20    schnam = schnam(1:clen)
      if (index(schnam,'.sch').eq.0) then
        schnam(clen+1:clen+4) = '.sch'
      endif
      inquire(file=schnam,exist=goahead)
c ... If the schedule file does not exist, ask for it again
      if (.not. goahead) then
        call message(' ')
        call message(' The schedule file could not be read.')
        call message(' Re-enter schedule file name or Q to quit):')
        read(*,*) schnam
        clen = getlen(schnam)
        if ((schnam .eq. 'q' .or. schnam .eq. 'Q') .and.
     &      clen .eq. 1) then
          call message(' ')
          STOP ' Abnormal termination.'
        else
          goto 20
        endif
      endif

c ... Obtain name of the binary output file
      call message(' ')
      call message('   Enter name of binary output file (no .bin):')
      read(*,*) newbin
      clen = getlen(newbin)
30    newbin = newbin(1:clen)
      if (index(newbin,'.bin').eq.0) then
        newbin(clen+1:clen+4) = '.bin'
      endif
      inquire(file=newbin,exist=goahead)
c ... If the binary file already exists, ask for it again
      if (goahead) then
        call message(' ')
        call message(' Binary file already exists.')
        call message(' Re-enter binary file name or Q to quit):')
        read(*,*) newbin
        clen = getlen(newbin)
        if ((newbin .eq. 'q' .or. newbin .eq. 'Q') .and.
     &      clen .eq. 1) then
          call message(' ')
          STOP ' Abnormal termination.'
        else
          goto 30
        endif
      endif

c ... Check if the user wishes to create an extended site file that can
c ... be used to initialize a subsequent simulation
      call message(' ')
      call message('   Do you want to create an extended site file?')
      call message('   This file can be used to initialize a')
      call message('   subsequent simulation.')
      call message('   Enter Y for Yes or N for No:')
      read(*,*) extsiteout
      clen = getlen(extsiteout)
      if ((extsiteout .eq. 'y' .or. extsiteout .eq. 'Y') .and.
     &    clen .eq. 1) then
        wrtsite = .true.
        call message(' ')
        call message('   Enter name of site file to create (no .100):')
        read(*,*) extsiteout
        clen = getlen(extsiteout)
50      extsiteout = extsiteout(1:clen)
        if (index(extsiteout,'.100').eq.0) then
          extsiteout(clen+1:clen+4) = '.100'
        endif
        inquire(file=extsiteout,exist=goahead)
c ..... If the extended site file already exists ask for it again
        if (goahead) then
          call message(' ')
          call message('   The extended site file already exists.')
          call message(' Reenter extend site file name or Q to quit):')
          read(*,*) extsiteout
          clen = getlen(extsiteout)
          if ((extsiteout .eq. 'q' .or. extsiteout .eq. 'Q') .and.
     &        clen .eq. 1) then
            call message(' ')
            STOP ' Abnormal termination.'
          else
            goto 50
          endif
        endif
      else
        wrtsite = .false.
      endif

c ... Check for extend simulation
      call message(' ')
      call message('   Is this an extend simulation?')
      call message('   Enter Y for Yes or N for No:')
      read(*,*) oldbin
      clen = getlen(oldbin)
      if ((oldbin .eq. 'y' .or. oldbin .eq. 'Y') .and.
     &    clen .eq. 1) then
c ..... Will the extend use a binary file or a site file
        call message(' ')
        call message('   Will you use a binary file or a site file')
        call message('   for the extend?')
70      call message('   Enter B for binary file or S for site file:')
        read(*,*) oldbin
        clen = getlen(oldbin)
c ..... Extend from an existing binary file
        if ((oldbin .eq. 'b' .or. oldbin .eq. 'B') .and.
     &      clen .eq. 1) then
          ext = .true.
          call message(' ')
          call message('   Enter name of old binary file (no .bin):')
          read(*,*) oldbin
          clen = getlen(oldbin)
40        oldbin = oldbin(1:clen)
          if (index(oldbin,'.bin').eq.0) then
            oldbin(clen+1:clen+4) = '.bin'
          endif
          inquire(file=oldbin,exist=goahead)
c ....... If the binary file does not exist, ask for it again
          if (.not. goahead) then
            call message(' ')
            call message(' Could not find extend binary file.')
            call message(' Re-enter extend binary file or Q to quit):')
            read(*,*) oldbin
            clen = getlen(oldbin)
            if ((oldbin .eq. 'q' .or. oldbin .eq. 'Q') .and.
     &          clen .eq. 1) then
              call message(' ')
              STOP ' Abnormal termination.'
            else
              goto 40
            endif
          endif
c ..... Extend from an existing extended site file
        elseif ((oldbin .eq. 's' .or. oldbin .eq. 'S') .and.
     &      clen .eq. 1) then
          ext_site = .true.
          call message(' ')
          call message('   Enter name of old site file (no .100):')
          read(*,*) extsitein
          clen = getlen(extsitein)
60        extsitein = extsitein(1:clen)
          if (index(extsitein,'.100').eq.0) then
            extsitein(clen+1:clen+4) = '.100'
          endif
          inquire(file=extsitein,exist=goahead)
          if (.not. goahead) then
            call message(' ')
            call message('   The extended site file could not be read.')
            call message(' Re-enter extended site file or Q to quit):')
            read(*,*) extsitein
            clen = getlen(extsitein)
            if ((extsitein .eq. 'q' .or. extsitein .eq. 'Q') .and.
     &          clen .eq. 1) then
              call message(' ')
              STOP ' Abnormal termination.'
            else
              goto 60
            endif
          endif
        else
          call message('   Invalid option.')
          goto 70
        endif
      else
        ext = .false.
        ext_site = .false.
      endif

      return
      end


c ... SUBROUTINE GETCMDLINE
      subroutine getcmdline(nargs, schnam, newbin, ext, oldbin,
     &                      ext_grid, wrtsite, extsiteout, ext_site,
     &                      extsitein, filpath)

      implicit none
      integer       nargs
      character*100 extflag, newbin, oldbin, schnam
      character*100 extsiteout, extsitein, filpath
      logical       ext, ext_grid, ext_site, wrtsite

c ... Function declarations
      integer getlen
 
c ... Local variables
      integer       clen, ii
      character*100 iname
      logical       goahead

c ... Process command line arguments
c ... Add a new argument to indicate reading an extended <site>.100 file
c ... which was created from a site.nc file from a Gridded Century run
c ... cak - 10/05/01
      ii = 1
10    if (ii .lt. nargs) then
        call getarg(ii, extflag)
        ii = ii + 1
        call getarg(ii, iname)
        ii = ii + 1
        clen = getlen(iname)

        if (extflag .eq. '-s') then
          schnam = iname(1:clen)
          if (index(schnam,'.sch').eq.0) schnam(clen+1:clen+4) = '.sch'
          inquire(file=schnam,exist=goahead)
          if (.not. goahead) then
            call message(' ')
            call message('   The schedule file could not be read.')
            call message(' ')
            STOP 'Execution error.'
          endif

        else if (extflag .eq. '-n') then
          newbin = iname(1:clen)
          if (index(newbin,'.bin').eq.0) newbin(clen+1:clen+4) = '.bin'
          inquire(file=newbin,exist=goahead)
          if (goahead) then
            call message(' ')
            call message('   The new binary file already exists.')
            call message(' ')
            STOP 'Execution error.'
          endif

        else if (extflag .eq. '-e') then
          if (ext_site) then
            call message(' ')
            call message('   Can not extend from both a binary and')
            call message('   a site file.')
            call message(' ')
            STOP 'Execution error.'
          else
            ext = .true.
            oldbin = iname
            if (index(oldbin,'.bin').eq.0) then
              oldbin(clen+1:clen+4) = '.bin'
            endif
            inquire(file=oldbin,exist=goahead)
            if (.not. goahead) then
              call message(' ')
              call message('   The old binary file could not be read.')
              call message(' ')
              STOP 'Execution error.'
            endif
          endif

        else if (extflag .eq. '-g') then
          ext_grid = .true.
          schnam = iname(1:clen)
          if (index(schnam,'.sch').eq.0) schnam(clen+1:clen+4) = '.sch'
          inquire(file=schnam,exist=goahead)
          if (.not. goahead) then
            call message(' ')
            call message('   The schedule file could not be read.')
            call message(' ')
            STOP 'Execution error.'
          endif

        else if (extflag .eq. '-c') then
          wrtsite = .true.
          extsiteout = iname(1:clen)
          if (index(extsiteout,'.100').eq.0)
     &      extsiteout(clen+1:clen+4) = '.100'
          inquire(file=extsiteout,exist=goahead)
          if (goahead) then
            call message(' ')
            call message('   The extended site file already exists.')
            call message(' ')
            STOP 'Execution error.'
          endif

        else if (extflag .eq. '-x') then
          if (ext) then
            call message(' ')
            call message('   Can not extend from both a site and')
            call message('   a binary file.')
            call message(' ')
            STOP 'Execution error.'
          else
            ext_site = .true.
            extsitein = iname(1:clen)
            if (index(extsitein,'.100').eq.0)
     &        extsitein(clen+1:clen+4) = '.100'
            inquire(file=extsitein,exist=goahead)
            if (.not. goahead) then
              call message(' ')
              call message('   The extended site file could not')
              call message('   be read.')
              call message(' ')
              STOP 'Execution error.'
            endif
          endif

        else if (extflag .eq. '-l') then
          filpath = iname(1:clen)
        else
          call message('   Unknown argument skipped.')
        endif
        goto 10
      endif

      return
      end
