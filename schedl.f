
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... SCHEDL.F

      subroutine schedl(curday)

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'schvar.inc'
      include 'seq.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Argument declarations
      integer curday

c ... Determine the next set of scheduling options from the 
c ... schedule file

c ... Function declarations
      real      line
      external  line

c ... Local variables
      integer     crtyr
      real        savedfert
      character   string*80
      character*8 curcult, curfert, curfire, curgraz, curharv, 
     &            curirri, curomad, curtrm

c ... Saved variables
      save        curcult, curfert, curfire, curgraz, curharv, 
     &            curirri, curomad, curtrm, savedfert

      data curcult / ' ' /
      data curfert / ' ' /
      data curfire / ' ' /
      data curgraz / ' ' /
      data curharv / ' ' /
      data curirri / ' ' /
      data curomad / ' ' /
      data curtrm / ' ' /
      data savedfert / 0.0 /

c ... Reset do variables to false
      docult = .false.
      doerod = .false.
      dofert = .false.
      doflst = .false.
      dofone = .false.
      dofrst = .false.
      dograz = .false.
      dohrvt = .false.
      doirri = .false.
      dolast = .false.
      doomad = .false.
      doplnt = .false.
      dosene = .false.
      dotrem = .false.
      dofire(CRPSYS) = .false.
      dofire(FORSYS) = .false.
      dofire(SAVSYS) = .false.
      aufert = 0.0
      harmth = 0

c ... Convert time to integer year
      crtyr = aint(time + .001)
      crtyr = mod((crtyr - strtyr + 1), rptyrs)

c ... Working with real numbers - inexact so check on small number -rm
      if (crtyr .lt. 0.1) then
        crtyr = rptyrs
      endif

10    continue

c ... Exit this subroutine when all of the events in the block have
c ... been handled
      if (evtptr .gt. ttlind) then
        if ((timary(evtptr-1,1) .eq. crtyr) .and.
     &      (timary(evtptr-1,2) .eq. curday)) then
          goto 999
        endif
      endif

c ... Determine if evtptr needs to go back to 1
      if (ttlind .ne. 1 .and. evtptr .gt. ttlind) then
        evtptr = 1
      endif

c ... Look for events in timary that match the current time
c ... If found, handle the event
      if ((timary(evtptr,1) .eq. crtyr) .and.
     &     timary(evtptr,2) .eq. curday) then

        if (cmdary(evtptr) .eq. 'CROP') then
          if (curcrp .ne. typary(evtptr)) then
            call cropin(typary(evtptr))
            call co2eff(time)
          endif

        elseif (cmdary(evtptr) .eq. 'PLTM') then
          doplnt = .true.
          plntday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'HARV') then
          dohrvt = .true.
          if (curharv .ne. typary(evtptr)) then
            call harvin(typary(evtptr),curharv)
          endif 
          hrvtday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'FRST') then
          dofrst = .true.
          frstday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'LAST') then
          dolast = .true.
          lastday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'SENM') then
          dosene = .true.
          seneday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'FERT') then
          dofert = .true.
          aufert = savedfert
          if (curfert .ne. typary(evtptr)) then
            call fertin(typary(evtptr),curfert,savedfert)
          endif
          fertday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'CULT') then
          docult = .true.
          if (curcult .ne. typary(evtptr)) then
            call cultin(typary(evtptr),curcult)
          endif
          cultday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'OMAD') then
          doomad = .true.
          if (curomad .ne. typary(evtptr)) then
            call omadin(typary(evtptr),curomad)
          endif
          omadday = timary(evtptr, 2)

        elseif (cmdary(evtptr).eq. 'IRRI') then
          doirri = .true.
          if (curirri .ne. typary(evtptr)) then
            call irrgin(typary(evtptr),curirri)
          endif
          irriday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'GRAZ') then
          dograz = .true.
          if (curgraz .ne. typary(evtptr)) then
            call grazin(typary(evtptr),curgraz)
          endif
          grazday = timary(evtptr, 2)

        elseif (cmdary(evtptr).eq. 'EROD') then
          doerod = .true.
          psloss = fltary(evtptr, 1)
          erodday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'FIRE') then
          dofire(cursys) = .true.
          if (curfire .ne. typary(evtptr)) then
            call firein(typary(evtptr),curfire)
          endif
          fireday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'TREE' .and. 
     &          curtre .ne. typary(evtptr)) then
          call treein(typary(evtptr))
          call co2eff(time)

        elseif (cmdary(evtptr) .eq. 'TREM') then
          dotrem = .true.
          if (curtrm .ne. typary(evtptr)) then
            call tremin(typary(evtptr),curtrm)
          endif
          tremday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'TFST') then
          dofone = .true.
c          forgrw = 1
          foneday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'TLST') then
          doflst = .true.
          flstday = timary(evtptr, 2)

        elseif (cmdary(evtptr) .eq. 'FLOD') then
          watertable = 1
          watrflag = intary(evtptr, 1)

        elseif (cmdary(evtptr) .eq. 'DRAN') then
          watertable = 0
          watrflag = 0

        endif

c ..... Check the next array 'record'
        evtptr = evtptr + 1
        goto 10
      else
        goto 999
      endif

      string = '   Type not found: ' // typary(evtptr)
      call message(string)
      STOP

999   return

      end
