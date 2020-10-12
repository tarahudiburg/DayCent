      program list100

c ... This routine opens the *.lis file and write the output
c ... records.

c ... VAX NOTE: $ used in format statements does NOT
c ... refrain from carriage return

c ... Modified:
c ...   Modified the output format to provide scientific notation
c ...   for output values that overflow or underflow the format
c ...   specification.  KLM - 04/16/02, implemented by CAK - 08/07/02

      implicit none
      include 'outval.inc'
      include 'table.inc'

      integer found(150,2), i, iargc, in, nargs, numchsn,
     &        out, recnum
      real time, etime, stime
      real outvar
      character*150 ascname, biname

c ... Initialization
      stime = -1
      etime = -1

c ... Determine command line arguments or interactive
      nargs = iargc()
      if (nargs .eq. 0) then
        call interactive(biname,ascname,numchsn,found,stime,etime)
      else
        call getcmdline(nargs,biname,ascname,numchsn,found,stime,etime)
      endif

c ... Assign file unit numbers
      in = 1
      out = 2

c ... Open binary file
      open(unit=in,file=biname,form='UNFORMATTED',status='OLD')

c ... Open ascii file
      open(unit=out,file=ascname,status='NEW')

c ... Write header labels to ascii file
      write(out,10) '  time'
c10    format(a6,$)
10    format(a11,$)
      do 40 i = 1, numchsn
        if (i .eq. 1) then
          write(out,20) table(found(i,1), found(i,2))
c20        format(10x,a10,$)
20        format(10x,a11,$)
        else
          write(out,30) table(found(i,1), found(i,2))
c30        format(5x,a10,$)
30        format(5x,a11,$)
        endif
40    continue
      write(out,'(/)')

c ... Read binary file, writing to ascii 
      recnum = 1
50    continue
      if (MAX .eq. 1) then
        read(unit=in,end=90) time,vals1
      else if (MAX .eq. 2) then
        read(unit=in,end=90) time,vals1,vals2
      else if (MAX .eq. 3) then
        read(unit=in,end=90) time,vals1,vals2,vals3
       else if (MAX .eq. 4) then
         read(unit=in,end=90) time,vals1,vals2,vals3,vals4
c      else if (MAX .eq. 5) then
c        read(unit=in,end=90) time,vals1,vals2,vals3,vals4,vals5
       else
         write(*,*) 'Error in list100.f. MAX parameter not found.'
      endif
c ... If the user didn't specify times or the user specified times 
c ... and it's within those times, write out the desired values
      if ((stime.eq.-1 .and. etime.eq.-1) .or.
     &    (stime.ne.-1 .and. time.ge.stime .and. etime.eq.-1) .or.
     &    (stime.eq.-1 .and. etime.ne.-1 .and. time.le.etime) .or.
     &    (stime.ne.-1 .and. time.ge.stime .and.
     &     etime.ne.-1 .and. time.le.etime)) then
        write(out,60) time
c60      format(f9.2,$)
60      format(f11.2,$)
        do 80 i = 1, numchsn
          if (found(i,1) .eq. 1) then
c            write(out,70) vals1(found(i,2))
            outvar = vals1(found(i,2))
          elseif (found(i,1) .eq. 2) then
c            write(out,70) vals2(found(i,2))
            outvar = vals2(found(i,2))
          elseif (found(i,1) .eq. 3) then
c            write(out,70) vals3(found(i,2))
            outvar = vals3(found(i,2))
          elseif (found(i,1) .eq. 4) then
c            write(out,70) vals4(found(i,2))
            outvar = vals4(found(i,2))
          endif
c70        format(5x,f10.4,$)
c70        format(5x,f11.4,$)
          if (abs(outvar) .lt. 1.0e-9) then
c ......... Zero the displayed output for very small values
            write(out, '(5x,f11.4,$)') 0.0
          else if ((abs(outvar) .ge. 1000000.0) .or.
     &        (abs(outvar) .lt. 0.0001)) then
c ......... Output in scientific notation values that won't display
c ......... in the f11.4 format
            write(out, '(5x,es11.4,$)') outvar
          else
c ......... Standard output format
            write(out, '(5x,f11.4,$)') outvar
          endif
80      continue
c        write(out,*)
        write(out,71)
71      format(a1)
      endif
c ... Update record numbers
      recnum = recnum + 1
c ... If past desired ending time, quit
      if (etime .ne. -1 .and. time .gt. etime) then
        goto 90
      endif
c ... Continue to next time
      goto 50

90    continue

c ... Close files
      close(unit=in)
      close(unit=out)

      write(*,*) 'Done.'

      stop
      end


c ... SUBROUTINE INTERACTIVE
      subroutine interactive(biname,ascname,numchsn,found,stime,etime)
      integer numchsn, found(150,2)
      character*(*) biname, ascname
      real stime, etime
 
      include 'table.inc'

c ... Local variables
      integer i, itemp
      logical existascii
 
c ... Print title
      write(*,*)
      write(*,*) '          DailydayCent_muvps'
      write(*,*) '       Binary to Ascii Utility'
      write(*,*) '             10/07/2018'
 
c ... Obtain name of binary input file
      write(*,*)
      write(*,*) '   Enter name of binary input file (no .bin):'
      read(*,*) biname
 
c ... Check that binary file exists
      call chkbin(biname)
 
c ... Obtain name of ascii output file
20    write(*,*)
      write(*,*) '   Enter name of ASCII output file (no .lis):'
      read(*,*) ascname
 
c ... If ascii file already exists, ask for it again
c      if (existascii(ascname) .eq. .TRUE.) then
      if (existascii(ascname)) then
        goto 20
      endif

c ... Obtain time interval
      write(*,*)
      write(*,*) '   Enter starting time,'
      write(*,*) '   <return> for time file begins:'
      read(*,'(i5)') itemp
      if (itemp .gt. 0) then
        stime = itemp
      endif
      write(*,*) '   Enter ending time,'
      write(*,*) '   <return> for time file ends:'
      read(*,'(i5)') itemp
      if (itemp .gt. 0) then
        etime = itemp
        call fixetime(etime)
      endif
 
c ... Obtain list of variables
      i = 1
      write(*,*)
      write(*,*) '   Enter variables, one per line,'
      write(*,*) '   <return> to quit:'
      call getlist(5,numchsn,found)

      return
      end


c ... SUBROUTINE GETCMDLINE
      subroutine getcmdline(nargs,biname,ascname,numchsn,
     &                      found,stime,etime)
      integer nargs, numchsn, found(150,2)
      character*(*) biname, ascname
      real stime, etime

c ... Local variables
      character*10 xstime, xetime
      character*150 vname
      logical existascii
      integer(2) n

c ... Get name of binary file
      n = 1
      call getarg(n, biname)
c ... Check that binary file exists
      call chkbin(biname)

c ... Get name of ascii file
      n = 2
      call getarg(n, ascname)
c ... Stop if ascii file exists
c      if (existascii(ascname) .eq. .TRUE.) then
      if (existascii(ascname)) then
        STOP
      endif

c ... Get name of input variables file
      n = 3
      call getarg(n, vname)
c ... Check that input variables file exists
      call chkfile(vname)
c ... Get the variables out of the file
      open(unit=9,file=vname)
      call getlist(9,numchsn,found)
      close(unit=9)

c ... Get starting time
      if (nargs .gt. 3) then
        n = 4
        call getarg(n, xstime)
c        read(xstime, '(f10.5)') stime
        read(xstime, *) stime
      endif

c ... Get ending time
      if (nargs .gt. 4) then
        n = 5
        call getarg(n, xetime)
c        read(xetime, '(f10.5)') etime
        read(xetime, *) etime
        call fixetime(etime)
      endif

      return
      end


c ... SUBROUTINE CHKBIN
      subroutine chkbin(biname)
      character*(*) biname

      integer clen, getlen
      logical goahead

      clen = getlen(biname)
      biname = biname(1:clen)
      biname(clen+1:clen+4) = '.bin'
      inquire(file=biname,exist=goahead)
      if (goahead .eqv. .false.) then
        write(*,*)
        write(*,*) '   The binary file could not be read.'
        write(*,*)
        STOP
      endif
         
      return
      end


c ... SUBROUTINE EXISTASCII
      logical function existascii(ascname)
      character*(*) ascname

      integer clen, getlen
      logical goahead

      clen = getlen(ascname)
      ascname = ascname(1:clen)
      ascname(clen+1:clen+4) = '.lis'

      inquire(file=ascname,exist=goahead)
      if (goahead) then
        write(*,*)
        write(*,*) "   ASCII file already exists."
        existascii = .TRUE.
      else
        existascii = .FALSE.
      endif

      return 
      end


c ... SUBROUTINE CHKFILE
      subroutine chkfile(filename)
      character*(*) filename

      logical goahead

      inquire(file=filename,exist=goahead)
      if (goahead .eqv. .false.) then
        write(*,*)
        write(*,*) '   The input variables file could not be read.'
        write(*,*)
        STOP
      endif

      return
      end


c ... SUBROUTINE FIXETIME
      subroutine fixetime(etime)
      real etime

c ... Add very small value to etime b/c of real numbers in binary
      etime = etime + 0.000999

      return
      end


c ... SUBROUTINE GETLIST
      subroutine getlist(unitnum,numchsn,found)
      integer unitnum, numchsn, found(150,2)

      include 'table.inc'

c ... Local variables
      integer i, incnum, search
      character*20 word
      logical isblank

c ... Obtain list of variables from input file
      i = 1
10    read(unitnum,'(a)',end=30) word
      if (.not. isblank(word)) then
        found(i,1) = 1
        found(i,2) = search(word,1)
        do 20 incnum = 2, MAX
          if (found(i,2) .eq. -1) then
            found(i,1) = incnum
            found(i,2) = search(word,incnum)
          endif
20      continue
        if (found(i,2) .eq. -1) then
          write(*,*)
          write(*,*) '   Variable not found, skipping: ', word
        else
          i = i + 1
        endif
        goto 10
      endif
          
30    numchsn = i - 1
      if (numchsn .eq. 0) then
        write(*,*)
        write(*,*) '   No variables chosen; no output file created.'
        STOP
      endif

      return
      end


c ... INTEGER FUNCTION GETLEN
      integer function getlen(name)
      character*(*) name

      integer rlen, max

      max = len(name)
      do 10 rlen = 1, max
        if (name(rlen:rlen) .eq. ' ') goto 20
10    continue
      rlen = 0

20    if (rlen .eq. 0) then
        getlen = max
      else
        getlen = rlen - 1
      endif

      return
      end


c ... FUNCTION ISBLANK
      logical function isblank(line)
      character*(*) line

c ... Local variables
      integer i, numchars
      logical blanko

      blanko = .true.
      numchars = len(line)
      do 10 i = 1, numchars
        if (line(i:i) .ne. ' ') then
          blanko = .false.
          goto 20
        endif
10    continue
20    isblank = blanko

      return
      end


c ... INTEGER FUNCTION SEARCH
      integer function search(word,incnum)
      integer incnum
      character*(*) word

      include 'table.inc'

c ... Local variables
      integer i

      do 10 i = 1, tvals(incnum)
        if (table(incnum,i) .eq. word) then
          search = i
          return
        endif
10    continue
      search = -1

      return
      end
