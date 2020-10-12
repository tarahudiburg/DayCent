
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine ckdata(routin,expect,found)

      implicit none

c ... Argument declarations
      character*(*) routin, expect, found

c ... Check a parameter name read by fixin or sitein.
c ...   routin is the name of the calling routine.
c ...   expect is the name that should have been read.
c ...   found  is the name that was read from the data file.

c ... Local variables
      integer ii, loc
      character par*8, string*80

c ... Extract junk from the parameter name

      loc = index(found,'(')
      if (loc .eq. 0) then
        loc = index(found,'|')
        if (loc .eq. 0) then
          loc = index(found,',')
          if (loc .eq. 0) then
            loc = index(found,' ')
c           write(*,*) 'loc =', loc
          endif
        endif
      endif
      if (loc .eq. 0) then
        par = found( :8)
      else
        par = found( :loc-1)
      endif

c ... Convert to lower case if needed (MACHINE DEPENDANT)
      do 100 ii = 1, 8
        if (par(ii:ii) .ge. 'A' .and. par(ii:ii) .le. 'Z') then
          par(ii:ii) = char( ichar(par(ii:ii)) + 32 )
        endif
100   continue

c ... Test name for expected name
      if (expect .ne. par) then
        if (routin .eq. 'fixin') then
          call message('   There is an error in your fix.100 file.')
        else if (routin .eq. 'sitein') then
          call message('   There is an error in your <site>.100 file.')
        else if (routin .eq. 'cropin') then
          call message('   There is an error in your crop.100 file.')
        else if (routin .eq. 'cultin') then
          call message('   There is an error in your cult.100 file.')
        else if (routin .eq. 'fertin') then
          call message('   There is an error in your fert.100 file.')
        else if (routin .eq. 'firein') then
          call message('   There is an error in your fire.100 file.')
        else if (routin .eq. 'harvin') then
          call message('   There is an error in your harv.100 file.')
        else if (routin .eq. 'grazin') then
          call message('   There is an error in your graz.100 file.')
        else if (routin .eq. 'irrgin') then
          call message('   There is an error in your irri.100 file.')
        else if (routin .eq. 'omadin') then
          call message('   There is an error in your omad.100 file.')
        else if (routin .eq. 'treein') then
          call message('   There is an error in your tree.100 file.')
        else if (routin .eq. 'tremin') then
          call message('   There is an error in your trem.100 file.')
        endif
        string = '   The data for ' // par // ' was read when the ' //
     &           'data for ' // expect // ' was expected.'
        call message(string)
        STOP
      endif

      return
      end
