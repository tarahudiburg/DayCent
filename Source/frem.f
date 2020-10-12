
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... FREM

      subroutine frem()

      implicit none
      include 'const.inc'
      include 'forrem.inc'

c ... Forest removal - fire or cutting (cutting includes storms)

c ... Called from:  simsom

c ... Local variables
      real     accum(ISOS)

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Removal for both CUT and FIRE events
      if ((evntyp .eq. 0) .or. (evntyp .eq. 1)) then

c ..... Live Removal
        call livrem(accum)

c ..... Standing Dead Removal. -mdh 9/19/2018
        call stdedrem(accum)

c ..... Live leaves and wood transferred to their dead standing 
c ..... counterparts. -mdh 9/19/2018
        call killiv(accum)

c ..... Death of Roots
        call killrt(accum)

      endif

      if (evntyp .eq. 0) then
c ..... Removal of dead wood occurs only during a cutting event,
c ..... removal of dead material due to burning is now done in
c ..... the grem subroutine during a FIRE event, cak - 01/02
        call dedrem(accum)
c ..... Returns from cutting event
        call cutrtn(accum)

      else if (evntyp .eq. 1) then
c ..... Returns from fire event
        call firrtn()

      endif

      return
      end
