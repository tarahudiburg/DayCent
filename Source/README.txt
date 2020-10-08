I converted the original Windows 32-bit code in PC_Orig_Dec_2014/ to code that can be
compiled for both Linux and Windows-64 bit.  I have not tested this code to see how
the results compare to the Windows-32 bit version. -mdh 9/15/2017
====================================================================================================

To build the DailyDayCent_muvp executable for Linux using gfortran:

make -f Makefile.txt DailyDayCent_muvp


To build the DailyDayCent_muvp.exe executable for Windows 64-bit:
You will need update makefile_Win64.txt with the path to the ming32-gfortran compile installed on 
your system; at NREL on rubel it is /data/paustian/Century/bin/MinGW/bin/x86_64-w64-mingw32-gfortran 

make -f Makefile_Win64.txt DailyDayCent_muvp.exe

====================================================================================================
Code updates:
Sept. 15, 2017

All *.c *.f *.h *.inc files must begin with small case letters to be consistent with other
DayCent Linux versions.

----------------------------------------------------------------------------------------------------
The soil surface temperature function (surftemp.c) was replaced since the previous version 
(surftemp_old.c) had a bug when calculating soil surface temperature under the snow. -mdh 9/15/2017

----------------------------------------------------------------------------------------------------
The file catanf.f was renamed to carctanf.f.  Also, calls to catanf wre renamed to carctanf in the 
following files:

carctanf.f
droot.f
litdec.f
prelim.f
pschem.f
somdec.f
tcalc.f
wdeath.f
woodec.f

We discovered a couple years ago that catanf function is an intrinisic gfortran compiler function
of some sort, and there is a conflict between this internal compiler function and DayCent's catanf
function.  Without this name substitution DayCent will generate randomly large numbers when the
function is called. -mdh 9/15/2017

====================================================================================================

Nov. 7, 2017

The indices in the agcmth, bgcjmth, and bgcmmth accumulators in csa_main.f were going
out of bounds when month=1.  I trapped for this condition.

Fixed all floating point comparisons to 0.0.  See file flteq0.txt.

Within calcdfac.c, tcalc (a FORTRAN function) was being called as a real function.
This worked for the old MSDOS FORTRAN compilier, but was causing a segmentation
violation in the executable compiled with gfortran.  I changed tcalc to a 
FORTRAN subroutine.

====================================================================================================

Nov. 8, 2017

There seemed to be an underflow problem with rleavc in leafa.f.  When rleavc went
slightly negative, DayCent put out and error message and quit.  This stopped the 
model from ever recovering from this underrflow condition. I trapped for the 
(rleavc < 0.0) condtion in treegrow, right before the call to leafa, and 
reset rleavc and its isotope components (rlvcis(*) to 0.0.

====================================================================================================
.
