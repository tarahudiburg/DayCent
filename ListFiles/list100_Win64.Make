#
#  Creates "list100_DailyDayCent_muvp.exe" executable for Windows 64-bit OS
#
#
#Source code:	list100.f
#               table.f
#Include files:	outval.inc
#		table.inc

#COMPGCC  = /data/paustian/Century/bin/MinGW/bin/x86_64-w64-mingw32-
#FCMD   = $(COMPGCC)gfortran

/data/paustian/Century/bin/MinGW/bin/x86_64-w64-mingw32-gfortran -u -g -c list100.f
/data/paustian/Century/bin/MinGW/bin/x86_64-w64-mingw32-gfortran -u -g -c table.f
/data/paustian/Century/bin/MinGW/bin/x86_64-w64-mingw32-gfortran -g -o list100_DailyDayCent_muvp.exe list100.o table.o -static

