#
#  Creates "list100_DailyDayCent_muvps" executable for Linux
#
#
#Source code:	list100.f
#               table.f
#Include files:	outval.inc
#		table.inc
gfortran -u -g -c list100.f
gfortran -u -g -c table.f
gfortran -g -o list100_DailyDayCent_muvps list100.o table.o

