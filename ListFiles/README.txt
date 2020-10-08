To build the list100_DayCent_Photosyn_UV executable for Linux:

chmod 755 list100.Make
./list100.Make

To build the list100_DailyDayCent_muvp.exe executable for Windows 64-bit:
You will need update makefile_Win64.txt with the path to the ming32-gfortran compile installed on your system;
at NREL on rubel it is /data/paustian/Century/bin/MinGW/bin/x86_64-w64-mingw32-gfortran 

chmod 755 list100_Win64.Make
./list100_Win64.Make
