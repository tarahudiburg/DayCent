
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine wrtextsite()

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'ligvar.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'wth.inc'

c ... Write an "extended" site file that can be used as initialization
c ... values for a new simulation

c ... Fortran to C prototype
      INTERFACE
        SUBROUTINE getswc(swcextend, numlyrs)
          !MS$ATTRIBUTES ALIAS:'_getswc' :: getswc
          REAL             swcextend(*)
          INTEGER          numlyrs
        END SUBROUTINE getswc
      END INTERFACE

c ... Local variables
      real          swcextend(SWMAXLYR)
      character*20  charval

      open(unit=100,file=extsiteout)

      write(100,15) 'EXSIT', 'Archived_site_file_for:', schnam

      write(100,25) '*** Climate parameters'
      call leftjust(precip(1), charval)
      write(100,35) charval, "'PRECIP(1)'"
      call leftjust(precip(2), charval)
      write(100,35) charval, "'PRECIP(2)'"
      call leftjust(precip(3), charval)
      write(100,35) charval, "'PRECIP(3)'"
      call leftjust(precip(4), charval)
      write(100,35) charval, "'PRECIP(4)'"
      call leftjust(precip(5), charval)
      write(100,35) charval, "'PRECIP(5)'"
      call leftjust(precip(6), charval)
      write(100,35) charval, "'PRECIP(6)'"
      call leftjust(precip(7), charval)
      write(100,35) charval, "'PRECIP(7)'"
      call leftjust(precip(8), charval)
      write(100,35) charval, "'PRECIP(8)'"
      call leftjust(precip(9), charval)
      write(100,35) charval, "'PRECIP(9)'"
      call leftjust(precip(10), charval)
      write(100,35) charval, "'PRECIP(10)'"
      call leftjust(precip(11), charval)
      write(100,35) charval, "'PRECIP(11)'"
      call leftjust(precip(12), charval)
      write(100,35) charval, "'PRECIP(12)'"
      call leftjust(prcstd(1), charval)
      write(100,35) charval, "'PRCSTD(1)'"
      call leftjust(prcstd(2), charval)
      write(100,35) charval, "'PRCSTD(2)'"
      call leftjust(prcstd(3), charval)
      write(100,35) charval, "'PRCSTD(3)'"
      call leftjust(prcstd(4), charval)
      write(100,35) charval, "'PRCSTD(4)'"
      call leftjust(prcstd(5), charval)
      write(100,35) charval, "'PRCSTD(5)'"
      call leftjust(prcstd(6), charval)
      write(100,35) charval, "'PRCSTD(6)'"
      call leftjust(prcstd(7), charval)
      write(100,35) charval, "'PRCSTD(7)'"
      call leftjust(prcstd(8), charval)
      write(100,35) charval, "'PRCSTD(8)'"
      call leftjust(prcstd(9), charval)
      write(100,35) charval, "'PRCSTD(9)'"
      call leftjust(prcstd(10), charval)
      write(100,35) charval, "'PRCSTD(10)'"
      call leftjust(prcstd(11), charval)
      write(100,35) charval, "'PRCSTD(11)'"
      call leftjust(prcstd(12), charval)
      write(100,35) charval, "'PRCSTD(12)'"
      call leftjust(prcskw(1), charval)
      write(100,35) charval, "'PRCSKW(1)'"
      call leftjust(prcskw(2), charval)
      write(100,35) charval, "'PRCSKW(2)'"
      call leftjust(prcskw(3), charval)
      write(100,35) charval, "'PRCSKW(3)'"
      call leftjust(prcskw(4), charval)
      write(100,35) charval, "'PRCSKW(4)'"
      call leftjust(prcskw(5), charval)
      write(100,35) charval, "'PRCSKW(5)'"
      call leftjust(prcskw(6), charval)
      write(100,35) charval, "'PRCSKW(6)'"
      call leftjust(prcskw(7), charval)
      write(100,35) charval, "'PRCSKW(7)'"
      call leftjust(prcskw(8), charval)
      write(100,35) charval, "'PRCSKW(8)'"
      call leftjust(prcskw(9), charval)
      write(100,35) charval, "'PRCSKW(9)'"
      call leftjust(prcskw(10), charval)
      write(100,35) charval, "'PRCSKW(10)'"
      call leftjust(prcskw(11), charval)
      write(100,35) charval, "'PRCSKW(11)'"
      call leftjust(prcskw(12), charval)
      write(100,35) charval, "'PRCSKW(12)'"
      call leftjust(tmn2m(1), charval)
      write(100,35) charval, "'TMN2M(1)'"
      call leftjust(tmn2m(2), charval)
      write(100,35) charval, "'TMN2M(2)'"
      call leftjust(tmn2m(3), charval)
      write(100,35) charval, "'TMN2M(3)'"
      call leftjust(tmn2m(4), charval)
      write(100,35) charval, "'TMN2M(4)'"
      call leftjust(tmn2m(5), charval)
      write(100,35) charval, "'TMN2M(5)'"
      call leftjust(tmn2m(6), charval)
      write(100,35) charval, "'TMN2M(6)'"
      call leftjust(tmn2m(7), charval)
      write(100,35) charval, "'TMN2M(7)'"
      call leftjust(tmn2m(8), charval)
      write(100,35) charval, "'TMN2M(8)'"
      call leftjust(tmn2m(9), charval)
      write(100,35) charval, "'TMN2M(9)'"
      call leftjust(tmn2m(10), charval)
      write(100,35) charval, "'TMN2M(10)'"
      call leftjust(tmn2m(11), charval)
      write(100,35) charval, "'TMN2M(11)'"
      call leftjust(tmn2m(12), charval)
      write(100,35) charval, "'TMN2M(12)'"
      call leftjust(tmx2m(1), charval)
      write(100,35) charval, "'TMX2M(1)'"
      call leftjust(tmx2m(2), charval)
      write(100,35) charval, "'TMX2M(2)'"
      call leftjust(tmx2m(3), charval)
      write(100,35) charval, "'TMX2M(3)'"
      call leftjust(tmx2m(4), charval)
      write(100,35) charval, "'TMX2M(4)'"
      call leftjust(tmx2m(5), charval)
      write(100,35) charval, "'TMX2M(5)'"
      call leftjust(tmx2m(6), charval)
      write(100,35) charval, "'TMX2M(6)'"
      call leftjust(tmx2m(7), charval)
      write(100,35) charval, "'TMX2M(7)'"
      call leftjust(tmx2m(8), charval)
      write(100,35) charval, "'TMX2M(8)'"
      call leftjust(tmx2m(9), charval)
      write(100,35) charval, "'TMX2M(9)'"
      call leftjust(tmx2m(10), charval)
      write(100,35) charval, "'TMX2M(10)'"
      call leftjust(tmx2m(11), charval)
      write(100,35) charval, "'TMX2M(11)'"
      call leftjust(tmx2m(12), charval)
      write(100,35) charval, "'TMX2M(12)'"

      write(100,25) '*** Site and control parameters'
      call leftjust(float(ivauto), charval)
      write(100,35) charval, "'IVAUTO'"
      call leftjust(float(nelem), charval)
      write(100,35) charval, "'NELEM'"
      call leftjust(sitlat, charval)
      write(100,35) charval, "'SITLAT'"
      call leftjust(sitlng, charval)
      write(100,35) charval, "'SITLNG'"
      call leftjust(sand, charval)
      write(100,35) charval, "'SAND'"
      call leftjust(silt, charval)
      write(100,35) charval, "'SILT'"
      call leftjust(clay, charval)
      write(100,35) charval, "'CLAY'"
      call leftjust(rock, charval)
      write(100,35) charval, "'ROCK'"
      call leftjust(bulkd, charval)
      write(100,35) charval, "'BULKD'"
      call leftjust(float(nlayer), charval)
      write(100,35) charval, "'NLAYER'"
      call leftjust(float(nlaypg), charval)
      write(100,35) charval, "'NLAYPG'"
      call leftjust(drain, charval)
      write(100,35) charval, "'DRAIN'"
      call leftjust(basef, charval)
      write(100,35) charval, "'BASEF'"
      call leftjust(stormf, charval)
      write(100,35) charval, "'STORMF'"
      call leftjust(precro, charval)
      write(100,35) charval, "'PRECRO'"
      call leftjust(fracro, charval)
      write(100,35) charval, "'FRACRO'"
      call leftjust(float(swflag), charval)
      write(100,35) charval, "'SWFLAG'"
      call leftjust(awilt(1), charval)
      write(100,35) charval, "'AWILT(1)'"
      call leftjust(awilt(2), charval)
      write(100,35) charval, "'AWILT(2)'"
      call leftjust(awilt(3), charval)
      write(100,35) charval, "'AWILT(3)'"
      call leftjust(awilt(4), charval)
      write(100,35) charval, "'AWILT(4)'"
      call leftjust(awilt(5), charval)
      write(100,35) charval, "'AWILT(5)'"
      call leftjust(awilt(6), charval)
      write(100,35) charval, "'AWILT(6)'"
      call leftjust(awilt(7), charval)
      write(100,35) charval, "'AWILT(7)'"
      call leftjust(awilt(9), charval)
      write(100,35) charval, "'AWILT(8)'"
      call leftjust(awilt(9), charval)
      write(100,35) charval, "'AWILT(9)'"
      call leftjust(awilt(10), charval)
      write(100,35) charval, "'AWILT(10)'"
      call leftjust(afiel(1), charval)
      write(100,35) charval, "'AFIEL(1)'"
      call leftjust(afiel(2), charval)
      write(100,35) charval, "'AFIEL(2)'"
      call leftjust(afiel(3), charval)
      write(100,35) charval, "'AFIEL(3)'"
      call leftjust(afiel(4), charval)
      write(100,35) charval, "'AFIEL(4)'"
      call leftjust(afiel(5), charval)
      write(100,35) charval, "'AFIEL(5)'"
      call leftjust(afiel(6), charval)
      write(100,35) charval, "'AFIEL(6)'"
      call leftjust(afiel(7), charval)
      write(100,35) charval, "'AFIEL(7)'"
      call leftjust(afiel(8), charval)
      write(100,35) charval, "'AFIEL(8)'"
      call leftjust(afiel(9), charval)
      write(100,35) charval, "'AFIEL(9)'"
      call leftjust(afiel(10), charval)
      write(100,35) charval, "'AFIEL(10)'"
      call leftjust(ph, charval)
      write(100,35) charval, "'PH'"
      call leftjust(pslsrb, charval)
      write(100,35) charval, "'PSLSRB'"
      call leftjust(sorpmx, charval)
      write(100,35) charval, "'SORPMX'"
 
      write(100,25) '*** External nutrient input parameters'
      call leftjust(epnfa(1), charval)
      write(100,35) charval, "'EPNFA(1)'"
      call leftjust(epnfa(2), charval)
      write(100,35) charval, "'EPNFA(2)'"
      call leftjust(epnfs(1), charval)
      write(100,35) charval, "'EPNFS(1)'"
      call leftjust(epnfs(2), charval)
      write(100,35) charval, "'EPNFS(2)'"
      call leftjust(satmos(1), charval)
      write(100,35) charval, "'SATMOS(1)'"
      call leftjust(satmos(2), charval)
      write(100,35) charval, "'SATMOS(2)'"
      call leftjust(sirri, charval)
      write(100,35) charval, "'SIRRI'"

      write(100,25) '*** Organic matter initial values'
      call leftjust(som1ci(1,1), charval)
      write(100,35) charval, "'SOM1CI(1,1)'"
      call leftjust(som1ci(1,2), charval)
      write(100,35) charval, "'SOM1CI(1,2)'"
      call leftjust(som1ci(2,1), charval)
      write(100,35) charval, "'SOM1CI(2,1)'"
      call leftjust(som1ci(2,2), charval)
      write(100,35) charval, "'SOM1CI(2,2)'"
      call leftjust(som2ci(1,1), charval)
      write(100,35) charval, "'SOM2CI(1,1)'"
      call leftjust(som2ci(1,2), charval)
      write(100,35) charval, "'SOM2CI(1,2)'"
      call leftjust(som2ci(2,1), charval)
      write(100,35) charval, "'SOM2CI(2,1)'"
      call leftjust(som2ci(2,2), charval)
      write(100,35) charval, "'SOM2CI(2,2)'"
      call leftjust(som3ci(1), charval)
      write(100,35) charval, "'SOM3CI(1)'"
      call leftjust(som3ci(2), charval)
      write(100,35) charval, "'SOM3CI(2)'"
      call leftjust(rces1(1,1), charval)
      write(100,35) charval, "'RCES1(1,1)'"
      call leftjust(rces1(1,2), charval)
      write(100,35) charval, "'RCES1(1,2)'"
      call leftjust(rces1(1,3), charval)
      write(100,35) charval, "'RCES1(1,3)'"
      call leftjust(rces1(2,1), charval)
      write(100,35) charval, "'RCES1(2,1)'"
      call leftjust(rces1(2,2), charval)
      write(100,35) charval, "'RCES1(2,2)'"
      call leftjust(rces1(2,3), charval)
      write(100,35) charval, "'RCES1(2,3)'"
      call leftjust(rces2(1,1), charval)
      write(100,35) charval, "'RCES2(1,1)'"
      call leftjust(rces2(1,2), charval)
      write(100,35) charval, "'RCES2(1,2)'"
      call leftjust(rces2(1,3), charval)
      write(100,35) charval, "'RCES2(1,3)'"
      call leftjust(rces2(2,1), charval)
      write(100,35) charval, "'RCES2(2,1)'"
      call leftjust(rces2(2,2), charval)
      write(100,35) charval, "'RCES2(2,2)'"
      call leftjust(rces2(2,3), charval)
      write(100,35) charval, "'RCES2(2,3)'"
      call leftjust(rces3(1), charval)
      write(100,35) charval, "'RCES3(1)'"
      call leftjust(rces3(2), charval)
      write(100,35) charval, "'RCES3(2)'"
      call leftjust(rces3(3), charval)
      write(100,35) charval, "'RCES3(3)'"
      call leftjust(clittr(1,1), charval)
      write(100,35) charval, "'CLITTR(1,1)'"
      call leftjust(clittr(1,2), charval)
      write(100,35) charval, "'CLITTR(1,2)'"
      call leftjust(clittr(2,1), charval)
      write(100,35) charval, "'CLITTR(2,1)'"
      call leftjust(clittr(2,2), charval)
      write(100,35) charval, "'CLITTR(2,2)'"
      call leftjust(rcelit(1,1), charval)
      write(100,35) charval, "'RCELIT(1,1)'"
      call leftjust(rcelit(1,2), charval)
      write(100,35) charval, "'RCELIT(1,2)'"
      call leftjust(rcelit(1,3), charval)
      write(100,35) charval, "'RCELIT(1,3)'"
      call leftjust(rcelit(2,1), charval)
      write(100,35) charval, "'RCELIT(2,1)'"
      call leftjust(rcelit(2,2), charval)
      write(100,35) charval, "'RCELIT(2,2)'"
      call leftjust(rcelit(2,3), charval)
      write(100,35) charval, "'RCELIT(2,3)'"
      call leftjust(aglcis(1), charval)
      write(100,35) charval, "'AGLCIS(1)'"
      call leftjust(aglcis(2), charval)
      write(100,35) charval, "'AGLCIS(2)'"
      call leftjust(aglive(1), charval)
      write(100,35) charval, "'AGLIVE(1)'"
      call leftjust(aglive(2), charval)
      write(100,35) charval, "'AGLIVE(2)'"
      call leftjust(aglive(3), charval)
      write(100,35) charval, "'AGLIVE(3)'"
      call leftjust(bglcisj(1), charval)
      write(100,35) charval, "'BGLCISJ(1)'"
      call leftjust(bglcisj(2), charval)
      write(100,35) charval, "'BGLCISJ(2)'"
      call leftjust(bglivej(1), charval)
      write(100,35) charval, "'BGLIVEJ(1)'"
      call leftjust(bglivej(2), charval)
      write(100,35) charval, "'BGLIVEJ(2)'"
      call leftjust(bglivej(3), charval)
      write(100,35) charval, "'BGLIVEJ(3)'"
      call leftjust(stdcis(1), charval)
      write(100,35) charval, "'STDCIS(1)'"
      call leftjust(stdcis(2), charval)
      write(100,35) charval, "'STDCIS(2)'"
      call leftjust(stdede(1), charval)
      write(100,35) charval, "'STDEDE(1)'"
      call leftjust(stdede(2), charval)
      write(100,35) charval, "'STDEDE(2)'"
      call leftjust(stdede(3), charval)
      write(100,35) charval, "'STDEDE(3)'"

      write(100,25) '*** Forest organic matter initial parameters'
      call leftjust(rlvcis(1), charval)
      write(100,35) charval, "'RLVCIS(1)'"
      call leftjust(rlvcis(2), charval)
      write(100,35) charval, "'RLVCIS(2)'"
      call leftjust(rleave(1), charval)
      write(100,35) charval, "'RLEAVE(1)'"
      call leftjust(rleave(2), charval)
      write(100,35) charval, "'RLEAVE(2)'"
      call leftjust(rleave(3), charval)
      write(100,35) charval, "'RLEAVE(3)'"
      call leftjust(fbrcis(1), charval)
      write(100,35) charval, "'FBRCIS(1)'"
      call leftjust(fbrcis(2), charval)
      write(100,35) charval, "'FBRCIS(2)'"
      call leftjust(fbrche(1), charval)
      write(100,35) charval, "'FBRCHE(1)'"
      call leftjust(fbrche(2), charval)
      write(100,35) charval, "'FBRCHE(2)'"
      call leftjust(fbrche(3), charval)
      write(100,35) charval, "'FBRCHE(3)'"
      call leftjust(rlwcis(1), charval)
      write(100,35) charval, "'RLWCIS(1)'"
      call leftjust(rlwcis(2), charval)
      write(100,35) charval, "'RLWCIS(2)'"
      call leftjust(rlwode(1), charval)
      write(100,35) charval, "'RLWODE(1)'"
      call leftjust(rlwode(2), charval)
      write(100,35) charval, "'RLWODE(2)'"
      call leftjust(rlwode(3), charval)
      write(100,35) charval, "'RLWODE(3)'"
      call leftjust(frtcisj(1), charval)
      write(100,35) charval, "'FRTCISJ(1)'"
      call leftjust(frtcisj(2), charval)
      write(100,35) charval, "'FRTCISJ(2)'"
      call leftjust(frootej(1), charval)
      write(100,35) charval, "'FROOTEJ(1)'"
      call leftjust(frootej(2), charval)
      write(100,35) charval, "'FROOTEJ(2)'"
      call leftjust(frootej(3), charval)
      write(100,35) charval, "'FROOTEJ(3)'"
      call leftjust(crtcis(1), charval)
      write(100,35) charval, "'CRTCIS(1)'"
      call leftjust(crtcis(2), charval)
      write(100,35) charval, "'CRTCIS(2)'"
      call leftjust(croote(1), charval)
      write(100,35) charval, "'CROOTE(1)'"
      call leftjust(croote(2), charval)
      write(100,35) charval, "'CROOTE(2)'"
      call leftjust(croote(3), charval)
      write(100,35) charval, "'CROOTE(3)'"
      call leftjust(wd1cis(1), charval)
      write(100,35) charval, "'WD1CIS(1)'"
      call leftjust(wd1cis(2), charval)
      write(100,35) charval, "'WD1CIS(2)'"
      call leftjust(wd2cis(1), charval)
      write(100,35) charval, "'WD2CIS(1)'"
      call leftjust(wd2cis(2), charval)
      write(100,35) charval, "'WD2CIS(2)'"
      call leftjust(wd3cis(1), charval)
      write(100,35) charval, "'WD3CIS(1)'"
      call leftjust(wd3cis(2), charval)
      write(100,35) charval, "'WD3CIS(2)'"

      write(100,25) '*** Mineral initial parameters'
      call leftjust(minerl(1,1), charval)
      write(100,35) charval, "'MINERL(1,1)'"
      call leftjust(minerl(2,1), charval)
      write(100,35) charval, "'MINERL(2,1)'"
      call leftjust(minerl(3,1), charval)
      write(100,35) charval, "'MINERL(3,1)'"
      call leftjust(minerl(4,1), charval)
      write(100,35) charval, "'MINERL(4,1)'"
      call leftjust(minerl(5,1), charval)
      write(100,35) charval, "'MINERL(5,1)'"
      call leftjust(minerl(6,1), charval)
      write(100,35) charval, "'MINERL(6,1)'"
      call leftjust(minerl(7,1), charval)
      write(100,35) charval, "'MINERL(7,1)'"
      call leftjust(minerl(8,1), charval)
      write(100,35) charval, "'MINERL(8,1)'"
      call leftjust(minerl(9,1), charval)
      write(100,35) charval, "'MINERL(9,1)'"
      call leftjust(minerl(10,1), charval)
      write(100,35) charval, "'MINERL(10,1)'"
      call leftjust(minerl(1,2), charval)
      write(100,35) charval, "'MINERL(1,2)'"
      call leftjust(minerl(2,2), charval)
      write(100,35) charval, "'MINERL(2,2)'"
      call leftjust(minerl(3,2), charval)
      write(100,35) charval, "'MINERL(3,2)'"
      call leftjust(minerl(4,2), charval)
      write(100,35) charval, "'MINERL(4,2)'"
      call leftjust(minerl(5,2), charval)
      write(100,35) charval, "'MINERL(5,2)'"
      call leftjust(minerl(6,2), charval)
      write(100,35) charval, "'MINERL(6,2)'"
      call leftjust(minerl(7,2), charval)
      write(100,35) charval, "'MINERL(7,2)'"
      call leftjust(minerl(8,2), charval)
      write(100,35) charval, "'MINERL(8,2)'"
      call leftjust(minerl(9,2), charval)
      write(100,35) charval, "'MINERL(9,2)'"
      call leftjust(minerl(10,2), charval)
      write(100,35) charval, "'MINERL(10,2)'"
      call leftjust(minerl(1,3), charval)
      write(100,35) charval, "'MINERL(1,3)'"
      call leftjust(minerl(2,3), charval)
      write(100,35) charval, "'MINERL(2,3)'"
      call leftjust(minerl(3,3), charval)
      write(100,35) charval, "'MINERL(3,3)'"
      call leftjust(minerl(4,3), charval)
      write(100,35) charval, "'MINERL(4,3)'"
      call leftjust(minerl(5,3), charval)
      write(100,35) charval, "'MINERL(5,3)'"
      call leftjust(minerl(6,3), charval)
      write(100,35) charval, "'MINERL(6,3)'"
      call leftjust(minerl(7,3), charval)
      write(100,35) charval, "'MINERL(7,3)'"
      call leftjust(minerl(8,3), charval)
      write(100,35) charval, "'MINERL(8,3)'"
      call leftjust(minerl(9,3), charval)
      write(100,35) charval, "'MINERL(9,3)'"
      call leftjust(minerl(10,3), charval)
      write(100,35) charval, "'MINERL(10,3)'"
      call leftjust(parent(1), charval)
      write(100,35) charval, "'PARENT(1)'"
      call leftjust(parent(2), charval)
      write(100,35) charval, "'PARENT(2)'"
      call leftjust(parent(3), charval)
      write(100,35) charval, "'PARENT(3)'"
      call leftjust(secndy(1), charval)
      write(100,35) charval, "'SECNDY(1)'"
      call leftjust(secndy(2), charval)
      write(100,35) charval, "'SECNDY(2)'"
      call leftjust(secndy(3), charval)
      write(100,35) charval, "'SECNDY(3)'"
      call leftjust(occlud, charval)
      write(100,35) charval, "'OCCLUD'"

      write(100,25) '*** Water initial parameters'
      call leftjust(rwcf(1), charval)
      write(100,35) charval, "'RWCF(1)'"
      call leftjust(rwcf(2), charval)
      write(100,35) charval, "'RWCF(2)'"
      call leftjust(rwcf(3), charval)
      write(100,35) charval, "'RWCF(3)'"
      call leftjust(rwcf(4), charval)
      write(100,35) charval, "'RWCF(4)'"
      call leftjust(rwcf(5), charval)
      write(100,35) charval, "'RWCF(5)'"
      call leftjust(rwcf(6), charval)
      write(100,35) charval, "'RWCF(6)'"
      call leftjust(rwcf(7), charval)
      write(100,35) charval, "'RWCF(7)'"
      call leftjust(rwcf(8), charval)
      write(100,35) charval, "'RWCF(8)'"
      call leftjust(rwcf(9), charval)
      write(100,35) charval, "'RWCF(9)'"
      call leftjust(rwcf(10), charval)
      write(100,35) charval, "'RWCF(10)'"
      call leftjust(snlq, charval)
      write(100,35) charval, "'SNLQ'"
      call leftjust(snow, charval)
      write(100,35) charval, "'SNOW'"

      write(100,25) '*** Enhanced extend parameters'
      call leftjust(metcis(1,1), charval)
      write(100,35) charval, "'METCIS(1,1)'"
      call leftjust(metcis(1,2), charval)
      write(100,35) charval, "'METCIS(1,2)'"
      call leftjust(metcis(2,1), charval)
      write(100,35) charval, "'METCIS(2,1)'"
      call leftjust(metcis(2,2), charval)
      write(100,35) charval, "'METCIS(2,2)'"
      call leftjust(metabe(1,1), charval)
      write(100,35) charval, "'METABE(1,1)'"
      call leftjust(metabe(1,2), charval)
      write(100,35) charval, "'METABE(1,2)'"
      call leftjust(metabe(1,3), charval)
      write(100,35) charval, "'METABE(1,3)'"
      call leftjust(metabe(2,1), charval)
      write(100,35) charval, "'METABE(2,1)'"
      call leftjust(metabe(2,2), charval)
      write(100,35) charval, "'METABE(2,2)'"
      call leftjust(metabe(2,3), charval)
      write(100,35) charval, "'METABE(2,3)'"
      call leftjust(wood1e(1), charval)
      write(100,35) charval, "'WOOD1E(1)'"
      call leftjust(wood1e(2), charval)
      write(100,35) charval, "'WOOD1E(2)'"
      call leftjust(wood1e(3), charval)
      write(100,35) charval, "'WOOD1E(3)'"
      call leftjust(wood2e(1), charval)
      write(100,35) charval, "'WOOD2E(1)'"
      call leftjust(wood2e(2), charval)
      write(100,35) charval, "'WOOD2E(2)'"
      call leftjust(wood2e(3), charval)
      write(100,35) charval, "'WOOD2E(3)'"
      call leftjust(wood3e(1), charval)
      write(100,35) charval, "'WOOD3E(1)'"
      call leftjust(wood3e(2), charval)
      write(100,35) charval, "'WOOD3E(2)'"
      call leftjust(wood3e(3), charval)
      write(100,35) charval, "'WOOD3E(3)'"
      call leftjust(crpstg(1), charval)
      write(100,35) charval, "'CRPSTG(1)'"
      call leftjust(crpstg(2), charval)
      write(100,35) charval, "'CRPSTG(2)'"
      call leftjust(crpstg(3), charval)
      write(100,35) charval, "'CRPSTG(3)'"
      call leftjust(strlig(1), charval)
      write(100,35) charval, "'STRLIG(1)'"
      call leftjust(strlig(2), charval)
      write(100,35) charval, "'STRLIG(2)'"

      call leftjust(som1e(1,1), charval)
      write(100,35) charval, "'SOM1E(1,1)'"
      call leftjust(som1e(1,2), charval)
      write(100,35) charval, "'SOM1E(1,2)'"
      call leftjust(som1e(1,3), charval)
      write(100,35) charval, "'SOM1E(1,3)'"
      call leftjust(som1e(2,1), charval)
      write(100,35) charval, "'SOM1E(2,1)'"
      call leftjust(som1e(2,2), charval)
      write(100,35) charval, "'SOM1E(2,2)'"
      call leftjust(som1e(2,3), charval)
      write(100,35) charval, "'SOM1E(2,3)'"
      call leftjust(som2e(1,1), charval)
      write(100,35) charval, "'SOM2E(1,1)'"
      call leftjust(som2e(1,2), charval)
      write(100,35) charval, "'SOM2E(1,2)'"
      call leftjust(som2e(1,3), charval)
      write(100,35) charval, "'SOM2E(1,3)'"
      call leftjust(som2e(2,1), charval)
      write(100,35) charval, "'SOM2E(2,1)'"
      call leftjust(som2e(2,2), charval)
      write(100,35) charval, "'SOM2E(2,2)'"
      call leftjust(som2e(2,3), charval)
      write(100,35) charval, "'SOM2E(2,3)'"
      call leftjust(som3e(1), charval)
      write(100,35) charval, "'SOM3E(1)'"
      call leftjust(som3e(2), charval)
      write(100,35) charval, "'SOM3E(2)'"
      call leftjust(som3e(3), charval)
      write(100,35) charval, "'SOM3E(3)'"
      call leftjust(bglcism(1), charval)
      write(100,35) charval, "'BGLCISM(1)'"
      call leftjust(bglcism(2), charval)
      write(100,35) charval, "'BGLCISM(2)'"
      call leftjust(bglivem(1), charval)
      write(100,35) charval, "'BGLIVEM(1)'"
      call leftjust(bglivem(2), charval)
      write(100,35) charval, "'BGLIVEM(2)'"
      call leftjust(bglivem(3), charval)
      write(100,35) charval, "'BGLIVEM(3)'"
      call leftjust(frtcism(1), charval)
      write(100,35) charval, "'FRTCISM(1)'"
      call leftjust(frtcism(2), charval)
      write(100,35) charval, "'FRTCISM(2)'"
      call leftjust(frootem(1), charval)
      write(100,35) charval, "'FROOTEM(1)'"
      call leftjust(frootem(2), charval)
      write(100,35) charval, "'FROOTEM(2)'"
      call leftjust(frootem(3), charval)
      write(100,35) charval, "'FROOTEM(3)'"
      call leftjust(strcis(1,1), charval)
      write(100,35) charval, "'STRCIS(1,1)'"
      call leftjust(strcis(1,2), charval)
      write(100,35) charval, "'STRCIS(1,2)'"
      call leftjust(strcis(2,1), charval)
      write(100,35) charval, "'STRCIS(2,1)'"
      call leftjust(strcis(2,2), charval)
      write(100,35) charval, "'STRCIS(2,2)'"
      call leftjust(struce(1,1), charval)
      write(100,35) charval, "'STRUCE(1,1)'"
      call leftjust(struce(1,2), charval)
      write(100,35) charval, "'STRUCE(1,2)'"
      call leftjust(struce(1,3), charval)
      write(100,35) charval, "'STRUCE(1,3)'"
      call leftjust(struce(2,1), charval)
      write(100,35) charval, "'STRUCE(2,1)'"
      call leftjust(struce(2,2), charval)
      write(100,35) charval, "'STRUCE(2,2)'"
      call leftjust(struce(2,3), charval)
      write(100,35) charval, "'STRUCE(2,3)'"
      call leftjust(carbostg(1,1), charval)
      write(100,35) charval, "'CARBOSTG(1,1)'"
      call leftjust(carbostg(1,2), charval)
      write(100,35) charval, "'CARBOSTG(1,2)'"
      call leftjust(carbostg(2,1), charval)
      write(100,35) charval, "'CARBOSTG(2,1)'"
      call leftjust(carbostg(2,2), charval)
      write(100,35) charval, "'CARBOSTG(2,2)'"
      call leftjust(csrsnk(1), charval)
      write(100,35) charval, "'CSRSNK(1)'"
      call leftjust(csrsnk(2), charval)
      write(100,35) charval, "'CSRSNK(2)'"
      call leftjust(esrsnk(1), charval)
      write(100,35) charval, "'ESRSNK(1)'"
      call leftjust(esrsnk(2), charval)
      write(100,35) charval, "'ESRSNK(2)'"
      call leftjust(esrsnk(3), charval)
      write(100,35) charval, "'ESRSNK(3)'"
      call leftjust(forstg(1), charval)
      write(100,35) charval, "'FORSTG(1)'"
      call leftjust(forstg(2), charval)
      write(100,35) charval, "'FORSTG(2)'"
      call leftjust(forstg(3), charval)
      write(100,35) charval, "'FORSTG(3)'"
      call leftjust(real(ammonium), charval)
      write(100,35) charval, "'AMMONIUM'"
      call leftjust(real(nitrate(1)), charval)
      write(100,35) charval, "'NITRATE(1)'"
      call leftjust(real(nitrate(2)), charval)
      write(100,35) charval, "'NITRATE(2)'"
      call leftjust(real(nitrate(3)), charval)
      write(100,35) charval, "'NITRATE(3)'"
      call leftjust(real(nitrate(4)), charval)
      write(100,35) charval, "'NITRATE(4)'"
      call leftjust(real(nitrate(5)), charval)
      write(100,35) charval, "'NITRATE(5)'"
      call leftjust(real(nitrate(6)), charval)
      write(100,35) charval, "'NITRATE(6)'"
      call leftjust(real(nitrate(7)), charval)
      write(100,35) charval, "'NITRATE(7)'"
      call leftjust(real(nitrate(8)), charval)
      write(100,35) charval, "'NITRATE(8)'"
      call leftjust(real(nitrate(9)), charval)
      write(100,35) charval, "'NITRATE(9)'"
      call leftjust(real(nitrate(10)), charval)
      write(100,35) charval, "'NITRATE(10)'"
      call leftjust(real(nitrate(11)), charval)
      write(100,35) charval, "'NITRATE(11)'"
      call leftjust(real(nitrate(12)), charval)
      write(100,35) charval, "'NITRATE(12)'"
      call leftjust(real(nitrate(13)), charval)
      write(100,35) charval, "'NITRATE(13)'"
      call leftjust(real(nitrate(14)), charval)
      write(100,35) charval, "'NITRATE(14)'"
      call leftjust(real(nitrate(15)), charval)
      write(100,35) charval, "'NITRATE(15)'"
      call leftjust(real(nitrate(16)), charval)
      write(100,35) charval, "'NITRATE(16)'"
      call leftjust(real(nitrate(17)), charval)
      write(100,35) charval, "'NITRATE(17)'"
      call leftjust(real(nitrate(18)), charval)
      write(100,35) charval, "'NITRATE(18)'"
      call leftjust(real(nitrate(19)), charval)
      write(100,35) charval, "'NITRATE(19)'"
      call leftjust(real(nitrate(20)), charval)
      write(100,35) charval, "'NITRATE(20)'"
      call leftjust(real(nitrate(21)), charval)
      write(100,35) charval, "'NITRATE(21)'"

      call getswc(swcextend, SWMAXLYR)
      call leftjust(swcextend(1), charval)
      write(100,35) charval, "'SWCEXTEND(1)'"
      call leftjust(swcextend(2), charval)
      write(100,35) charval, "'SWCEXTEND(2)'"
      call leftjust(swcextend(3), charval)
      write(100,35) charval, "'SWCEXTEND(3)'"
      call leftjust(swcextend(4), charval)
      write(100,35) charval, "'SWCEXTEND(4)'"
      call leftjust(swcextend(5), charval)
      write(100,35) charval, "'SWCEXTEND(5)'"
      call leftjust(swcextend(6), charval)
      write(100,35) charval, "'SWCEXTEND(6)'"
      call leftjust(swcextend(7), charval)
      write(100,35) charval, "'SWCEXTEND(7)'"
      call leftjust(swcextend(8), charval)
      write(100,35) charval, "'SWCEXTEND(8)'"
      call leftjust(swcextend(9), charval)
      write(100,35) charval, "'SWCEXTEND(9)'"
      call leftjust(swcextend(10), charval)
      write(100,35) charval, "'SWCEXTEND(10)'"
      call leftjust(swcextend(11), charval)
      write(100,35) charval, "'SWCEXTEND(11)'"
      call leftjust(swcextend(12), charval)
      write(100,35) charval, "'SWCEXTEND(12)'"
      call leftjust(swcextend(13), charval)
      write(100,35) charval, "'SWCEXTEND(13)'"
      call leftjust(swcextend(14), charval)
      write(100,35) charval, "'SWCEXTEND(14)'"
      call leftjust(swcextend(15), charval)
      write(100,35) charval, "'SWCEXTEND(15)'"
      call leftjust(swcextend(16), charval)
      write(100,35) charval, "'SWCEXTEND(16)'"
      call leftjust(swcextend(17), charval)
      write(100,35) charval, "'SWCEXTEND(17)'"
      call leftjust(swcextend(18), charval)
      write(100,35) charval, "'SWCEXTEND(18)'"
      call leftjust(swcextend(19), charval)
      write(100,35) charval, "'SWCEXTEND(19)'"
      call leftjust(swcextend(20), charval)
      write(100,35) charval, "'SWCEXTEND(20)'"
      call leftjust(swcextend(21), charval)
      write(100,35) charval, "'SWCEXTEND(21)'"

15    format(a,t7,a,a)
25    format(a)
35    format(a,t19,a) 

      close(unit=100)

      return
      end


      subroutine leftjust(value, charval)

      implicit none

c ... Reformat numeric values as character strings so that they can be
c ... printed left justified in the output file

c ... Argument declarations
      real value
      character*20 charval

c ... Internal write
      write(charval, '(F12.4)') value

c ... Strip off leading blanks as necessary
      do while (charval(1:1) .eq. ' ')
        charval = charval(2:20)
      end do

c ... Add leading zero if necessary
      if (charval(1:1) .eq. '.') then
        charval = '0' // charval(1:19)
      endif

      return
      end
