
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... SITEIN_EXT.F

      subroutine sitein_ext(swcextend)

      implicit none
      include 'const.inc'
      include 'doubles.inc'
      include 'ligvar.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'wth.inc'

c ... Argument declarations
      real swcextend(SWMAXLYR)

c ... Read the parameter file in the form needed by time0.

c ... Local variables
      integer     ii, jj
      character*6 name
      real        temp

      read(7,*)
      read(7,*)

      do 10 ii = 1, MONTHS
        read(7,*) precip(ii), name
        call ckdata('sitein_ext','precip',name)
10    continue

      do 20 ii = 1, MONTHS
        read(7,*) prcstd(ii), name
        call ckdata('sitein_ext','prcstd',name)
20    continue

      do 30 ii = 1, MONTHS
        read(7,*) prcskw(ii), name
        call ckdata('sitein_ext','prcskw',name)
30    continue

      do 40 ii = 1, MONTHS
        read(7,*) tmn2m(ii), name
        call ckdata('sitein_ext','tmn2m',name)
40    continue

c ... maxt calculation for nitrify added 4/03/98 - mdh
      maxt = -99.0
      do 50 ii = 1, MONTHS
        read(7,*) tmx2m(ii), name
        if (tmx2m(ii) .gt. maxt) maxt = tmx2m(ii)
        call ckdata('sitein_ext','tmx2m',name)
50    continue

c ... Check to be sure that the tmax value read from the <site>.100
c ... file is greater than the tmin value read from the <site>.100
c ... file, cak - 09/17/03
      do 55 ii = 1, MONTHS
        if (tmx2m(ii) .lt. tmn2m(ii)) then
          write(*,*) 'ERROR: Invalid weather data in site file, ',
     &               'tmx2m < tmn2m for month ', ii
          write(*,*) 'tmx2m(', ii, ') = ', tmx2m(ii)
          write(*,*) 'tmn2m(', ii, ') = ', tmn2m(ii)
          STOP
        endif
55    continue

      read(7,*)
      read(7,*) temp, name
      ivauto = int(temp)
      call ckdata('sitein_ext','ivauto',name)
      read(7,*) temp, name
      nelem = int(temp)
      call ckdata('sitein_ext','nelem',name)
 
      read(7,*) sitlat, name
      call ckdata('sitein_ext','sitlat',name)
      read(7,*) sitlng, name
      call ckdata('sitein_ext','sitlng',name)

      read(7,*) sand, name
      call ckdata('sitein_ext','sand',name)
      read(7,*) silt, name
      call ckdata('sitein_ext','silt',name)
      read(7,*) clay, name
      call ckdata('sitein_ext','clay',name)
c ... Add rock fraction, cak - 05/27/03
      read(7,*) rock, name
      call ckdata('sitein_ext','rock',name)
      if (rock .gt. 0.90) then 
        write(*,*) 'Rock fraction too large, rock = ', rock
        STOP
      endif
      read(7,*) bulkd, name
      call ckdata('sitein_ext','bulkd',name)

      read(7,*) temp, name
      nlayer = int(temp)
      if (nlayer .gt. 9) then
        nlayer = 9
        call message('   Warning: nlayer value too large, reset to 9')
      endif
      call ckdata('sitein_ext','nlayer',name)
      read(7,*) temp, name
      nlaypg = int(temp)
      call ckdata('sitein_ext','nlaypg',name)
      read(7,*) drain, name
      call ckdata('sitein_ext','drain',name)
      read(7,*) basef, name
      call ckdata('sitein_ext','basef', name)
      read(7,*) stormf, name
      call ckdata('sitein_ext','stormf', name)
c ... Add precipitation amount required for runoff and fraction of
c ... precipitation above amount required for runoff which is lost
c ... via runoff, cak - 05/27/03
      read(7,*) precro, name
      call ckdata('sitein_ext','precro', name)
      read(7,*) fracro, name
      call ckdata('sitein_ext','fracro', name)
      read(7,*) temp, name
      swflag = int(temp)
      call ckdata('sitein_ext','swflag', name)

      do 60 ii = 1, MAXLYR
        read(7,*) awilt(ii), name
        call ckdata('sitein_ext','awilt', name)
60    continue

      do 70 ii = 1, MAXLYR
        read(7,*) afiel(ii), name
        call ckdata('sitein_ext','afiel', name)
70    continue

      read(7,*) ph, name
      call ckdata('sitein_ext','ph',name)
c ... New phstart variable added for pH shift, cak - 08/02/02
      phstart = ph
      read(7,*) pslsrb, name
      call ckdata('sitein_ext','pslsrb',name)
      read(7,*) sorpmx, name
      call ckdata('sitein_ext','sorpmx',name)
      read(7,*)
      read(7,*) epnfa(INTCPT), name
      call ckdata('sitein_ext','epnfa',name)
      read(7,*) epnfa(SLOPE), name
      call ckdata('sitein_ext','epnfa',name)
      read(7,*) epnfs(INTCPT), name
      call ckdata('sitein_ext','epnfs',name)
      read(7,*) epnfs(SLOPE), name
      call ckdata('sitein_ext','epnfs',name)
      read(7,*) satmos(INTCPT), name
      call ckdata('sitein_ext','satmos',name)
      read(7,*) satmos(SLOPE), name
      call ckdata('sitein_ext','satmos',name)
      read(7,*) sirri, name
      call ckdata('sitein_ext','sirri',name)

      read(7,*)
      read(7,*) som1ci(SRFC,UNLABL), name
      call ckdata('sitein_ext','som1ci',name)
      read(7,*) som1ci(SRFC,LABELD), name
      call ckdata('sitein_ext','som1ci',name)
      read(7,*) som1ci(SOIL,UNLABL), name
      call ckdata('sitein_ext','som1ci',name)
      read(7,*) som1ci(SOIL,LABELD), name
      call ckdata('sitein_ext','som1ci',name)

      read(7,*) som2ci(SRFC,UNLABL), name
      call ckdata('sitein_ext','som2ci',name)
      read(7,*) som2ci(SRFC,LABELD), name
      call ckdata('sitein_ext','som2ci',name)
      read(7,*) som2ci(SOIL,UNLABL), name
      call ckdata('sitein_ext','som2ci',name)
      read(7,*) som2ci(SOIL,LABELD), name
      call ckdata('sitein_ext','som2ci',name)

      read(7,*) som3ci(UNLABL), name
      call ckdata('sitein_ext','som3ci',name)
      read(7,*) som3ci(LABELD), name
      call ckdata('sitein_ext','som3ci',name)

      do 90 ii = SRFC, SOIL
        do 80 jj = 1, MAXIEL
          read(7,*) rces1(ii,jj), name
          call ckdata('sitein_ext','rces1',name)
80      continue
90    continue

      do 105 ii = SRFC, SOIL
        do 100 jj = 1, MAXIEL
          read(7,*) rces2(ii,jj), name
          call ckdata('sitein_ext','rces2',name)
100     continue
105   continue

      do 110 ii = 1, MAXIEL
        read(7,*) rces3(ii), name
        call ckdata('sitein_ext','rces3',name)
110   continue

      read(7,*) clittr(SRFC,UNLABL), name
      call ckdata('sitein_ext','clittr',name)
      read(7,*) clittr(SRFC,LABELD), name
      call ckdata('sitein_ext','clittr',name)
      read(7,*) clittr(SOIL,UNLABL), name
      call ckdata('sitein_ext','clittr',name)
      read(7,*) clittr(SOIL,LABELD), name
      call ckdata('sitein_ext','clittr',name)
c ... Add check for initial litter values to prevent divide by zero error
c ... in calciv subroutine
      if ((ivauto .eq. 0) .and. 
     &    (clittr(SRFC,UNLABL) + clittr(SRFC,LABELD) .le. 0.0)) then
C    &    (clittr(SRFC,UNLABL) + clittr(SRFC,LABELD) .eq. 0.0)) then
        call message('   Initial surface litter values, clittr(1,*),')
        call message('   must be greater than zero.')
        STOP
      endif
      if ((ivauto .eq. 0) .and. 
     &    (clittr(SOIL,UNLABL) + clittr(SOIL,LABELD) .le. 0.0)) then
C    &    (clittr(SOIL,UNLABL) + clittr(SOIL,LABELD) .eq. 0.0)) then
        call message('   Initial soil litter values, clittr(2,*),')
        call message('   must be greater than zero.')
        STOP
      endif

      do 130 ii = SRFC, SOIL
        do 120 jj = 1, MAXIEL
          read(7,*) rcelit(ii,jj), name
          call ckdata('sitein_ext','rcelit',name)
120     continue
130   continue

      read(7,*) aglcis(UNLABL), name
      call ckdata('sitein_ext','aglcis',name)
      read(7,*) aglcis(LABELD), name
      call ckdata('sitein_ext','aglcis',name)

      do 140 ii = 1, MAXIEL
        read(7,*) aglive(ii), name
        call ckdata('sitein_ext','aglive',name)
140   continue

c ... Read the fine root values into the juvenile fine root pool
      read(7,*) bglcisj(UNLABL), name
      call ckdata('sitein_ext','bglcis',name)
      read(7,*) bglcisj(LABELD), name
      call ckdata('sitein_ext','bglcis',name)

      do 150 ii = 1, MAXIEL
        read(7,*) bglivej(ii), name
        call ckdata('sitein_ext','bglive',name)
150   continue

      read(7,*) stdcis(UNLABL), name
      call ckdata('sitein_ext','stdcis',name)
      read(7,*) stdcis(LABELD), name
      call ckdata('sitein_ext','stdcis',name)

      do 160 ii = 1, MAXIEL
        read(7,*) stdede(ii), name
        call ckdata('sitein_ext','stdede',name)
160   continue

      read(7,*)
      read(7,*) rlvcis(UNLABL), name
      call ckdata('sitein_ext','rlvcis',name)
      read(7,*) rlvcis(LABELD), name
      call ckdata('sitein_ext','rlvcis',name)

      do 170 ii = 1, MAXIEL
        read(7,*) rleave(ii), name
        call ckdata('sitein_ext','rleave',name)
170   continue

      read(7,*) fbrcis(UNLABL), name
      call ckdata('sitein_ext','fbrcis',name)
      read(7,*) fbrcis(LABELD), name
      call ckdata('sitein_ext','fbrcis',name)

      do 180 ii = 1, MAXIEL
        read(7,*) fbrche(ii), name
        call ckdata('sitein_ext','fbrche',name)
180   continue

      read(7,*) rlwcis(UNLABL), name
      call ckdata('sitein_ext','rlwcis',name)
      read(7,*) rlwcis(LABELD), name
      call ckdata('sitein_ext','rlwcis',name)

      do 190 ii = 1, MAXIEL
        read(7,*) rlwode(ii), name
        call ckdata('sitein_ext','rlwode',name)
190   continue

c ... Read the fine root values into the juvenile fine root pools
      read(7,*) frtcisj(UNLABL), name
      call ckdata('sitein_ext','frtcis',name)
      read(7,*) frtcisj(LABELD), name
      call ckdata('sitein_ext','frtcis',name)

      do 200 ii = 1, MAXIEL
        read(7,*) frootej(ii), name
        call ckdata('sitein_ext','froote',name)
200   continue

      read(7,*) crtcis(UNLABL), name
      call ckdata('sitein_ext','crtcis',name)
      read(7,*) crtcis(LABELD), name
      call ckdata('sitein_ext','crtcis',name)

      do 210 ii = 1, MAXIEL
        read(7,*) croote(ii), name
        call ckdata('sitein_ext','croote',name)
210   continue

      read(7,*) wd1cis(UNLABL), name
      call ckdata('sitein_ext','wd1cis',name)
      read(7,*) wd1cis(LABELD), name
      call ckdata('sitein_ext','wd1cis',name)
      read(7,*) wd2cis(UNLABL), name
      call ckdata('sitein_ext','wd2cis',name)
      read(7,*) wd2cis(LABELD), name
      call ckdata('sitein_ext','wd2cis',name)
      read(7,*) wd3cis(UNLABL), name
      call ckdata('sitein_ext','wd3cis',name)
      read(7,*) wd3cis(LABELD), name
      call ckdata('sitein_ext','wd3cis',name)

      read(7,*)

      do 230 ii = 1, MAXIEL
        do 220 jj = 1, MAXLYR
          read(7,*) minerl(jj,ii), name
          call ckdata('sitein_ext','minerl',name)
220     continue
230   continue

      do 240 ii = 1, MAXIEL
        read(7,*) parent(ii), name
        call ckdata('sitein_ext','parent',name)
240   continue

c ... The secndy and occlud input values can now be double precision,
c ... cak - 03/20/02
      do 250 ii = 1, MAXIEL
        read(7,*) secndy_double(ii), name
        call ckdata('sitein_ext','secndy',name)
250   continue

      read(7,*) occlud_double, name
      call ckdata('sitein_ext','occlud',name)
      read(7,*)

c ... Save the double precision secndy and occlud variables read into their
c ... single precision counterparts, cak - 03/20/02
      secndy(1) = real(secndy_double(1))
      secndy(2) = real(secndy_double(2))
      secndy(3) = real(secndy_double(3))
      occlud = real(occlud_double)

      do 260 ii = 1, MAXLYR
        read(7,*) rwcf(ii), name
        call ckdata('sitein_ext','rwcf',name)
260   continue

      read(7,*) snlq, name
      call ckdata('sitein_ext','snlq',name)
      read(7,*) snow, name
      call ckdata('sitein_ext','snow', name)

c ... Enhanced site parameters for an extend simulation
      read(7,*)

      read(7,*) metcis(SRFC,UNLABL), name
      call ckdata('sitein_ext','metcis',name)
      read(7,*) metcis(SRFC,LABELD), name
      call ckdata('sitein_ext','metcis',name)
      read(7,*) metcis(SOIL,UNLABL), name
      call ckdata('sitein_ext','metcis',name)
      read(7,*) metcis(SOIL,LABELD), name
      call ckdata('sitein_ext','metcis',name)
      do 265 ii = 1, MAXIEL
        read(7,*) metabe(SRFC,ii), name
        call ckdata('sitein_ext','metabe',name)
265   continue
      do 270 ii = 1, MAXIEL
        read(7,*) metabe(SOIL,ii), name
        call ckdata('sitein_ext','metabe',name)
270   continue

      do 275 ii = 1, MAXIEL
        read(7,*) wood1e(ii), name
        call ckdata('sitein_ext','wood1e',name)
275   continue
      do 280 ii = 1, MAXIEL
        read(7,*) wood2e(ii), name
        call ckdata('sitein_ext','wood2e',name)
280   continue
      do 285 ii = 1, MAXIEL
        read(7,*) wood3e(ii), name
        call ckdata('sitein_ext','wood3e',name)
285   continue

      do 290 ii = 1, MAXIEL
        read(7,*) crpstg(ii), name
        call ckdata('sitein_ext','crpstg',name)
290   continue

      read(7,*) strlig(SRFC), name
      call ckdata('sitein_ext','strlig',name)
      read(7,*) strlig(SOIL), name
      call ckdata('sitein_ext','strlig',name)


      do 295 ii = 1, MAXIEL
        read(7,*) som1e(SRFC,ii), name
        call ckdata('sitein_ext','som1e',name)
295   continue
      do 300 ii = 1, MAXIEL
        read(7,*) som1e(SOIL,ii), name
        call ckdata('sitein_ext','som1e',name)
300   continue
      do 305 ii = 1, MAXIEL
        read(7,*) som2e(SRFC,ii), name
        call ckdata('sitein_ext','som2e',name)
305   continue
      do 310 ii = 1, MAXIEL
        read(7,*) som2e(SOIL,ii), name
        call ckdata('sitein_ext','som2e',name)
310   continue
      do 315 ii = 1, MAXIEL
        read(7,*) som3e(ii), name
        call ckdata('sitein_ext','som3e',name)
315   continue

c ... Read the fine root values into the mature fine root pools
      read(7,*) bglcism(UNLABL), name
      call ckdata('sitein_ext','bglcis',name)
      read(7,*) bglcism(LABELD), name
      call ckdata('sitein_ext','bglcis',name)

      do 320 ii = 1, MAXIEL
        read(7,*) bglivem(ii), name
        call ckdata('sitein_ext','bglive',name)
320   continue

      read(7,*) frtcism(UNLABL), name
      call ckdata('sitein_ext','frtcis',name)
      read(7,*) frtcism(LABELD), name
      call ckdata('sitein_ext','frtcis',name)

      do 325 ii = 1, MAXIEL
        read(7,*) frootem(ii), name
        call ckdata('sitein_ext','froote',name)
325   continue

      read(7,*) strcis(SRFC,UNLABL), name
      call ckdata('sitein_ext','strcis',name)
      read(7,*) strcis(SRFC,LABELD), name
      call ckdata('sitein_ext','strcis',name)
      read(7,*) strcis(SOIL,UNLABL), name
      call ckdata('sitein_ext','strcis',name)
      read(7,*) strcis(SOIL,LABELD), name
      call ckdata('sitein_ext','strcis',name)
      do 330 ii = 1, MAXIEL
        read(7,*) struce(SRFC,ii), name
        call ckdata('sitein_ext','struce',name)
330   continue
      do 335 ii = 1, MAXIEL
        read(7,*) struce(SOIL,ii), name
        call ckdata('sitein_ext','struce',name)
335   continue

      read(7,*) carbostg(CRPSYS,UNLABL), name
      call ckdata('sitein_ext','carbos',name)
      read(7,*) carbostg(CRPSYS,LABELD), name
      call ckdata('sitein_ext','carbos',name)
      read(7,*) carbostg(FORSYS,UNLABL), name
      call ckdata('sitein_ext','carbos',name)
      read(7,*) carbostg(FORSYS,LABELD), name
      call ckdata('sitein_ext','carbos',name)

      read(7,*) csrsnk(UNLABL), name
      call ckdata('sitein_ext','csrsnk',name)
      read(7,*) csrsnk(LABELD), name
      call ckdata('sitein_ext','csrsnk',name)
      do 340 ii = 1, MAXIEL
        read(7,*) esrsnk(ii), name
        call ckdata('sitein_ext','esrsnk',name)
340   continue

      do 345 ii = 1, MAXIEL
        read(7,*) forstg(ii), name
        call ckdata('sitein_ext','forstg',name)
345   continue

      read(7,*) ammonium, name
      call ckdata('sitein_ext','ammoni',name)
      do 355 ii = 1, SWMAXLYR
        read(7,*) nitrate(ii), name
        call ckdata('sitein_ext','nitrat',name)
355   continue

      do 360 ii = 1, SWMAXLYR
        read(7,*) swcextend(ii), name
        call ckdata('sitein_ext','swcext',name)
360   continue

      close(unit=7)

999   continue

      return
      end
