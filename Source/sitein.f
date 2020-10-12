
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... SITEIN.F

      subroutine sitein(ivopt)

      implicit none
      include 'const.inc'
      include 'doubles.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'site.inc'
      include 'wth.inc'

c ... Argument declarations
      logical   ivopt

c ... Read the parameter file in the form needed by time0.

c ... Local variables
      integer     ii, jj
      character*8 name
      real        temp

      read(7,*)
      read(7,*)

      do 10 ii = 1, MONTHS
        read(7,*) precip(ii), name
        call ckdata('sitein','precip',name)
10    continue

      do 20 ii = 1, MONTHS
        read(7,*) prcstd(ii), name
        call ckdata('sitein','prcstd',name)
20    continue

      do 30 ii = 1, MONTHS
        read(7,*) prcskw(ii), name
        call ckdata('sitein','prcskw',name)
30    continue

      do 40 ii = 1, MONTHS
        read(7,*) tmn2m(ii), name
        call ckdata('sitein','tmn2m',name)
40    continue

c ... maxt calculation for nitrify added 4/03/98 - mdh
      maxt = -99.0
      do 50 ii = 1, MONTHS
        read(7,*) tmx2m(ii), name
        if (tmx2m(ii) .gt. maxt) maxt = tmx2m(ii)
        call ckdata('sitein','tmx2m',name)
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
      call ckdata('sitein','ivauto',name)
      read(7,*) temp, name
      nelem = int(temp)
      call ckdata('sitein','nelem',name)
 
      read(7,*) sitlat, name
      call ckdata('sitein','sitlat',name)
      read(7,*) sitlng, name
      call ckdata('sitein','sitlng',name)

      read(7,*) sand, name
      call ckdata('sitein','sand',name)
      read(7,*) silt, name
      call ckdata('sitein','silt',name)
      read(7,*) clay, name
      call ckdata('sitein','clay',name)
c ... Add rock fraction, cak - 05/27/03
      read(7,*) rock, name
      call ckdata('sitein','rock',name)
      if (rock .gt. 0.90) then 
        write(*,*) 'Rock fraction too large, rock = ', rock
        STOP
      endif
      read(7,*) bulkd, name
      call ckdata('sitein','bulkd',name)

      read(7,*) temp, name
      nlayer = int(temp)
      if (nlayer .gt. 9) then
        nlayer = 9
        call message('   Warning: nlayer value too large, reset to 9')
      endif
      call ckdata('sitein','nlayer',name)
      read(7,*) temp, name
      nlaypg = int(temp)
      call ckdata('sitein','nlaypg',name)
      read(7,*) drain, name
      call ckdata('sitein','drain',name)
      read(7,*) basef, name
      call ckdata('sitein','basef', name)
      read(7,*) stormf, name
      call ckdata('sitein','stormf', name)
c ... Add precipitation amount required for runoff and fraction of
c ... precipitation above amount required for runoff which is lost
c ... via runoff, cak - 05/27/03
      read(7,*) precro, name
      call ckdata('sitein','precro', name)
      read(7,*) fracro, name
      call ckdata('sitein','fracro', name)
      read(7,*) temp, name
      swflag = int(temp)
      call ckdata('sitein','swflag', name)

      do 60 ii = 1, CMXLYR
        read(7,*) awilt(ii), name
        call ckdata('sitein','awilt', name)
60    continue

      do 70 ii = 1, CMXLYR
        read(7,*) afiel(ii), name
        call ckdata('sitein','afiel', name)
70    continue

      read(7,*) ph, name
      call ckdata('sitein','ph',name)
c ... New phstart variable added for pH shift, cak - 08/02/02
      phstart = ph
      read(7,*) pslsrb, name
      call ckdata('sitein','pslsrb',name)
      read(7,*) sorpmx, name
      call ckdata('sitein','sorpmx',name)
      read(7,*)
      read(7,*) epnfa(INTCPT), name
      call ckdata('sitein','epnfa',name)
      read(7,*) epnfa(SLOPE), name
      call ckdata('sitein','epnfa',name)
      read(7,*) epnfs(INTCPT), name
      call ckdata('sitein','epnfs',name)
      read(7,*) epnfs(SLOPE), name
      call ckdata('sitein','epnfs',name)
      read(7,*) satmos(INTCPT), name
      call ckdata('sitein','satmos',name)
      read(7,*) satmos(SLOPE), name
      call ckdata('sitein','satmos',name)
      read(7,*) sirri, name
      call ckdata('sitein','sirri',name)

c ... If extending, do not read in initial conditions
      if (ivopt) then
        goto 999
      endif

      read(7,*)
      read(7,*) som1ci(SRFC,UNLABL), name
      call ckdata('sitein','som1ci',name)
      read(7,*) som1ci(SRFC,LABELD), name
      call ckdata('sitein','som1ci',name)
      read(7,*) som1ci(SOIL,UNLABL), name
      call ckdata('sitein','som1ci',name)
      read(7,*) som1ci(SOIL,LABELD), name
      call ckdata('sitein','som1ci',name)

      read(7,*) som2ci(SRFC,UNLABL), name
      call ckdata('sitein','som2ci',name)
      read(7,*) som2ci(SRFC,LABELD), name
      call ckdata('sitein','som2ci',name)
      read(7,*) som2ci(SOIL,UNLABL), name
      call ckdata('sitein','som2ci',name)
      read(7,*) som2ci(SOIL,LABELD), name
      call ckdata('sitein','som2ci',name)

      read(7,*) som3ci(UNLABL), name
      call ckdata('sitein','som3ci',name)
      read(7,*) som3ci(LABELD), name
      call ckdata('sitein','som3ci',name)

      do 90 ii = SRFC, SOIL
        do 80 jj = 1, MAXIEL
          read(7,*) rces1(ii,jj), name
          call ckdata('sitein','rces1',name)
80      continue
90    continue

      do 105 ii = SRFC, SOIL
        do 100 jj = 1, MAXIEL
          read(7,*) rces2(ii,jj), name
          call ckdata('sitein','rces2',name)
100     continue
105   continue

      do 110 ii = 1, MAXIEL
        read(7,*) rces3(ii), name
        call ckdata('sitein','rces3',name)
110   continue

      read(7,*) clittr(SRFC,UNLABL), name
      call ckdata('sitein','clittr',name)
      read(7,*) clittr(SRFC,LABELD), name
      call ckdata('sitein','clittr',name)
      read(7,*) clittr(SOIL,UNLABL), name
      call ckdata('sitein','clittr',name)
      read(7,*) clittr(SOIL,LABELD), name
      call ckdata('sitein','clittr',name)
c ... Add check for initial litter values to prevent divide by zero error
c ... in calciv subroutine
      if ((ivauto .eq. 0) .and. 
     &    (clittr(SRFC,UNLABL) + clittr(SRFC,LABELD) .le. 0.0)) then
        call message('   Initial surface litter values, clittr(1,*),')
        call message('   must be greater than zero.')
        STOP
      endif
      if ((ivauto .eq. 0) .and. 
     &    (clittr(SOIL,UNLABL) + clittr(SOIL,LABELD) .le. 0.0)) then
        call message('   Initial soil litter values, clittr(2,*),')
        call message('   must be greater than zero.')
        STOP
      endif

      do 130 ii = SRFC, SOIL
        do 120 jj = 1, MAXIEL
          read(7,*) rcelit(ii,jj), name
          call ckdata('sitein','rcelit',name)
120     continue
130   continue

      read(7,*) aglcis(UNLABL), name
      call ckdata('sitein','aglcis',name)
      read(7,*) aglcis(LABELD), name
      call ckdata('sitein','aglcis',name)

      do 140 ii = 1, MAXIEL
        read(7,*) aglive(ii), name
        call ckdata('sitein','aglive',name)
140   continue

c ... Read the fine root values into the juvenile fine root pool, these
c ... values will be split between the fine root and mature root pool
c ... based on the fraction of fine root production that goes to mature
c ... roots parameter value as read from the crop.100 file in the detiv
c ... subroutine, cak - 05/24/2007
      read(7,*) bglcisj(UNLABL), name
      call ckdata('sitein','bglcis',name)
      read(7,*) bglcisj(LABELD), name
      call ckdata('sitein','bglcis',name)

      do 150 ii = 1, MAXIEL
        read(7,*) bglivej(ii), name
        call ckdata('sitein','bglive',name)
150   continue

      read(7,*) stdcis(UNLABL), name
      call ckdata('sitein','stdcis',name)
      read(7,*) stdcis(LABELD), name
      call ckdata('sitein','stdcis',name)

      do 160 ii = 1, MAXIEL
        read(7,*) stdede(ii), name
        call ckdata('sitein','stdede',name)
160   continue

      read(7,*)
      read(7,*) rlvcis(UNLABL), name
      call ckdata('sitein','rlvcis',name)
      read(7,*) rlvcis(LABELD), name
      call ckdata('sitein','rlvcis',name)

      do 170 ii = 1, MAXIEL
        read(7,*) rleave(ii), name
        call ckdata('sitein','rleave',name)
170   continue

      read(7,*) fbrcis(UNLABL), name
      call ckdata('sitein','fbrcis',name)
      read(7,*) fbrcis(LABELD), name
      call ckdata('sitein','fbrcis',name)

      do 180 ii = 1, MAXIEL
        read(7,*) fbrche(ii), name
        call ckdata('sitein','fbrche',name)
180   continue

      read(7,*) rlwcis(UNLABL), name
      call ckdata('sitein','rlwcis',name)
      read(7,*) rlwcis(LABELD), name
      call ckdata('sitein','rlwcis',name)

      do 190 ii = 1, MAXIEL
        read(7,*) rlwode(ii), name
        call ckdata('sitein','rlwode',name)
190   continue

c ... Read the fine root values into the juvenile fine root pool, these
c ... values will be split between the fine root and mature root pool
c ... based on the fraction of fine root production that goes to mature
c ... roots parameter value as read from the tree.100 file in the detiv
c ... subroutine, cak - 05/24/2007
      read(7,*) frtcisj(UNLABL), name
      call ckdata('sitein','frtcis',name)
      read(7,*) frtcisj(LABELD), name
      call ckdata('sitein','frtcis',name)

      do 200 ii = 1, MAXIEL
        read(7,*) frootej(ii), name
        call ckdata('sitein','froote',name)
200   continue

      read(7,*) crtcis(UNLABL), name
      call ckdata('sitein','crtcis',name)
      read(7,*) crtcis(LABELD), name
      call ckdata('sitein','crtcis',name)

      do 210 ii = 1, MAXIEL
        read(7,*) croote(ii), name
        call ckdata('sitein','croote',name)
210   continue

      read(7,*) wd1cis(UNLABL), name
      call ckdata('sitein','wd1cis',name)
      read(7,*) wd1cis(LABELD), name
      call ckdata('sitein','wd1cis',name)
      read(7,*) wd2cis(UNLABL), name
      call ckdata('sitein','wd2cis',name)
      read(7,*) wd2cis(LABELD), name
      call ckdata('sitein','wd2cis',name)
      read(7,*) wd3cis(UNLABL), name
      call ckdata('sitein','wd3cis',name)
      read(7,*) wd3cis(LABELD), name
      call ckdata('sitein','wd3cis',name)

c ... Six new pools for standing dead tree biomass.
c ... These are the carbon pools. Elemental
c ... pools are initialize in calciv. 
c ... -mdh 9/21/2018
      read(7,*) dlvcis(UNLABL), name
      call ckdata('sitein','dlvcis',name)
      read(7,*) dlvcis(LABELD), name
      call ckdata('sitein','dlvcis',name)
      read(7,*) dfbrcis(UNLABL), name
      call ckdata('sitein','dfbrcis',name)
      read(7,*) dfbrcis(LABELD), name
      call ckdata('sitein','dfbrcis',name)
      read(7,*) dlwcis(UNLABL), name
      call ckdata('sitein','dlwcis',name)
      read(7,*) dlwcis(LABELD), name
      call ckdata('sitein','dlwcis',name)

c      read(7,*) w1lig, name
c      call ckdata('sitein','w1lig',name)
c      read(7,*) w2lig, name
c      call ckdata('sitein','w2lig',name)
c      read(7,*) w3lig, name
c      call ckdata('sitein','w3lig',name)

      read(7,*)

      do 230 ii = 1, MAXIEL
        do 220 jj = 1, CMXLYR
          read(7,*) minerl(jj,ii), name
          call ckdata('sitein','minerl',name)
220     continue
230   continue

      do 240 ii = 1, MAXIEL
        read(7,*) parent(ii), name
        call ckdata('sitein','parent',name)
240   continue

c ... The secndy and occlud input values can now be double precision,
c ... cak - 03/20/02
      do 250 ii = 1, MAXIEL
        read(7,*) secndy_double(ii), name
        call ckdata('sitein','secndy',name)
250   continue

      read(7,*) occlud_double, name
      call ckdata('sitein','occlud',name)
      read(7,*)

c ... Save the double precision secndy and occlud variables read into their
c ... single precision counterparts, cak - 03/20/02
      secndy(1) = real(secndy_double(1))
      secndy(2) = real(secndy_double(2))
      secndy(3) = real(secndy_double(3))
      occlud = real(occlud_double)

      do 260 ii = 1, CMXLYR
        read(7,*) rwcf(ii), name
        call ckdata('sitein','rwcf',name)
260   continue

      read(7,*) snlq, name
      call ckdata('sitein','snlq',name)
      read(7,*) snow, name
      call ckdata('sitein','snow', name)

      close(unit=7)

999   continue

      return
      end
