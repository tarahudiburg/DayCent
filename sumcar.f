
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine sumcar

      implicit none
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'

c ... Sum unlabeled and labeled carbon to get totals.

      strucc(1)=strcis(1,1)+strcis(1,2)
      strucc(2)=strcis(2,1)+strcis(2,2) 
      metabc(1)=metcis(1,1)+metcis(1,2)
      metabc(2)=metcis(2,1)+metcis(2,2)

c ... Sum som1 surface and soil isotopes separately.  vek 08-91
      som1c(1)=som1ci(1,1)+som1ci(1,2)
      som1c(2)=som1ci(2,1)+som1ci(2,2)
      som2c(1)=som2ci(1,1)+som2ci(1,2)
      som2c(2)=som2ci(2,1)+som2ci(2,2)
      som3c=som3ci(1)+som3ci(2)

      wood1c = wd1cis(1) + wd1cis(2)
      wood2c = wd2cis(1) + wd2cis(2)
      wood3c = wd3cis(1) + wd3cis(2)

      somtc=som1c(2) + som2c(2) + som3c + strucc(2) + metabc(2)
      aglivc=aglcis(1)+aglcis(2)
      stdedc=stdcis(1)+stdcis(2)
      bglivcj=bglcisj(1)+bglcisj(2)
      bglivcm=bglcism(1)+bglcism(2)
      agcacc=agcisa(1)+agcisa(2)
      bgcjacc=bgcisja(1)+bgcisja(2)
      bgcmacc=bgcisma(1)+bgcisma(2)

      rleavc = rlvcis(1) + rlvcis(2)
      frootcj = frtcisj(1) + frtcisj(2)
      frootcm = frtcism(1) + frtcism(2)
      fbrchc = fbrcis(1) + fbrcis(2)
      rlwodc = rlwcis(1) + rlwcis(2)
      crootc = crtcis(1) + crtcis(2)

c ... New standing dead forest components. -mdh 9/20/2018
      dleavc = dlvcis(1) + dlvcis(2)
      dfbrchc = dfbrcis(1) + dfbrcis(2)
      dlwodc = dlwcis(1) + dlwcis(2)

      rlvacc = alvcis(1) + alvcis(2)
      frtjacc = afrcisj(1) + afrcisj(2)
      frtmacc = afrcism(1) + afrcism(2)
      fbracc = afbcis(1) + afbcis(2)
      rlwacc = alwcis(1) + alwcis(2)
      crtacc = acrcis(1) + acrcis(2)
      fcacc  = rlvacc + frtjacc + frtmacc + fbracc + rlwacc + crtacc

      return
      end
