
c               Copyright 1993 Colorado State University
c                       All Rights Reserved
c
c Added parameters for standing dead tree biomass. -mdh 9/20/2018
c   dlvfalrt - fall rate of dead attached leaves (fraction per month, in absence of disturbance) 
c   dfbfalrt - fall rate of dead attached branches (fraction per month, in absence of disturbance) 
c   dlwfalrt - fall rate of dead standing large wood (fraction per month, in absence of disturbance) 
c   decw4 - intrinsic decomposition rate constant for attached dead fine branches (yr-1)
c   decw5 - intrinsic decomposition rate constant for standing dead (yr-1) 

      common/parfs/basfct,btolai,ccefor(2,5,3),cerfor(3,5,3),
     &    decid,decw1,decw2,decw3,decw4,decw5,
     &    dfbfalrt,dlvfalrt,dlwfalrt,
     &    fcfrac(5,2),forgrw,forrtf(3),
     &    klai,laitop,ldrmlt,leafdr(12),
     &    maxlai,maxldr,maxnp,sapk,swold,
     &    tfrtcn(2),tfrtcw(2),tmxturn,
     &    wdlig(6),wmrtfrac,woodb,wooddr(6),wrdsrfc

      integer decid,forgrw

      real basfct,btolai,ccefor,cerfor,decw1,decw2,decw3,
     &     decw4,decw5,dfbfalrt,dlvfalrt,dlwfalrt,
     &     fcfrac,forrtf,klai,laitop,ldrmlt,leafdr,
     &     maxlai,maxldr,maxnp,sapk,swold,
     &     tfrtcn,tfrtcw,tmxturn,
     &     wdlig,wmrtfrac,woodb,wooddr,wrdsrfc

      save /parfs/
