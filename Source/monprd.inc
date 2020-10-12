
c               Copyright 1996 Colorado State University
c                       All Rights Reserved


c ... Monthly production variables. 

      common/monprd/nit_amt_year, nit_amt_month,
     &              mcprd(3), mfprd(6),
     &              mrspdyflux(2), cmrspdyflux(3), fmrspdyflux(6),
     &              grspdyflux(2), cgrspdyflux(3), fgrspdyflux(6),
     &              mrspTempEffect(2,2), mrspWaterEffect(2),
     &              N2O_year, NO_year, N2_year,
     &              ch4_oxid_year, stempmth, annppt,
     &              N2O_month, NO_month, N2_month, ch4_oxid_month,
     &              pptmonth, cwstress, gwstress,
     &              twstress

      double precision nit_amt_year, nit_amt_month

      real mcprd, mfprd,
     &     mrspdyflux, cmrspdyflux, fmrspdyflux,
     &     grspdyflux, cgrspdyflux, fgrspdyflux,
     &     mrspTempEffect, mrspWaterEffect,
     &     N2O_year, NO_year, N2_year, ch4_oxid_year, stempmth, annppt,
     &     N2O_month, NO_month, N2_month, ch4_oxid_month,
     &     pptmonth, cwstress, gwstress, twstress

      save /monprd/
