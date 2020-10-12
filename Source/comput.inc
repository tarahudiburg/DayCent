
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Variables used internally which must be shared among routines.
c ...
c ... Change p1co2 to p1co2(2) to represent surface and soil layers.
c ... Change ratnew to ratnew1 and ratnew2 for surface and soil
c ... decomposition
c ...
c ... baseNdep - base amount of atmospheric N deposition for the year
c ...            based on average annual precipitation only
c ... cercrp(mx/mn, nparts, nelem)
c ... fps1s3 - the fraction of som1 decomposing to som3
c ... fps2s3 - the fraction of som2 decomposing to som3
c ... lhzci(pool,iso) - the lower horizon carbon used as an input 
c ...                   to the surface soil in an erosion event
c ... lhze(pool,iel) - the lower horizon N,P,S used as an input 
c ...                  to the surface soil in an erosion event
c ... orglch - the fraction of organics that leach from soil som1 when
c ...          there is sufficient water flow.  Computed in predec.
c ... ratnew1(iel,1) - the C/E ratio for new material created when a
c ...                  lignin component decomposes to SRFC som1.
c ... ratnew2(iel,1) - the C/E ratio for new material created when a
c ...                  lignin component decomposes to SOIL som1.
c ... ratnew1(iel,2) - the C/E ratio for new material created when a
c ...                  lignin component decomposes to SRFC som2.
c ... ratnew2(iel,2) - the C/E ratio for new material created when a
c ...                  lignin component decomposes to SOIL som2.
c ... wc - the water capacity for the top layer afiel(1)-awilt(1)

      common/comput/agdefacm(12),bgdefacm(12),baseNdep,
     &  cercrp(2,2,3),eftext,fps1s3,fps2s3,
     &  lhzci(3,2),lhze(3,3),orglch,p1co2(2),h2ogef(3),
     &  ratnew1(3,2),ratnew2(3,2),wc

      real agdefacm,bgdefacm,baseNdep,cercrp,
     &  eftext,fps1s3,fps2s3,lhzci,lhze,orglch,p1co2,h2ogef,
     &  ratnew1,ratnew2,wc

      save /comput/
