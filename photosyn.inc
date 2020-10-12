
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      common/photosyn/aMax(2), aMaxFrac(2),
     &                aMaxScalar1(2), aMaxScalar2(2), aMaxScalar3(2),
     &                aMaxScalar4(2), attenuation(2),
     &                baseFolRespFrac(2), cFracLeaf(2), dVpdExp(2),
     &                dVpdSlope(2), growthDays1(2), growthDays2(2),
     &                growthDays3(2), growthDays4(2), halfSatPar(2),
     &                leafCSpWt(2), psnTMin(2), psnTOpt(2)

c ... aMax            - maximum photosynthesis
c ...                   (nmol CO2 * g^-1 leaf * sec^-1)
c ...                   assuming maximum possible par, all intercepted,
c ...                   no temperatre, water or vpd stress
c ... aMaxFrac        - average daily aMax as fraction of instantaneous
c ... aMaxScalar1     - multiplier used to adjust aMax based on
c ...                   growthDays1 days since germination
c ... aMaxScalar2     - multiplier used to adjust aMax based on
c ...                   growthDays2 days since germination
c ... aMaxScalar3     - multiplier used to adjust aMax based on
c ...                   growthDays3 days since germination
c ... aMaxScalar4     - multiplier used to adjust aMax based on
c ...                   growthDays4 days since germination
c ... attenuation     - light attenuation coefficient
c ... baseFolRespFrac - basal foliage respiration rate, as % of maximum
c ...                   net photosynthesis rate
c ... cFracLeaf       - factor for converting leaf biomass to carbon
c ...                   (leaf biomass * cFracLeaf = leaf carbon)
c ... dVpdExp         - exponential value in dVpd equation
c ... dVpdSlope       - slope value in dVpd equation
c ... growthDays1     - number of days after germination to start using
c ...                   aMaxScalar1
c ... growthDays2     - number of days after germination to start using
c ...                   aMaxScalar2
c ... growthDays3     - number of days after germination to start using
c ...                   aMaxScalar3
c ... growthDays4     - number of days after germination to start using
c ...                   aMaxScalar4
c ... halfSatPar      - par at which photosynthesis occurs at 1/2
c ...                   theoretical maximum
c ...                   (Einsteins * m^-2 ground area * day^-1)
c ... leafCSpWt       - grams of carbon in a square meter of leaf area
c ... psnTMin         - minimum tempature at which net photosynthesis
c ...                   occurs (degrees C)
c ... psnTOpt         - optimal temperature at which net photosynthesis
c ....                  occurs (degrees C)

      double precision aMax, aMaxFrac, aMaxScalar1, aMaxScalar2,
     &                 aMaxScalar3, aMaxScalar4, attenuation,
     &                 baseFolRespFrac, cFracLeaf, dVpdExp, dVpdSlope,
     &                 growthDays1, growthDays2, growthDays3,
     &                 growthDays4, halfSatPar, leafCSpWt,
     &                 psnTMin, psnTOpt

      save /photosyn/
