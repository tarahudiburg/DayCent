
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Variables used to track carbon flows from decomposition
c ... identifying source and destination.
c ...
c ... ametc1tosom11  - annual accumulator for carbon flow from surface
c ...                  metabolic pool to fast surface organic matter
c ...                  pool
c ... ametc2tosom12  - annual accumulator for  carbon flow from soil
c ...                  metabolic pool to fast soil organic matter pool
c ... asom11tosom21  - annual accumulator for carbon flow from fast
c ...                  surface organic matter pool to intermediate
c ...                  surface organic matter pool
c ... asom12tosom22  - annual accumulator for carbon flow from fast soil
c ...                  organic matter pool to intermediate soil organic
c ...                  matter pool
c ... asom12tosom3   - annual accumulator for carbon flow from fast soil
c ...                  organic matter pool to slow soil organic matter
c ...                  pool
c ... asom21tosom11  - annual accumulator for carbon flow from
c ...                  intermediate surface organic matter pool to fast
c ...                  surface organic matter pool
c ... asom21tosom22  - annual accumulator for carbon flow from
c ...                  intermediate surface organic matter pool to
c ...                  intermediate soil organic matter pool
c ... asom22tosom12  - annual accumulator for carbon flow from
c ...                  intermediate soil organic matter pool to fast
c ...                  soil organic matter pool
c ... asom22tosom3   - annual accumulator for carbon flow from
c ...                  intermediate soil organic matter pool to slow
c ...                  soil organic matter pool
c ... asom3tosom12   - annual accumulator for carbon flow from slow soil
c ...                  organic matter pool to fast soil organic matter
c ...                  pool
c ... astruc1tosom11 - annual accumulator for carbon flow from surface
c ...                  structural pool to fast surface organic matter
c ...                  pool
c ... astruc1tosom21 - annual accumulator for carbon flow from surface
c ...                  structural pool to intermediate surface organic
c ...                  matter pool
c ... astruc2tosom12 - annual accumulator for carbon flow from soil
c ...                  structural pool to fast soil organic matter pool
c ... astruc2tosom22 - annual accumulator for carbon flow from soil
c ...                  structural pool to intermediate soil organic
c ...                  matter pool
c ... awood1tosom11  - annual accumulator for carbon flow from dead
c ...                  fine branch pool to fast surface organic matter pool
c ... awood1tosom21  - annual accumulator for carbon flow from dead fine
c ...                  branch pool pool to intermediate surface organic matter
c ...                  pool
c ... awood2tosom11  - annual accumulator for carbon flow from dead large wood
c ...                  pool to fast surface organic matter pool
c ... awood2tosom21  - annual accumulator for carbon flow from dead large wood
c ...                  pool to intermediate surface organic matter pool
c ... awood3tosom12  - annual accumulator for carbon flow from dead coarse
c ...                  root pool to fast soil organic matter pool
c ... awood3tosom22  - annual accumulator for carbon flow from dead coarse
c ...                  root pool to intermediate soil organic matter pool
c ... metc1tosom11   - carbon flow from surface metabolic pool to fast
c ...                  surface organic matter pool
c ... metc2tosom12   - carbon flow from soil metabolic pool to fast soil
c ...                  organic matter pool
c ... som11tosom21   - carbon flow from fast surface organic matter pool
c ...                  to intermediate surface organic matter pool
c ... som12tosom22   - carbon flow from fast soil organic matter pool
c ...                  to intermediate soil organic matter pool
c ... som12tosom3    - carbon flow from fast soil organic matter pool to
c ...                  slow soil organic matter pool
c ... som21tosom11   - carbon flow from intermediate surface organic
c ...                  matter pool to fast surface organic matter pool
c ... som21tosom22   - carbon flow from intermediate surface organic
c ...                  matter pool to intermediate soil organic matter
c ...                  pool
c ... som22tosom12   - carbon flow from intermediate soil organic matter
c ...                  pool to fast soil organic matter pool
c ... som22tosom3    - carbon flow from intermediate soil organic matter
c ...                  pool to slow soil organic matter pool
c ... som3tosom12    - carbon flow from slow soil organic matter pool to
c ...                  fast soil organic matter pool
c ... struc1tosom11  - carbon flow from surface structural pool to fast
c ...                  surface organic matter pool
c ... struc1tosom21  - carbon flow from surface structural pool to
c ...                  intermediate surface organic matter pool
c ... struc2tosom12  - carbon flow from soil structural pool to fast
c ...                  soil organic matter pool
c ... struc2tosom22  - carbon flow from soil structural pool to
c ...                  intermediate soil organic matter pool
c ... wood1tosom11   - carbon flow from dead fine branch pool to fast
c ...                  surface organic matter pool
c ... wood1tosom21   - carbon flow from dead fine branch pool to
c ...                  intermediate surface organic matter pool
c ... wood2tosom11   - carbon flow from dead large wood pool to fast
c ...                  surface organic matter pool
c ... wood2tosom21   - carbon flow from dead large wood pool pool to
c ...                  intermediate surface organic matter pool
c ... wood3tosom12   - carbon flow from dead coarse root pool to fast
c ...                  soil organic matter pool
c ... wood3tosom22   - carbon flow from dead coarse root pool to
c ...                  intermediate soil organic matter pool

      common/cflows/ametc1tosom11, ametc2tosom12,
     &  asom11tosom21, asom12tosom22, asom12tosom3, asom21tosom11,
     &  asom21tosom22, asom22tosom12, asom22tosom3, asom3tosom12,
     &  astruc1tosom11, astruc1tosom21, astruc2tosom12, astruc2tosom22,
     &  awood1tosom11, awood1tosom21, awood2tosom11, awood2tosom21,
     &  awood3tosom12, awood3tosom22,
     &  metc1tosom11, metc2tosom12,
     &  som11tosom21, som12tosom22, som12tosom3, som21tosom11,
     &  som21tosom22, som22tosom12, som22tosom3, som3tosom12,
     &  struc1tosom11, struc1tosom21, struc2tosom12, struc2tosom22,
     &  wood1tosom11, wood1tosom21, wood2tosom11, wood2tosom21,
     &  wood3tosom12, wood3tosom22

      real ametc1tosom11, ametc2tosom12,
     &  asom11tosom21, asom12tosom22, asom12tosom3, asom21tosom11,
     &  asom21tosom22, asom22tosom12, asom22tosom3, asom3tosom12,
     &  astruc1tosom11, astruc1tosom21, astruc2tosom12, astruc2tosom22,
     &  awood1tosom11, awood1tosom21, awood2tosom11, awood2tosom21,
     &  awood3tosom12, awood3tosom22,
     &  metc1tosom11, metc2tosom12,
     &  som11tosom21, som12tosom22, som12tosom3, som21tosom11,
     &  som21tosom22, som22tosom12, som22tosom3, som3tosom12,
     &  struc1tosom11, struc1tosom21, struc2tosom12, struc2tosom22,
     &  wood1tosom11, wood1tosom21, wood2tosom11, wood2tosom21,
     &  wood3tosom12, wood3tosom22

      save /cflows/
