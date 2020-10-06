testfiles_DAYCENT/
  * Original files from Jeff Stenzel with a few modifications.
  * I spilt the schedule file into a 1976-year spinup then another for 1977-2160.
  * I updated outvars.txt
  * Uses DailyDayCent_muvp.exe (unmodified from the version I recieved)
  * Represents TREM events at evntyp=0 (cutting)

testfiles_DAYCENT_snag_no_fire/
  * See R script compareLVresults.R that plots results from testfiles_DAYCENT and testfiles_DAYCENT_snag_no_fire.
    If you don't have perl installed on your computer you won't be able to convert the .lis file to a .csv file 
    so you will need to change the read.csv() to read.table().
  * I updated outvars.txt to include the new output variables for the snag model
  * Uses DailyDayCent_muvps.exe (muvp + snags)
  * Represents TREM events at evntyp=0 (cutting)
  * All things being equal, when large wood that remains standing vs. falling to the
    ground, one should expect to see N immobilization, at least initially.

testfiles_DAYCENT_snag/
  * See R script compareLVresults.R that plots results from testfiles_DAYCENT_snag_no_fire and testfiles_DAYCENT_snag (with fire).
    If you don't have perl installed on your computer you won't be able to convert the .lis file to a .csv file 
    so you will need to change the read.csv() to read.table().
  * I updated outvars.txt to include the new output variables for the snag model
  * Uses DailyDayCent_muvps.exe (muvp + snags)
  * Represents TREM events at evntyp=1 (FIRE). The fire.100 file was updated to burn wood1c and wood2c.
  * All things being equal, fire returns some C as charcoal, increasing som3c, and 
    elemental returns increase minerl N increasing productivity.

See comments in the fire.100, tree.100, and trem.100 files.

The examples above use an evergreen tree (DECID=0).  I tested the model for the same tree as a decidous type.
The big difference between evergreen and deciduous is that when deciduous leaves die, they
immdeiately drop to the ground instead of becoming dead attached leaves which may fall to the ground
at a slower rate.  This had feedbacks on soil litter, mineral N, and production. Dead attached leaves lose 
some mass (C and N) from photodegration but this amount is small.
