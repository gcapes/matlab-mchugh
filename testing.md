## Running on the CSF (i.e. compiled)
### From /scratch
- The output file is written using parfor for the smaller problem size only. Sample size=2; population=50. 
This works for R2014a, R2015a, R2016a.
- The output file is not written for the larger problem size: samplesize=100; population=1600. 
Same MATLAB version numbers tested.

### From /home
- R2014a and full size problem doesn't write the output `.mat` file

## Running on desktop
### Using GUI
- File is written for full size problem, using R2016a.

### Running compiled on desktop
- Full size problem using R2016a writes out the `.mat` file.

## Summary of Desktop and CSF
The full size problem runs correctly on my desktop whether using GUI or compiled. 
It doesn't run on the CSF for the full problem size, only the smaller size.
File system and MATLAB version don't look to affect the degree of success on the CSF.

## DPSF
Running full size problem using R2015a, from scratch, does yield the output file.
DPSF has a newer OS than CSF - this is the main difference.
