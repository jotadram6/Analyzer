# Setting up

To make the code run, simply type
```
setenv SCRAM_ARCH slc6_amd64_gcc491
source /cvmfs/cms.cern.ch/cmsset_default.csh 
cmsrel CMSSW_7_4_15
cd CMSSW_7_4_15/src
cmsenv
git clone https://github.com/dteague/Analyzer
cd Analyzer
git branch TNT74x
make
``` 
This will create the file ```Analyzer``` 

# Running the code
Running is simply typing
```
./Analyzer <INFILE> <OUTFILE>
```
