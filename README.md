# Setting up

To make the code run, simply type
```
setenv SCRAM_ARCH slc6_amd64_gcc530
source /cvmfs/cms.cern.ch/cmsset_default.csh 
cmsrel CMSSW_8_0_10
cd CMSSW_8_0_10/src
cmsenv
git clone https://github.com/dteague/Analyzer
cd Analyzer
git branch TNT80x
make
``` 
This will create the file ```Analyzer``` 

# Running the code
Running is simply typing
```
./Analyzer <INFILE> <OUTFILE>
```
