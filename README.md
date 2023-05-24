# Histogram / Plotter for MuonHLT

## Making Histograms from CRAB Outputs
```
# Outputs from https://github.com/wonpoint4/MuonHLTTool

cd Analyzers
python submit_batch.py
condor_submit condor_submit.sub
```

## Hadd Hisograms
```
python submit_batch.py hadd
condor_submit condor_hadd.sub
```

## Draw Plots
```
# Need to set histograms's directory, names correctly

root -l -b -q drawtnpCompEffL3wrtL1.C
root -l -b -q drawtnpCompEffL3wrtOff.C
```

## Old Commands (but could be useful)
```
condor_hold msoh

condor_release msoh

python submit_batch.py recover

ll -Sr condor/job.*.err

ll -S Outputs_v*/*/*/*.root

find ./Outputs_v*/* \
-type f -name '*--*.root' \
-size -10k \
-delete

scp msoh@cms.knu.ac.kr:/u/user/msoh/MuonHLT/Run3/Local/test_20211201_Run3Review/Outputs_vRun3Review*/Eff/hist-*.root .

cdq 397147 -long | grep -v Environment | grep CPU
```
