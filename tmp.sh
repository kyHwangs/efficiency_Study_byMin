

python submit_batch.py

condor_submit condor_submit.sub

#condor_hold msoh

#condor_release msoh

#python submit_batch.py recover

python submit_batch.py hadd

condor_submit condor_hadd.sub



ll -Sr condor/job.*.err

ll -S Outputs_v*/*/*/*.root

find ./Outputs_v*/* \
-type f -name '*--*.root' \
-size -10k \
-delete

scp msoh@cms.knu.ac.kr:/u/user/msoh/MuonHLT/Run3/Local/test_20211201_Run3Review/Outputs_vRun3Review*/Eff/hist-*.root .

cdq 397147 -long | grep -v Environment | grep CPU
