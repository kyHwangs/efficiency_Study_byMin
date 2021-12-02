

python submit_batch.py

condor_submit condor_submit.sub

condor_hold msoh

condor_release msoh

python submit_batch.py recover

python submit_batch.py hadd >&hadd.log&


ll -Sr condor/job.*.err

ll -S Outputs_v51*/*/*/*.root

find ./Outputs_v51*/* \
-type f -name '*--*.root' \
-size -10k \
-delete




