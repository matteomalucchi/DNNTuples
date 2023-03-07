# DNNTuples for AK4 jets

## Setup
```
# use CMSSW_11_1_0_pre8 which has Puppi V14
cmsrel CMSSW_11_1_0
cd CMSSW_11_1_0/src
cmsenv

# clone this repo into "DeepNTuples" directory
git clone git@github.com:hqucms/DNNTuples.git DeepNTuples -b ak4/UL/11_1_X

scram b -j8
```


## Run
```bash
cd $CMSSW_BASE/src/
cmsenv
voms-proxy-init -rfc -voms cms --valid 168:00
scram b -j24 # if change something the code must be compiled
cd DeepNTuples/Ntupler/test/
cmsRun DeepNtuplizerAK4Puppi.py
```

## Submit jobs via CRAB

**Step 0**: switch to the crab production directory and set up grid proxy, CRAB environment, etc.

```bash
cd $CMSSW_BASE/src/DeepNTuples/Ntupler/run
# set up grid proxy
voms-proxy-init -rfc -voms cms --valid 168:00
# set up CRAB env (must be done after cmsenv)
source /cvmfs/cms.cern.ch/common/crab-setup.sh
```

**Step 1**: use the `crab.py` script to submit the CRAB jobs:

`python crab.py --set-input-dataset -p ../test/DeepNtuplizerAK4[CHS|Puppi].py --site T2_CH_CERN -o /store/user/$USER/DeepNtuples/[version] -t DeepNtuplesAK4-[version] --no-publication -i [ABC].conf -s FileBased -n 1 --work-area crab_projects_[ABC] --send-external [--input_files JEC.db] --dryrun`

python crab.py -p ../test/DeepNtuplizerAK4Puppi.py -o /store/user/mmalucch/dnntuple_output/all/ -i ../test/datasets.conf --site T2_IT_Pisa/T1_IT_CNAF_Disk --work-area crab_projects_[] -n 10000000 --splitting Automatic

python crab.py -p ../test/DeepNtuplizerAK4Puppi.py -o /store/user/mmalucch/dnntuple_output/big/ -i ../test/datasets.conf --site T1_IT_CNAF_Disk --work-area crab_projects_big -n 10000000 --splitting Automatic

python crab.py -p ../test/DeepNtuplizerAK4Puppi.py -o /store/user/mmalucch/dnntuple_output/big/ -i ../test/datasets.conf --site T1_IT_CNAF_Disk --work-area crab_projects_big -n 300 --splitting Automatic --max-units 10000000


These command will perform a "dryrun" to print out the CRAB configuration files. Please check everything is correct (e.g., the output path, version number, requested number of cores, etc.) before submitting the actual jobs. To actually submit the jobs to CRAB, just remove the `--dryrun` option at the end.

**[Note] For the QCD samples use `-n 1 --max-units 50` to run one file per job, and limit the total files per job to 50.**


**Step 2**: check job status

The status of the CRAB jobs can be checked with:

```bash
./crab.py --status --work-area crab_projects_[ABC]
```

Note that this will also resubmit failed jobs automatically.

The crab dashboard can also be used to get a quick overview of the job status:
`https://dashb-cms-job.cern.ch/dashboard/templates/task-analysis`

More options of this `crab.py` script can be found with:

```bash
./crab.py -h
```
