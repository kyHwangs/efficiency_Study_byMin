#!/usr/bin/env python3
import sys
import os
import time
import glob
import subprocess
from shutil import copyfile
import gc

def jobSpritting( path, nfiles, prefix = "" ):
    str_dcap  = ""
    out = []

    lines = glob.glob(path+"/"+prefix+"*.root")
    lines.sort(key=os.path.getmtime)

    n = 0
    files = ''
    for i, line in enumerate(lines):
        if not line.endswith(".root"):
            continue

        files += '\\"%s%s\\",' % (str_dcap, line)

        if i == nfiles*(n+1) - 1:
            filesout = '\'{'+files+'}\''
            filesout = filesout.replace(',}', '}')
            out.append( ( n, filesout) )
            n = n+1
            files = ''

        if i == len(lines)-1 and files != '':
            filesout = '\'{'+files+'}\''
            filesout = filesout.replace(',}', '}')
            out.append( ( n, filesout) )

    return out

samples = [
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240521/Muon0/crab_Muon0_Run2024C_hlt_muon_data_20240521/240521_153146/0000/",
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240521/Muon0/crab_Muon0_Run2024C_hlt_muon_data_chaining_20240521/240521_153333/0000/",
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240521/Muon0/crab_Muon0_Run2024C_hlt_muon_data_chaining_v2_20240521/240521_180934/0000/",

    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240528/Muon0/crab_Muon0_Run2024C_hlt_muon_data_20240528/240528_210936/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240528/Muon0/crab_Muon0_Run2024C_hlt_muon_data_CCC_20240528/240528_210106/0000/",

    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1409/20240627/Muon0/crab_Muon0_Run2024F_hlt_muon_data_20240627/240627_160740/0000",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1409/20240627/Muon0/crab_Muon0_Run2024F_hlt_muon_data_CSC_20240627/240627_160753/0000",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1409/20240628/Muon0/crab_Muon0_Run2024F_hlt_muon_data_CSC_RPC_20240628/240628_100651/0000/",
]

dates = [
    '_20240521',
    '20240521',
    '_20240528',
    '20240528',
    '_20240627',
    '20240627',
    '_20240628',
    '20240628',
]

menus = [
    #'_hlt_muon_data_Run2018_',
    '_hlt_muon_mc_',
    '_hlt_muon_data_',
]

analyzers = {
    'Muon0_Run2024B': ('Eff'),
    'Muon1_Run2024B': ('Eff'),
    'Muon0_Run2024C': ('Eff'),
    'Muon1_Run2024C': ('Eff'),
    'Muon0_Run2024D': ('Eff'),
    'Muon1_Run2024D': ('Eff'),
    'Muon0_Run2024Ev1': ('Eff'),
    'Muon1_Run2024Ev1': ('Eff'),
    'Muon0_Run2024Ev2': ('Eff'),
    'Muon1_Run2024Ev2': ('Eff'),
    'Muon0_Run2024F': ('Eff'),
    'Muon1_Run2024F': ('Eff'),

    'Muon0_Run2023B': ('Eff'),
    'Muon1_Run2023B': ('Eff'),
    'Muon0_Run2023Cv1': ('Eff'),
    'Muon1_Run2023Cv1': ('Eff'),
    'Muon0_Run2023Cv2': ('Eff'),
    'Muon1_Run2023Cv2': ('Eff'),
    'Muon0_Run2023Cv3': ('Eff'),
    'Muon1_Run2023Cv3': ('Eff'),
    'Muon0_Run2023Cv4': ('Eff'),
    'Muon1_Run2023Cv4': ('Eff'),
    'Muon0_Run2023Dv1': ('Eff'),
    'Muon1_Run2023Dv1': ('Eff'),
    'Muon0_Run2023Dv2': ('Eff'),
    'Muon1_Run2023Dv2': ('Eff'),

    'SingleMuon_RunUL2018D': ('Eff'),
    'SingleMuon_Run2022B': ('Eff'),
    'SingleMuon_Run2022C': ('Eff'),
    'Muon_Run2022C': ('Eff'),
    'Muon_Run2022Dv1': ('Eff'),
    'Muon_Run2022Dv2': ('Eff'),
    'Muon_Run2022E': ('Eff'),
    'Muon_Run2022F': ('Eff'),
    'Muon_Run2022G': ('Eff'),

    'DYToLL_M50_133X': ('Eff', 'Effgen'),
    'JPsi_133X': ('Eff', 'Effgen'),

    'DYToLL_M50_130X': ('Eff'),

    'Bs_126X': ('Eff', 'Effgen'),
    'JPsi_126X': ('Eff'),
    'DYToLL_M50_126X': ('Eff'),
    'Zprime_126X': ('Eff'),
}

# python3 submit_batch.py
if __name__ == '__main__':
    VER_base = 'vRun3_08'
    tag_prefix = 'crab_'

    doHadd = False
    if len(sys.argv) > 1 and 'hadd' == sys.argv[1]:
        doHadd = True

    doRecover = False
    if len(sys.argv) > 1 and 'recover' == sys.argv[1]:
        doRecover = True

    doTest = False
    if len(sys.argv) > 1 and 'test' == sys.argv[1]:
        doTest = True

    PWD = os.getcwd()

    condor_dir = PWD+'/condor/'
    if not os.path.isdir(condor_dir):
        os.makedirs(condor_dir)

    if doHadd:
        haddlist = open('haddlist.txt', 'w')
    else:
        joblist = open('joblist.txt', 'w')

    for i, path in enumerate(samples):
        info = path.split('/')

        VER = None
        TAG = None
        for x in info:
            if tag_prefix in x:
                info = x.replace(tag_prefix, '')
                tag_split = ''
                for menu in menus:
                    if menu in info:
                        tag_split = menu
                TAG = info.split(tag_split)[0]
                extra = info.split(tag_split)[-1]
                for da in dates:
                    extra = extra.replace(da, '')
                if '_mc' in extra:
                    extra = extra.replace('_mc', '')
                if extra != '':
                    extra = '-' + extra
                VER = VER_base + extra

        if type(analyzers[TAG]) != tuple:
            analyzers[TAG] = (analyzers[TAG], )

        for an in analyzers[TAG]:
            output_dir_base = PWD+'/Outputs_'+VER+'/'+an+'/'
            output_dir = output_dir_base+TAG+'/'
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)

            if doHadd:
                source_files = output_dir+'hist-%(VER)s-%(TAG)s--*-%(an)s.root' % locals()
                target_file = output_dir_base+'hist-%(VER)s-%(TAG)s-%(an)s.root' % locals()
                log_file = output_dir+'hadd-%(VER)s-%(TAG)s.log' % locals()
                if os.path.isfile(target_file):
                    continue

                cmd = "\"hadd -j 15 %(target_file)s %(source_files)s\"\n" % locals()
                haddlist.write(cmd)
                sys.stdout.flush()
            else:
                nfiles = 2

                doDimuon = "false"
                if "DYToLL_M" in TAG or "Zprime" in TAG:
                    doDimuon = "true"

                jobid_files = jobSpritting(path, nfiles)

                MACRO = 'HLT%sAnalyzer' % an # HLTEffAnalyzer
                #MACRO = 'HLTHLTgenAnalyzer' % - efficiency w.r.t gen muons
                #MACRO = 'trackQualAnalyzer'  # trackQualAnalyzer - investigate trk, muon's qualities

                for jobid, files in jobid_files:
                    strjobid = "Job"+str(jobid)
                    arg_str = '{},{},{},{},{},{},{}\n'.format(
                        MACRO,
                        VER,
                        TAG,
                        strjobid,
                        output_dir,
                        doDimuon,
                        files
                    )

                    if doRecover:
                        outfile = output_dir+'hist-%(VER)s-%(TAG)s--%(strjobid)s-%(an)s.root' % locals()
                        if not os.path.isfile(outfile):
                            joblist.write(arg_str)
                            sys.stdout.flush()
                    else:
                        joblist.write(arg_str)
                        sys.stdout.flush()

                    if doTest:
                        break
