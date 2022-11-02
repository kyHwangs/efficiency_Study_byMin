#!/usr/bin/env python
import sys
import os
import time
import glob
import subprocess
from shutil import copyfile
import gc

def jobSpritting( path, nfiles, prefix = "" ):
    str_dcap  = "dcap://cluster142.knu.ac.kr/"
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
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/SingleMuon/crab_SingleMuon_Run2022B_hlt_muon_data_full_Run2022_20220927/220927_155858/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/SingleMuon/crab_SingleMuon_Run2022B_hlt_muon_data_full_Run2022_NoWP_20220927/220927_161133/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/SingleMuon/crab_SingleMuon_Run2022C_hlt_muon_data_full_Run2022_20220927/220927_160108/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/SingleMuon/crab_SingleMuon_Run2022C_hlt_muon_data_full_Run2022_NoWP_20220927/220927_161233/0000/",

    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/Muon/crab_Muon_Run2022C_hlt_muon_data_full_Run2022_20220927/220927_160335/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/Muon/crab_Muon_Run2022C_hlt_muon_data_full_Run2022_NoWP_20220927/220927_161316/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/Muon/crab_Muon_Run2022Dv1_hlt_muon_data_full_Run2022_20220927/220927_160648/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/Muon/crab_Muon_Run2022Dv1_hlt_muon_data_full_Run2022_NoWP_20220927/220927_161425/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/Muon/crab_Muon_Run2022Dv2_hlt_muon_data_full_Run2022_20220927/220927_160845/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220927/Muon/crab_Muon_Run2022Dv2_hlt_muon_data_full_Run2022_NoWP_20220927/220927_161528/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20221017/Muon/crab_Muon_Run2022E_hlt_muon_data_full_Run2022_20221017/221017_144420/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20221017/Muon/crab_Muon_Run2022E_hlt_muon_data_full_Run2022_NoWP_20221017/221017_144544/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20221019/Muon/crab_Muon_Run2022F_hlt_muon_data_full_Run2022_20221019/221019_141625/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20221019/Muon/crab_Muon_Run2022F_hlt_muon_data_full_Run2022_NoWP_20221019/221019_141643/0000/",

    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220930/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_120X_hlt_muon_mc_full_Run3_20220930/220930_173539/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20220930/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_120X_hlt_muon_mc_full_Run3_NoWP_20220930/220930_173555/0000/",

    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20221017/SingleMuon/crab_SingleMuon_Run2018D_hlt_muon_data_full_Run2018_20221017/221017_160415/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1248/20221017/SingleMuon/crab_SingleMuon_RunUL2018D_hlt_muon_data_full_Run2018_20221017/221017_160439/0000/",
]

dates = [
    '_20220927',
    '20220927',
    '_20220930',
    '20220930',
    '_20221017',
    '20221017',
    '_20221019',
    '20221019',
]

analyzers = {
    'SingleMuon_Run2022B': ('Eff'),
    'SingleMuon_Run2022C': ('Eff'),
    'Muon_Run2022C': ('Eff'),
    'Muon_Run2022Dv1': ('Eff'),
    'Muon_Run2022Dv2': ('Eff'),
    'Muon_Run2022E': ('Eff'),
    'Muon_Run2022F': ('Eff'),
    'DYToLL_M50_120X': ('Eff'),

    'SingleMuon_Run2018D': ('Eff'),
    'SingleMuon_RunUL2018D': ('Eff'),
    # 'JPsi': ('EffSim', 'Eff'),
    # 'Bs': ('EffSim', 'Eff'),
    # 'Zprime': ('EffSim', 'Eff'),

    # 'Wprime': ('EffSim', 'Eff'),
    # 'MuGun': ('EffSim'),
}


# python submit_batch.py
if __name__ == '__main__':
    VER_base = 'vRun3_02'
    tag_prefix = 'crab_'
    #tag_split = '_hlt_muon_data_full_Run2022_'
    tag_split = '_hlt_muon_data_full_Run2018_'
    #tag_split = '_hlt_muon_mc_full_Run3_'

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

                cmd = "\"hadd -j 6 %(target_file)s %(source_files)s\"\n" % locals()
                haddlist.write(cmd)
                sys.stdout.flush()
            else:
                nfiles = 3

                doDimuon = "false"
                if "DYToLL_M" in TAG or "Zprime" in TAG:
                    doDimuon = "true"

                jobid_files = jobSpritting(path, nfiles)

                MACRO = 'HLT%sAnalyzer' % an

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



