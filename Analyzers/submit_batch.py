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
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3_mc_wp00_20211128/211128_160946/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3_mc_wp00_OI_20211128/211128_160454/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3_mc_wp01_20211128/211128_161215/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3_mc_wp01_OI_20211128/211128_160721/0000/",

    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3_mc_wp00_20211128/211128_153549/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3_mc_wp00_OI_20211128/211128_153155/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3_mc_wp01_20211128/211128_153732/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3_mc_wp01_OI_20211128/211128_153430/0000/",

    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3_mc_wp00_20211128/211128_160835/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3_mc_wp00_OI_20211128/211128_160342/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3_mc_wp01_20211128/211128_161100/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3_mc_wp01_OI_20211128/211128_160609/0000/",

    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3_mc_wp00_20211128/211128_/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3_mc_wp00_OI_20211128/211128_/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3_mc_wp01_20211128/211128_/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3_mc_wp01_OI_20211128/211128_/0000/",

    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3_mc_wp00_20211128/211128_162026/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3_mc_wp00_OI_20211128/211128_161758/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3_mc_wp01_20211128/211128_162141/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211128/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3_mc_wp01_OI_20211128/211128_161911/0000/",
]

dates = [
    '_20211128',
]

analyzers = {
    'DYToLL_M50': ('Eff'),
    #'WJets': ('Eff'),
    'Zprime_M6000': ('Eff'),
    'MuGunPU': ('Eff'),
    'Jpsi': ('Eff'),
}


# python submit_batch.py > joblist.txt
if __name__ == '__main__':
    VER_base = 'v1'
    tag_prefix = 'crab_'
    tag_split = '_121X_hlt_muon_Run3_mc_'

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
                VER = VER_base + '-' + extra

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

                cmd = "hadd -j 6 %(target_file)s %(source_files)s >%(log_file)s" % locals()

                print("")
                print(cmd)
                sys.stdout.flush()

                os.system(cmd)
                sys.stdout.flush()
            else:
                nfiles = 10
                if 'Zprime_M6000' in TAG:
                    nfiles = 10

                doDimuon = "false"
                if "DYToLL_M" in TAG or "Zprime_M6000" in TAG:
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



