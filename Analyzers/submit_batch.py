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
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_default_20220218/220218_072531/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_upgradeIO_20220218/220218_072613/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_upgradeIO_ROI1n5_20220218/220218_072655/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_default_20220218/220218_113112/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_upgradeIO_20220218/220218_113244/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_upgradeIO_ROI1n5_20220218/220218_113412/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_120X_hlt_muon_mc_default_20220218/220218_151523/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_120X_hlt_muon_mc_upgradeIO_20220218/220218_151636/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_120X_hlt_muon_mc_upgradeIO_ROI1n5_20220218/220218_151734/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_120X_hlt_muon_mc_default_20220218/220218_113029/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_120X_hlt_muon_mc_upgradeIO_20220218/220218_113200/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_120X_hlt_muon_mc_upgradeIO_ROI1n5_20220218/220218_113329/0000/",

    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_default_propa_20220218/220218_073004/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_upgradeIO_propa_20220218/220218_073053/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220218/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_upgradeIO_ROI1n5_propa_20220218/220218_073153/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220220/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_default_propa_20220220/220220_151127/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220220/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_upgradeIO_propa_20220220/220220_151250/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220220/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_upgradeIO_ROI1n5_propa_20220220/220220_151413/0000/",

    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220220/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_120X_hlt_muon_mc_default_propa_20220220/220220_150546/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220220/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_120X_hlt_muon_mc_upgradeIO_propa_20220220/220220_150630/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220220/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_120X_hlt_muon_mc_upgradeIO_ROI1n5_propa_20220220/220220_150713/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220220/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_120X_hlt_muon_mc_default_propa_20220220/220220_151046/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220220/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_120X_hlt_muon_mc_upgradeIO_propa_20220220/220220_151208/0000/",
    #"/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre5/20220220/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_120X_hlt_muon_mc_upgradeIO_ROI1n5_propa_20220220/220220_151332/0000/",

]

dates = [
    '_20220218',
    '_20220220',
]

analyzers = {
    'JPsi': ('Eff'),
    'Bs': ('Eff'),
    'DYToLL_M50': ('Eff'),
    'Zprime': ('Eff'),
    #'MuGunPU': ('Eff'),
    #'WJets': ('Eff'),
}


# python submit_batch.py
if __name__ == '__main__':
    VER_base = 'v4'
    tag_prefix = 'crab_'
    tag_split = '_120X_hlt_muon_mc_'

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

                cmd = "\"hadd -j 6 %(target_file)s %(source_files)s\"\n" % locals()
                haddlist.write(cmd)
                sys.stdout.flush()
            else:
                nfiles = 3

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



