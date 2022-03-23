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
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_default_20220321/220321_161053/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_ROIL1_2_20220321/220321_161245/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_ROIL1_2_ROIL2_3_20220321/220321_161340/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_ROIL1_2_ROIL2_5_20220321/220321_161312/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_ROIL1_3_20220321/220321_161217/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_ROIL1_4_20220321/220321_161150/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/crab_Bs_120X_hlt_muon_mc_ROIL1_5_20220321/220321_161121/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_default_20220321/220321_161040/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_ROIL1_2_20220321/220321_161232/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_ROIL1_2_ROIL2_3_20220321/220321_161327/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_ROIL1_2_ROIL2_5_20220321/220321_161259/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_ROIL1_3_20220321/220321_161204/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_ROIL1_4_20220321/220321_161137/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1230pre6/20220321/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_JPsi_120X_hlt_muon_mc_ROIL1_5_20220321/220321_161108/0000/",

]

dates = [
    '_20220321',
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
    VER_base = 'v6'
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
                nfiles = 5

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



