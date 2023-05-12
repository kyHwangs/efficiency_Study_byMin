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
    #str_dcap  = ""
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
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230416/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/crab_Bs_126X_hlt_muon_mc_Run3_20230416/230416_202019/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230410/JPsiTo2Mu_Pt-0To100_pythia8-gun/crab_JPsi_126X_hlt_muon_mc_Run3_20230410/230410_213927/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230410/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/crab_DYToLL_M50_126X_hlt_muon_mc_Run3_20230410/230410_213653/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230410/ZprimeToMuMu_M-6000_TuneCP5_13p6TeV_pythia8/crab_Zprime_126X_hlt_muon_mc_Run3_20230410/230410_213807/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230412/Muon/crab_Muon_Run2022G_hlt_muon_data_Run2022G_Full_20230412/230412_054909/0000/",

    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230413/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/crab_Bs_126X_hlt_muon_mc_Run3_wp02_20230413/230413_163853/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230414/JPsiTo2Mu_Pt-0To100_pythia8-gun/crab_JPsi_126X_hlt_muon_mc_Run3_wp02_20230414/230414_091718/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230415/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/crab_DYToLL_M50_126X_hlt_muon_mc_Run3_wp02_20230415/230415_125456/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230415/ZprimeToMuMu_M-6000_TuneCP5_13p6TeV_pythia8/crab_Zprime_126X_hlt_muon_mc_Run3_wp02_20230415/230415_130128/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230413/Muon/crab_Muon_Run2022G_hlt_muon_data_Run2022G_Full_wp02_20230413/230413_163307/0000/",

    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230413/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/crab_Bs_126X_hlt_muon_mc_Run3_wp04_20230413/230413_164024/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230414/JPsiTo2Mu_Pt-0To100_pythia8-gun/crab_JPsi_126X_hlt_muon_mc_Run3_wp04_20230414/230414_174842/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230415/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/crab_DYToLL_M50_126X_hlt_muon_mc_Run3_wp04_20230415/230415_125533/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230415/ZprimeToMuMu_M-6000_TuneCP5_13p6TeV_pythia8/crab_Zprime_126X_hlt_muon_mc_Run3_wp04_20230415/230415_141434/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230413/Muon/crab_Muon_Run2022G_hlt_muon_data_Run2022G_Full_wp04_20230413/230413_163534/0000/",

    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230414/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/crab_Bs_126X_hlt_muon_mc_Run3_wp10_20230414/230414_083333/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230415/JPsiTo2Mu_Pt-0To100_pythia8-gun/crab_JPsi_126X_hlt_muon_mc_Run3_wp10_20230415/230415_124301/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230415/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/crab_DYToLL_M50_126X_hlt_muon_mc_Run3_wp10_20230415/230415_125610/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230415/ZprimeToMuMu_M-6000_TuneCP5_13p6TeV_pythia8/crab_Zprime_126X_hlt_muon_mc_Run3_wp10_20230415/230415_141601/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1303/20230413/Muon/crab_Muon_Run2022G_hlt_muon_data_Run2022G_Full_wp10_20230413/230413_163720/0000/",
]

dates = [
    '_20230410',
    '20230410',
    '_20230412',
    '20230412',
    '_20230413',
    '20230413',
    '_20230414',
    '20230414',
    '_20230415',
    '20230415',
    '_20230416',
    '20230416',
]

menus = [
    '_hlt_muon_data_Run2018_',
    '_hlt_muon_data_Run2022_',
    '_hlt_muon_mc_Run3_',
    '_hlt_muon_data_Run2022G_',
]

analyzers = {
    'SingleMuon_RunUL2018D': ('Eff'),
    'SingleMuon_Run2022B': ('Eff'),
    'SingleMuon_Run2022C': ('Eff'),
    'Muon_Run2022C': ('Eff'),
    'Muon_Run2022Dv1': ('Eff'),
    'Muon_Run2022Dv2': ('Eff'),
    'Muon_Run2022E': ('Eff'),
    'Muon_Run2022F': ('Eff'),
    'Muon_Run2022G': ('Eff'),

    'Bs_126X': ('Eff'),
    'JPsi_126X': ('Eff'),
    'DYToLL_M50_126X': ('Eff'),
    'Zprime_126X': ('Eff'),

    # 'Wprime': ('EffSim', 'Eff'),
    # 'MuGun': ('EffSim'),
}


# python submit_batch.py
if __name__ == '__main__':
    VER_base = 'vRun3_07'
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
                nfiles = 1

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



