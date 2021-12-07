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
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_mc_run2_20211130/211130_170045/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_mc_Run3_20211130/211130_170425/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_mc_Run3_OI_20211130/211130_170656/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_mc_Run3_OI_wp01_20211130/211130_170849/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_mc_Run3_OI_wp01_Iso_20211130/211130_171046/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_mc_Run3_wp01_Iso_20211130/211130_171235/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_mc_run2_20211130/211130_142105/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_mc_Run3_20211130/211130_142218/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_mc_Run3_OI_20211130/211130_142326/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_mc_Run3_OI_wp01_20211130/211130_142434/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_mc_Run3_OI_wp01_Iso_20211130/211130_164609/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_mc_Run3_wp01_Iso_20211130/211130_164739/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_mc_run2_20211130/211130_200916/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_mc_Run3_20211130/211130_201226/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_mc_Run3_OI_20211130/211130_201638/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_mc_Run3_OI_wp01_20211130/211130_202057/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_mc_Run3_OI_wp01_Iso_20211130/211130_202516/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_mc_Run3_wp01_Iso_20211130/211130_202856/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_mc_run2_20211130/211130_201010/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_mc_Run3_20211130/211130_201408/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_mc_Run3_OI_20211130/211130_201806/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_mc_Run3_OI_wp01_20211130/211130_202215/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_mc_Run3_OI_wp01_Iso_20211130/211130_202637/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_mc_Run3_wp01_Iso_20211130/211130_202958/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_mc_run2_20211130/211130_200823/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_mc_Run3_20211130/211130_201120/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_mc_Run3_OI_20211130/211130_201525/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_mc_Run3_OI_wp01_20211130/211130_201931/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_mc_Run3_OI_wp01_Iso_20211130/211130_202341/0000/",
    "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211130/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_mc_Run3_wp01_Iso_20211130/211130_202753/0000/",



    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3v8_mc_wp00_ROI1p5_20211111/211111_084027/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3v8_mc_wp00_ROI2p0_20211111/211111_084225/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3v8_mc_wp01_ROI1p5_20211111/211111_084421/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3v8_mc_wp01_ROI2p0_20211111/211111_084645/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3v8_mc_wp00_ROI1p5_20211111/211111_082110/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3v8_mc_wp00_ROI2p0_20211111/211111_082149/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3v8_mc_wp01_ROI1p5_20211111/211111_082227/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3v8_mc_wp01_ROI2p0_20211111/211111_082304/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3v8_mc_wp00_ROI1p5_20211111/211111_083938/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3v8_mc_wp00_ROI2p0_20211111/211111_084146/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3v8_mc_wp01_ROI1p5_20211111/211111_084342/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3v8_mc_wp01_ROI2p0_20211111/211111_084604/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3v8_mc_wp00_ROI1p5_20211111/211111_084107/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3v8_mc_wp00_ROI2p0_20211111/211111_084303/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3v8_mc_wp01_ROI1p5_20211111/211111_084504/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3v8_mc_wp01_ROI2p0_20211111/211111_084722/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3v8_mc_wp00_ROI1p5_20211111/211111_083606/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3v8_mc_wp00_ROI2p0_20211111/211111_083645/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3v8_mc_wp01_ROI1p5_20211111/211111_083728/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211111/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3v8_mc_wp01_ROI2p0_20211111/211111_083809/0000/",


    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211110/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3Base_mc_20211110/211110_090152/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211110/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3Base_mc_20211110/211110_085556/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211110/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3Base_mc_20211110/211110_090112/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211110/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3Base_mc_20211110/211110_090233/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211110/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3Base_mc_20211110/211110_085717/0000/",

    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211109/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3v8_mc_wp00_20211109/211109_164700/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211109/DoubleMuon_Pt-1To1000-gun/crab_MuGunPU_121X_hlt_muon_Run3v8_mc_wp01_20211109/211109_164856/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211109/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3v8_mc_wp00_20211109/211109_160249/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211109/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_121X_hlt_muon_Run3v8_mc_wp01_20211109/211109_160327/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211110/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3v8_mc_wp00_20211110/211110_023139/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211109/JPsiToMuMu_Pt-0To100-pythia8-gun/crab_Jpsi_121X_hlt_muon_Run3v8_mc_wp01_20211109/211109_164818/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211109/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3v8_mc_wp00_20211109/211109_164740/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211109/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/crab_WJets_121X_hlt_muon_Run3v8_mc_wp01_20211109/211109_164935/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211109/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3v8_mc_wp00_20211109/211109_164430/0000/",
    # "/pnfs/knu.ac.kr/data/cms/store/user/wjun/MuonHLTRun3_cmssw1210/20211109/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/crab_Zprime_M6000_121X_hlt_muon_Run3v8_mc_wp01_20211109/211109_164509/0000/",
]

dates = [
    '_20211109',
    '_20211110',
    '_20211111',
    '_20211130',
]

analyzers = {
    'DYToLL_M50': ('Eff'),
    'WJets': ('Eff'),
    'Zprime_M6000': ('Eff'),
    # 'MuGunPU': ('Eff'),
    # 'Jpsi': ('Eff'),
}


# python submit_batch.py
if __name__ == '__main__':
    VER_base = 'vRun3Review'
    tag_prefix = 'crab_'
    tag_split = '_121X_hlt_muon_mc_'

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



