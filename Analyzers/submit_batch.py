#!/usr/bin/env python
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
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1306/20230526/Muon/crab_Muon_Run2022G_hlt_muon_data_20230526/230526_181751/0000/",
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1306/20230526/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/crab_DYToLL_M50_126X_hlt_muon_mc_Run3_20230526/230526_183156/0000/",

    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1326/20231107/Muon0/crab_Muon0_Run2023Dv2_hlt_muon_data_20231107/231107_180843/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1326/20231107/Muon0/crab_Muon0_Run2023Dv2_hlt_muon_data_Doublet_20231107/231107_180901/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1326/20231114/Muon0/crab_Muon0_Run2023Dv2_hlt_muon_data_Doublet_L1L2_20231114/231114_201631/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1326/20231113/Muon0/crab_Muon0_Run2023Dv2_hlt_muon_data_Doublet_L2noROI2_20231113/231113_154427/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1326/20231113/Muon0/crab_Muon0_Run2023Dv2_hlt_muon_data_Doublet_L2noROI_20231113/231113_114124/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1326/20231113/Muon0/crab_Muon0_Run2023Dv2_hlt_muon_data_Doublet_L2noROI_fullPix_20231113/231113_114549/0000/",

]

dates = [
    '_20230526',
    '20230526',
    '_20231107',
    '20231107',
    '_20231113',
    '20231113',
    '_20231114',
    '20231114',
]

menus = [
    #'_hlt_muon_data_Run2018_',
    '_hlt_muon_mc_Run3_',
    '_hlt_muon_data_',
]

analyzers = {
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

    'Bs_126X': ('Eff'),
    'JPsi_126X': ('Eff'),
    'DYToLL_M50_126X': ('Eff'),
    'Zprime_126X': ('Eff'),

    # 'Wprime': ('EffSim', 'Eff'),
    # 'MuGun': ('EffSim'),
}

# python submit_batch.py
if __name__ == '__main__':
    VER_base = 'vRun3_03'
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
                nfiles = 3

                doDimuon = "false"
                if "DYToLL_M" in TAG or "Zprime" in TAG:
                    doDimuon = "true"

                jobid_files = jobSpritting(path, nfiles)

                MACRO = 'HLT%sAnalyzer' % an # HLTEffAnalyzer
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
