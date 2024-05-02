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
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240308/Muon0/crab_Muon0_Run2023Dv1_hlt_muon_data_20240308/240308_110624/0000/",
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240308/Muon0/crab_Muon0_Run2023Dv1_hlt_muon_data_chaining_20240308/240308_110644/0000/",
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240311/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/crab_DYToLL_M50_133X_hlt_muon_mc_20240311/240311_210512/0000/",
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240311/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/crab_DYToLL_M50_133X_hlt_muon_mc_chaining_20240311/240311_210542/0000/",
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240311/JPsiToMuMu_PT-0to100_pythia8-gun/crab_JPsi_133X_hlt_muon_mc_20240311/240311_210526/0000/",
    #"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240311/JPsiToMuMu_PT-0to100_pythia8-gun/crab_JPsi_133X_hlt_muon_mc_chaining_20240311/240311_210555/0000/",

    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1405/20240501/Muon0/crab_Muon0_Run2024B_hlt_muon_data_20240501/240501_214423/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1405/20240501/Muon0/crab_Muon0_Run2024C_hlt_muon_data_20240501/240501_214450/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1405/20240501/Muon1/crab_Muon1_Run2024B_hlt_muon_data_20240501/240501_214437/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1405/20240501/Muon1/crab_Muon1_Run2024C_hlt_muon_data_20240501/240501_214502/0000/",

    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon0/crab_Muon0_Run2024B_hlt_muon_data_20240430/240429_222438/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon0/crab_Muon0_Run2024B_hlt_muon_data_BDT_wp04_20240430/240429_222633/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon0/crab_Muon0_Run2024B_hlt_muon_data_BDT_wp10_20240430/240429_222830/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon0/crab_Muon0_Run2024B_hlt_muon_data_BDT_wp20_20240430/240429_223026/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon0/crab_Muon0_Run2024C_hlt_muon_data_20240430/240429_222534/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon0/crab_Muon0_Run2024C_hlt_muon_data_BDT_wp04_20240430/240429_222730/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon0/crab_Muon0_Run2024C_hlt_muon_data_BDT_wp10_20240430/240429_222927/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon0/crab_Muon0_Run2024C_hlt_muon_data_BDT_wp20_20240430/240429_223124/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon1/crab_Muon1_Run2024B_hlt_muon_data_20240430/240429_222506/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon1/crab_Muon1_Run2024B_hlt_muon_data_BDT_wp04_20240430/240429_222702/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon1/crab_Muon1_Run2024B_hlt_muon_data_BDT_wp10_20240430/240429_222858/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon1/crab_Muon1_Run2024B_hlt_muon_data_BDT_wp20_20240430/240429_223055/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon1/crab_Muon1_Run2024C_hlt_muon_data_20240430/240429_222603/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon1/crab_Muon1_Run2024C_hlt_muon_data_BDT_wp04_20240430/240429_222759/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon1/crab_Muon1_Run2024C_hlt_muon_data_BDT_wp10_20240430/240429_222955/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1406/20240430/Muon1/crab_Muon1_Run2024C_hlt_muon_data_BDT_wp20_20240430/240429_223153/0000/",
]

dates = [
    '_20240308',
    '20240308',
    '_20240311',
    '20240311',
    '_20240319',
    '20240319',
    '_20240329',
    '20240329',
    '_20240331',
    '20240331',
    '_20240430',
    '20240430',
    '_20240501',
    '20240501',
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
    VER_base = 'vRun3_06'
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
