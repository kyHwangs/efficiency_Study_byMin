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
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240308/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/crab_DYToLL_M50_133X_hlt_muon_mc_20240308/240308_101938/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240308/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/crab_DYToLL_M50_133X_hlt_muon_mc_chaining_20240308/240308_102234/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240308/JPsiToMuMu_PT-0to100_pythia8-gun/crab_JPsi_133X_hlt_muon_mc_20240308/240308_213445/0000/",
    "/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTRun3_cmssw1401/20240308/JPsiToMuMu_PT-0to100_pythia8-gun/crab_JPsi_133X_hlt_muon_mc_chaining_20240308/240308_213516/0000/",
]

dates = [
    '_20240308',
    '20240308',
]

menus = [
    '_hlt_muon_mc_',
]

analyzers = {
    'DYToLL_M50_133X': ('Eff'),
    'JPsi_133X': ('Eff'),

    'JPsi_126X': ('Eff'),
    'Bs_126X': ('Eff'),
    'DYToLL_M50_126X': ('Eff'),
    'Zprime_126X': ('Eff'),
    #'MuGunPU': ('Eff'),
    #'WJets': ('Eff'),
}

# python3 submit_batch.py
if __name__ == '__main__':
    VER_base = 'vRun3_05'
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
