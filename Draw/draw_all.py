import os, sys
from time import sleep

VER = 'v6'

samples = [
    ('Jpsi gun', 'JPsi'),
    ('Bs', 'Bs'),
    ('DY', 'DYToLL_M50'),
    ('Zprime 6 TeV', 'Zprime'),
    #('WJets', 'WJets'),
    #('Muon gun', 'MuGunPU'),
]

l1types = [
    #('L1SQ22', 'L1 qual > 11, p_{T}^{L1} > 22 GeV'),
    #('L1DQ22', 'L1 qual > 7, p_{T}^{L1} > 22 GeV'),
    #('L1SQ8', 'L1 qual > 11, p_{T}^{L1} > 8 GeV'),
    #('L1DQ8', 'L1 qual > 7, p_{T}^{L1} > 8 GeV'),
    ('L1SQ0', 'L1 qual > 11, no p_{T}^{L1} cut'),
    ('L1DQ0', 'L1 qual > 7, no p_{T}^{L1} cut'),
]

macros = [
    #'drawEffPixelTrackwrtL1.C',
    #'drawEffCompIOFromL1wrtL1.C',
    #'drawEffOIIOFromL2wrtL1.C',
    #'drawEffCompL3TrackwrtL1.C',
    #'drawEffCompL3MuonIDwrtL1.C',
    #'drawEffCompL3MuonwrtL1.C',
    #'drawEffOIwrtL1.C',
    #'drawEffL2wrtL1.C',
    #'drawEffL3FromL2wrtL1.C',

    #'drawPurCompIOFromL1wrtL1.C',
    #'drawPurCompL3MuonwrtL1.C',
    'drawPurCompL3MuonIDwrtL1.C',
]

l3types = [
    #("hltOI", "OI"),
    #("hltL3FromL2Merged", "OI + IO(L2)"),
    #("hltL3Merged", "OI + IO(L2) + IO(L1)"),
    #("hltIterL3MuonNoID", "L3 muon before ID"),
    #("hltIterL3Muon", "L3 muon after ID")
]

for SAMPLE, TAG in samples:
    for MACRO in macros:
        if 'Pur' in MACRO:
            cmd = 'root -l -b -q \'{}("{}", "{}", "{}")\''.format(
                MACRO,
                VER,
                SAMPLE,
                TAG
            )
            print ''
            print cmd
            sys.stdout.flush()
            os.system(cmd)
            sys.stdout.flush()
            sleep(0.1)

        else:
            for L1SQTAG, L1SQTAGSTR in l1types:
                cmd = 'root -l -b -q \'{}("{}", "{}", "{}", "{}", "{}")\''.format(
                    MACRO,
                    VER,
                    SAMPLE,
                    TAG,
                    L1SQTAG,
                    L1SQTAGSTR
                )
                print ''
                print cmd
                sys.stdout.flush()
                os.system(cmd)
                sys.stdout.flush()
                sleep(0.1)

    if TAG in ["Jpsi", "WJets"]:
        continue
    for L3type, L3typestr in l3types:
        cmd = 'root -l -b -q \'drawRes.C("{}", "{}", "{}", "{}", "{}")\' >res.log'.format(
            VER,
            SAMPLE,
            TAG,
            L3type,
            L3typestr
        )
        #print ''
        #print cmd
        #sys.stdout.flush()
        #os.system(cmd)
        #sys.stdout.flush()
        #sleep(0.1)



