#!/usr/bin/env python
import sys,os

MACRO=sys.argv[1]
VER=sys.argv[2]
TAG=sys.argv[3]
JobID=sys.argv[4]
outputDir=sys.argv[5]
doDimuon=sys.argv[6]
DATASET=sys.argv[7]

cmd = 'root -l -b -q \'%(MACRO)s.C("%(VER)s","%(TAG)s",%(DATASET)s,"%(JobID)s","%(outputDir)s",%(doDimuon)s)\'' % locals()

print(cmd)
sys.stdout.flush()

os.system(cmd)
sys.stdout.flush()