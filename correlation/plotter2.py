#!/usr/bin/env python
import sys, getopt
import array
import datetime
import csv
import math
import numpy as np
import time
import re
import matplotlib.pyplot as plt

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'y:h')

for opt, arg in options:
    if opt == '-y':
        maxY = int(arg)
    elif opt == '-h':
        print "\nUsage: "+sys.argv[0]+" <output basename> <root files>"
        print "Arguments: "
        print "\t-y: max-Y cut to apply to data"
        print "\n"
        sys.exit(0)

from ROOT import gROOT
gROOT.SetBatch(True)
from ROOT import gStyle, TFile, TTree, TChain, TCanvas, gDirectory, TH1, TH1F, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex, TEfficiency, TMath, TVectorD

def decode(filename,i=17):
    ints = [int(c) for c in re.split(r'(\d+)',filename) if c.isdigit()]
    #print(ints)
    return ints[-i]

gStyle.SetOptStat(110011)
gStyle.SetOptFit(1)

infiles = remainder

#extracting data from trees
halos=[]
for i in infiles:
    halos.append(decode(i))
#    print("Adding "+str(i))

halos=list(dict.fromkeys(halos))
halos.sort()
infiles.sort(key=decode)

data=[TChain("pvalues") for i in halos]
counter=0
for i in infiles:
	if halos[counter]!=decode(i):
		counter+=1
	data[counter].Add(i)
	print("Adding "+str(i))



pvalues=[[],[],[],[]]
pvaluesErr=[[],[],[],[]]
fig=plt.figure()
for hdu in [1,2,3,4]:
    for i in range(0,len(halos)):
        # pvalues[hdu-1].append(data[i].GetEntries("pvalues<0.1 && ohdu=="+str(hdu))/float(data[i].GetEntries("ohdu=="+str(hdu))))
        # pvaluesErr[hdu-1].append((np.sqrt(data[i].GetEntries("ohdu=="+str(hdu))*(3.0/10))/data[i].GetEntries("ohdu=="+str(hdu))))
        pvalues[hdu-1].append(data[i].GetEntries("pvalues2<0.1 && ohdu=="+str(hdu))/float(data[i].GetEntries("ohdu=="+str(hdu))))
        pvaluesErr[hdu-1].append((np.sqrt(data[i].GetEntries("ohdu=="+str(hdu))*(np.sqrt(0.1*(1-0.1))/10))/data[i].GetEntries("ohdu=="+str(hdu))))
    plt.errorbar(halos,pvalues[hdu-1],yerr=pvaluesErr[hdu-1],label="ohdu="+str(hdu))
#	plt.plot(halos,pvalues[hdu-1],'.-',label="ohdu="+str(hdu))

plt.title("module"+str(decode(infiles[0],2)))
plt.legend(loc="lower left")
plt.show()
