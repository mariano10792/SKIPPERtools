#!/usr/bin/env python
import sys, getopt, os
import array
import re
from ROOT import gROOT, gStyle, TFile, TTree, TChain, TMVA, TCut, TCanvas, gDirectory, TH1, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex, TList
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import datetime
import pandas as pd
from openpyxl import load_workbook
import path
from scipy.stats import chi2

gStyle.SetStatW(0.25)
gStyle.SetStatH(0.25)
gROOT.SetBatch(True)
gStyle.SetOptStat(100011)
gStyle.SetOptFit(10001)
gStyle.SetPalette(57)



NROWS = 3072 #number of rows read out, per file
gain = 1.0
halo1eOpt=30

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'H:h')

for opt, arg in options:
    if opt == '-H':
        halo1eOpt = int(arg)
    elif opt == '-h':
        print "\nUsage: "+sys.argv[0]+" <output basename> <root files>"
        print "Arguments: "
        print "\t-H: set 1ehalo"
        sys.exit(0)

#some alarm help
if (len(remainder)<2):
    print sys.argv[0]+' <output basename> <root files>'
    sys.exit()

	
########################################################################################################################################################################################################
########################################################################################################################################################################################################
#from "options, remainder" until here this code decodifies what instructions you've gave to it. The if up here helps you with the remainders: output basename (name we want for the output) and the rootfile name (input file to analyze).
########################################################################################################################################################################################################
########################################################################################################################################################################################################

#this guy gets the value of a header.
def getHeaderValue(tree,name):
    tree.GetEntry(0)
    stringdata = tree.GetBranch(name).GetLeaf("string").GetValuePointer()
    b = bytearray()
    for iword in range(0,len(stringdata)):
        word = stringdata[iword]
        #print(word)
        for ibyte in range(0,8):
            b.append(word & 0xFF)
            word >>= 8
    return float(b)

def getHeaderString(tree,name):
    tree.GetEntry(0)
    return getattr(tree,name).split(b'\0',1)[0]

def decodeRunnum(filename,i=1):
    ints = [int(c) for c in re.split(r'(\d+)',filename) if c.isdigit()]
    #print(ints)
    return ints[-i]


###########################################################################################################################################################
################################################ Those were the functions we use. ################################################
###########################################################################################################################################################

#default dimensions - for prototype CCD
ccdPrescan = 7
ccdCol = 886
ccdRow = 6144
nRow = 700
nCol = 450

#organize input files.
outfilename = remainder[0]
infiles = remainder[1:]
infiles.sort(key=decodeRunnum) #oh, okey, this sorts in files by the -2 number.

files=[[],[]]
# nmods=len(files)
nmods=array.array('i', [len(files)])
modules=array.array('i', range(len(files)+1)[1:])

for i in range(0,len(infiles)):
    thefile = TFile(infiles[i])
    header = thefile.Get("headerTree_1")
    header.GetEntry(0)
    if i==0:
        cddPrescan = int(getattr(header, "CCDNPRES"))
        ccdCol = int(getattr(header, "CCDNCOL"))
        ccdRow = int(getattr(header, "CCDNROW"))
        nRow = int(getattr(header, "NROW"))
        nCol = int(getattr(header, "NCOL"))
        try:
            binRows = int(getattr(header, "NBINROW"))
        except:
            binRows = 1
        try:
            binCols = int(getattr(header, "NBINCOL"))
        except:
            binCols = 1
    files[int(getattr(header, "LTANAME"))-1].append(infiles[i])
    print("Adding "+str(infiles[i])+" to module = "+str(int(getattr(header, "LTANAME"))))


#starters
hitSumm = TChain("hitSumm")

#enable root branch.
hitSumm.SetBranchStatus("*",1)		

#create auxiliar branches
hitSummAux1 = TTree()
hitSummAux2 = TTree()

#set geometry
prescanEdge = ccdPrescan + 1 #X of first active pixel
activeEdgeX = ccdPrescan + ccdCol/2 #X of last active pixel
activeEdgeY = ccdRow/2 #Y of last active pixel

#set cuts
hdus=['1','2','3','4'] #set hdus to analyze
hdus=['ohdu=='+hdu for hdu in hdus]
hduscut='('+' || '.join(hdus)+')'
halo=20
halo1e=halo1eOpt
edge=halo
masks='(4+16+128+256+512+1024+4096)'


"ohdu=="," + ".join(hdus)



#### calculate expected 2e- events. As a rough estimation I'll take all quads from each module to calculate the expected 2e- events per pix.
#range
pixmax = 4.0*gain; pixmin = -1.3*gain; pixnbin = 250; pixbinning = "({0},{1},{2})".format(pixnbin,pixmin,pixmax); binwidth = (pixmax-pixmin)/pixnbin; aduZero = 0; aduPerElectron = 1; noise = 0.15;
#fit function
ipeak = 0; funcformulas = [];
for ipeak in range(0,2):
    funcformulas.append("[0]*(TMath::Gaus(x,[1]+{0}*[2],[2]*[3],1)*TMath::Poisson({0},[4]) )".format(ipeak))
fitfunc = TF1("fitfunc"," + ".join(funcformulas)) #funcformulas is a list so "join" joins its values.
#canvas and rootfile creation
c = TCanvas("c","c",1200,900); c.Clear();
outfile = TFile(outfilename+".root","RECREATE");
hists = [] #need it to store hists after plot them


edgedistX = "min(x*{3}-{0}, {1}-((x+1)*{3}-1))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
edgedistY = "min((y-1)*{4}, {2}-(y*{4}))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area


exp2e=[[],[]]
exp2eErr=[[],[]]
meas2e=[[],[]]

for halo1e in [5,10,20,30]:
    calPixTree = [TChain("calPixTree") for i in range(0,len(files))]
    print(halo1e)
    edgedistX = "min(x*{3}-{0}, {1}-((x+1)*{3}-1))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
    edgedistY = "min((y-1)*{4}, {2}-(y*{4}))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
    data = [array.array('f',[0,0]) for i in range(0,7)] #entries, 1e- rate, 1e- rateErr, expected 1pix2e- rate, expected 1pix2e- rateErr, expected 1pix2e- events, expected 1pix2e- events uncertainty
    data[0]=array.array('i',[0,0]) #first one is integer!
    for I in range(0,len(files)):
        for iFile in files[I]:
            calPixTree[I].Add(iFile)
        regioncuts="min(distance,edgedist)>="+str(halo)+" && distance1e>"+str(halo1e)+" && (mask & "+str(masks)+")==0 &&"+hduscut
        hname='h'
        calPixTree[I].SetAlias("edgedist","min({0},{1})".format(edgedistX, edgedistY))
        calPixTree[I].Draw("ePix>>"+hname+"Module"+str(I+1)+pixbinning,regioncuts)
        gPad.SetLogy(1)
        h = gDirectory.Get(hname+"Module"+str(I+1))
        hists.append(h)
        fitfunc.SetParameters(h.GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.00001)
        s = h.Fit(fitfunc,"QSL","",aduZero-1*aduPerElectron,aduZero+1.5)
        data[0][I]=int(h.GetEntries())
        data[1][I]=float(s.Parameter(4))
        data[2][I]=float(s.Error(4))
        data[3][I]=float((s.Parameter(4)**2)/2)
        data[4][I]=float(s.Error(4)*s.Parameter(4))
        data[5][I]=float(h.GetEntries()*(s.Parameter(4)**2)/2)
        data[6][I]=float(h.GetEntries()*s.Error(4)*s.Parameter(4))
        exp2e[I].append(data[5][I])
        exp2eErr[I].append(data[6][I])

    #### At this point I will count the 2e- events ####
    edgedistX = "min(xMin*{3}-{0}, {1}-((xMax+1)*{3}-1))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
    edgedistY = "min((yMin-1)*{4}, {2}-(yMax*{4}))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
    hitSumm.SetAlias("edgedist","min({0},{1})".format(edgedistX, edgedistY))
    for I in range(0,len(files)):
        hitSumm.Reset()
        for iFile in files[I]:
            hitSumm.Add(iFile)
        if I==0: #first module
            hitSummAux1 = hitSumm.CopyTree("min(distance,edgedist)>="+str(halo)+" && distance1e>"+str(halo1e)+" && e==2 && n==2 && (flag & "+str(masks)+")==0 &&"+hduscut);
            meas2e[I].append(hitSummAux1.GetEntries())
        else: #second module
            hitSummAux2 = hitSumm.CopyTree("min(distance,edgedist)>="+str(halo)+" && distance1e>"+str(halo1e)+" && e==2 && n==2 && (flag & "+str(masks)+")==0 &&"+hduscut);
            meas2e[I].append(hitSummAux2.GetEntries())

        

print(hitSummAux1.GetEntries())
print(hitSummAux2.GetEntries())


lists = TList()
lists.Add(hitSummAux1)
lists.Add(hitSummAux2)

del hitSummAux1
del hitSummAux2
del hitSumm

outfile.Write()
outfile.Close()

#python -i /home/mariano/tools/SENSEI_Analysis/2e_events/2pix2e_events_extractor.py salida hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run1**7200*t