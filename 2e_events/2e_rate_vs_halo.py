#!/usr/bin/env python
import pickle
import sys, getopt, os
import array
import re
from ROOT import gROOT, gStyle, TFile, TTree, TChain, TMVA, TCut, TCanvas, gDirectory, TH1, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex ,TH1F
from decimal import *
getcontext().prec=4
import numpy as np
import matplotlib.pyplot as plt
# sphinx_gallery_thumbnail_number = 3
from matplotlib.backends.backend_pdf import PdfPages
import datetime
import pandas as pd
from openpyxl import load_workbook
import path
from scipy.stats import chi2
from decimal import Decimal
import time, sys

sys.path.insert(1, '/mnt/mariano/sho/AnalysisTools/pyroot')
import skipper_utils

gStyle.SetStatW(0.25)
gStyle.SetStatH(0.25)
gROOT.SetBatch(True)
gStyle.SetOptStat(100011)
gStyle.SetOptFit(10001)
gStyle.SetPalette(57)



NROWS = 3072 #number of rows read out, per file
gain = 1.0


HALO=False
ZONES = False

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'e:r:g:o:h:H:f', ['exposure','readout','gain','ohdu','help','halo','fit'])



#some options, actually outdated
for opt, arg in options:
    if opt in ('-e', '--exposure'):
        EXPOSURE = float(arg)
    if opt in ('-H', '--halo'):
        HALO = True 
        n = int(arg)
    if opt in ('-r', '--readout'):
        READOUT = float(arg)
    if opt in ('-g', '--gain'):
        gain = float(arg)
    if opt in ('-o', '--ohdu'):
        OHDU = int(arg)
    if opt in ('-f', '--fit'):
        FIT = True



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
    # return getattr(tree,name).split(b'\0',1)[0]
    return getattr(tree,name).split(b'\0',1)[0]

def decodeRunnum(filename):
    ints = [int(c) for c in re.split(r'(\d+)',filename) if c.isdigit()]
    #print(ints)
    return ints[-1]
#This thing splits our filename using doubles as delimiters so if our filename is "a1b2c3", we do "strs = [str(c) for c in re.split(r'(\d+)',filename)" and then strs=['a','1','b','2','c','3']. This allows us to know what it's in each part of our filename. On the other hand, if we use c.isdigit() now strs=['1','2','3']. In fact that's why sho used int(c), cause he knows that we're dealing with integers

def decodeRunnum2(filename,i):
    ints = [int(c) for c in re.split(r'(\d+)',filename) if c.isdigit()]
    #print(ints)
    return ints[-i]


########################################################################################################################################################################################################
########################################################################################################################################################################################################
#Those were the functions we use. 
########################################################################################################################################################################################################
########################################################################################################################################################################################################    
    

#organize input files.
outfilename = remainder[0]
infiles = remainder[1:]
infiles.sort(key=decodeRunnum) #oh, okey, this sorts in files by the -2 number.

mods=[]
for i in infiles:
    mods.append(decodeRunnum2(i,2))
mods = list(dict.fromkeys(mods))
nmods=len(mods)

#default dimensions - for prototype CCD
ccdPrescan = 7
ccdCol = 724
ccdRow = 1248
nRow = 700
nCol = 450
binRows = 1
binCols = 1
exposure = 0

#default dimensions - for prototype CCD
ccdPrescan = 7
ccdCol = 886
ccdRow = 6144
nRow = 700
nCol = 450


pix = [TChain("calPixTree") for i in range(0,nmods)]
hit = [TChain("hitSumm") for i in range(0,nmods)]

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
    pix[int(getattr(header, "LTANAME"))-1].Add(infiles[i])
    hit[int(getattr(header, "LTANAME"))-1].Add(infiles[i])
    print("Adding "+str(infiles[i])+" | module = "+str(int(getattr(header, "LTANAME"))))

#set geometry
prescanEdge = ccdPrescan + 1 #X of first active pixel
activeEdgeX = ccdPrescan + ccdCol/2 #X of last active pixel
activeEdgeY = ccdRow/2 #Y of last active pixel
edgedistXpix = "min(x*{3}-{0}, {1}-((x+1)*{3}-1))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
edgedistYpix = "min((y-1)*{4}, {2}-(y*{4}))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
edgedistXhit = "min(xMin*{3}-{0}, {1}-((xMax+1)*{3}-1))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
edgedistYhit = "min((yMin-1)*{4}, {2}-(yMax*{4}))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area

#create auxiliar branches
calPixTreeAux = TTree()
hitSummAux = TTree()

#fit function
ipeak = 0
funcformulas = []
for ipeak in range(0,2):
    funcformulas.append("[0]*(TMath::Gaus(x,[1]+{0}*[2],[2]*[3],1)*TMath::Poisson({0},[4]) )".format(ipeak)) #{0} is the ipeak number, [0] not really sure what it is, [1] is mean of gaussian, [2] is gain of calibration, [3] is sigma of gaussian,[4] is mu of poisson signal.
fitfunc = TF1("fitfunc"," + ".join(funcformulas)) #funcformulas is a list so "join" joins its values.



#range
pixmax = 4.0*gain
pixmin = -1.3*gain
pixnbin = 250
pixbinning = "({0},{1},{2})".format(pixnbin,pixmin,pixmax)
binwidth = (pixmax-pixmin)/pixnbin



#canvas creation
c = TCanvas("c","c",1200,900);

outfile = TFile(outfilename+".root","RECREATE")
latex = TLatex()
latex.SetNDC(True) #no idea,  something with the font surely

#regions of interest
htypes = ["unmasked fraction"]
#starters
numRuns = 0

print("Processing...")	
hists = [] #need it to store hists after plot them
c.Clear()
latex.SetTextSize(0.03)

# ##################################################################################################################################################################
# ##################################################################################################################################################################
# ##################################################################################################################################################################

#I will fit the first two hdu's to get rate and eff.

aduZero = 0 #0 of first gaussian
aduPerElectron = gain #gain
noise = 0.15 #noise



hdus=[1,2,3,4]
halos=[5,10]

c.Print(outfilename+".pdf[") #opens pdf and it's multiple pages. Each page has a canvas.
rates=np.zeros([len(pix),len(halos),len(hdus)])
ratesErr=np.zeros([len(pix),len(halos),len(hdus)])
start = time.time()
for I in range(0,len(pix)):
    for i,halo in enumerate(halos):
        latex.DrawLatex(0.1,0.9,"MOD={0},HALO={1}".format(I+1,halo))    
        c.Divide(2,2)
        for j,hdu in enumerate(hdus):
            c.cd(hdu)
            regioncuts = "min(edgedist,distance)>="+str(halo)+" && min(bleedX,bleedY)>50 && !(mask & (1+16+32+128+256+512+1024+2048+4096)) && ohdu=="+str(hdu)
            print(regioncuts)

            hname='h1'
            gPad.SetLogy(1)
            pix[I].SetAlias("edgedist","min({0},{1})".format(edgedistXpix, edgedistYpix))
            pix[I].Draw("ePix>>"+hname+"Module"+str(I+1)+"Halo"+str(halo)+"Quad"+str(hdu)+pixbinning,regioncuts,"colz")
            h = gDirectory.Get(hname+"Module"+str(I+1)+"Halo"+str(halo)+"Quad"+str(hdu))
            hists.append(h)

            # pix[I].SetAlias("edgedist","min({0},{1})".format(edgedistXpix, edgedistYpix))
            # calPixTreeAux.Reset()
            # calPixTreeAux = pix[I].CopyTree(regioncuts)

            hit[I].SetAlias("edgedist","min({0},{1})".format(edgedistXhit, edgedistYhit))
            hitSummAux.Reset()
            hitSummAux = hit[I].CopyTree("ohdu=="+str(hdu)+" && min(edgedist,distance)>="+str(halo)+" && min(bleedX,bleedY)>50 && e==2 && n==1 && (flag & (16+32+128+256+512+1024+2048+4096))==0");
                        
            # rates[I,i,j]=float(hitSummAux.GetEntries())/calPixTreeAux.GetEntries() #mods, halos, hdus
            # ratesErr[I,i,j]=float(hitSummAux.GetEntries())/calPixTreeAux.GetEntries() #mods, halos, hdus
            
            rates[I,i,j]=float(hitSummAux.GetEntries())/h.GetEntries() #mods, halos, hdus
            ratesErr[I,i,j]=float(hitSummAux.GetEntries())/h.GetEntries() #mods, halos, hdus
            
                        
        c.cd()
        c.Print(outfilename+".pdf");
        c.Clear()

c.cd()
c.Print(outfilename+".pdf]");

end = time.time()
print(str(end-start)+" seconds!\n")

with PdfPages(outfilename+'_2e_halos.pdf') as pdf:
    for I in range(0,len(pix)):
        for j,hdu in enumerate(hdus):
            plt.subplot2grid((2,2), (int(hdu-1)/2, int(hdu-1)%2)  , rowspan=1, colspan=1)
            plt.errorbar(halos,rates[I,:,j],yerr=ratesErr[I,:,j],fmt='r.')
            plt.annotate('mod'+str(I+1)+'ohdu '+str(hdu), xy=(-2, -5), xycoords='axes points', size=10, ha='right', va='top', bbox=dict(boxstyle='round', fc='w'))


        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        pdf.savefig()
        plt.close()

