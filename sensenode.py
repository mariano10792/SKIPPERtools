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

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'e:r:g:o:h:n:f', ['exposure','readout','gain','ohdu','help','nocheck','fit'])

#some options, actually outdated
for opt, arg in options:
    if opt in ('-h', '--help'):
        print "\nUsage: "+sys.argv[0]+" <output basename> <root files>"
        print "\n"
        sys.exit(0)



def decodeRunnum(filename,i=1):
    ints = [int(c) for c in re.split(r'(\d+)',filename) if c.isdigit()]
    #print(ints)
    return ints[-i]


#gains
gain=225 #hdu2

#organize input files.
outfilename = remainder[0]
infiles = remainder[1:]
infiles.sort(key=decodeRunnum) #oh, okey, this sorts in files by the -2 number.

thefile = TFile(infiles[0])
header = thefile.Get("headerTree_1")
header.GetEntry(0)
nrow=int(getattr(header, "NROW"))
ncol=int(getattr(header, "NCOL"))
vdd=int(float(getattr(header, "VDD")))

if vdd==-17:
    gain=gain*0.8

splPixTree = TChain("splPixTree")

#canvas creation
c = TCanvas("c","c",1200,900);
outfile = TFile(outfilename+".root","RECREATE")
hists = [] #need it to store hists after plot them
c.Clear()

c.Print(outfilename+".pdf[") #opens pdf and it's multiple pages. Each page has a canvas.

stds=[]
stdsErr=[]
Events=0
Denom=0
for hdu in [2]: #[2,3]
    # for i in infiles:
    for i,file in enumerate(infiles):
        splPixTree.Reset()
        for x in range(ncol):
            for y in range(nrow):
                if x+y==0:
                    continue
                splPixTree.Add(file)
                hname='h1'+str(i)+''+str(x)+''+str(y)
                splPixTree.Draw("pix/"+str(gain)+":Iteration$>>"+str(hname)+"(50,0,5000000,,,)","ohdu=="+str(hdu)+" && x=="+str(x)+" && y=="+str(y),"prof")
                hists.append(gDirectory.Get(hname))
                bins=[hists[-1].GetBinContent(i) for i in range(1,hists[-1].GetXaxis().GetNbins()+1)]
                
                if round(bins[-1]-bins[0])>0 and round(bins[-1]-bins[0])<10:
                    Events+=int(round(bins[-1]-bins[0]))
                elif round(bins[-1]-bins[0])>10:
                    continue
                Denom+=1

                print(hists[-1].GetStdDev(2))
                stds.append(hists[-1].GetStdDev(2))
                stdsErr.append(hists[-1].GetStdDevError(2))
                print("Done")
                c.cd()
                c.Print(outfilename+".pdf");
                c.Clear()
c.cd()
c.Print(outfilename+".pdf]");
