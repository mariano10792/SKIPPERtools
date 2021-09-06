#!/usr/bin/env python
import sys, getopt
import array
import datetime
import json
import numpy as np
import xml.etree.ElementTree as ET

#these parameters will need to be changed as appropriate for the data you're looking at - they are not auto-detected
#HDU_LIST = [1,2,3,4]
#HDU_LIST = [3]
#EXPOSURE_DAYS = 0.0
#HDU_LIST = [1,2,3]
#EXPOSURE_DAYS = 0.5
#EXPOSURE_DAYS = 1.0/24
HDU_LIST = [1]

outfilename=""

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'o:h')

for opt, arg in options:
    if opt == '-o':
        outfilename = arg
    elif opt == '-h':
        print "\nUsage: "+sys.argv[0]+" <output basename> <root files>"
        print "Arguments: "
        print "\t-o: basename for output file"
        print "\n"
        sys.exit(0)

#import ROOT now, since otherwise PyROOT eats the command-line options
#from ROOT import gROOT, gStyle, TFile, TTree, TChain, TCanvas, gDirectory, TH1, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex, TEfficiency, TMath
from ROOT import gROOT
gROOT.SetBatch(True)
from ROOT import gStyle
from ROOT import TFile
from ROOT import TTree
from ROOT import TChain
from ROOT import TCanvas
from ROOT import gDirectory
from ROOT import TH1, TH1F
from ROOT import gPad
from ROOT import THStack
from ROOT import TLegend
from ROOT import TLatex
from ROOT import TEfficiency
from ROOT import TFeldmanCousins
from ROOT import TF1



gStyle.SetOptStat(110011)
gStyle.SetOptFit(1)
#gStyle.SetPalette(57)

if (len(remainder)<1):
    print sys.argv[0]+' <root files>'
    sys.exit()

infiles = remainder
#if len(infiles)>1:
#    infiles.sort(key=skipper_utils.decodeRunnum)

#if (len(sys.argv)<3):
#    print("not enough input args")
#    sys.exit(1)
if outfilename=="":
    filename = remainder[0].split('/')[-1] #strip off directory path
    outfilename = filename.rpartition(".")[0] #strip off file extension, if any
    if outfilename=="":
        outfilename = filename
    outfilename = "hotcol_"+outfilename
    print "no output filename supplied, using default: "+outfilename

data = TChain("hitSumm")
pix = TChain("calPixTree")

#default dimensions - for prototype CCD
ccdPrescan = 7
ccdCol = 886
ccdRow = 6144
nRow = 700
nCol = 450

numRuns = 0
readoutDays = 1.0
for i in range(0,len(infiles)):
    print(infiles[i])
    data.Add(infiles[i])
    pix.Add(infiles[i])
    numRuns += 1
    thefile = TFile(infiles[i])
    header = thefile.Get("headerTree_0")
    header.GetEntry(0)
    if i==0: #get CCD dimensions from the first file
        ccdPrescan = int(getattr(header, "CCDNPRES"))
        ccdCol = int(getattr(header, "CCDNCOL"))
        ccdRow = int(getattr(header, "CCDNROW"))
        nRow = int(getattr(header, "NROW"))
        nCol = int(getattr(header, "NCOL"))
    #ogl = skipper_utils.getHeaderValue(header,"OGAL")
    #swl = skipper_utils.getHeaderValue(header,"SWAL")
    #print("OGL={0}, SWL={1}".format(ogl,swl))



#x=0 and y=0 look weird (negative or extreme values)
#for prototype CCD:
#x=[370,450] is overscan
#x=[1,7] is prescan
#x=[8,369] is active area - matches dimension in paper (362 pixels)
#y=[625,700] is vertical overscan
#y=[1,624] is active area - matches dimension in paper (624 pixels)
#y=624 (top row of active area) seems to have excess charge, so we exclude it
#x=8 (first column of active area) seems to have excess charge, so we exclude it

prescanEdge = ccdPrescan + 1 #X of first active pixel
activeEdgeX = ccdPrescan + ccdCol/2 #X of last active pixel
#if nRow>ccdRow:
#    activeEdgeY = ccdRow
#else:
#    activeEdgeY = ccdRow/2 #Y of last active pixel
activeEdgeY = ccdRow/2


pix.SetAlias("edgedist","min(min(x-{0}, {1}-x), min(y-1, {2}-y))".format(prescanEdge, activeEdgeX, activeEdgeY)) #distance from closest edge, <0 for pixels outside the active area
print pix.GetAlias("edgedist")

epixThr = 0.5

edgeMask = 0
activecutX = "xBary>10"
activecutY = "yBary>0 && yBary<={0}".format(activeEdgeY)
activecutpix = "edgedist>=0"
edgecutpix = "edgedist>{0}".format(edgeMask)
edgecut20 = "edgedist>={0} && (ohdu!=3 || x>=100)".format(20)
flagcut = "!(flag & (4+16+128+256+512+4096))" #4=inBleedZone, 16=crossTalk, 128=serial register hit, 256=low-E cluster, 512=bad pixel, 4096=extended bleed
maskcut = "!(mask & (4+16+128+256+512+4096))" #4=inBleedZone, 16=crossTalk, 128=serial register hit, 256=low-E cluster, 512=bad pixel, 4096=extended bleed
halocutloose = "min(distance, edgedist)>=10"
halocut2e = "min(distance, edgedist)>=20 && bleedY>=50"
halocut1e = "min(distance, edgedist)>=60"
pixcut = "ePix>0.6"
bigpixcut = "ePix>1.6"
pix2ecut = "ePix>1.6 && ePix<2.6"
pix1ecut = "ePix>0.6 && ePix<1.6"


c = TCanvas("c","c",1200,900);
c.Print(outfilename+".pdf[")

outfile = TFile(outfilename+".root","RECREATE")

latex = TLatex()
latex.SetNDC(True)

c.SetLogy(1)

hduCut = "(" + "||".join(["ohdu=={0}".format(x) for x in HDU_LIST]) + ")"
nHdu = len(HDU_LIST)


histList = [] #not used, but needed so Python doesn't delete the histograms


gStyle.SetOptStat(0)

pol0 = TF1("pol0","[0]",0,500)

if True:
    hotcoldict = {}
    c.Divide(3,len(HDU_LIST))
    for i,hdu in enumerate(HDU_LIST):
        hotcoldict[hdu] = []
        c.cd(3*i+1)
        gPad.SetLogz(1)

        #we subtract epixThr-0.5 from the pixel value, so 0.5 is the boundary between 0e- and 1e-
        #this shifts some 0e- pixels into the [-1.5,-0.5] bin
        pix.Draw("ePix-{1}:x>>hotcols2d{0}(500,-0.5,499.5,7,-1.5,5.5)".format(hdu,epixThr-0.5),"&&".join([edgecut20, maskcut, "ohdu=={0}".format(hdu)]),"colz")
        hotcols2d = gDirectory.Get("hotcols2d{0}".format(hdu))
        histList.append(hotcols2d)

        hotcolsDenom = hotcols2d.ProjectionX("hotcolsDenom{0}".format(hdu))
        histList.append(hotcolsDenom)

        bleedRates = TH1F("bleedRates{0}".format(hdu), "bleed rates, HDU {0}".format(hdu), 100, 0.0, 0.02)
        histList.append(bleedRates)

        hotcols = hotcols2d.ProjectionX("hotcols{0}".format(hdu))
        hotcols.Reset()
        histList.append(hotcols)
        for ibinX in range(1,hotcols2d.GetXaxis().GetNbins()+1):#fill with the number of electrons
            x = hotcols2d.GetXaxis().GetBinCenter(ibinX)
            n_electrons = 0.0

            for ibinY in range(1,hotcols2d.GetYaxis().GetNbins()+1):
                y = hotcols2d.GetYaxis().GetBinCenter(ibinY)
                bincount = hotcols2d.GetBinContent(ibinX, ibinY)
                if y>0:
                    n_electrons += bincount*y
            hotcols.SetBinContent(ibinX,n_electrons)
            hotcols.SetBinError(ibinX,np.sqrt(n_electrons))
            if hotcolsDenom.GetBinContent(ibinX)>0:
                bleedRates.Fill(hotcols.GetBinContent(ibinX)/hotcolsDenom.GetBinContent(ibinX))

        c.cd(3*i+2)
        bleedRates.Draw()

        c.cd(3*i+3)
        meanrate = hotcols.Integral()/hotcolsDenom.Integral()
        hotcols.Divide(hotcolsDenom)
        hotcols.Draw()

        rates = []
        for ibin in range(1,hotcols.GetXaxis().GetNbins()+1):
            x = int(hotcols.GetXaxis().GetBinCenter(ibin))
            binrate = hotcols.GetBinContent(ibin)
            rates.append(binrate)
        rates = np.array(rates)
        rate_median = np.median(rates)
        rate_MAD = np.median(np.abs(rates - rate_median))
        #rate_cut = rate_median + 2.5*1.4826*rate_MAD
        rate_cut = 2.0*meanrate
        print "meanrate {0}, median {1}, MAD {2}, MAD cut {3}".format(meanrate, rate_median, rate_MAD, rate_cut)
        pol0.SetParameter(0,rate_cut)
        pol0.DrawCopy("same")

        for ibin in range(1,hotcols.GetXaxis().GetNbins()+1):
            x = int(hotcols.GetXaxis().GetBinCenter(ibin))
            binrate = hotcols.GetBinContent(ibin)
            #if (binrate > 2.0*meanrate):
            if (binrate > rate_cut):
                print hdu,x,binrate
                hotcoldict[hdu].append((x,binrate))

        for ibin in range(1,hotcols2d.GetXaxis().GetNbins()+1):
            x = int(hotcols2d.GetXaxis().GetBinCenter(ibin))
            bincount = hotcols2d.GetBinContent(ibin,4) #bin 4 is the 2e- bin
            if (bincount >= 2):
                print hdu,x,bincount
                hotcoldict[hdu].append((x,bincount))
    c.cd()
    c.cd()
    c.Print(outfilename+".pdf");
    c.Clear()

    with open(outfilename+'.xml','w') as f:
        hotcolXML = ET.Element('badCols')
        hotcolXML.tail='\n'
        for hdu in hotcoldict:
            for x,n in hotcoldict[hdu]:
                colEle = ET.SubElement(hotcolXML,'column')
                colEle.set('hdu',str(hdu))
                colEle.set('x',str(x))
                colEle.tail='\n'
        hotcoltree = ET.ElementTree(hotcolXML)
        hotcoltree.write(f)
    
    #with open(outfilename+"_hotcol.json", 'w') as hotcolfile:
    #    json.dump(hotcoldict, hotcolfile)





c.Print(outfilename+".pdf]");
outfile.Write()
outfile.Close()

sys.exit(0)
