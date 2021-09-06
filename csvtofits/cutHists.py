#!/usr/bin/env python
import sys, getopt
import array

#these parameters will need to be changed as appropriate for the data you're looking at - they are not auto-detected
HDU_LIST = [1]
drawPix = False #set to True if you want to see the pixel histos - they make the PDF big

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
from ROOT import TH1, TH1F, TH2F, TH3F
from ROOT import gPad
from ROOT import TLatex

gStyle.SetOptStat(110011)

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
    outfilename = "cutHists_"+outfilename
    print "no output filename supplied, using default: "+outfilename

pix = TChain("calPixTree")

#default dimensions - for prototype CCD
ccdPrescan = 7
ccdCol = 886
ccdRow = 6144
nRow = 700
nCol = 450
maxBleedX = 100
maxHalo = 100

numRuns = 0
for i in range(0,len(infiles)):
    print(infiles[i])
    pix.Add(infiles[i])
    numRuns += 1
    thefile = TFile(infiles[i])
    if i==0: #get CCD dimensions from the first file
        config = thefile.Get("config")
        config.GetEntry(0)
        maxBleedX = config.bleedX
        maxHalo = config.halo
        header = thefile.Get("headerTree_0")
        header.GetEntry(0)
        ccdPrescan = int(getattr(header, "CCDNPRES"))
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
#Y of first active pixel is 1
activeEdgeY = ccdRow/2 #Y of last active pixel

minCol = (prescanEdge)//binCols #image X of first active pixel
maxCol = (activeEdgeX)//binCols #image X of last active pixel
maxRow = (activeEdgeY-1)//binRows + 1 #image Y of last active pixel (first is always 1)
pixdims = "({0},{1},{2},{3},{4},{5})".format(maxCol-minCol+1, minCol-0.5, maxCol+0.5, maxRow, 0.5, maxRow+0.5)

#image pix x is physical [x*bin,(x+1)*bin-1]
#image pix y is physical [(y-1)*bin+1,y*bin]
edgedistX = "min(x*{3}-{0}, {1}-((x+1)*{3}-1))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
edgedistY = "min((y-1)*{4}, {2}-(y*{4}))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
pix.SetAlias("edgedistX",edgedistX)
pix.SetAlias("edgedistY",edgedistY)
pix.SetAlias("edgedist","min({0},{1})".format(edgedistX, edgedistY))
print(pix.GetAlias("edgedistX"), pix.GetAlias("edgedistY"))


epixThr = 0.5
nPixBins = 6
pixBins = array.array('d')
pixBins.append(-1.0)
pixBins.extend([i+epixThr for i in range(nPixBins)])

#bins for bleed histogram
rowBins = array.array('d', [x-0.5 for x in range(1, maxRow+1)])
colBins = array.array('d', [x-0.5 for x in range(minCol, maxCol+2)])
bleedBins = array.array('d', [x-0.5 for x in range(maxBleedX+3)])

pixcut = "ePix>{0}".format(epixThr)
pix1ecut = "ePix>{0} && ePix<{1}".format(epixThr,epixThr+1.0)
pix2ecut = "ePix>{0}".format(epixThr+1.0)

#edgecut20 = "edgedist>={0} && (ohdu!=3 || x>=100)".format(20)
edgecutY20 = "edgedistY>={0}".format(20)
edgecut20 = "edgedist>={0}".format(20)
edgecut0 = "edgedist>=0"

maskcutCommon = "!(mask & (16+128+256))" #16=crossTalk, 128=serial register hit, 256=low-E cluster, 16384=full well
maskcutHotCol = "!(mask & (4+16+128+256+512+4096))" #4=inBleedZone, 16=crossTalk, 128=serial register hit, 256=low-E cluster, 512=bad pixel, 4096=extended bleed
maskcutHotPix = "!(mask & (16+128+256))" #16=crossTalk, 128=serial register hit, 256=low-E cluster
maskcutBleed = "!(mask & (16+128+256+512+1024))" #16=crossTalk, 128=serial register hit, 256=low-E cluster, 512=bad pixel, 1024=bad col
#maskcutHalo = "!(mask & (4+16+128+256+512+1024))" #4=bleed, 16=crossTalk, 128=serial register hit, 256=low-E cluster, 512=bad pixel, 1024=bad col
maskcutHalo = "!(mask & (4+16+128+256+512+1024+2048+4096+8192))" #4=bleed, 16=crossTalk, 128=serial register hit, 256=low-E cluster, 512=bad pixel, 1024=bad col
#maskcutHotPix = "!(mask & (16+128+256+1024))" #16=crossTalk, 128=serial register hit, 256=low-E cluster, 1024=bad column

c = TCanvas("c","c",1200,900);
c.Print(outfilename+".pdf[")

outfile = TFile(outfilename+".root","RECREATE")

latex = TLatex()
latex.SetNDC(True)

histList = [] #not used, but needed so Python doesn't delete the histograms

gStyle.SetOptStat(0)

if drawPix:
    drawPixOpt = "goff"
    nplots = 5
else:
    drawPixOpt = "colz"
    nplots = 2

c.Divide(nplots,len(HDU_LIST))
for i,hdu in enumerate(HDU_LIST):
    #set an event list so we only loop over the current HDU
    pix.SetEventList(0)
    
    elistname = "e{0}".format(hdu)
    pix.Draw(">>"+elistname,"&&".join(["ohdu=={0}".format(hdu), maskcutCommon, edgecut0]))
    pix.SetEventList(gDirectory.Get(elistname))
    print(pix.GetEntries())

    c.cd(nplots*i+1)
    gPad.SetLogz(1)
    hotcols2d = TH2F("hotcols2d_{0}".format(hdu), "e vs. col", maxCol-minCol+1, minCol-0.5, maxCol+0.5, nPixBins, pixBins)
    pix.Draw("ePix:x>>+hotcols2d_{0}".format(hdu),"&&".join([maskcutHotCol, "ohdu=={0}".format(hdu)]),"colz")
    print(hotcols2d.GetEntries())
    print("&&".join([edgecutY20, maskcutHotCol, "ohdu=={0}".format(hdu)]))
    histList.append(hotcols2d)

    c.cd(nplots*i+2)
    gPad.SetLogz(1)
    halo = TH2F("halo_{0}".format(hdu), "e vs. halo distance", maxHalo+2, -0.5, maxHalo+1.5, nPixBins, pixBins)
    pix.Draw("ePix:min(distance, edgedist)>>+halo_{0}".format(hdu),"&&".join([edgecut0, maskcutHalo, "ohdu=={0}".format(hdu)]), "colz")
    histList.append(halo)

    if i==0:
        halo_all = halo.Clone("halo_all")
    else:
        halo_all.Add(halo)

    hotpix = TH3F("hotpix_{0}".format(hdu), "e vs. pix", len(colBins)-1, colBins, len(rowBins)-1, rowBins, nPixBins, pixBins)
    pix.Draw("ePix:y:x>>+hotpix_{0}".format(hdu),"&&".join([maskcutHotPix, "ohdu=={0}".format(hdu)]), "goff")
    histList.append(hotpix)

    #pix.Draw("y:x>>hotpixDenom_{0}{1}".format(hdu,pixdims),"&&".join([maskcutHotPix, "ohdu=={0}".format(hdu)]),drawPixOpt)
    #hotpixDenom = gDirectory.Get("hotpixDenom{0}".format(hdu))
    #histList.append(hotpixDenom)
#
#    pix.Draw("y:x>>hotpix1_{0}{1}".format(hdu,pixdims),"&&".join([maskcutHotPix, pix1ecut, "ohdu=={0}".format(hdu)]),drawPixOpt)
#    hotpix1 = gDirectory.Get("hotpix1{0}".format(hdu))
#    histList.append(hotpix1)
#
#    pix.Draw("y:x>>hotpix2_{0}{1}".format(hdu,pixdims),"&&".join([maskcutHotPix, pix2ecut, "ohdu=={0}".format(hdu)]),drawPixOpt)
#    hotpix2 = gDirectory.Get("hotpix2{0}".format(hdu))
#    histList.append(hotpix2)

    #c.cd(nplots*i+6)
    #gPad.SetLogz(1)
    #bleedX = TH3F("bleedX_{0}".format(hdu), "e vs. X, bleed pixels", len(colBins)-1, colBins, len(bleedBins)-1, bleedBins, nPixBins, pixBins)
    #pix.Draw("ePix:bleedX:x>>+bleedX_{0}".format(hdu),"&&".join([edgecut0, maskcutBleed, "ohdu=={0}".format(hdu)]), "goff")
    #histList.append(bleedX)

    bleedX_hitX = TH3F("bleedX_hitX_{0}".format(hdu), "e vs. X, bleed pixels", len(colBins)-1, colBins, len(bleedBins)-1, bleedBins, nPixBins, pixBins)
    pix.Draw("ePix:bleedX:x-bleedX>>+bleedX_hitX_{0}".format(hdu),"&&".join([edgecut0, maskcutBleed, "ohdu=={0}".format(hdu)]), "goff")
    histList.append(bleedX_hitX)


c.cd()
c.Print(outfilename+".pdf");
c.Clear()

c.Print(outfilename+".pdf]");
outfile.Write()
outfile.Close()

sys.exit(0)
