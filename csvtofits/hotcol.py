#!/usr/bin/env python
import sys, getopt
import array
import datetime
import json
import numpy as np
import xml.etree.ElementTree as ET
from scipy.stats import poisson, combine_pvalues
#from scipy.special import lambertw #for calculating mu from the 1e count
import struct

#these parameters will need to be changed as appropriate for the data you're looking at - they are not auto-detected
HDU_LIST = [1]

#number of CCDs in the dataset: this is used to set p-values
N_CCD = 1

outfilename=""
findPixels = False

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'o:ph')

for opt, arg in options:
    if opt == '-o':
        outfilename = arg
    elif opt == '-p':
        findPixels = True
    elif opt == '-h':
        print "\nUsage: "+sys.argv[0]+" <output basename> <root files>"
        print "Arguments: "
        print "\t-o: basename for output file"
        print "\t-p: find hot pixels"
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

def stripFilename(name):
    name = name.split('/')[-1] #strip off directory path
    name = name.rpartition(".")[0] #strip off file extension, if any
    return name

#if (len(sys.argv)<3):
#    print("not enough input args")
#    sys.exit(1)
if outfilename=="":
    outfilename = stripFilename(infiles[0])
    if outfilename=="":
        outfilename = filename
    outfilename = "hotcol_"+outfilename
    print "no output filename supplied, using default: "+outfilename


c = TCanvas("c","c",1200,900);
c.Print(outfilename+".pdf[")

outfile = TFile(outfilename+".root","RECREATE")

latex = TLatex()
latex.SetNDC(True)

histList = [] #not used, but needed so Python doesn't delete the histograms

gStyle.SetOptStat(0)

pol0 = TF1("pol0","[0]",0,10000)

np.set_printoptions(precision=2, linewidth=300)

def combinePvals(pvalList): #a few ways to do this - we do the simplest for now
    return np.min(pvalList, axis=0)
    #return np.prod(pvalList,0)
    #return np.apply_along_axis(combine_pvalues, 0, pvalList)[1,...]

def calculateChunkPvals(goodCells, rateList, cellArrayList, nChunks):
    chunkPlist = np.zeros((len(cellArrayList),nChunks))
    #if nDim==1:
    for i, cellArray in enumerate(cellArrayList):
        chunks = zip(np.array_split(cellArray,nChunks), np.array_split(goodCells,nChunks))
        chunkPlist[i] = [poisson.sf(cData[cGood,1].sum()-0.5, cData[cGood,0].sum()*2.0*rateList[i]) for cData, cGood in chunks]
    return chunkPlist

def calculatePvals(goodCells, cellArrayList):
    cellShape = cellArrayList[0].shape[:-1]
    nDim = len(cellShape)
    pvalList = np.zeros((len(cellArrayList),)+cellShape)
    rateList = np.zeros(len(cellArrayList))
    for i, cellArray in enumerate(cellArrayList):
        goodDenom = np.sum(cellArray[goodCells,0])
        goodHits = np.sum(cellArray[goodCells,1])
        rateGood = goodHits/goodDenom
        rateList[i] = rateGood

        #update p-values
        #survival function = probability of getting nEle (or greater) given mu=rate*nPix
        #we double the rate to account for variations (e.g. the SC at the right edge is double the mean value)
        pvalList[i] = poisson.sf(cellArray[...,1]-0.5, cellArray[...,0]*2.0*rateGood)
    return rateList, pvalList

def addNeighbors(pvals, goodCells, badX, addCut):
    for iseed in badX:
        for iadd in range(iseed-1, -1, -1): #go backward
            if not goodCells[iadd] or pvals[iadd]>addCut: #already found/not a bad cell
                break
            #print("adding {0}, p={1}".format(iadd,pvals[iadd]))
            goodCells[iadd] = False
        for iadd in range(iseed+1, len(goodCells)): #go forward
            if not goodCells[iadd] or pvals[iadd]>addCut: #already found/not a bad cell
                break
            #print("adding {0}, p={1}".format(iadd,pvals[iadd]))
            goodCells[iadd] = False

def findBadCells(data, nCells):
    cellArrayList, nameList = zip(*data)

    #"cells" can be pixels, columns, whatever
    #iterate:
    # 1. mark cells as bad "seeds" if they have an improbably high number of hits given the estimated rate
    # 2. add neighbors to the seeds if they also have elevated rates
    # 3. the average of the good cells is the new estimated rate
    #we start by assuming all cells are good
    #as we mark more cells as bad, the estimated rate goes down
    #when the rate stops changing, we have a stable set of good cells

    pCut = 0.2/(nCells*4*N_CCD) #expect to remove 0.2 cells over the CCD
    addCut = 5e-2 #p-value threshold for also removing a neighbor
    #print("pCut=",pCut)
    nChunks = 8
    chunkCut = 2.7e-3/(nChunks*4*N_CCD) #3 sigma global

    cellShape = cellArrayList[0].shape[:-1]
    goodCells = np.full(cellShape, True)

    nDim = len(cellShape)

    while True:
        while True: #run until we have a stable set of good cells
            oldGood = goodCells

            rateList, pvalList = calculatePvals(goodCells, cellArrayList)
            pvals = combinePvals(pvalList)
            goodCells = pvals>pCut

            badX = np.transpose(np.logical_not(goodCells).nonzero())
            if nDim==1: #check neighbors of bad cells
                addNeighbors(pvals, goodCells, badX, addCut)
                #TODO: check blocks of good cells

            goodCells = np.logical_and(oldGood, goodCells) #since the rate monotonically decreases, goodCells should also be a strictly shrinking set, but we have to do this so we don't lose any cells that were marked bad by the chunk check
            nNew = np.count_nonzero(goodCells!=oldGood)
            if nNew==0: break
            print("removed {0} cells".format(nNew))
        
        if nDim==1:
            chunkPlist = calculateChunkPvals(goodCells, rateList, cellArrayList, nChunks)
            chunkPvals = combinePvals(chunkPlist)
            #print(chunkPvals)
            badChunks = (chunkPvals<chunkCut).nonzero()[0]
            if len(badChunks)==0: break #no bad chunks, we are done
            for iChunk in badChunks:
                np.copyto(np.array_split(goodCells,nChunks)[iChunk], False)
            print("removed {0} chunks".format(len(badChunks)))
        else:
            break


    badX = np.transpose(np.logical_not(goodCells).nonzero())
    #print(badX)
    #badX = badCells[:,2:].astype(int)
    if badX.shape[1]==1: #reduce to 1-D if we can
        badX = badX.ravel()
    else: #convert to list of tuples
        badX = map(tuple, badX)

    if nDim==1:
        cellType = "cols"
    else:
        cellType = "pix"
    print("found {0} bad {1}".format(len(badX), cellType))
    for i, name in enumerate(nameList):
        medianN = np.mean(cellArrayList[i][goodCells,0])
        #cellArray, name in data:
        rateGood = rateList[i]
        print("{0}: rate {1}, median denominator {2}, expected count {3}, seed >{4}, add >{5}".format(name, rateGood, medianN, rateGood*medianN, poisson.isf(pCut, rateGood*medianN), poisson.isf(addCut, rateGood*medianN)))

    return (badX, goodCells, pvalList)

def plotPvals(goodCells, pvals, hname, stackGood, stackBad, minP, color):
    #make a histogram of the p-values, with logarithmic binning
    #add the histogram to the stack
    pedges = array.array('d')
    #minP = np.floor(np.log10(np.min(pvals)))
    #minP = -5.0
    maxP = 0.0
    nPbins = int(-10*minP)
    for ix in range(0,nPbins+1):
        pedges.append(10**(minP+(ix-0.5)*(maxP-minP)/(nPbins-1)))
    
    pvals_toplot = np.maximum(pvals,10**minP)
    hPvals_good = TH1F(hname+"_good", "p-values", nPbins, pedges)
    hPvals_good.SetLineColor(color)
    for p in pvals_toplot[goodCells]:
        hPvals_good.Fill(p)
    histList.append(hPvals_good)
    stackGood.Add(hPvals_good)
    
    hPvals_bad = TH1F(hname+"_bad", "p-values", nPbins, pedges)
    hPvals_bad.SetLineColor(color)
    for p in pvals_toplot[np.logical_not(goodCells)]:
        hPvals_bad.Fill(p)
    histList.append(hPvals_bad)
    stackBad.Add(hPvals_bad)

def readTH1D(h): #get the histogram's data array, discard over/underflow
    #note that TH2.Projection always returns a TH1D, even if you started with a TH2F
    #arrTemp = struct.unpack_from(str(h.GetSize())+'d', h.GetArray())
    arrTemp = np.frombuffer(h.GetArray(), dtype=np.float64, count=h.GetSize())
    return arrTemp[1:-1]
def readTH2F(h, nX, nY): #get the histogram's data array, reshape to 2-D and discard over/underflow
    #arrTemp = struct.unpack_from(str(h.GetSize())+'f', h.GetArray())
    arrTemp = np.frombuffer(h.GetArray(), dtype=np.float32, count=h.GetSize())
    return arrTemp.reshape(((nY+2),(nX+2)))[1:-1,1:-1].T
def readTH2D(h, nX, nY): #get the histogram's data array, reshape to 2-D and discard over/underflow
    #arrTemp = struct.unpack_from(str(h.GetSize())+'f', h.GetArray())
    arrTemp = np.frombuffer(h.GetArray(), dtype=np.float64, count=h.GetSize())
    return arrTemp.reshape(((nY+2),(nX+2)))[1:-1,1:-1].T

nplots = 3

coldatadict = {}
pixdatadict = {}
hotcoldict = {}
hotpixdict = {}

for iFile, filename in enumerate(infiles):
    thefile = TFile(filename)
    outfile.cd()
    c.Divide(nplots,len(HDU_LIST))
    for i,hdu in enumerate(HDU_LIST):
        print("reading {0}, HDU {1}".format(filename, hdu))
        if hdu not in hotcoldict: #initialize dictionaries
            coldatadict[hdu] = []
            pixdatadict[hdu] = []
            hotcoldict[hdu] = set()
            hotpixdict[hdu] = set()
        
        c.cd(nplots*i+1)
        gPad.SetLogz(1)

        hotcols2d = thefile.Get("hotcols2d_{0}".format(hdu))
        histList.append(hotcols2d)
        hotcols2d.Draw("colz")

        #number of unmasked pixels per col
        hotcolsDenom = hotcols2d.ProjectionX("hotcols_all_{0}".format(hdu))
        histList.append(hotcolsDenom)

        #technically, we should use binomial statistics for these? but mu is small so it doesn't matter
        #number of pixels with >=1e per col
        hotcols1e = hotcols2d.ProjectionX("hotcols_1_{0}".format(hdu),2,2)
        histList.append(hotcols1e)
        #number of pixels with >=2e per col
        hotcols2e = hotcols2d.ProjectionX("hotcols_2_{0}".format(hdu),3,-1)
        histList.append(hotcols2e)

        #number of electrons per col
        hotcols = hotcols2d.ProjectionX("hotcols_{0}".format(hdu))
        hotcols.Reset()
        histList.append(hotcols)
        nX = hotcols2d.GetXaxis().GetNbins()
        nY = hotcols2d.GetYaxis().GetNbins()
        hotcolsData = readTH2F(hotcols2d, nX, nY)
        n_electrons = np.zeros(nX)
        for ibinY in range(1,nY+1):
            y = max(0, round(hotcols2d.GetYaxis().GetBinCenter(ibinY)))
            if y>0:
                n_electrons += y*hotcolsData[:,ibinY-1]
        for ibinX in range(1,nX+1):
            x = hotcols2d.GetXaxis().GetBinCenter(ibinX)
            hotcols.SetBinContent(ibinX,n_electrons[ibinX-1])

        c.cd(nplots*i+2)
        
        gPad.SetLogy(1)
        hs = THStack("hs1d_{0}".format(hdu),"pixel counts")
        histList.append(hs)
        hs.Add(hotcolsDenom)
        hs.Add(hotcols.Clone())
        hs.Draw("nostack")

        c.cd(nplots*i+3)
        gPad.SetLogy(1)
        colrates = hotcols.Clone("colrates_{0}".format(hdu))
        histList.append(colrates)
        meanrate = colrates.Integral()/hotcolsDenom.Integral()
        colrates.Sumw2()
        colrates.Divide(hotcolsDenom)
        colrates.Draw()

        rates = readTH1D(colrates)
        rate_median = np.median(rates)
        rate_MAD = np.median(np.abs(rates - rate_median))
        rate_cut = rate_median + 2.5*1.4826*rate_MAD
        #rate_cut = 2.0*meanrate
        print "old method: meanrate {0}, median {1}, MAD {2}, MAD cut {3}".format(meanrate, rate_median, rate_MAD, rate_cut)
        #pol0.SetParameter(0,rateGoodCols)
        #pol0.DrawCopy("same")


        nCols = hotcolsDenom.GetXaxis().GetNbins()
        minX = int(hotcolsDenom.GetXaxis().GetBinCenter(1)) #center of first bin
        maxX = int(hotcolsDenom.GetXaxis().GetBinCenter(nCols)) #center of last bin
        dataDenom = readTH1D(hotcolsDenom)
        data1e = readTH1D(hotcols1e)
        data2e = readTH1D(hotcols2e)

        col1eArray = np.zeros((maxX+1,2))
        col2eArray = np.zeros((maxX+1,2))
        col1eArray[minX:,0] = dataDenom
        col1eArray[minX:,1] = data1e
        col2eArray[minX:,0] = dataDenom
        col2eArray[minX:,1] = data2e

        coldatadict[hdu].append((col1eArray, stripFilename(filename)+" cols (1 e)"))
        coldatadict[hdu].append((col2eArray, stripFilename(filename)+" cols (2+ e)"))

        if findPixels:
            hotpix = thefile.Get("hotpix_{0}".format(hdu))
            hotpixDenom = hotpix.Project3D("all_yx")
            hotpix.GetZaxis().SetRange(2,2)
            hotpix1 = hotpix.Project3D("1_yx")
            hotpix.GetZaxis().SetRange(3,hotpix.GetZaxis().GetNbins()+1)
            hotpix2 = hotpix.Project3D("2_yx")

            nX = hotpix1.GetXaxis().GetNbins()
            nY = hotpix1.GetYaxis().GetNbins()
            nPix = nX*nY
            dataDenom = readTH2D(hotpixDenom, nX, nY)
            data1e = readTH2D(hotpix1, nX, nY)
            data2e = readTH2D(hotpix2, nX, nY)
            minX = int(hotpix1.GetXaxis().GetBinCenter(1)) #center of first bin
            minY = int(hotpix1.GetYaxis().GetBinCenter(1)) #center of first bin
            maxX = int(hotpix1.GetXaxis().GetBinCenter(nX)) #center of last bin
            maxY = int(hotpix1.GetYaxis().GetBinCenter(nY)) #center of last bin

            pixArray1 = np.zeros((maxX+1, maxY+1, 2))
            pixArray2 = np.zeros((maxX+1, maxY+1, 2))
            pixArray1[minX:,minY:,0] = dataDenom
            pixArray1[minX:,minY:,1] = data1e
            pixArray2[minX:,minY:,0] = dataDenom
            pixArray2[minX:,minY:,1] = data2e
            
            pixdatadict[hdu].append((pixArray1, stripFilename(filename)+" pix (1 e)"))
            pixdatadict[hdu].append((pixArray2, stripFilename(filename)+" pix (2+ e)"))

    c.cd()
    c.Print(outfilename+".pdf");
    c.Clear()

nplots = 4
c.Divide(nplots,len(HDU_LIST))
for iHdu,hdu in enumerate(HDU_LIST):
    print("masking cols, HDU {0}".format(hdu))

    badX, goodCols, pvalList = findBadCells(coldatadict[hdu], nCols)
    hotcoldict[hdu].update(badX) #add these to the bad cols

    hsPcol1_good = THStack("hs_colpvals1e_{0}_good".format(hdu),"good col p-values, 1e")
    hsPcol2_good = THStack("hs_colpvals2e_{0}_good".format(hdu),"good col p-values, 2+ e")
    histList.append(hsPcol1_good)
    histList.append(hsPcol2_good)
    hsPcol1_bad = THStack("hs_colpvals1e_{0}_bad".format(hdu),"bad col p-values, 1e")
    hsPcol2_bad = THStack("hs_colpvals2e_{0}_bad".format(hdu),"bad col p-values, 2+ e")
    histList.append(hsPcol1_bad)
    histList.append(hsPcol2_bad)

    for iFile in range(len(pvalList)/2):
        plotPvals(goodCols, pvalList[2*iFile], "h_colpvals1e_{0}_{1}".format(iFile,hdu), hsPcol1_good, hsPcol1_bad, -5.0, iFile+1)
        plotPvals(goodCols, pvalList[2*iFile+1], "h_colpvals2e_{0}_{1}".format(iFile,hdu), hsPcol2_good, hsPcol2_bad, -5.0, iFile+1)

    c.cd(nplots*iHdu+1)
    gPad.SetLogy(1)
    gPad.SetLogx(1)
    hsPcol1_good.Draw("nostack")
    c.cd(nplots*iHdu+2)
    gPad.SetLogy(1)
    gPad.SetLogx(1)
    hsPcol2_good.Draw("nostack")
    c.cd(nplots*iHdu+3)
    gPad.SetLogy(1)
    gPad.SetLogx(1)
    hsPcol1_bad.Draw("nostack")
    c.cd(nplots*iHdu+4)
    gPad.SetLogy(1)
    gPad.SetLogx(1)
    hsPcol2_bad.Draw("nostack")


c.cd()
c.Print(outfilename+".pdf");
c.Clear()

if findPixels:
    nplots = 4
    c.Divide(nplots,len(HDU_LIST))
    for iHdu,hdu in enumerate(HDU_LIST):
        print("masking pix, HDU {0}".format(hdu))

        for cellArray, name in pixdatadict[hdu]: #mask bad cols
            cellArray[list(hotcoldict[hdu]),:,:] = 0

        badPix, goodPix, pvalList = findBadCells(pixdatadict[hdu], nPix)
        hotpixdict[hdu].update(badPix)

        hsPpix1_good = THStack("hs_pixpvals1e_{0}_good".format(hdu),"good pix p-values, 1e")
        hsPpix2_good = THStack("hs_pixpvals2e_{0}_good".format(hdu),"good pix p-values, 2+ e")
        histList.append(hsPpix1_good)
        histList.append(hsPpix2_good)
        hsPpix1_bad = THStack("hs_pixpvals1e_{0}_bad".format(hdu),"bad pix p-values, 1e")
        hsPpix2_bad = THStack("hs_pixpvals2e_{0}_bad".format(hdu),"bad pix p-values, 2+ e")
        histList.append(hsPpix1_bad)
        histList.append(hsPpix2_bad)

        for iFile in range(len(pvalList)/2):
            plotPvals(goodPix, pvalList[2*iFile], "h_pixpvals1e_{0}_{1}".format(iFile,hdu), hsPpix1_good, hsPpix1_bad, -8.0, iFile+1)
            plotPvals(goodPix, pvalList[2*iFile+1], "h_pixpvals2e_{0}_{1}".format(iFile,hdu), hsPpix2_good, hsPpix2_bad, -8.0, iFile+1)
        #print(badPix)

        c.cd(nplots*iHdu+1)
        gPad.SetLogy(1)
        gPad.SetLogx(1)
        hsPpix1_good.Draw("nostack")
        c.cd(nplots*iHdu+2)
        gPad.SetLogy(1)
        gPad.SetLogx(1)
        hsPpix2_good.Draw("nostack")
        c.cd(nplots*iHdu+3)
        gPad.SetLogy(1)
        gPad.SetLogx(1)
        hsPpix1_bad.Draw("nostack")
        c.cd(nplots*iHdu+4)
        gPad.SetLogy(1)
        gPad.SetLogx(1)
        hsPpix2_bad.Draw("nostack")
        print("heu")

        print("{0} hot pixels before merging".format(len(hotpixdict[hdu])))
        pix2col = {}
        for x,y in hotpixdict[hdu]:
            if x not in pix2col:
                pix2col[x] = set()
            pix2col[x].add(y)
        for x,ySet in pix2col.items():
            yList = sorted(ySet)
            if (len(yList)>2 and yList[-1]-yList[0]>10): #at least 3 pixels in the same col, separated sufficiently in Y (so it's not just some local hot spot)
                print("found {0} pixels with x={1}: merging into a bad column".format(len(yList), x))
                hotcoldict[hdu].add(x)
                for y in yList:
                    hotpixdict[hdu].remove((x,y))

    c.cd()
    c.Print(outfilename+".pdf");
    c.Clear()

for hdu in hotcoldict:
    hotcolsmerged = []
    firstX = -99
    lastX = -99
    for x in sorted(hotcoldict[hdu]):
        if x==lastX+1:
            lastX = x
        else:
            if firstX>=0:
                if firstX==lastX:
                    hotcolsmerged.append(firstX)
                else:
                    hotcolsmerged.append((firstX, lastX))
            firstX = x
            lastX = x

    if firstX>=0:
        if firstX==lastX:
            hotcolsmerged.append(firstX)
        else:
            hotcolsmerged.append((firstX, lastX))
    #print("bad cols:", sorted(hotcoldict[hdu]))
    print("bad cols:", hotcolsmerged)
    hotcoldict[hdu] = hotcolsmerged
    if findPixels: print("bad pix:", sorted(hotpixdict[hdu]))

with open(outfilename+'.xml','w') as f:
    hotcolXML = ET.Element('badCols')
    hotcolXML.text='\n'
    hotcolXML.tail='\n'
    for hdu in hotcoldict:
        firstX = -99
        lastX = -99
        for x in hotcoldict[hdu]:
            if isinstance(x,tuple):
                colEle = ET.SubElement(hotcolXML,'cRange')
                colEle.set('hdu',str(hdu))
                colEle.set('x1',str(x[0]))
                colEle.set('x2',str(x[1]))
            else:
                colEle = ET.SubElement(hotcolXML,'column')
                colEle.set('hdu',str(hdu))
                colEle.set('x',str(x))
            colEle.tail='\n'
    hotcoltree = ET.ElementTree(hotcolXML)
    hotcoltree.write(f)

    if findPixels:
        hotpixXML = ET.Element('badPixels')
        hotpixXML.text='\n'
        hotpixXML.tail='\n'
        for hdu in hotpixdict:
            for x,y in sorted(hotpixdict[hdu]):
                pixEle = ET.SubElement(hotpixXML,'pixel')
                pixEle.set('hdu',str(hdu))
                pixEle.set('x',str(x))
                pixEle.set('y',str(y))
                pixEle.tail='\n'
        hotpixtree = ET.ElementTree(hotpixXML)
        hotpixtree.write(f)

#with open(outfilename+"_hotcol.json", 'w') as hotcolfile:
#    json.dump(hotcoldict, hotcolfile)


c.Print(outfilename+".pdf]");
outfile.Write()
outfile.Close()

sys.exit(0)


