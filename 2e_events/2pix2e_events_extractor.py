#!/usr/bin/env python
import sys, getopt, os
import array
import re
from ROOT import gROOT, gStyle, TFile, TTree, TChain, TMVA, TCut, TCanvas, gDirectory, TH1, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex, TList, TEfficiency
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import datetime
import pandas as pd
from openpyxl import load_workbook
import path
from scipy.stats import chi2
import time

gStyle.SetStatW(0.25)
gStyle.SetStatH(0.25)
gROOT.SetBatch(True)
gStyle.SetOptStat(100011)
gStyle.SetOptFit(10001)
gStyle.SetPalette(57)

seconds = time.time()


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

def decodeRunnum(filename,i=1):
    ints = [int(c) for c in re.split(r'(\d+)',filename) if c.isdigit()]
    #print(ints)
    return ints[-i]

def readconfig(infiles,ccdPrescan,ccdCol,ccdRow,nRow,nCol,binRows,binCols,files):
    for i in range(0,len(infiles)):
        thefile = TFile(infiles[i])
        header = thefile.Get("headerTree_0")
        header.GetEntry(0)
        if i==0:
            try:
                cddPrescan = int(getattr(header, "CCDNPRES"))
            except:
                cddPrescan = 7
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
            try:
                start = datetime.datetime.strptime(getattr(header,"DATESTART").split(b'\0',1)[0],"%Y-%m-%dT%H:%M:%S")
                end = datetime.datetime.strptime(getattr(header,"DATEEND").split(b'\0',1)[0],"%Y-%m-%dT%H:%M:%S")
                readoutTime = end-start
            except:
                readoutTime = 0
            #get exposure time
            try:
                exposure = int(getattr(header, "EXPOSURE"))
            except:
                exposure = 0
            try:
                readoutDays = readoutTime.days + (float(readoutTime.seconds) - exposure)/(24*60*60)
            except:
                readoutDays = 0
        try:
            files[int(getattr(header, "LTANAME"))-1].append(infiles[i])
        except:
            files[0].append(infiles[i])
        # print("Adding "+str(infiles[i])+" to module = "+str(int(getattr(header, "LTANAME"))))

    return ccdPrescan,ccdCol,ccdRow,nRow,nCol,binRows,binCols,files,readoutDays,exposure

def writeConfig(nmods,modules,data,outfile):
    outfile.cd()
    #write into configuration of rootfile
    config = TTree("","")
    config.SetName("config")
    config.Branch('nmods',nmods,'nmods/I')
    config.Branch('module',modules,'module[nmods]/I')
    config.Branch('entries',data[0],'entries[nmods]/I')
    config.Branch('effexposure',data[1],'effexposure[nmods]/F')
    config.Branch('1erate',data[2],'1erate[nmods]/F')
    config.Branch('1erateErr',data[3],'1erateErr[nmods]/F')
    config.Branch('1pix2erateUL',data[4],'1pix2erateUL[nmods]/F')
    config.Branch('noise',data[5],'noise[nmods]/F')
    # config.Branch('1pix2erate',data[4],'1pix2erate[nmods]/F')
    # config.Branch('1pix2erateErr',data[5],'1pix2erateErr[nmods]/F')
    # config.Branch('1pix2events',data[6],'1pix2events[nmods]/F')
    # config.Branch('1pix2eventsErr',data[7],'1pix2eventsErr[nmods]/F')
    # config.Fill()
    return config

def fit1e(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows, halo, xcut, halo1e, masks, hduscut ,files, f, outfile, outfilename, PLOTS, readoutDays, exposure):


    
    #### calculate expected 2e- events. As a rough estimation I'll take all quads from each module to calculate the expected 2e- events per pix.
    #range
    pixmax = 4.0*gain; pixmin = -1.3*gain; pixnbin = 250; pixbinning = "({0},{1},{2})".format(pixnbin,pixmin,pixmax); binwidth = (pixmax-pixmin)/pixnbin; aduZero = 0; aduPerElectron = 1; noise = 0.164;

    #fit function
    funcformulas=["[0]*(TMath::Gaus(x,[1]+{0}*[2],[2]*[3],1)*TMath::Poisson({0},[4]) )".format(ipeak) for ipeak in range(0,2)]
    fitfunc = TF1("fitfunc"," + ".join(funcformulas)) #funcformulas is a list so "join" joins its values.

    #canvas and rootfile creation
    c = TCanvas("c","c",1200,900); c.Clear();
    c.Print(outfilename+"_1epeaks.pdf[")

    hists = [] #need it to store hists after plot them

    edgedistX = "min(x*{3}-{0}, {1}-((x+1)*{3}-1))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
    edgedistY = "min((y-1)*{4}, {2}-(y*{4}))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area

    calPixTree = [TChain("calPixTree") for i in range(0,len(files))]
    # data = [array.array('f',[0,0]) for i in range(0,8)] #entries, exposure in gram days,1e- rate, 1e- rateErr, expected 1pix2e- rate, expected 1pix2e- rateErr, expected 1pix2e- events, expected 1pix2e- events uncertainty
    data = [array.array('f',[0,0]) for i in range(0,6)] #entries, exposure in gram days,1e- rate, 1e- rateErr, 1pix2e- events UL (90%)
    data[0]=array.array('i',[0,0]) #first one is integer!

    entries = np.empty((1), dtype="int32"); entries1e = np.empty((1), dtype="float32"); mu = np.empty((1), dtype="float32"); muErr = np.empty((1), dtype="float32"); modules = np.empty((1), dtype="int32"); runs = np.empty((1), dtype="int32"); ohdu = np.empty((1), dtype="int32");
    calPixTree1e = TTree("","")
    calPixTree1e.SetName("calPixTree1e")
    calPixTree1e.Branch('entries',entries,'entries/I')
    calPixTree1e.Branch('entries1e',entries1e,'entries1e/F')
    calPixTree1e.Branch('mu',mu,'mu/F')
    calPixTree1e.Branch('muErr',muErr,'muErr/F')
    calPixTree1e.Branch('LTANAME',modules,'LTANAME/I')
    calPixTree1e.Branch('runID',runs,'runID/I')
    calPixTree1e.Branch('ohdu',ohdu,'ohdu/I')

    f.write("regioncuts \n")
    for I in range(0,len(files)):
        #starters for plot
        # regioncuts="x<"+str(xcut)+" && min(distance,edgedist)>="+str(halo)+" && distance1e>"+str(halo1e)+" && (mask & "+str(masks)+")==0"# 
        # regioncuts="x<"+str(xcut)+" && min(distance,edgedist)>="+str(halo)+" && (mask & "+str(masks)+")==0"
        regioncuts="(mask & "+str(masks)+")==0"

        for iFile in files[I]:
            if True:
                print(iFile)
                thefile = TFile(iFile)
                calPixTreeAux = thefile.Get("calPixTree")
                calPixTreeAux.SetAlias("edgedist","min({0},{1})".format(edgedistX, edgedistY))
                outfile.cd()
                for hdu in range(1,5):
                    hname='h_temp'
                    calPixTreeAux.Draw("ePix>>"+hname+"Module"+str(I+1)+pixbinning,regioncuts+" && ohdu=="+str(hdu))
                    gPad.SetLogy(1)
                    h = gDirectory.Get(hname+"Module"+str(I+1))
                    fitfunc.SetParameters(h.GetEntries()*binwidth,aduZero,aduPerElectron,noise,4e-3)
                    s = h.Fit(fitfunc,"QSL","",aduZero-1*aduPerElectron,aduZero+1.4)
                    entries[0]=int(h.GetEntries())
                    mu[0]=float(s.Parameter(4))
                    muErr[0]=float(s.Error(4))
                    entries1e[0]=mu*entries
                    modules[0]=I+1
                    runs[0]=decodeRunnum(iFile)
                    ohdu[0]=hdu
                    calPixTree1e.Fill()
        outfile.cd()
        for iFile in files[I]:
            calPixTree[I].Add(iFile)
        
        
        calPixTree[I].SetAlias("edgedist","min({0},{1})".format(edgedistX, edgedistY))
        
        #write regioncuts into txt
        f.write(str(regioncuts+" && "+hduscut[I])+"\n")
        
        #draw and fit
        hname='h_ePix'
        calPixTree[I].Draw("ePix>>"+hname+"Module"+str(I+1)+pixbinning,regioncuts+" && "+hduscut[I])
        gPad.SetLogy(1)
        h = gDirectory.Get(hname+"Module"+str(I+1))
        print(h.GetEntries())
        hists.append(h)
        fitfunc.SetParameters(h.GetEntries()*binwidth,aduZero,aduPerElectron,noise,0.00001)
        s = h.Fit(fitfunc,"QSL","",aduZero-1*aduPerElectron,aduZero+1.5)

        #plots into outfile
        if PLOTS==True:
            plot(calPixTree, regioncuts, edgedistX, edgedistY,hists, I)
        
        pixelmass=3.53716875e-06 #grams per pixel
        
        #to write into config tree of outfile 
        data[0][I]=int(h.GetEntries())
        data[1][I]=float(h.GetEntries()*pixelmass*(exposure/86400.0+readoutDays/2))
        data[2][I]=float(s.Parameter(4))
        data[3][I]=float(s.Error(4))
        fPoiss_cdf = TF1("f","ROOT::Math::poisson_cdf("+str(h.GetEntries()*(s.Parameter(4)**2)/2)+", x)",0,4)
        data[4][I]=float(fPoiss_cdf.GetX(.1))
        data[5][I]=float(s.Parameter(3))
        # data[4][I]=float((s.Parameter(4)**2)/2)
        # data[5][I]=float(s.Error(4)*s.Parameter(4))
        # data[6][I]=float(h.GetEntries()*(s.Parameter(4)**2)/2)
        # data[7][I]=float(h.GetEntries()*s.Error(4)*s.Parameter(4))

        c.cd()
        c.Print(outfilename+"_1epeaks.pdf")
        c.Clear()


    # entries = array.array('i',entries)
    # entries1e = array.array('f',entries1e)
    # mu = array.array('f',mu)
    # muErr = array.array('f',muErr)
    # modules = array.array('i',modules)
    # runs = array.array('i',runs)

    print(entries)
    print(modules)
    print(runs)



    c.cd()
    c.Print(outfilename+"_1epeaks.pdf]")

    return calPixTree1e, data

def plot(calPixTree, regioncuts, edgedistX, edgedistY, hists, I):
    calPixTree1e=calPixTree[I].CopyTree(regioncuts+" && "+hduscut[I]+" && ePix>0.7") #need to correct this number at some point
    calPixTree1e.SetAlias("edgedist","min({0},{1})".format(edgedistX, edgedistY))
    for hdu in [1,2,3,4]:
        hname='h_x_ohdu'+str(hdu)
        # calPixTree[I].Draw("x>>"+hname+"Module"+str(I+1)+"(6,0,3200)",regioncuts+" && "+hduscut[I])
        calPixTree[I].Draw("x>>"+hname+"Module"+str(I+1)+"(6,0,3200)",regioncuts+" && ohdu=="+str(hdu))
        h = gDirectory.Get(hname+"Module"+str(I+1))
        hists.append(h)
        hname='h_y_ohdu'+str(hdu)
        calPixTree[I].Draw("y>>"+hname+"Module"+str(I+1)+"(4,0,520)",regioncuts+" && ohdu=="+str(hdu))
        h = gDirectory.Get(hname+"Module"+str(I+1))
        hists.append(h)
        # hname='h_x_y_ohdu'+str(hdu)
        # calPixTree[I].Draw("x:y>>"+hname+"Module"+str(I+1)+"(4,0,520,6,0,3200)",regioncuts+" && "+hduscut[I],"colz")
        # h = gDirectory.Get(hname+"Module"+str(I+1))
        # hists.append(h)
        hname='h_x_1e_ohdu'+str(hdu)
        calPixTree1e.Draw("x>>"+hname+"Module"+str(I+1)+"(6,0,3200)",regioncuts+" && ohdu=="+str(hdu))
        h = gDirectory.Get(hname+"Module"+str(I+1))
        hists.append(h)
        hname='h_y_1e_ohdu'+str(hdu)
        calPixTree1e.Draw("y>>"+hname+"Module"+str(I+1)+"(4,0,520)",regioncuts+" && ohdu=="+str(hdu))
        h = gDirectory.Get(hname+"Module"+str(I+1))
        hists.append(h)
        # hname='h_x_y_1e_ohdu'+str(hdu)
        # calPixTree1e.Draw("x:y>>"+hname+"Module"+str(I+1)+"(4,0,520,6,0,3200)",regioncuts+" && "+hduscut[I],"colz")
        # h = gDirectory.Get(hname+"Module"+str(I+1))
        # hists.append(h)
    return hists

def count2e(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows, halo, xcut, halo1e, masks, hduscut ,modules,files, f):

    #### At this point I will count the 2e- events ####
    edgedistX = "min(xMin*{3}-{0}, {1}-((xMax+1)*{3}-1))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area
    edgedistY = "min((yMin-1)*{4}, {2}-(yMax*{4}))".format(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows) #distance from closest edge, <0 for pixels outside the active area

    #starters
    hitSumm = TChain("hitSumm")
    hitSumm.SetBranchStatus("*",1)		
    hitSummAux=[TTree() for i in modules]
    hitSumm.SetAlias("edgedist","min({0},{1})".format(edgedistX, edgedistY))

    f.write("regioncuts2e \n")
    for I in range(0,len(modules)):
        hitSumm.Reset()
        for iFile in files[I]:
            hitSumm.Add(iFile)
        # regioncuts2e="min(distance,edgedist)>="+str(halo)+" && xMin<"+str(xcut)+" && distance1e>"+str(halo1e)+" && e==2 && n<3 && (flag & "+str(masks)+")==0 &&"+hduscut[I]
        regioncuts2e="min(distance,edgedist)>="+str(halo)+" && xMin<"+str(xcut)+" && e>=2 && e<=4 && n>=1 && n<=4 && (flag & "+str(masks)+")==0 &&"+hduscut[I]
        # regioncuts2e="min(distance,edgedist)>="+str(halo)+" && xMin<"+str(xcut)+" && e==2 && n==1 && ePix[0]>1.7 && (flag & "+str(masks)+")==0 &&"+hduscut[I]
        print(regioncuts2e)
        hitSummAux[I] = hitSumm.CopyTree(regioncuts2e)
        f.write(str(regioncuts2e)+"\n")
            


    lists = TList()
    for i,tree in enumerate(hitSummAux):
        print("2e- events in module "+str(i+1)+": "+str(tree.GetEntries()))
        lists.Add(tree)

    # del hitSummAux
    del hitSumm
    return lists




###########################################################################################################################################################
################################################ Those were the functions we use. ################################################
###########################################################################################################################################################

#default dimensions - for prototype CCD
ccdPrescan = 7
ccdCol = 886
ccdRow = 6144
nRow = 700
nCol = 450
binRows=1
binCols=1

#organize input files.
outfilename = remainder[0]
infiles = remainder[1:]
infiles.sort(key=decodeRunnum) #oh, okey, this sorts in files by the -2 number.
# files=[[],[]]
files=[[]]
nmods=array.array('i', [len(files)])
modules=array.array('i', range(len(files)+1)[1:])

#txt of configuration
f = open(outfilename+".txt", "w")

#read config
ccdPrescan,ccdCol,ccdRow,nRow,nCol,binRows,binCols,files,readoutDays,exposure = readconfig(infiles,ccdPrescan,ccdCol,ccdRow,nRow,nCol,binRows,binCols,files)


#set geometry
prescanEdge = ccdPrescan + 1 #X of first active pixel
activeEdgeX = ccdPrescan + ccdCol/2 #X of last active pixel
activeEdgeY = ccdRow/2 #Y of last active pixel

#set cuts
hdus=[[],[]]
# hdus[0]=['ohdu=='+hdu for hdu in ['1','2','3','4']] #set hdus to analyze
# # hdus[1]=['ohdu=='+hdu for hdu in ['1','2','3','4']] #set hdus to analyze
# hdus[1]=['ohdu=='+hdu for hdu in ['3','4']] #set hdus to analyze
# hduscut=['('+' || '.join(i)+')' for i in hdus]
hdus=[[]]
hdus[0]=['ohdu=='+hdu for hdu in ['1']] #set hdus to analyze
# hdus[1]=['ohdu=='+hdu for hdu in ['1','2','3','4']] #set hdus to analyze
# hdus[1]=['ohdu=='+hdu for hdu in ['3','4']] #set hdus to analyze
hduscut=['('+' || '.join(i)+')' for i in hdus]

#set other cuts and options
halo=00
halo1e=halo1eOpt
masks='(4+16+128+256+512+1024+4096+8192+16384)' #2048==looseCluster
masks='(512)' #2048==looseCluster
xcut=555555
PLOTS=False

file0 = TFile(infiles[0])
#create output file
outfile = TFile(outfilename+".root","RECREATE");

#copy header TTree
header = file0.Get("headerTree_0")
header.CloneTree().Write()

# #analyse 1e- events
# calPixTree1e, data = fit1e(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows, halo, xcut, halo1e, masks, hduscut,files, f, outfile, outfilename, PLOTS, readoutDays, exposure)
# calPixTree1e.Print()
# calPixTree1e.Write()


# #write config branch
# config = writeConfig(nmods,modules,data,outfile)
# config.Fill()

#count 2e events
lists = count2e(prescanEdge, activeEdgeX, activeEdgeY, binCols, binRows, halo, xcut, halo1e, masks, hduscut ,modules,files, f)

#merge lists to create hitSumm2e
hitSummSel = TTree.MergeTrees(lists)
hitSummSel.SetName("hitSumm2e")

outfile.Write()
outfile.Close()

f.close()

print("\n Seconds for first loop =", time.time()-seconds)	
#python -i /home/mariano/tools/SENSEI_Analysis/2e_events/2pix2e_events_extractor.py salida hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run1**7200*t -H 0