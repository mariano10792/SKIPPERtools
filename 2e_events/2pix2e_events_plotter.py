#!/usr/bin/env python
import sys, getopt, os
import array
import re
from ROOT import gROOT, gStyle, TFile, TTree, TChain, TMVA, TCut, TCanvas, gDirectory, TH1, TH1F, TGraph, gPad, TF1, THStack, TLegend, TGraphErrors, TLatex, TList, TEfficiency
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import datetime
import pandas as pd
from openpyxl import load_workbook
import path
from scipy.stats import chi2

gStyle.SetOptStat(100011)
gStyle.SetOptFit(10001)

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'H:h')

for opt, arg in options:
    if opt == '-H':
        halo1eOpt = int(arg)
    elif opt == '-h':
        print "\nUsage: "+sys.argv[0]+" <output basename> <root files>"
        print "Arguments: "
        print "\t-H: set 1ehalo"
        sys.exit(0)


outfilename = remainder[0]
infiles = remainder[1:]


thefile = TFile(infiles[0])
hitSumm2e = thefile.Get("hitSumm2e")
calPixTree1e = thefile.Get("calPixTree1e")

nmods=2
eff=[]; 
hists=[]
outfile = TFile(outfilename+".root","RECREATE");
c = TCanvas("c","c",1200,900); c.Clear();
outfile.cd()



#canvas and rootfile creation
c.Print("plots1.pdf[")



##################################
############## xBary #############
##################################
binning="(6,0,3200)"

for I in range(0,nmods):
    #1pix1e- events
    c.Divide(2,2)
    for hdu in [1,2,3,4]:
        hname='h_1pix1e_xBary_ohdu'+str(hdu)
        c.cd(hdu)
        eff.append(TEfficiency(thefile.Get("h_x_1e_ohdu"+str(hdu)+"Module"+str(I+1)),thefile.Get("h_x_ohdu"+str(hdu)+"Module"+str(I+1))))
        eff[-1].Draw()
        eff[-1].SetTitle(hname+"Module"+str(I+1))
        eff[-1].SetName("eff"+hname+"LTANAME"+str(I+1))
        eff[-1].Write()
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

    hadd=eff[-4]; hadd.Add(eff[-3]); hadd.Add(eff[-2]); hadd.Add(eff[-1]); hadd.Draw();
    c.cd(); c.Print("plots1.pdf"); c.Clear();

for I in range(0,nmods):
    #1pix2e- events
    c.Divide(2,2)
    for hdu in [1,2,3,4]:
        hname='h_1pix2e_xBary_ohdu'+str(hdu)
        c.cd(hdu)
        hitSumm2e.Draw("xBary>>"+hname+"Module"+str(I+1)+binning,"n==1 && ohdu=="+str(hdu)+" && LTANAME=="+str(I+1))
        h = gDirectory.Get(hname+"Module"+str(I+1))
        hists.append(h)
        eff.append(TEfficiency(gDirectory.Get(hname+"Module"+str(I+1)),thefile.Get("h_x_ohdu"+str(hdu)+"Module"+str(I+1))))
        eff[-1].Draw()
        eff[-1].SetTitle(hname+"Module"+str(I+1))
        eff[-1].SetName("eff"+hname+"LTANAME"+str(I+1))
        eff[-1].Write()
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

    hadd=eff[-4]; hadd.Add(eff[-3]); hadd.Add(eff[-2]); hadd.Add(eff[-1]); hadd.Draw();
    c.cd(); c.Print("plots1.pdf"); c.Clear();

for I in range(0,nmods):
    #2pix2e- events
    c.Divide(2,2)
    for hdu in [1,2,3,4]:
        hname='h_2pix2e_xBary_ohdu'+str(hdu)
        c.cd(hdu)
        hitSumm2e.Draw("xBary>>"+hname+"Module"+str(I+1)+binning,"n==2 && ohdu=="+str(hdu)+" && LTANAME=="+str(I+1))
        h = gDirectory.Get(hname+"Module"+str(I+1))
        hists.append(h)
        eff.append(TEfficiency(gDirectory.Get(hname+"Module"+str(I+1)),thefile.Get("h_x_ohdu"+str(hdu)+"Module"+str(I+1))))
        eff[-1].Draw()
        eff[-1].SetTitle(hname+"Module"+str(I+1))
        eff[-1].SetName("eff"+hname+"LTANAME"+str(I+1))
        eff[-1].Write()
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

    hadd=eff[-4]; hadd.Add(eff[-3]); hadd.Add(eff[-2]); hadd.Add(eff[-1]); hadd.Draw();
    c.cd(); c.Print("plots1.pdf"); c.Clear();

##################################
############## yBary #############
##################################
binning="(4,0,520)"

for I in range(0,nmods):
    #1pix1e- events
    c.Divide(2,2)
    for hdu in [1,2,3,4]:
        hname='h_1pix1e_yBary_ohdu'+str(hdu)
        c.cd(hdu)
        eff.append(TEfficiency(thefile.Get("h_y_1e_ohdu"+str(hdu)+"Module"+str(I+1)),thefile.Get("h_y_ohdu"+str(hdu)+"Module"+str(I+1))))
        eff[-1].Draw()
        eff[-1].SetTitle(hname+"Module"+str(I+1))
        eff[-1].SetName("eff"+hname+"LTANAME"+str(I+1))
        eff[-1].Write()
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

    hadd=eff[-4]; hadd.Add(eff[-3]); hadd.Add(eff[-2]); hadd.Add(eff[-1]); hadd.Draw();
    gPad.Update(); 
    graph = hadd.GetPaintedGraph(); 
    graph.SetMinimum(0);
    graph.SetMaximum(1.2e-3); 
    gPad.Update(); 
    c.cd(); c.Print("plots1.pdf"); c.Clear();

for I in range(0,nmods):
    #1pix2e- events
    c.Divide(2,2)
    for hdu in [1,2,3,4]:
        hname='h_1pix2e_yBary_ohdu'+str(hdu)
        c.cd(hdu)
        hitSumm2e.Draw("yBary>>"+hname+"Module"+str(I+1)+binning,"n==1 && ohdu=="+str(hdu)+" && LTANAME=="+str(I+1))
        h = gDirectory.Get(hname+"Module"+str(I+1))
        hists.append(h)
        eff.append(TEfficiency(gDirectory.Get(hname+"Module"+str(I+1)),thefile.Get("h_y_ohdu"+str(hdu)+"Module"+str(I+1))))
        eff[-1].Draw()
        eff[-1].SetTitle(hname+"Module"+str(I+1))
        eff[-1].SetName("eff"+hname+"LTANAME"+str(I+1))
        eff[-1].Write()
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

    hadd=eff[-4]; hadd.Add(eff[-3]); hadd.Add(eff[-2]); hadd.Add(eff[-1]); hadd.Draw();
    c.cd(); c.Print("plots1.pdf"); c.Clear();

for I in range(0,nmods):
    #2pix2e- events
    c.Divide(2,2)
    for hdu in [1,2,3,4]:
        hname='h_2pix2e_yBary_ohdu'+str(hdu)
        c.cd(hdu)
        hitSumm2e.Draw("yBary>>"+hname+"Module"+str(I+1)+binning,"n==2 && ohdu=="+str(hdu)+" && LTANAME=="+str(I+1))
        h = gDirectory.Get(hname+"Module"+str(I+1))
        hists.append(h)
        eff.append(TEfficiency(gDirectory.Get(hname+"Module"+str(I+1)),thefile.Get("h_y_ohdu"+str(hdu)+"Module"+str(I+1))))
        eff[-1].Draw()
        eff[-1].SetTitle(hname+"Module"+str(I+1))
        eff[-1].SetName("eff"+hname+"LTANAME"+str(I+1))
        eff[-1].Write()
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

    hadd=eff[-4]; hadd.Add(eff[-3]); hadd.Add(eff[-2]); hadd.Add(eff[-1]); hadd.Draw();
    gPad.Update(); 
    graph = hadd.GetPaintedGraph(); 
    graph.SetMinimum(0);
    graph.SetMaximum(1e-5); 
    gPad.Update(); 
    c.cd(); c.Print("plots1.pdf"); c.Clear();


##################################
########### xBary:yBary ##########
##################################
# binning="(4,0,520,6,0,3200)"
# for I in range(0,nmods):
    #1pix2e- events
    # hname='h_1pix2e_xBary_yBary'
    # hitSumm2e.Draw("xBary:yBary>>"+hname+"Module"+str(I+1)+binning,"n==1 && LTANAME=="+str(I+1))
    # h = gDirectory.Get(hname+"Module"+str(I+1))
    # hists.append(h)
    # eff.append(TEfficiency(gDirectory.Get(hname+"Module"+str(I+1)),thefile.Get("h_x_y"+"Module"+str(I+1))))
    # eff[-1].SetName("eff"+hname+"LTANAME"+str(I+1))
    # eff[-1].Write()

# for I in range(0,nmods):
    #2pix2e- events
    # hname='h_2pix2e_xBary_yBary'
    # hitSumm2e.Draw("xBary:yBary>>"+hname+"Module"+str(I+1)+binning,"n==2 && LTANAME=="+str(I+1))
    # h = gDirectory.Get(hname+"Module"+str(I+1))
    # hists.append(h)
    # eff.append(TEfficiency(gDirectory.Get(hname+"Module"+str(I+1)),thefile.Get("h_x_y"+"Module"+str(I+1))))
    # eff[-1].SetName("eff"+hname+"LTANAME"+str(I+1))
    # eff[-1].Write()
    # c.cd()
    # c.Print("plots1.pdf")
    # c.Clear()

##################################
############# 1pix2e #############
##################################
binning="(4,1,5)"

for I in range(0,nmods):
    hname='h_1pix2e_ohdu'
    hitSumm2e.Draw("ohdu>>"+hname+"Module"+str(I+1)+binning,"n==1 && LTANAME=="+str(I+1))
    h = gDirectory.Get(hname+"Module"+str(I+1))
    hists.append(h)
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

##################################
############# 2pix2e #############
##################################

for I in range(0,nmods):
    hname='h_2pix2e_ohdu'
    hitSumm2e.Draw("ohdu>>"+hname+"Module"+str(I+1)+binning,"n==2 && LTANAME=="+str(I+1))
    h = gDirectory.Get(hname+"Module"+str(I+1))
    hists.append(h)
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

##################################
###### 2pix2e (horizontal) #######
##################################

for I in range(0,nmods):
    hname='h_2pix2e_ohdu'
    hitSumm2e.Draw("ohdu>>"+hname+"Module"+str(I+1)+binning,"n==2 && yVar==0 && LTANAME=="+str(I+1))
    h = gDirectory.Get(hname+"Module"+str(I+1))
    hists.append(h)
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

##################################
######## 2pix2e (vertical) #######
##################################

for I in range(0,nmods):
    hname='h_2pix2e_ohdu'
    hitSumm2e.Draw("ohdu>>"+hname+"Module"+str(I+1)+binning,"n==2 && xVar==0 && LTANAME=="+str(I+1))
    h = gDirectory.Get(hname+"Module"+str(I+1))
    hists.append(h)
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()

##################################
######## 2pix2e (diagonal) #######
##################################

for I in range(0,nmods):
    hname='h_2pix2e_ohdu'
    hitSumm2e.Draw("ohdu>>"+hname+"Module"+str(I+1)+binning,"n==2 && (xVar!=0 && yVar!=0) && LTANAME=="+str(I+1))
    h = gDirectory.Get(hname+"Module"+str(I+1))
    hists.append(h)
    c.cd()
    c.Print("plots1.pdf")
    c.Clear()



c.cd()
c.Print("plots1.pdf]")

#range
pixmax = 0.005
pixmin = 0
pixnbin = 400
pixbinning = "({0},{1},{2})".format(pixnbin,pixmin,pixmax)
binwidth = (pixmax-pixmin)/pixnbin

#hdus
hdus=[1,2,3,4]

#canvas and rootfile creation
c.Print("plots2.pdf[")

#funtion
fitfunc = TF1("fitfunc","[0]*(TMath::Gaus(x,[1],[2],1))") # [0]:area / [1]:mean / [2]:sigma 

for j in [1,2]:

    plotfunc = [TF1("fitfunc","[0]*(TMath::Gaus(x,[1],[2],1))",pixmin,pixmax) for i in hdus] # [0]:area / [1]:mean / [2]:sigma 

    for i, hdu in enumerate(hdus):

        calPixTree1e.Draw("mu>>hMod"+str(j)+"hdu"+str(hdu)+pixbinning,"LTANAME=="+str(j)+" && ohdu=="+str(hdu))
        h = gDirectory.Get("hMod"+str(j)+"hdu"+str(hdu))
        hists.append(h)
        fitfunc.SetParameters(h.GetEntries()*binwidth,h.GetMean(),h.GetStdDev())
        s = h.Fit(fitfunc,"QSL","",pixmin,pixmax)

        plotfunc[i].SetParameters(s.Parameter(0),s.Parameter(1),s.Parameter(2))
        plotfunc[i].SetNpx(10000)
        c.cd()
        c.Clear()


    legend = TLegend(0.60,0.25,0.85,0.40); 
    for i, hdu in enumerate(hdus):
        plotfunc[i].SetLineColor(i+1)
        legend.AddEntry(plotfunc[i],"Mod"+str(j)+"Hdu"+str(hdu),"f");
        if i==0:
            plotfunc[i].Draw("")
            plotfunc[i].GetYaxis().SetRangeUser(0,5)
            plotfunc[i].GetXaxis().SetTitle("mu");
            plotfunc[i].SetTitle("1e- events. Module"+str(j))
            plotfunc[i].Draw("")
        else:
            plotfunc[i].Draw("SAME")
        legend.Draw();

    c.cd()
    c.Print("plots2.pdf")
    c.Clear()

c.cd()
c.Print("plots2.pdf]")

outfile.Write()
outfile.Close()

