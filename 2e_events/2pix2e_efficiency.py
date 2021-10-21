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
from scipy.stats import mvn
from scipy.stats import multivariate_normal as MVN


gStyle.SetStatW(0.25)
gStyle.SetStatH(0.25)
gROOT.SetBatch(True)
gStyle.SetOptStat(100011)
gStyle.SetOptFit(10001)
gStyle.SetPalette(57)

options, remainder = getopt.gnu_getopt(sys.argv[1:], 'N:h')

for opt, arg in options:
    if opt == '-N': # number of electrons
        N = int(arg)
    elif opt == '-h':
        print "\nUsage: "+sys.argv[0]+" <output basename> <root files>"
        print "Arguments: "
        print "\t-H: set 1ehalo"
        sys.exit(0)

#some alarm help
if (len(remainder)<2):
    print sys.argv[0]+' <output basename> <root files>'
    sys.exit()


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
            start = datetime.datetime.strptime(getattr(header,"DATESTART").split(b'\0',1)[0],"%Y-%m-%dT%H:%M:%S")
            end = datetime.datetime.strptime(getattr(header,"DATEEND").split(b'\0',1)[0],"%Y-%m-%dT%H:%M:%S")
            readoutTime = end-start
            #get exposure time
            try:
                exposure = int(getattr(header, "EXPOSURE"))
            except:
                exposure = 0
            readoutDays = readoutTime.days + (float(readoutTime.seconds) - exposure)/(24*60*60)
        files[int(getattr(header, "LTANAME"))-1].append(infiles[i])
        print("Adding "+str(infiles[i])+" to module = "+str(int(getattr(header, "LTANAME"))))

    return ccdPrescan,ccdCol,ccdRow,nRow,nCol,binRows,binCols,files,readoutDays,exposure

def count(clustercut,electroncut,threshold,threshold2=None):
    aducut=' && '.join(['ePix['+str(i)+']>'+str(thr) for i,thr in enumerate([threshold]*clustercut)])
    if threshold2!=None:
        aducut2=' || '.join(['(ePix['+str(i)+']>'+str(thr)+')' for i,thr in enumerate([threshold2]*clustercut)])
        entries=hitSumm2e.GetEntries("e=={0} && n=={1}".format(electroncut,clustercut)+" && "+aducut+" && "+aducut2)
        print("Cuts: e=={0} && n=={1}".format(electroncut,clustercut)+" && "+aducut+" && "+aducut2)
    else:
        entries=hitSumm2e.GetEntries("e=={0} && n=={1}".format(electroncut,clustercut)+" && "+aducut)
        print("Cuts: e=={0} && n=={1}".format(electroncut,clustercut)+" && "+aducut)
    configuration="{0}pix{1}e".format(clustercut,electroncut)
    print("Found {0} events for {1} using a threshold of {2:.2} electrons".format(entries,configuration,threshold))
    
    
    


# #default dimensions - for prototype CCD
# ccdPrescan = 7
# ccdCol = 886
# ccdRow = 6144
# nRow = 700
# nCol = 450
# binRows=1
# binCols=1

#organize input files.
outfilename = remainder[0]
infiles = remainder[1:]

# #read config
# ccdPrescan,ccdCol,ccdRow,nRow,nCol,binRows,binCols,files,readoutDays,exposure = readconfig(infiles,ccdPrescan,ccdCol,ccdRow,nRow,nCol,binRows,binCols,files)

#open file
file = TFile(infiles[0])
#open hitSumm2e
hitSumm2e = file.Get("hitSumm2e")
# #open config
# config = file.Get("config")
# config.GetEntry(0)
# entries = config.GetLeaf("entries").GetValue(0)+config.GetLeaf("entries").GetValue(1)
# entries = config.GetLeaf("entries").GetValue(0)
# mu = np.mean([config.GetLeaf("1erate").GetValue(0),config.GetLeaf("1erate").GetValue(1)])
# mu=4e-3
# print("mu: "+str(mu))
# noise = np.mean([config.GetLeaf("noise").GetValue(0),config.GetLeaf("noise").GetValue(1)])
# noise=.164
entries=1160*11600*100.0*10
noise=.164
mu=4e-4
print("noise: "+str(noise))
#1% contamination
#create gaussian TF1
fGaussian_cdf = TF1("f","ROOT::Math::normal_cdf_c(x,[0],[1])",-1,5) #[0] is noise, [1] is mean
fGaussian_cdf.SetParameter(0,noise)

#thresholds
thr1e=0.7
thr2e=1.7
thr3e=2.7

print("------------------------")
####################################################################################################
####################################################################################################
####################################################################################################
count(1,1,thr1e)
#1pix2e
thresholds=[thr2e,thr3e]
count(1,2,thresholds[0])
# efficiency

fGaussian_cdf.SetParameter(1,2)
efficiency=fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1])
print("Efficiency = {0:.2f}%".format(efficiency*100))
# contamination (misclassification)
contaminationvector=[]
#first term of contamination
fGaussian_cdf.SetParameter(1,0)
contaminationvector.append(fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1]))
#second term of contamination
fGaussian_cdf.SetParameter(1,1)
contaminationvector.append((fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1]))*mu)
fGaussian_cdf.SetParameter(1,3)
contaminationvector.append((fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1]))*mu*mu*mu/6)
# fGaussian_cdf.SetParameter(1,2)
# contaminationvector.append(-(fGaussian_cdf.Eval(-555)-fGaussian_cdf.Eval(thr2e))*mu*mu/2)
# fGaussian_cdf.SetParameter(1,2)
# contaminationvector.append(-(fGaussian_cdf.Eval(thr3e)-fGaussian_cdf.Eval(555))*mu*mu/2)
false_events=np.sum(contaminationvector)
true_events=(mu*mu/2)*efficiency
print("Contamination = {0:.2f}%".format(
    100*(false_events/(true_events+false_events))
    ))
print("------------------------")
####################################################################################################
####################################################################################################
####################################################################################################
#2pix2e
thresholds=[thr1e,thr2e]
count(2,2,thresholds[0])
# efficiency
center = np.array([1,1])
covariance = np.array([[noise, 0],[0, noise]])
dist = MVN(mean=center, cov=covariance)
efficiency = dist.cdf(np.array([thresholds[1]]*2))-dist.cdf(np.array([thresholds[0]]*2))
# center = np.array([1.0, 1.0])
# S = np.array([[noise]*2,[noise]*2])
# low = np.array([thresholds[0]]*2)
# upp = np.array([thresholds[1]]*2)
# p,i = mvn.mvnun(low,upp,center,S)
print("Efficiency = {0:.2f}%".format(efficiency*100))
# contamination (misclassification)
contaminationvector=[]
# first term of contamination. 0- events that are classified as 1e- and are next to a 1e- event
fGaussian_cdf.SetParameter(1,0)
contaminationvector.append(mu*8*(fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1])))
# second term of contamination. 0- events that are classified as 1e- and are next to a 0- events that is classified as 1e-
center = np.array([0,0])
covariance = np.array([[noise, 0],[0, noise]])
dist = MVN(mean=center, cov=covariance)
contaminationvector.append((mu*mu/2)*8*(dist.cdf(np.array([thresholds[1]]*2))-dist.cdf(np.array([thresholds[0]]*2))))
#third term of contamination. 2- events that are classified as 1e- and are next to a 1e- event
fGaussian_cdf.SetParameter(1,2)
contaminationvector.append((mu**3)/2*8*(fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1])))
# #fourth term of contamination. 1- events that are classified as 0e- and are next to a 1e- event
# fGaussian_cdf.SetParameter(1,1)
# contaminationvector.append(-(mu**2)*2*8*(fGaussian_cdf.Eval(-555)-fGaussian_cdf.Eval(thresholds[0])))
# #fifth term of contamination. 1- events that are classified as 2e- and are next to a 1e- event
# fGaussian_cdf.SetParameter(1,1)
# contaminationvector.append(-(mu**2)*2*8*(fGaussian_cdf.Eval(thr2e)-fGaussian_cdf.Eval(thr3e)))
false_events=np.sum(contaminationvector)
true_events=8*(mu*mu/2)*efficiency
print("Contamination = {0:.2f}%".format(100*false_events/(true_events+false_events)))
print("------------------------")
####################################################################################################
####################################################################################################
####################################################################################################
#1pix3e
count(1,3,thr3e) #note above 2 thresholds are at .5
# efficiency
thresholds=[thr3e,3.7]
fGaussian_cdf.SetParameter(1,3)
efficiency=fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1])
print("Efficiency = {0:.2f}%".format(efficiency*100))
# contamination (misclassification)
contaminationvector=[]
#first term of contamination
fGaussian_cdf.SetParameter(1,0) # 0e- events that are classified as 3e-
contaminationvector.append(fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1]))
#second term of contamination
fGaussian_cdf.SetParameter(1,1) # 1e- events that are classified as 3e-
contaminationvector.append((fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1]))*mu)
#third term of contamination
fGaussian_cdf.SetParameter(1,2) # 2e- events that are classified as 3e-
contaminationvector.append((fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1]))*(mu*mu/2))
# #fourth term of contamination
# fGaussian_cdf.SetParameter(1,3) # 3e- events that are classified as 0,1 or 2e-
# contaminationvector.append(-(fGaussian_cdf.Eval(-555)-fGaussian_cdf.Eval(thresholds[0]))*(mu*mu*mu/8))
# #fifth term of contamination
# fGaussian_cdf.SetParameter(1,3) # 3e- events that are classified as 4 or more e-
# contaminationvector.append(-(fGaussian_cdf.Eval(thresholds[1])-fGaussian_cdf.Eval(555))*(mu*mu*mu/8))
# print(contaminationvector)
false_events=np.sum(contaminationvector)
# print(false_events)
true_events=(mu*mu*mu/6)*efficiency
# print(true_events)
print("Contamination = {0:.2f}%".format(100*false_events/(true_events+false_events)))
print("------------------------")
####################################################################################################
####################################################################################################
####################################################################################################
#2pix3e
thresholds=[[thr1e,thr2e],[thr2e,thr3e]]
count(2,3,thresholds[0][0],thresholds[0][1]) #note above 2 thresholds are at .5
# efficiency
center = np.array([1,2])
covariance = np.array([[noise, 0],[0, noise]])
dist = MVN(mean=center, cov=covariance)
efficiency = dist.cdf(np.array(thresholds[1]))-dist.cdf(np.array(thresholds[0]))
print("Efficiency = {0:.2f}%".format(efficiency*100))
# contamination (misclassification)
contaminationvector=[]
#first term of contamination
fGaussian_cdf.SetParameter(1,0)
contaminationvector.append((fGaussian_cdf.Eval(thr1e)-fGaussian_cdf.Eval(thr2e))*(mu*mu/2))
#second term of contamination
fGaussian_cdf.SetParameter(1,0)
contaminationvector.append(fGaussian_cdf.Eval(thr2e)-fGaussian_cdf.Eval(thr3e)*mu)
#third term of contamination
fGaussian_cdf.SetParameter(1,1)
contaminationvector.append((fGaussian_cdf.Eval(thr2e)-fGaussian_cdf.Eval(thr3e))*(mu*mu/2))
false_events=np.sum(contaminationvector)
# #fourth term of contamination
# fGaussian_cdf.SetParameter(1,2)
# contaminationvector.append(-(fGaussian_cdf.Eval(-555)-fGaussian_cdf.Eval(thr2e))*(mu*mu*mu/2))
# #fifth term of contamination
# fGaussian_cdf.SetParameter(1,2)
# contaminationvector.append(-(fGaussian_cdf.Eval(thr3e)-fGaussian_cdf.Eval(555))*(mu*mu*mu/2))
# #sixth term of contamination
# fGaussian_cdf.SetParameter(1,1)
# contaminationvector.append(-(fGaussian_cdf.Eval(-555)-fGaussian_cdf.Eval(thr1e))*(mu*mu*mu/2))
# #septh term of contamination
# fGaussian_cdf.SetParameter(1,1)
# contaminationvector.append(-(fGaussian_cdf.Eval(thr1e)-fGaussian_cdf.Eval(thr2e))*(mu*mu*mu/2))
false_events=np.sum(contaminationvector)*8
true_events=(mu*mu*mu/2)*8*efficiency
print("Contamination = {0:.2f}%".format(100*false_events/(true_events+false_events)))
# print("Contamination = {0:.2f}%".format(100*entries*mu*8*np.sum(contaminationvector)))
print("------------------------")
####################################################################################################
####################################################################################################
####################################################################################################
#3pix3e
count(3,3,thr1e)
thresholds=[thr1e,thr2e]
# efficiency
center = np.array([1,1,1])
covariance = np.array([[noise, 0, 0],[0, noise, 0],[0, 0, noise]])
dist = MVN(mean=center, cov=covariance)
p = dist.cdf(np.array([thresholds[1]]*3))-dist.cdf(np.array([thresholds[0]]*3))
print("Efficiency = {0:.2f}%".format(p*100))
# contamination (misclassification)
contaminationvector=[]
fGaussian_cdf.SetParameter(1,0)
contaminationvector.append(mu*mu*44*(fGaussian_cdf.Eval(thresholds[0])-fGaussian_cdf.Eval(thresholds[1])))
false_events=np.sum(contaminationvector)
true_events=(mu*mu*mu/6)*120*efficiency
print("Contamination = {0:.2f}%".format(100*false_events/(true_events+false_events)))
print("------------------------")
####################################################################################################
####################################################################################################
####################################################################################################
#1pix4e
count(1,4,4.5)