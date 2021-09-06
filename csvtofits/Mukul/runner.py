#!/usr/bin/python2.7
import os,sys,glob,re

#os.system("for i in /home/mariano/tools/SKIPPERtools/csvtofits/Mukul/image*fits; do skipper2root -irBd $i -g1; done;")

procs=glob.glob1(os.popen("pwd").read()[0:-1],"proc*.fits") #procs availables
iterator=len(procs)

for i in range(0,iterator):
	print("~/tools/skextract/skExtract.exe -c extractConfig_1.xml -C cal.xml "+str(procs[i])+" -a")
	os.system("~/tools/skextract/skExtract.exe -c extractConfig_1.xml -C cal.xml "+str(procs[i])+" -a")


