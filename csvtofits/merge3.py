from astropy.io import fits
import numpy as np
import csv
import sys

def csvtofits(data,filename,electrons):
    results = []
    with open(filename) as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
        for row in reader: # each row is a list
            results.append(row)
    for i,j in results:
        print("x = "+str(i))
        print("y = "+str(j))
        data[int(i)][int(j)+443]+=electrons


data=np.zeros([3072,443*2+1])

csvtofits(data,"./Mukul/oneelectronpixels"+str(sys.argv[1])+".csv",1)
csvtofits(data,"./Mukul/twoelectronpixels"+str(sys.argv[1])+".csv",2)
csvtofits(data,"./Mukul/trackpixels"+str(sys.argv[1])+".csv",1000)

# data=np.transpose(data)
cards=['NSAMP','SSAMP','NROW','NCOL','CCDNROW','CCDNCOL','CCDNPRES','NBINCOL','NBINROW','RUNID']
cardsval=[1,1,3072,443*2+1,3072*2,(443*2+1)*2,0,1,1,sys.argv[1]]

hduinit = fits.PrimaryHDU(data)
hduinitl = fits.HDUList([hduinit])
for i in range(0,len(cards)):
    hduinit.header[cards[i]] = cardsval[i]
hduinitl.writeto('./Mukul/image_1_'+str(sys.argv[1])+'.fits')
