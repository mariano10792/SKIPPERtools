import numpy as np
from scipy import stats
import csv
import matplotlib.pyplot as plt

with open('counts.txt','r') as fd:
    reader = csv.reader(fd)
    for row in reader:
        row=row[0:-1]

row=[float(i) for i in row]



rowSim=np.sort(np.random.poisson(np.mean(row),len(row)))
#rowSim2=np.sort(np.random.poisson(np.mean(row),len(row)))
rowSimCum=[np.sum(rowSim[0:i])/float(np.sum(rowSim)) for i in range(1,len(rowSim)+1)]
#rowSimCum2=[np.sum(rowSim2[0:i])/float(np.sum(rowSim2)) for i in range(1,len(rowSim2)+1)]
#kolmo2=stats.ks_2samp(rowSim2,rowSim)
#print(kolmo2)

#plt.plot(rowSimCum,rowSim,'.b')
#plt.plot(rowSimCum2,rowSim2,'.r')

#plt.show()

kolmo=stats.ks_2samp(row,rowSim)
print(kolmo)
rowCum=[np.sum(row[0:i])/np.sum(row) for i in range(1,len(row)+1)]
#rowSimCum=[np.sum(rowSim[0:i])/float(np.sum(rowSim)) for i in range(1,len(rowSim)+1)]

plt.plot(row,rowCum,'.r')
plt.plot(rowSim,rowSimCum,'.b')
#plt.savefig('counts100.pdf')  

plt.show()


