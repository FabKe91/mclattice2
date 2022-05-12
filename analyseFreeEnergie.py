import numpy as np
import matplotlib as ma
ma.use("Agg")
import matplotlib.pylab as plt
plt.rc('text', usetex=True)


data=np.zeros((5,151))
counter=np.zeros((5,151))


makeIm=400
with open("freeEnergieOut.txt") as f:
    for line in f:
        line=[float(x) for x in line.split()]
        if line[0]==makeIm:
            print(makeIm)
            makeIm+=100
            meanData=data/counter
            #print(data)
            #print(counter)
            fig=plt.figure(figsize=(4.980614173228346,3.2))
            ax=plt.subplot(111)

            for i in range(5):
                ax.plot(np.arange(-0.5,1.01,0.01),meanData[i]-meanData[i,20],".-",markersize=2,label=str(i))
                #ax.plot(np.arange(0,151,1),meanData[i]-meanData[i,0],"-",label=str(i))
            plt.legend()
            
            ax.set_xlabel(r'$S_{\textrm{\tiny CD}$')
            ax.set_ylabel(r'$G$ [kJ/mol]')
            
            plt.savefig("freeEnergie.png",bbox_inches="tight",dpi=300)
            plt.close(fig)
            fig=None
            #exit()
            
        if line[0]>300:
            data[line[2],line[1]]+=line[3]+300
            counter[line[2],line[1]]+=1

#meanData=data/counter
#print(counter)
#print(data)

#for i in range(5):
    #plt.plot(np.arange(0,151,1),data[i],label=str(i))

#plt.savefig("freeEnergie.png")
#plt.show()
