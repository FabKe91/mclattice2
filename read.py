#!/usr/bin/python3
'''
class data
    .getStep   <- dataSetName, step -> returns the dataSet of name <dataSetName> of <step>
    .getOrderDistr <- -> Normalized orderparameter distribution for frames 20 to end
    .getSingleOrderDistr <- imageNumber -> Normalized order parameter distribution of frame <imageNumber> 
    .getMeanNN <- imageNumber, Type -> Average number of like-neighbor in frame <imageNumber> for type <Type>

orderDistrAni <- instance of data -> creates an animation axs[0]=orderparadistr axs[1]=typedistr
snap

'''
import h5py
import numpy as np
import pandas as pd 
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.animation as animation
from matplotlib.colors import ListedColormap

############


cm = 1/2.54
a169 = 1/(16/9)
singlex=16*cm



DATAFILE = "out.h5"


###############

def set_context():
    sns.set_context("poster")
    #sns.set_context("paper", font_scale=1.2, rc={"line.linewidth":1, **rc})

set_context()


class data(object):

    def __init__(self,filename):
        self.f = h5py.File(filename,  "r")
        
        self.images=self.f["/orderPara"].shape[0]
        self.paras={}
        
        for item in self.f.attrs.keys():
            self.paras[item]=self.f.attrs[item][0]

        print("Initialized parameters:\n", *["{}:{}\n".format(str(i), str(j)) for i,j in self.paras.items()])
        
        
    def getStep(self,dataSetName,step):
        return np.array(self.f["/"+dataSetName][step,:,:])

    def getOrderDistr(self):
        counter=np.zeros( int(self.paras["maxOrderIndex"]+1) )
        for image in range(20,self.images):
            for order in np.nditer(self.getStep("orderPara",image)):
                counter[order]+=1
        return counter/np.sum(counter)/self.paras["DeltaOrder"]

    def getSingleOrderDistr(self,imageNumber):
            counter=np.zeros( int(self.paras["maxOrderIndex"]+1) )
            for order in np.nditer(self.getStep("orderPara",imageNumber)):
                counter[order]+=1
            return counter/np.sum(counter)/self.paras["DeltaOrder"]
    
    def getMeanNN(self,imageNumber,Type):
        '''
            Looks for pairs in each direction for numpy roll: axis=0/1 and step 1/-1 
            All equal pairs are summed in each direction

        '''
        step = np.array( self.f["/Type"][imageNumber,:,:] ) == Type
        sameN = 0
        sameN += np.sum(step[np.roll(step,1,axis=0)])
        sameN += np.sum(step[np.roll(step,-1,axis=0)])
        sameN += np.sum(step[np.roll(step,1,axis=1)])
        sameN += np.sum(step[np.roll(step,-1,axis=1)])
        return sameN/np.sum(step)        
    

    
def orderDistrAni(data):
    fig=plt.figure()
    ax=plt.subplot(111)
    plot,=ax.plot(np.arange(self.paras["minOrder"],self.paras["minOrder"]+self.paras["DeltaOrder"],self.paras["DeltaOrder"]),data.getSingleOrderDistr(0),"k-")
    ax.set_ylim([0,5])

    def update (i):
        plot.set_ydata(data.getSingleOrderDistr(i))
        print(i)

    ani=animation.FuncAnimation(fig, update, blit=False,frames=300, interval=100, repeat_delay=600)
    mywriter = animation.FFMpegWriter(fps=10)
    ani.save('orderDistrAni.avi',writer=mywriter,dpi=400)
    #plt.show()


    
def orderDistrAndType(data):
    fig=plt.figure(figsize=(singlex*2,singlex*a169))
    
    ax2=plt.subplot(121)
    empty=np.zeros( ( int(data.paras["width"]), int(data.paras["height"]) ) )
    im=ax2.imshow(empty,vmax=1,vmin=-0.3,interpolation='None',cmap='gnuplot')
    plt.colorbar(im,ax=ax2)
    ax2.set_title('axes title')

    ax3=plt.subplot(122)
    empty=np.zeros( ( int(data.paras["width"]), int(data.paras["height"]) ) )
    my_cmap=ListedColormap([(1.0, 0.6, 0.6),(0.1, 0.7, 0.8)])
    im2=ax3.imshow(empty,vmax=1,vmin=0,interpolation='None',cmap=my_cmap)
    plt.colorbar(im2,ax=ax3)
    
    def update (i):
        j=i*100
        ax2.set_title('%s'%j)
        im.set_array((data.getStep("orderPara",j)-50)*0.01)
        im2.set_array(data.getStep("Type",j))
        print("step:",i,"sameN DPPC",data.getMeanNN(j,0))

    ani=animation.FuncAnimation(fig, update, blit=False,frames=data.images//100, interval=20,repeat_delay=600)
    mywriter = animation.FFMpegWriter(fps=1)
    print("Saving animation...")
    ani.save('Ani.avi', writer=mywriter, dpi=150)
    #plt.show()
    
def snap(data,imageNumber):
    fig, axs = plt.subplots(ncols=2, figsize=(singlex*2,singlex) )

    axs[0].set_axis_off()
    im=axs[0].imshow((data.getStep("orderPara",imageNumber)-50)*0.01,vmax=1,vmin=-0.3,interpolation='None',cmap='gnuplot')
    plt.colorbar(im,ax=axs[0])
    
    axs[1].set_axis_off()
    my_cmap=ListedColormap([(1.0, 0.6, 0.6),(0.1, 0.7, 0.8)])
    im2=axs[1].imshow(data.getStep("Type",imageNumber),vmax=1,vmin=0,interpolation='None',cmap=my_cmap)
    plt.colorbar(im2,ax=axs[1])
    axs[1].axis('off')



    plt.tight_layout()
    fig.savefig("snap_image%s.png"%(imageNumber), dpi=300)




data = data(DATAFILE)

temp = int(data.paras["T"])

print("Saving snapshot of last frame...")
snap(data, -1)

print("Animate simulation...")
orderDistrAndType(data)

svals = np.arange(
    data.paras["minOrder"], data.paras["maxOrder"]+0.00001, data.paras["DeltaOrder"]
	)
print("Calculating order distribution...")
orderdistr=data.getOrderDistr()
dat = pd.DataFrame({"order": svals, "freq": orderdistr})
print("Avg Order:",  (dat["order"] * dat["freq"]).sum() / dat.freq.sum())
dat.to_csv("orderdistr.csv", index=False)






#plt.plot(np.arange(-0.5,1.01,0.01),orderdistr,"k-")
###plt.show()




############### Plot system evolution  ###################
#fig, axs = plt.subplots(ncols=2, figsize=(singlex,singlex*a43) )
#
#ax = sns.lineplot(x="step", y="order", data=dt_s,  ax=axs[0])
#ax = sns.lineplot(x="step", y="NN",    data=dt_NN, ax=axs[1])
#
#fig.savefig("system_evolution_T{temp}.png")








############## Plot meanNN evolution ###################

