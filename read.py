#!/usr/bin/python3

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
from matplotlib.colors import ListedColormap
import matplotlib.pylab as plt
import matplotlib.animation as animation
import h5py
from sys import argv



class data(object):
    def __init__(self,filename):
        self.f = h5py.File(filename,  "r")
        
        self.images=self.f["/orderPara"].shape[0]
        self.paras={}
        
        for item in self.f.attrs.keys():
            self.paras[item]=self.f.attrs[item][0]
        
        
    def getStep(self,dataSetName,step):
        return np.array(self.f["/"+dataSetName][step,:,:])

    def getOrderDestr(self):
        counter=np.zeros((self.paras["maxOrderIndex"]+1))
        for image in range(20,self.images):
            for order in np.nditer(self.getStep("orderPara",image)):
                counter[order]+=1
        return counter/np.sum(counter)/self.paras["DeltaOrder"]

    def getSingleOrderDestr(self,imageNumber):
            counter=np.zeros((self.paras["maxOrderIndex"]+1))
            for order in np.nditer(self.getStep("orderPara",imageNumber)):
                counter[order]+=1
            return counter/np.sum(counter)/self.paras["DeltaOrder"]
    
    def getMeanNN(self,imageNumber,Type):
        step=np.array(self.f["/Type"][imageNumber,:,:])==Type
        sameN=0
        sameN+=np.sum(step[np.roll(step,1,axis=0)])
        sameN+=np.sum(step[np.roll(step,-1,axis=0)])
        sameN+=np.sum(step[np.roll(step,1,axis=1)])
        sameN+=np.sum(step[np.roll(step,-1,axis=1)])
        return sameN/np.sum(step)        
    

    
def orderDestrAni(data):
    fig=plt.figure()
    ax=plt.subplot(111)
    plot,=ax.plot(np.arange(self.paras["minOrder"],self.paras["minOrder"]+self.paras["DeltaOrder"],self.paras["DeltaOrder"]),data.getSingleOrderDestr(0),"k-")
    ax.set_ylim([0,5])

    def update (i):
        plot.set_ydata(data.getSingleOrderDestr(i))
        print(i)

    ani=animation.FuncAnimation(fig, update, blit=False,frames=300, interval=100, repeat_delay=600)
    #mywriter = animation.FFMpegWriter(fps=10)
    #ani.save('orderDestrAni.avi',writer=mywriter,dpi=400)
    plt.show()


    
def orderDestrAndType(data):
    fig=plt.figure(figsize=(12,4))
    
    ax2=plt.subplot(121)
    empty=np.zeros((data.paras["width"],data.paras["height"]))
    im=ax2.imshow(empty,vmax=1,vmin=-0.3,interpolation='None',cmap='gnuplot')
    plt.colorbar(im,ax=ax2)
    ax2.set_title('axes title')

    ax3=plt.subplot(122)
    empty=np.zeros((data.paras["width"],data.paras["height"]))
    my_cmap=ListedColormap([(1.0, 0.6, 0.6),(0.1, 0.7, 0.8)])
    im2=ax3.imshow(empty,vmax=1,vmin=0,interpolation='None',cmap=my_cmap)
    plt.colorbar(im2,ax=ax3)
    
    def update (i):
        j=i*100
        ax2.set_title('%s'%j)
        im.set_array((data.getStep("orderPara",j)-50)*0.01)
        im2.set_array(data.getStep("Type",j))
        print("step:",i,"sameN DPPC",data.getMeanNN(j,0))

    ani=animation.FuncAnimation(fig, update, blit=False,frames=data.images, interval=20,repeat_delay=600)
    #mywriter = animation.FFMpegWriter(fps=10)
    #ani.save('Ani.avi',writer=mywriter,dpi=150)
    plt.show()
    
def snap(data,imageNumber):
    #fig=plt.figure(frameon=False)
    #fig.set_size_inches(12,12)
    #ax1 = plt.Axes(fig, [0., 0., 1., 1.])
    #ax1.set_axis_off()
    #fig.add_axes(ax1)
    #im=ax1.imshow((data.getStep("orderPara",imageNumber)-50)*0.01,vmax=1,vmin=-0.3,interpolation='None',cmap='gnuplot')
    #plt.colorbar(im,ax=ax1)
    #plt.savefig("snap_order_image%s.png"%imageNumber,bbox_inches="tight",dpi=300,pad_inches=0)
    #fig=None
    
    fig=plt.figure(frameon=False)
    fig.set_size_inches(12,12)
    ax1 = plt.Axes(fig, [0., 0., 1., 1.])
    ax1.set_axis_off()
    fig.add_axes(ax1)
    my_cmap=ListedColormap([(1.0, 0.6, 0.6),(0.1, 0.7, 0.8)])
    im2=ax1.imshow(data.getStep("Type",imageNumber),vmax=1,vmin=0,interpolation='None',cmap=my_cmap)
    ax1.axis('off')
    plt.savefig("n%s_snap_type_image%s.png"%(data.shape[0],imageNumber),bbox_inches="tight",dpi=300,pad_inches=0)
    fig=None
    



data1=data("../MCLattice2/out.h5")
#data1=data("out.h5")



orderDestrAndType(data1)


#snap(data1,int(argv[1]))
#orderdestr=data1.getOrderDestr()
#np.save("orderdestr.npy",orderdestr)

#plt.plot(np.arange(-0.5,1.01,0.01),orderdestr,"k-")
###plt.show()

#print(argv[1],np.sum(orderdestr*np.arange(-0.5,1.01,0.01)*0.01))
#plt.savefig("T_%s.png"%argv[1],bbox_inches="tight",dpi=400)


