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
            if self.f.attrs[item]==int(self.f.attrs[item]):
                self.paras[item]=int(self.f.attrs[item])
            else:
                self.paras[item]=self.f.attrs[item]

        
    def getStep(self,dataSetName,step):
        return np.array(self.f["/"+dataSetName][step,:,:])

    def getOrderDistr(self):
        counter=np.zeros((self.paras["maxOrderIndex"]+1))
        for image in range(self.findEqui(),self.images):
            for order in np.nditer(self.getStep("orderPara",image)):
                counter[order]+=1
        return counter/np.sum(counter)/self.paras["DeltaOrder"]

    def getSingleOrderDistr(self,imageNumber):
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
    
    def getCholMeanNN(self,imageNumber):
        step=np.array(self.f["/Chol"][imageNumber,:,:],dtype=bool)
        hit=0
        hit+=np.sum(step[np.roll(step,1,axis=0)])
        hit+=np.sum(step[np.roll(step,-1,axis=0)])
        hit+=np.sum(step[np.roll(step,1,axis=1)])
        hit+=np.sum(step[np.roll(step,-1,axis=1)])
        return hit/np.sum(step)**2/4*self.paras["width"]*self.paras["height"]
    
    def findEqui(self):
        lastMeanOrder=np.mean(self.getStep("orderPara",-1))
        i=self.images-1
        while abs(lastMeanOrder-np.mean(self.getStep("orderPara",i)))<2.5 and i>0.1*self.images-1:
            i-=1
        #print(" equi at",i,"of",self.images,"images")
        else:
            return i

    
    def getMeanCholMeanNN(self):
        n=0
        mean=0
        for image in range(self.findEqui(),self.images):
            mean+=self.getCholMeanNN(image)
            n+=1
        return mean/n

        
    

    
def orderDistrAni(data):
    fig=plt.figure()
    ax=plt.subplot(111)
    plot,=ax.plot(np.arange(data.paras["minOrder"],data.paras["maxOrder"]+data.paras["DeltaOrder"],data.paras["DeltaOrder"]),data.getSingleOrderDistr(0),"k-")
    ax.set_ylim([0,5])

    def update (i):
        plot.set_ydata(data.getSingleOrderDistr(i))
        print(i)

    ani=animation.FuncAnimation(fig, update, blit=False,frames=300, interval=100, repeat_delay=600)
    #mywriter = animation.FFMpegWriter(fps=10)
    #ani.save('orderDistrAni.avi',writer=mywriter,dpi=400)
    plt.show()


    
def orderDistrAndType(data):
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
        j=i
        ax2.set_title('%s'%j)
        im.set_array(data.getStep("orderPara",j)*data.paras["DeltaOrder"]+data.paras["minOrder"])
        im2.set_array(data.getStep("Type",j))
        print("step:",i,"sameN DPPC",data.getMeanNN(j,0))

    ani=animation.FuncAnimation(fig, update, blit=False,frames=data.images, interval=20,repeat_delay=600)
    #mywriter = animation.FFMpegWriter(fps=10)
    #ani.save('Ani.avi',writer=mywriter,dpi=150)
    plt.show()
    
def cholAni(data):
    fig=plt.figure(figsize=(12,4))
    empty=np.zeros((data.paras["width"],data.paras["height"]))

    ax1=plt.subplot(121)
    im1=ax1.imshow(empty,vmax=1,vmin=-0.3,interpolation='None',cmap='gnuplot')
    plt.colorbar(im1,ax=ax1)
    ax1.set_title('axes title')

    ax3=plt.subplot(122)
    #my_cmap=ListedColormap([(1.0, 0.6, 0.6),(0.1, 0.7, 0.8)])
    im3=ax3.imshow(empty,vmax=1,vmin=0,interpolation='None',cmap="Greys")
    plt.colorbar(im3,ax=ax3)
    
    def update (i):
        j=i
        ax1.set_title('%s'%j)
        im1.set_array(data.getStep("orderPara",j)*data.paras["DeltaOrder"]+data.paras["minOrder"])
        im3.set_array(data.getStep("Chol",j))
        print("step:",j,"chol neighbours",data.getCholMeanNN(j),"mean order",np.mean(data.getStep("orderPara",j)))

    ani=animation.FuncAnimation(fig, update, blit=False,frames=data.images, interval=10,repeat_delay=600)
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
    

def subgridPlot(data,imageNumber):
    fig=plt.figure()
    ax=plt.subplot(111)
    ax.set_axis_off()
    im=ax.imshow(data.getStep("orderPara",imageNumber)*data.paras["DeltaOrder"]+data.paras["minOrder"],vmax=1,vmin=-0.3,interpolation='None',cmap='gnuplot')
    plt.colorbar(im,ax=ax)

    
    chol=data.getStep("Chol",imageNumber)
    chol=np.array(chol,dtype=bool)

    occ=np.where(chol)
    empt=np.where(np.invert(chol))
    scat=ax.scatter(occ[0]-0.5,occ[1]-0.5,c="green",s=100)#,edgecolor='k')
    #ax.scatter(empt[0]-0.5,empt[1]-0.5,facecolors='none',edgecolor='k')

    plt.show()
    
def subgridAni(data):
    fig=plt.figure()
    ax=plt.subplot(111)
    ax.set_axis_off()
    
    empty=np.zeros((data.paras["width"],data.paras["height"]))
    im=ax.imshow(empty,vmax=1,vmin=-0.3,interpolation='None',cmap='gnuplot')
    plt.colorbar(im,ax=ax)
    ax.set_title('axes title')

    occ=[np.array([-0.5]),np.array([-0.5])]
    scat=ax.scatter(occ[0]-0.5,occ[1]-0.5,c="white",s=5)#,edgecolor='k')

    def update (i):
        j=i
        ax.set_title('%s'%j)
        im.set_array(data.getStep("orderPara",j)*data.paras["DeltaOrder"]+data.paras["minOrder"])
        occ=np.where(data.getStep("Chol",j))
        occ=np.array([occ[1],occ[0]]).T-0.5
        scat.set_offsets(occ)

        print("step:",i)

    ani=animation.FuncAnimation(fig, update, blit=False,frames=data.images, interval=500,repeat_delay=600)
    #mywriter = animation.FFMpegWriter(fps=10)
    #ani.save('Ani.avi',writer=mywriter,dpi=150)
    plt.show()
    
    
#data1=data("../MCLattice2/out.h5")
data1=data("out.h5")



#subgridAni(data1)
#orderDistrAni(data1)
#orderDistrAndType(data1)
#cholAni(data1)

#from matplotlib.collections import LineCollection
#from matplotlib.colors import ListedColormap, BoundaryNorm

#n=0
#fig=plt.figure()
#ax=plt.subplot(111)

#ids=np.arange(0,10000)
#np.random.shuffle(ids)
#for j in ids :
    #if data1.getStep("Chol",0)[np.where(data1.getStep("CholIDs",0)==j)]:
        #n+=1
        #print(j)
        #x=[]
        #y=[]
        #last=np.where(data1.getStep("CholIDs",0)==j)

        #for i in range(data1.images):
            #now=np.where(data1.getStep("CholIDs",i)==j)
            #if now[0]-last[0]!=0 or now[1]-last[1]!=0:
                    #x.append(now[0][0])
                    #y.append(now[1][0])
            #last=now

        #counter=np.arange(len(x)-1)

        ##print(x)
        ##print(y)
        #points = np.array([x, y]).T.reshape(-1, 1, 2)
        #segments_raw = np.concatenate([points[:-1], points[1:]], axis=1)
        #segments=[]
        #for seg in segments_raw:
            ##print(abs(np.sum(seg[1]-seg[0])))
            #if np.sum(np.abs(seg[1]-seg[0])) <3:
                ##print(seg)
                #segments.append(seg)
        ##print(segments)
        #norm = plt.Normalize(counter.min(), counter.max())
        #lc = LineCollection(segments, cmap='rainbow', norm=norm)
        ## Set the values used for colormapping
        #lc.set_array(counter)
        ##lc.set_linewidth(2)
        #line = ax.add_collection(lc)
        ##ax.plot(x,y)

        #ax.set_xlim(0,100)
        #ax.set_ylim(0,100)
        #if n==10:
            #break
#fig.colorbar(line, ax=ax)

#plt.show()
    
##exit()








#last=np.where(data1.getStep("CholIDs",0)==45)
#direction=np.zeros((4))

#for j in range(0,10000):
    #if data1.getStep("Chol",0)[np.where(data1.getStep("CholIDs",0)==j)]:
        #for i in range(data1.images):
            #now=np.where(data1.getStep("CholIDs",i)==j)
            ##print(i, now)
            #if now[0]-last[0]!=0 or now[1]-last[1]!=0:
                #if now[0]-last[0]==0 and (now[1]-last[1])%10==1:
                    #direction[0]+=1
                #elif now[0]-last[0]==0 and (now[1]-last[1])%10==9:
                    #direction[1]+=1
                #elif (now[0]-last[0])%10==1 and now[1]-last[1]==0:
                    #direction[2]+=1
                #elif (now[0]-last[0])%10==9 and now[1]-last[1]==0:
                    #direction[3]+=1
                ##else:
                    ##print( now[0]-last[0],now[1]-last[1])
            #last=now

        #print(j, direction/np.sum(direction))
        
        
        

#print("",data1.getMeanCholMeanNN())
#snap(data1,int(argv[1]))
orderdistr=data1.getOrderDistr()
np.save("orderdistr.npy",orderdistr)

#plt.plot(np.arange(-0.5,1.01,0.01),orderdistr,"k-")
#plt.show()

#print(argv[1],np.sum(orderdistr*np.arange(-0.5,1.01,0.01)*0.01))
#plt.savefig("T_%s.png"%argv[1],bbox_inches="tight",dpi=400)


