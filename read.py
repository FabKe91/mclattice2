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
from matplotlib.colors import ListedColormap
import matplotlib.animation as animation

np.set_printoptions(precision=6,linewidth=150)

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
            # following lines important?
            #if isinstance(self.f.attrs[item], int):
            #    self.paras[item] = int(self.f.attrs[item])
            #else:
            #    self.paras[item] = self.f.attrs[item]

            self.paras[item] = self.f.attrs[item][0]

        print("Initialized parameters:\n", *["{}:{}\n".format(str(i), str(j)) for i,j in self.paras.items()])
        
        
    def getStep(self,dataSetName,step):
        return np.array(self.f["/"+dataSetName][step,:,:])

    def getOrderDistr(self):
        counter=np.zeros(int(self.paras["maxOrderIndex"]+1))
        for image in range(self.findEqui(),self.images):
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
    
    def getCholMeanNN(self,imageNumber):
        step=np.array(self.f["/Chol"][imageNumber,:,:],dtype=bool)
        hit=0
        hit+=np.sum(step[np.roll(step,1,axis=0)])
        hit+=np.sum(step[np.roll(step,-1,axis=0)])
        hit+=np.sum(step[np.roll(step,1,axis=1)])
        hit+=np.sum(step[np.roll(step,-1,axis=1)])
        return hit/np.sum(step)**2/4*self.paras["width"]*self.paras["height"]
    
    def getLipidCholNeig(self,imageNumber):
        chol=np.array(self.f["/Chol"][imageNumber,:,:],dtype=bool)
        neigh=np.zeros(chol.shape)
        neigh+=chol
        neigh+=np.roll(chol,1,axis=0)
        neigh+=np.roll(chol,1,axis=1)
        neigh+=np.roll(np.roll(chol,1,axis=0),1,axis=1)
        #print(neigh)
        return neigh

    
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
    
    def getOrderOfNeigh(self):
        orders=np.zeros((5))
        counter=np.zeros((5))
        #print("images: ",self.images)
        for image in range(self.findEqui(),self.images):
            #print(image,orders,counter)
            neighs=self.getLipidCholNeig(image)
            orderParas=self.getStep("orderPara",image)
            for i in range(self.paras["width"]):
                for j in range (self.paras["height"]):
                    orders[neighs[i,j]]+=orderParas[i,j]
                    counter[neighs[i,j]]+=1
                
        print(orders/counter*self.paras["DeltaOrder"]+self.paras["minOrder"],counter/np.sum(counter))

        
def orderDistrAni(data):
    fig=plt.figure()
    ax=plt.subplot(111)
    plot,=ax.plot(np.arange(data.paras["minOrder"],data.paras["maxOrder"]+data.paras["DeltaOrder"],data.paras["DeltaOrder"]),data.getSingleOrderDistr(0),"k-")
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
        j=i
        ax2.set_title('%s'%j)
        im.set_array(data.getStep("orderPara",j)*data.paras["DeltaOrder"]+data.paras["minOrder"])
        im2.set_array(data.getStep("Type",j))
        print("step:",i,"sameN DPPC",data.getMeanNN(j,0))

    ani=animation.FuncAnimation(fig, update, blit=False,frames=data.images//100, interval=20,repeat_delay=600)
    mywriter = animation.FFMpegWriter(fps=1)
    print("Saving animation...")
    ani.save('Ani.avi', writer=mywriter, dpi=150)
    #plt.show()
    
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
    fig, axs = plt.subplots(ncols=2, figsize=(singlex*2,singlex) )

    axs[0].set_axis_off()
    im=axs[0].imshow((data.getStep("orderPara",imageNumber)-50)*0.01,vmax=1,vmin=-0.3,interpolation='None',cmap='gnuplot')
    plt.colorbar(im,ax=axs[0])
    
    axs[1].set_axis_off()
    my_cmap=ListedColormap([(1.0, 0.6, 0.6),(0.1, 0.7, 0.8)])
    im2=axs[1].imshow(data.getStep("Type",imageNumber),vmax=1,vmin=0,interpolation='None',cmap=my_cmap)
    plt.colorbar(im2,ax=axs[1])
    axs[1].axis('off')

def subgridPlot(data,imageNumber):
    plt.rc('text', usetex=True)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig=plt.figure(figsize=(4.980614173228346,3.2))
    ax=plt.subplot(111)
    ax.set_axis_off()
    ax.set_xlim([0,50])
    ax.set_ylim([0,50])
    im=ax.imshow(data.getStep("orderPara",imageNumber)[:50,:50]*data.paras["DeltaOrder"]+data.paras["minOrder"],vmax=1,vmin=-0.3,interpolation='None',cmap='gnuplot')
    plt.colorbar(im,ax=ax,pad=0.005)
    #plt.colorbar(im, cax=make_axes_locatable(ax).append_axes("right", size="5%", pad=0))

    
    chol=data.getStep("Chol",imageNumber)[:51,:51]
    chol=np.array(chol,dtype=bool)

    occ=np.where(chol)
    empt=np.where(np.invert(chol))
    scat=ax.scatter(occ[0]-0.5,occ[1]-0.5,c="w",s=6,edgecolor='k',linewidth=0.2,)
    #ax.scatter(empt[0]-0.5,empt[1]-0.5,facecolors='none',edgecolor='k')

    plt.savefig("snap_chol_image%s.png"%imageNumber,bbox_inches="tight",dpi=300,pad_inches=0)

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

    ani=animation.FuncAnimation(fig, update, blit=False,frames=data.images, interval=1000,repeat_delay=600)
    #mywriter = animation.FFMpegWriter(fps=10)
    #ani.save('Ani.avi',writer=mywriter,dpi=150)
    plt.show()
    


    
#data1=data("../MCLattice2/out.h5")
data1=data("out.h5")


#data1.getOrderOfNeigh()

#subgridAni(data1)
subgridPlot(data1,int(argv[1]))
#orderDistrAni(data1)
##orderDistrAndType(data1)
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
#orderdistr=data1.getOrderDistr()
#np.save("orderdistr.npy",orderdistr)

#plt.plot(np.arange(-0.5,1.01,0.01),orderdistr,"k-")
#plt.show()




############### Plot system evolution  ###################
#fig, axs = plt.subplots(ncols=2, figsize=(singlex,singlex*a43) )
#
#ax = sns.lineplot(x="step", y="order", data=dt_s,  ax=axs[0])
#ax = sns.lineplot(x="step", y="NN",    data=dt_NN, ax=axs[1])
#
#fig.savefig("system_evolution_T{temp}.png")








############## Plot meanNN evolution ###################

