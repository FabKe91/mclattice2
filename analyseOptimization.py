import numpy as np
import matplotlib.pylab as plt

scale=1
fig=plt.figure(figsize=(6.4*scale,4.2*scale))
def polynom(paras,x):
    y=0
    for n in range(len(paras)):
        y+=paras[n]*x**n
    return y


distrOut=np.loadtxt("../MCLattice2/DistrOut.txt")

#print(distrOut)
#print(distrOut.shape)

#MD_paras=np.array([-0.14122, 7.51277, -9.36903, -4.43679, -97.86418, 192.92704, 19.37517, -168.20577])#dupc 330
MD_paras=np.array([-0.9767356, 8.69286553, -12.7808724, 12.12000201, -21.41776641, 7.14478559])#dppc 330



order=np.arange(-0.5,1.01,0.01)
MD_distr=np.exp(polynom(MD_paras,order))

ax1=plt.subplot(211)
ax1.plot(order[::2],MD_distr[::2],"rx",label="MD",zorder=100)

for i in range (distrOut.shape[0]):
    ax1.plot(order,distrOut[i,:],label="run %s"%i)
    
ax1.legend()

OptimzeOut=np.loadtxt("../MCLattice2/OptimzeOut.txt")


order=np.arange(-0.5,1.01,0.01)


ax3=plt.subplot(212)


for i in range (OptimzeOut.shape[0]):
    ax3.plot(order,OptimzeOut[i,:]-max(OptimzeOut[i,:]),label="run %s"%i)
#ax3.set_ylim([min(davit_enth),max(davit_enth)])



fit_paras=np.polyfit(order[20:-20],OptimzeOut[-1,20:-20]-max(OptimzeOut[-1,:]),13)[::-1]
for fp in fit_paras:
    print(fp,end=" ")
print()
print(fit_paras)

ax3.plot(order,polynom(fit_paras,order),"r-",label="final",zorder=100)

ax3.legend()


plt.show()
#plt.savefig("optimize_copy_and3N_64steps.png",dpi=500,bbox_inches="tight")
