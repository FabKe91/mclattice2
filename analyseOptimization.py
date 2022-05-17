import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from sizes import *

def polynom(paras,x):
    y=0
    for n in range(len(paras)):
        y+=paras[n]*x**n
    return y


set_context()

#############################

OUTNAME = "optimize_omegaDPPC_test.png"

#MD_DIST_PARAM = np.array([-0.14122, 7.51277, -9.36903, -4.43679, -97.86418, 192.92704, 19.37517, -168.20577])#dupc 330
MD_DIST_PARAM = np.array([-0.9767356, 8.69286553, -12.7808724, 12.12000201, -21.41776641, 7.14478559])#dppc 330

order = np.arange(-0.5,1.01,0.01)
MD_distr = np.exp(polynom(MD_DIST_PARAM,order))


###################


fig, axs = plt.subplots(nrows=2, figsize=(singlex,singlex*1.5*a43))



################################### Visualize Distr out #####################################
distrOut=np.loadtxt("./DistrOut.txt")
print(distrOut)

#pal = sns.color_palette("colorblind", n_colors=distrOut.shape[0]+1)
palb = sns.color_palette("Blues", n_colors=distrOut.shape[0])
palr = sns.color_palette("Reds", n_colors=1)


### MC distributions ###
for i in range (distrOut.shape[0]):
    ax = sns.lineplot(x=order,y=distrOut[i,:], label="run %s"%i, color=palb[i], ax=axs[0])


### MD distributions ###
ax = sns.lineplot(x=order[::2], y=MD_distr[::2], label="MD", zorder=100, ax=axs[0], color=palr[0], lw=1.0)


    
ax.set_ylim(0, 4.01)
ax.set_xlim(-0.55, 1.05)
ax.set_ylabel("p(S)")
axs[0].get_legend().remove()

handles, labels = axs[0].get_legend_handles_labels()
axs[0].legend(handles=[handles[-1]], labels=[labels[-1]], title="", loc="best")


################################################################################################




################################### Visualize Optimize out #####################################
OptimizeOut=np.loadtxt("./OptimizeOut.txt")


order=np.arange(-0.5,1.01,0.01)

palb = sns.color_palette("Blues", n_colors=OptimizeOut.shape[0])
palr = sns.color_palette("Greys", n_colors=1)


#### Visualize convergence of lnOmega ###
for i in range (OptimizeOut.shape[0]):
    ax = sns.lineplot(x=order, y=OptimizeOut[i,:]-max(OptimizeOut[i,:]), label="run %s"%i, ax=axs[1], color=palb[i])


#### Visualize fit ###
fit_paras = np.polyfit(
        order[20:-20], 
        OptimizeOut[-1,20:-20] - max(OptimizeOut[-1,:]),
        12
        )[::-1]

for fp in fit_paras:
    print(fp,end=" ")
print()
print(fit_paras)
ax = sns.lineplot(x=order, y=polynom(fit_paras,order), label="fit", zorder=100, palette=palr, ax=axs[1], color=palr[0], lw=1.0)

ax.set_ylim(-50, 5)
ax.set_xlim(-0.55, 1.05)
ax.set_xlabel("ln$\Omega$")
ax.set_xlabel("Order parameter S")
axs[1].get_legend().remove()

handles, labels = axs[1].get_legend_handles_labels()
axs[1].legend(handles=[handles[-1]], labels=[labels[-1]], title="", loc="best")

################################################################################################

#plt.show()
plt.tight_layout()
plt.savefig(OUTNAME)

