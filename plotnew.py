import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
import pandas as pd

################################################################################
def histedges_equalN(x, nbin):
    npt=len(x)
    return np.interp(np.linspace(0, npt, nbin + 1), np.arange(npt), np.sort(x))
    
#################################################################################
nam=[r"$\log_{10}[v_{\rm{r}}(\rm{km}/s)]$",r"$\log_{10}[\epsilon_{v_{\rm{r}}}(\rm{km}/s)]$",
r"$\log_{10}[v_{\rm{p}}(\rm{km}/s)]$", 
     r"$\log_{10}[\epsilon_{v_{\rm{p}}}(\rm{km}/s)]$", 
     r"$\log_{10}[v_{\rm{DEC}}(\rm{km}/s)]$", 
     r"$\log_{10}[\epsilon_{v_{\rm{DEC}}}(\rm{km}/s)]$", 
     r"$\log_{10}[|z|(\rm{pc})]$", r"$\log_{10}[\epsilon_{\rm{z}}(\rm{pc})]$",
     r"$\tau(\rm{Gyr})$", r"$\epsilon_{\tau}(\rm{Gyr})$"]#, r"$\sigma_{2,~\tau}(\rm{Gyrs})$"]
###=============================================================================
ndim=11
num= 536967
data=np.zeros((num,ndim))
data=np.loadtxt("./files/Gaia_filtered_raw.txt") 
################################################################################
fid=open("./tab1.txt","w")
fid.close()
nd=6
for i  in range(5):
    b=int(i*2)   
    dat=np.zeros((nd+1,4))
    n0, bins, patches = plt.hist( data[:,b],histedges_equalN( data[:,b],nd))
    dat[:,0]=bins 
    for j in range(nd):
        for k in range(num):
            if( data[k,b]>=dat[j,0] and data[k,b]<dat[j+1,0] and data[k,b+1]>0.0 and data[k,b+1] is not None and 
                    data[k,b+1]<float(2.0*data[k,b]) ):   
                dat[j,1]+=data[k,b]
                dat[j,2]+=data[k,b+1]
                dat[j,3]+=1.0
                break
        dat[j,1]=float(dat[j,1]/(dat[j,3]+0.00000045))
        dat[j,2]=float(dat[j,2]/(dat[j,3]+0.00000045))
        dat[j,0]=float(dat[j,0]+dat[j+1,0])*0.5 
    fid=open("./tab1.txt","a+")
    #np.savetxt(fid, dat[:9,0].reshape(-1,9),fmt="$%.2f$ &")
    np.savetxt(fid, dat[:nd,1].reshape(-1,nd),fmt="$%.2f$ &")
    np.savetxt(fid, dat[:nd,2].reshape(-1,nd),fmt="$%.2f$ &")
    fid.write("\n***************************************************\n")
    fid.close()
################################################################################
x1=[-1.0,-2.5,-2.7,-2.0, 0.0]
x2=[2.2,2.5,2.2, 2.5 ,13.7]
for i in range(5): 
    ll=0
    data1=np.zeros((num,2))
    for k in range(num): 
        if(data[k,i*2] is not None and data[k,i*2+1] is not None and data[k,i*2+1]>0.000000111):  
            if(i<4): data1[ll,:]= np.log10(data[k,i*2]), np.log10(data[k, i*2+1])
            else:    data1[ll,:]= data[k,i*2] , data[k,i*2+1]
            ll+=1
    plt.clf()
    plt.cla()
    fig, ax1= plt.subplots(figsize=(8,6))
    plt.hist(data1[:ll,0],40,histtype='bar',ec='k',facecolor='green',alpha=0.45, rwidth=1.5, label=str(nam[2*i]))
    plt.hist(data1[:ll,1],25,histtype='bar',ec='k',facecolor='red',alpha=0.55, rwidth=1.5, label=str(nam[2*i+1]))
    ax1.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    y_vals = ax1.get_yticks()
    ax1.set_yticklabels(['{:.2f}'.format(x/num) for x in y_vals]) 
    y_vals = ax1.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    plt.xlim([x1[i], x2[i]])
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.legend()
    plt.legend(prop={"size":17.0}, loc='best')
    plt.subplots_adjust(hspace=.0)
    fig3= plt.gcf()
    fig3.savefig("./Histos/histplotd{0:d}.jpg".format(i) , dpi=200)
################################################################################  
''' 
for i in range(5):
    plt.clf()
    plt.cla()
    fig, ax1= plt.subplots(figsize=(8,6))
    plt.scatter(data[:,int(2*i)], data[:, int(2*i+1) ],marker='o',c='b', s=9.0, alpha=1.0)
    ax1.set_xlabel(str(nam[2*i]),fontsize=19, labelpad=0.1)
    ax1.set_ylabel(str(nam[2*i+1]),fontsize=19, labelpad=0.1)
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig3= plt.gcf()
    fig3.savefig("./Histos/scatterplot{0:d}.jpg".format(i) , dpi=200)
'''      
    
