import numpy as np 
import matplotlib.pyplot as plt 
from pylab import *
import pandas as pd 
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import iqr
cmap=plt.get_cmap('viridis')
head0=['age', 'z', 'SigmaR',  'SigmaT',  'SigmaZ', 'Sigmatot' ]
lab=[r"$Age$",r"$z(\rm{pc})$",r"$\sigma_{\rm{R}}(\rm{km/s})$",r"$\sigma_{\rm{T}}(\rm{km/s})$",r"$\sigma_{\rm{Z}}(\rm{km/s})$",r"$\sigma_{\rm{tot}}(\rm{km/s})$"]
lab0=[r"$\rm{Age}$", r"$z$", r"$\sigma_{\rm{R}}$", r"$\sigma_{\rm{T}}$", r"$\sigma_{\rm{Z}}$", r"$\sigma_{\rm{tot}}$"]

filw=open("./files/Fitt.txt","w")
filw.close()

filv=open("./files/Fitt_coef.txt","w")
filv.close()
#######################################################################
def histedges_equalN(x, nbin):
    npt=len(x)
    return np.interp(np.linspace(0, npt, nbin + 1), np.arange(npt), np.sort(x))
#######################################################################    
def discrete_cmap(N, base_cmap=None):
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return(base.from_list(cmap_name, color_list, N) )    
#######################################################################
def FunFitt(x,c1,c2,c3):
    return(c1 + c2*x + c3*x**2.0)##x**0.1 + c2*x**0.3 + c3*x**0.5)
#######################################################################
head=['VR', 'VT', 'VZ',  'Vtot', 'b', 'l', 'age',  'R', 'zb']
df= pd.read_csv("./files/Gaia_filtered.txt", sep=" ",  skipinitialspace = True, header=None, usecols=[0,1,2,3,4,5,6,7],names=head)
print("describe:     ",  df.describe())
print("Columns:  ",  df.columns, "len(columns):  ",  len(df.columns))
print("******************************************************")
f1=open("./files/Gaia_filtered.txt","r")
nm= sum(1 for line in f1) 
par=np.zeros((nm,9)) 
par= np.loadtxt("./files/Gaia_filtered.txt")
par[:,4]=abs(par[:,4])
par[:,8]=abs(par[:,8]*1000.0)##parsec

print("min_z,  max_z:  ", np.min(par[:,8]) ,  np.max(par[:,8]) )
#input("Enter a number ")
#######################################################################
N=25;
ng=N-1
R2S= np.zeros((ng, 4))
MAPE=np.zeros((ng, 4))
MSE= np.zeros((ng, 4))
RMSE=np.zeros((ng, 4))
stds= np.zeros((ng*ng, 6))
iqrs= np.zeros((ng*ng, 6))

n0, bins0, patches = plt.hist( par[:,6], histedges_equalN(par[:,6],N))
print("Age_Bins:    ",  np.min(par[:,6]), np.max(par[:,6]),  n0,   bins0)
age= bins0
cfit=np.zeros((ng,4,3))             
count=0;
 
for i in range(ng):
    l=0
    arry=np.zeros(( int(n0[i]+1) ,9))   
    for j in range(nm): 
        if(float((par[j,6]-age[i])*(par[j,6]-age[i+1]))<0.0  or par[j,6]==age[i] ): 
            arry[l,:]=par[j,:]
            l+=1
    if(l!=int(n0[i])): 
        print("Error:  ",  l,  int(n0[i]) )
        input("Enter a number ")   
    n1, bins1, patches = plt.hist(arry[:l,8],histedges_equalN(arry[:l,8],N) ) 
    print("***************************************************************") 
    print("AGE:  ",  age[i],   n0[i] )   
    print("Height_bins:   ",np.min(arry[:l,8]), np.max(arry[:l,8]),  n1,   bins1)
    zz= bins1
    ########################################
    test=np.zeros((ng,5))
    
    for k in range(ng):
        m=0
        grid=np.zeros(( int(n1[k]+1) , 9))                
        for j in range(l): 
            if(float((arry[j,8]-zz[k])*(arry[j,8]-zz[k+1]))<0.0 or arry[j,8]==zz[k] ): 
                grid[m,:]=arry[j,:]
                m+=1
        if(m!=int(n1[k])): 
            print("Error: ",  m,  int(n1[k]) )
            input("Enter a number ")      
        iq1= iqr(grid[:m,0],axis=None, rng=(25, 75), scale=1.0, nan_policy='propagate', interpolation='linear', keepdims=False) 
        iq2= iqr(grid[:m,1],axis=None, rng=(25, 75), scale=1.0, nan_policy='propagate', interpolation='linear', keepdims=False)   
        iq3= iqr(grid[:m,2],axis=None, rng=(25, 75), scale=1.0, nan_policy='propagate', interpolation='linear', keepdims=False)   
        iq4= iqr(grid[:m,3],axis=None, rng=(25, 75), scale=1.0, nan_policy='propagate', interpolation='linear', keepdims=False)     
        st1= np.std(grid[:m,0])
        st2= np.std(grid[:m,1])
        st3= np.std(grid[:m,2])
        st4= np.std(grid[:m,3])
        Age=   np.mean(grid[:m,6])
        Height=np.mean(grid[:m,8])
        
        
        
            
        stds[count,:]=np.array([ Age, Height, st1, st2, st3, st4 ])  
        iqrs[count,:]=np.array([ Age, Height, iq1, iq2, iq3, iq4 ])
        test[k,:]=  Height, iq1, iq2, iq3, iq4
        count+=1         
    ##################################################
    for j in range(4):    
          
        #m=float(np.max(test[:,j+1]) - np.min(test[:,j+1]))/float(np.max(test[:,0]) - np.min(test[:,0]))
        #ini1 = np.array([np.min(test[:,j+1]) , m*0.9 , 0.1])
        
        ini1=np.array([14.0, 1.0, 1.0])
        fit,pcov=curve_fit(FunFitt, test[:,0], test[:,j+1], ini1 )#,  maxfev = 80000 )
        
        ytest=test[:,j+1]
        ypred=FunFitt(test[:,0] , fit[0] , fit[1] , fit[2] )
        
        R2S[ i,j]= np.abs(1.0-np.sum((ytest-ypred)**2.0)/np.sum((ytest-np.mean(ytest))**2.0)) 
        MAPE[i,j]= np.abs(np.mean(np.abs((ytest-ypred)/ytest)))*100.0
        MSE[ i,j]=np.mean(np.power(ytest-ypred,2.0))  
        RMSE[i,j]=np.sqrt(MSE[i,j])
        #print ("**********************************************" )
        #print ("Age,  Velicity_index:   " ,       Age , j ) 
        #print ("f0:  ", fit[0],  "+/-:  ",  np.sqrt(pcov[0,0])  )
        #print ("f1:  ", fit[1],  "+/-:  ",  np.sqrt(pcov[1,1])  )
        #print ("f2:  ", fit[2],  "+/-:  ",  np.sqrt(pcov[2,2])  )
        #print ("r2s, mape, mse, rmse:  ",  R2S[i,j],  MAPE[i,j],  MSE[i,j],   RMSE[i,j] )    
        #print ("**********************************************" ) 
        cfit[i,j,:]=fit[0], fit[1], fit[2]
        filw=open("./files/Fitt.txt","a")
        filv=open("./files/Fitt_coef.txt","a")
        arr=np.array([np.mean(arry[:l,6]),j,fit[0],np.sqrt(pcov[0,0]),fit[1]*100.0,np.sqrt(pcov[1,1])*100.0,fit[2]*10000.0,np.sqrt(pcov[2,2])*10000.0]) 
        np.savetxt(filw, arr.reshape(-1,8),fmt="$%-3.1f$  &  $%d$  &  $%-4.2f\pm %-3.2lf$  &  $%-4.2f\pm %-3.2lf$  &  $%-4.2f\pm %-3.2lf$")
        filw.close()
        if(j<3): 
            arr2=np.array([np.mean(arry[:l,6]) , j  ,  fit[0]  ,  fit[1] ,  fit[2]   ])
            np.savetxt(filv, arr2.reshape(-1,5),fmt="%-4.2f    %d    %-8.6f    %-8.6f     %-8.6f ")
            filv.close()
        
        '''
        plt.cla()
        plt.clf()
        plt.subplots(figsize=(8,6))
        plt.scatter(test[:,0] , test[:,j+1] ,s=32.0, alpha=1.0)
        plt.plot(zz,FunFitt(zz,*fit),  "r--", lw=2.0)
        plt.xticks(fontsize=17,rotation=0)
        plt.yticks(fontsize=17,rotation=0)
        plt.xlabel(str(lab[1]),  fontsize=18,labelpad=0.0)
        plt.ylabel(str(lab[2+j]),fontsize=18,labelpad=0.0)
        plt.grid("True")
        plt.grid(linestyle='dashed')
        plt.show()
        '''

for i in range(4): 
    filw=open("./files/Fitt.txt","a") 
    arr=np.array([ i, np.mean(R2S[:ng,i]) , np.mean(MAPE[:ng,i]) , np.mean(MSE[:ng,i]) , np.mean(RMSE[:ng,i]) ])
    np.savetxt(filw, arr.reshape(-1,5),fmt="$%d$  & $%.3f$ &  $%.3f$  &  $%.3f$  &  $%.3f$")
    filw.close()
    
##############################################################

fild=open("./files/stds.txt","w")
np.savetxt(fild,stds[:count,:].reshape(-1,6),fmt="%-9.5f    %-9.5lf       %9.5f    %9.5f     %9.5f     %9.5f ")
fild.close()
filf=open("./files/iqrs.txt","w")
np.savetxt(filf,iqrs[:count,:].reshape(-1,6),fmt="%-9.5f    %-9.5lf       %9.5f    %9.5f     %9.5f     %9.5f ")
filf.close()

##############################################################  
colors = plt.cm.jet(np.linspace(0, 1, ng)) 
#print("colors:  ",  colors)
hz=np.arange(0.0, np.max( iqrs[:count,1]) ,  0.005,dtype=float)    
print("*************** HZ:  ", hz,  np.max(par[:,8]))
y1=[29.0, 19.0, 12.0, 17.5]
y2=[59.0, 41.0, 28.0, 37.8]  
for j in range(4):
    plt.cla()
    plt.clf()
    plt.subplots(figsize=(8,6))
    plo=plt.scatter(iqrs[:count,1], iqrs[:count,2+j],s=28.0,c=iqrs[:count,0],cmap=discrete_cmap(24, 'jet'), alpha=1.0)
    for i in range(ng): 
        plt.plot(hz,FunFitt(hz,cfit[i,j,0],cfit[i,j,1],cfit[i,j,2]), color=colors[i], lw=0.5, ls='--')
    cbar = plt.colorbar(plo)
    cbar.set_label(r"$\rm{Age}(\rm{Gyr})$", fontsize=18.0)
    cbar.ax.tick_params(labelsize=16.)
    plt.rcParams["figure.autolayout"]= True
    plt.xticks(fontsize=17,rotation=0)
    plt.yticks(fontsize=17,rotation=0)
    plt.xlim([0.0 , np.max(iqrs[:count,1]) ])
    plt.ylim([ y1[j], y2[j] ])
    plt.xlabel(str(lab[1]),  fontsize=18,labelpad=0.0)
    plt.ylabel(str(lab[2+j]),fontsize=18,labelpad=0.0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig3= plt.gcf()
    fig3.savefig("./Histos/SVIQR_{0:d}.jpg".format(j),dpi=200)
    ###############################################
    plt.cla()
    plt.clf()
    plt.subplots(figsize=(8,6))
    plo= plt.scatter(stds[:count,1] , stds[:count,2+j],s=32.0 , c=stds[:count,0] , cmap=discrete_cmap(24, 'jet'), alpha=1.0)
    cbar = plt.colorbar(plo)
    cbar.set_label(r"$\rm{Age}(\rm{Gyr})$", fontsize=18.0)
    cbar.ax.tick_params(labelsize=16.)
    plt.rcParams["figure.autolayout"]= True
    plt.xticks(fontsize=17,rotation=0)
    plt.yticks(fontsize=17,rotation=0)
    plt.xlim([0.0 , np.max(stds[:count,1]) ])
    #plt.ylim([ y1[j], y2[j] ])
    plt.xlabel(str(lab[1]),  fontsize=18,labelpad=0.0)
    plt.ylabel(str(lab[2+j]),fontsize=18,labelpad=0.0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig3= plt.gcf()
    fig3.savefig("./Histos/SVSTD_{0:d}.jpg".format(j),dpi=200)
#############################################################    
plt.cla()
plt.clf()
plt.subplots(figsize=(8,6))
plt.scatter(iqrs[:count,0] , iqrs[:count,1],s=32.0,c=iqrs[:count,0], alpha=1.0)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlabel(str(lab[0]),fontsize=18,labelpad=0.0)
plt.ylabel(str(lab[1]),fontsize=18,labelpad=0.0)
plt.grid("True")
plt.grid(linestyle='dashed')
fig3= plt.gcf()
fig3.savefig("./Histos/AGE_Z.jpg".format(i),dpi=200)
        
##############################################################
df= pd.DataFrame(stds[:count,:], columns = head0)
corrM = df.corr()
fig, ax = plt.subplots(figsize=(10, 8))
corrM.style.background_gradient(cmap='coolwarm').set_precision(2)
ax= sns.heatmap(corrM, annot=True, xticklabels=lab0, yticklabels=lab0 , annot_kws={"size": 19}, square=True, 
linewidth=1.0, cbar_kws={"shrink": .99}, linecolor="k",fmt=".3f", cbar=True, vmax=1, vmin=0, center=0.0, ax=None, robust=True)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)
plt.xticks(rotation=0,horizontalalignment='right',fontweight='light', fontsize=18)
plt.yticks(rotation=0, horizontalalignment='right',fontweight='light', fontsize=18)
plt.title(r"$\rm{Correlation}~\rm{Matrix}$", fontsize=19)
fig.tight_layout()
plt.savefig("./Histos/corrMatrix_std.jpg", dpi=200)
print("**** Correlation matrix was calculated ******** ")
##############################################################



df= pd.DataFrame(iqrs[:count,:], columns = head0)
corrM = df.corr()
#mask = np.triu(np.ones_like(corrM, dtype=bool))
fig, ax = plt.subplots(figsize=(10, 8))
corrM.style.background_gradient(cmap='coolwarm').set_precision(2)
ax= sns.heatmap(corrM,  annot=True, xticklabels=lab0, yticklabels=lab0,annot_kws={"size": 19}, square=True, 
linewidth=1.0, cbar_kws={"shrink": .99}, linecolor="k",fmt=".3f", cbar=True, vmax=1, vmin=0, center=0.0, ax=None, robust=True)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)
plt.xticks(rotation=0,horizontalalignment='right',fontweight='light', fontsize=18)
plt.yticks(rotation=0, horizontalalignment='right',fontweight='light', fontsize=18)
plt.title(r"$\rm{Correlation}~\rm{Matrix}$", fontsize=19)
fig.tight_layout()
plt.savefig("./Histos/corrMatrix_iqr.jpg", dpi=200)
print("**** Correlation matrix was calculated ******** ")
##############################################################












