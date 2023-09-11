import numpy as np 
import matplotlib.pyplot as plt 
from pylab import *
import pandas as pd 
from matplotlib import cm 
from scipy.ndimage.filters import gaussian_filter
from matplotlib import rcParams
from matplotlib import colors
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"

################################################################################
ng=int(70)
nh=int(70)

age1=0.1
age2=12.99
dg=float(age2-age1)/ng/1.00

h1=0.0
h2=200.0
dh=float(h2-h1)/nh/1.00
################################################################################
head0=['VR', 'VT', 'VZ',  'Vtot', 'b', 'l', 'age',  'R', 'zb', 'EVR', 'EVT', 'EVZ', 'Evtot' ,'Eage', 'Ezb']
df= pd.read_csv("./files/Gaia_filtered.txt", sep=" ",  skipinitialspace = True, header=None,names=head0)
print("describe:     ",  df.describe())
print("Columns:  ",  df.columns, "len(columns):  ",  len(df.columns))
print("******************************************************")
f1=open("./files/Gaia_filtered.txt","r")
nm= sum(1 for line in f1) 
par=np.zeros((nm,15)) 
par= np.loadtxt("./files/Gaia_filtered.txt")



head=['age', 'z', 'VR',  'VT',  'VZ', 'Vtot', 'Eage', 'Ez', 'EVR', 'EVT', 'EVZ', 'EVTOT']
f5=open("./files/iqrs.txt","r")
nq=sum(1 for line in f5) 
iqrs=np.zeros(( nq , 12 )) 
iqrs=np.loadtxt("./files/iqrs.txt")


'''

f1=open("./files/Fitt_coef.txt","r")
nf= sum(1 for line in f1) 
fitt=np.zeros((nf,8)) 
fitt= np.loadtxt("./files/Fitt_coef.txt")
fitt[:,4]=fitt[:,4]*0.01
fitt[:,5]=fitt[:,5]*0.01
fitt[:,6]=fitt[:,6]*0.0001
fitt[:,7]=fitt[:,7]*0.0001
Fit=np.zeros((int(nf/3),3,8))
kl=0;  i=0;
while(i<nf):
    Fit[kl,0,:]= fitt[i,  :];   
    Fit[kl,1,:]= fitt[i+1,:];  
    Fit[kl,2,:]= fitt[i+2,:];   
    kl+=1;  
    i+=3
print ("Fit[0]:  ",  Fit[:,0,:])
print ("Fit[1]:  ",  Fit[:,1,:])
print ("Fit[2]:  ", Fit[:,2,:])
print (kl,    int(nf/3))
#input("Enter a number ")
'''
################################################################################
def prior1(age, dg):
    ngood=0.0
    for i in range(nm): 
        if(float((age-dg*0.5-par[i,6])*(age+dg*0.5-par[i,6]))<0.0  or age==par[i,6]): 
            ngood+=1.0
    return(float(ngood/nm/1.0))
    
########################################################## 
def prior2(height, dh): 
    ngood=0.0
    for i in range(nm): 
        if(float((height-dh*0.5-par[i,8])*(height+dh*0.5-par[i,8]))<0.0  or height==par[i,8] ): 
            ngood+=1.0
    return( float(ngood/nm/1.0) )
    
#########################################################
def prior3(velocity, dv, nv): 
    ngood=0.0
    for i in range(nm): 
        if(float((velocity-dv*0.5-par[i,nv])*(velocity+dv*0.5-par[i,nv]))<0.0  or velocity==par[i,nv] ): 
            ngood+=1.0
    return( float(ngood/nm/1.0) )    
    
################################################################################   
def likelihood(age, height, velocity, nv):
    co=-1
    if(age<Fit[0,nv,0]): co=0
    elif(age>Fit[kl-1,nv,0] or age==Fit[kl-1,nv,0] ):  co=kl-1
    else:   
        for i in range(kl-1):  
            if( float((age-Fit[i,nv,0])*(age-Fit[i+1,nv,0]))<0.0   or  age==Fit[i,nv,0]  ): 
                 co=i
                 break
    if(co>=0 and co<kl): 
        vmean=float(Fit[co,nv,2]+Fit[co,nv,4]*height+Fit[co,nv,6]*height**2.0)
        dv1=abs(vmean-float(Fit[co,nv,2]+Fit[co,nv,3]+(Fit[co,nv,4]+Fit[co,nv,5])*height+(Fit[co,nv,6]+Fit[co,nv,7])*height**2.0))
        dv2=abs(vmean-float(Fit[co,nv,2]-Fit[co,nv,3]+(Fit[co,nv,4]-Fit[co,nv,5])*height+(Fit[co,nv,6]-Fit[co,nv,7])*height**2.0))
    sigma=float(dv1+dv2)*0.5 
    #print ("param:  ", age, height,    vmean,    sigma)
    #input("Enter a number ")    
    pro= np.exp(-(velocity-vmean)**2.0/(2.0*sigma**2.0))            
    if(sigma>vmean  or pro>1.0  or vmean<0.0  or sigma<0.0 or age<0.0 or age>age2 or height<h1 or height>h2): 
        print("Errors:  ",  sigma,  vmean,   pro,    age,   age2,    age1,  height,   h1,  h2)
        input("Enter a number ")
    return(pro)
################################################################################

f1=open("./files/postb.txt","r")
nk= sum(1 for line in f1) 
post0=np.zeros((nk , 4)) 
post0= np.loadtxt("./files/postb.txt")
cont=np.zeros((ng,nh,2))
k=0
for i in range(ng): 
    age=float(age1+i*dg) 
    print ("Age:  ", age,  i )
    #prob1=np.zeros((nh))   
    #prob2=np.zeros((nh))   
    for j in range(nh): 
        height=float(h1+j*dh)
        cont[i,j,0]=post0[k,2]
        cont[i,j,1]=post0[k,3]
        k+=1
'''

fild=open("./files/post2.txt","w")
filc=open("./files/post3.txt","w")
fild.close();  
filc.close()
vr, vt, vz= 45.0 , 28.5, 17.0
post=np.zeros((ng*nh, 4))
cont=np.zeros((ng,nh,2))
k=0
for i in range(ng): 
    age=float(age1+i*dg) 
    print ("Age:  ", age,  i )
    prob1=np.zeros((nh))   
    prob2=np.zeros((nh))   
    for j in range(nh): 
        height=float(h1+j*dh)
        prior= prior1(age,dg)*prior2(height,dh)
        pro1= abs(prior*likelihood(age,height,vr,0)/prior3(vr,2.0,0)) *100.0
        pro2= abs(prior*likelihood(age,height,vt,1)/prior3(vt,2.0,1)) *100.0
        prob1[j], prob2[j]=pro1, pro2
        post[k,:]=age, height, pro1, pro2 
        cont[i,j,0]=pro1
        cont[i,j,1]=pro2
        print ("height:  ", height,  post[k,:]   )
        k+=1
        
    fild= open("./files/post2.txt","a")
    np.savetxt(fild,prob1.reshape(1,nh),fmt="%.9f ")
    fild.close()  
    filc= open("./files/post3.txt","a")
    np.savetxt(filc,prob2.reshape(1,nh),fmt="%.9f ")
    filc.close()     
    print("************************************") 
fild=open("./files/postb.txt","w")
np.savetxt(fild,post.reshape(-1,4),fmt="%.4f   %.4f   %.10f   %.10f")
fild.close()        
'''
################################################################################
xx=np.zeros((ng))
xx=[age1+dg*i for i in range(ng)]#np.arange(age1, age2, dg, dtype=float)
yy=np.zeros((ng))
yy=[h1+dh*i for i in range(ng)]
print (xx, yy)
cont[:,:,0] = gaussian_filter(cont[:,:,0],0.95)
cont[:,:,1] = gaussian_filter(cont[:,:,1],0.95)
mm= np.max(cont[:,:,0]+cont[:,:,1])
m1= np.max(cont[:,:,0])
m2= np.max(cont[:,:,1])

################################################################################


plt.cla()       
plt.clf()
fig,ax=plt.subplots(figsize=(8,6))

im =plt.imshow((cont[:,:,1]+cont[:,:,0])/mm, origin='lower',cmap='YlGnBu', extent=(age1,age2, h1, h2),interpolation='nearest',aspect='auto',alpha=0.8)

plo= plt.contour( xx, yy , cont[:,:,0]/m1 ,levels=[0.7, 0.8, 0.9], 
origin='lower',extend='both',colors='k',linewidths=np.arange(0.5,3.5,0.9), extent=(age1, age2, h1, h2))

plt.contour( xx, yy , cont[:,:,1]/m2 ,levels=[0.7, 0.8, 0.9], origin='lower',extend='both',colors='m',linestyles='dashed',linewidths=np.arange(0.5,3.5,0.9), extent=(age1, age2, h1, h2))


plt.scatter(4.5,50.0, marker='*',s=35, c='c')

cbar = plt.colorbar(im, shrink=1.0,pad=0.03)
cbar.set_label(r"$\rm{Posterior}~\rm{Probability}$", fontsize=18.0,labelpad=0.0)
cbar.ax.tick_params(labelsize=16.)
cbar.add_lines(plo)

plt.xticks(fontsize=18,rotation=0)
plt.yticks(fontsize=18,rotation=0)
plt.xlim([age1, age2 ])
plt.ylim([h1,   h2 ])
plt.xlabel(r"$\tau(\rm{Gyr})$",fontsize=19,labelpad=0.0)
plt.ylabel(r"$|z|(\rm{pc})$",fontsize=19,labelpad=0.0)
plt.title(r"$v_{\rm{R}}(\rm{km}/s)=45.0\pm 2.0,~v_{\rm{T}}(\rm{km}/s)=28.5\pm 2.0$", fontsize=19.0)
plt.subplots_adjust(hspace=.0)
fig.tight_layout(pad=0.1)
fig= plt.gcf()
fig.savefig("./mapb.jpg",dpi=200)


################################################################################

plt.cla()       
plt.clf()
fig,ax=plt.subplots(figsize=(8,6))
im = ax.imshow(cont[:,:,0], interpolation='bilinear', origin='lower',cmap='viridis')#, extent=(age1,age2, h1, h2))
fig.colorbar(im) # Add a colorbar to a plot
ax.set_title('Filled Contours Plot')
fig= plt.gcf()
fig.savefig("./map2b.jpg",dpi=200)









        
