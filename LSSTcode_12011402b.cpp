/// checked 12/01/1402
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>

using std::cout;
using std::endl;
using std::cin;
///**************************************************************
const int Num=2000;
const double MaxD=20.0;///kpc
const double RA=180.0/M_PI;
const double step=MaxD/(double)Num/1.0;///step in kpc
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity=299792458.0;//velosity of light
const double M_sun=1.98892*pow(10.,30); //in [kg].
const double vro_sun=226.0;
const double logT_sun=log10(5778.0);
const double AU=1.495978707*pow(10.0,11.0);
const double year=365.2425;//days
const double Rsun=6.9634*pow(10.0,8.0); ///solar radius [meter]
const double ddeg=0.25;
const int NM=20000;///number of rows in file 'Monte.txt'
const int M=6;///No. of filter  ugrizy
const int LL=5;
const double binary_fraction=double(2.0/3.0);
const int    GG=900;///number of bins of tE distribution
const double tE_min=0.0;
const double tE_max=300.0;  
const double stept=double((tE_max-tE_min)/GG); 

const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;
const int Nag=24; 

///============================ Besancon constant ==============///
const double Dsun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]=  {3.1,2.5,3.1,3.1};
const int N1=27004, N2=36558, N3=2612, N4=492;///CMD_BESANCON, thinD, bulge, thickD, halo

///=========================LSST constant====================================///
const double Tobs=10.0*year;///LSST observational time 10 years
const double cadence=3.0;//in days 
const double texp=30.0;///in second
const double Delta2=0.005;///systematic errors
const double gama[M]={0.037,0.038,0.039,0.039,0.040,0.040};///a function of sky brightness and airmass in different wavelengths.
const double seeing[M]={0.77,0.73,0.70,0.67,0.65,0.63};///seeing
const double msky[M]={22.9,22.3,21.2,20.5,19.6,18.6};
const double Cm[M]=  {22.92,24.29,24.33,24.20,24.07,23.69};
const double Dci[M]= {0.67,0.21,0.11,0.08,0.05,0.04};
const double km[M]=  {0.451,0.163,0.087,0.065,0.043,0.138};///sky extinction
const double wav[M+2]= {0.3543,0.477,0.6231,0.7625,0.9134,1.004,2.1739,0.7865};///in[micrometer] u g r i z y k Ic
const double sigma[M+2]={0.022,0.02,0.017,0.017,0.027,0.027,0.04,0.017};
const double detect[M]={23.4,24.6,24.3,23.6,22.9,21.7};///depth of single visit in ugriz
const double satu[M]=  {15.2,16.3,16.0,15.3,14.6,13.4};///saturation limit of single visit in ugriz
const double M50[M]={23.68,24.89,24.43,24.00,24.45,22.60};
const double FWHM[M+1]={1.22087,1.10136,0.993103,0.967076,0.951766,0.936578,0.967076*1.65};//LSST [arcsec]
const int YZ=3578;//// No.yzma.txt rows
///======================================================
struct source{
    int nums,struc,cl;
    int Nevent,Nstar;
    double Ds,TET,FI,ratios,u0;
    double Mab[M],Map[M];
    double lon,lat,col;
    double logl,logt;
    double od_disk,od_ThD,od_bulge,od_halo,opd;///optical  depth
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart,nstarti;
    double nsdis[Num],nddis[Num];
    double ro_star,Rs,magni;
    double Fluxb[M],magb[M]; 
    double blend[M],nsbl[M]; 
    double FIsource,Isource,MIsource;
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
    double blendI, deltao;
};
struct lens{
    double Ml,Dl,ratiol,vl,vs,Vt,xls;
    double rhomaxl,tE,RE;
    int numl,struc;
};
struct detection{
   int det,peak;
   double deltam[M],Dcm[M],dchi;
   double t,tmin,tmax,t0;
   double ave_re,ave_vl,ave_vs;
   double ave_opt,ave_opt2,ave_te,ave_u0,ave_Ds,ave_Dl;
   double ave_vt,ave_ml,ave_col,aveI; 
   double ave_aps[M],ave_ex[M],ave_dmag[M],neve[M],ave_apbd[M],ave_bl[M],ave_apb[M];
   double ratio,nsdet[M],nstar,nevent,nfin; 
   double ave_npsf[M]; 
   double m5c[NM]; int filter[NM];
   double ave_tei;///<effi/tE>
   double te[GG+1],effs[GG+1],effd[GG+1]; 
   double effs_com[GG+1],effd_com[GG+1];
   double ave_fred;  
   double weight[M],ave_w;
   int N_sim,N_det;
};
struct CMD{
    double logt_d[N1],logl_d[N1],Mab_d[M][N1],col_d[N1],MI_d[N1]; int cl_d[N1];  ///thin disk
    double logt_b[N2],logl_b[N2],Mab_b[M][N2],col_b[N2],MI_b[N2]; int cl_b[N2];  /// bulge
    double logt_t[N3],logl_t[N3],Mab_t[M][N3],col_t[N3],MI_t[N3]; int cl_t[N3];  ///thick disk
    double logt_h[N4],logl_h[N4],Mab_h[M][N4],col_h[N4],MI_h[N4]; int cl_h[N4];  /// halo
};
struct extinc{
   double x[M+2]; 
   double dis[100];///distance 
   double Extks[100];///ks-band extinction
   double Ai[M],Av,Aks,AI;
   double a[M+2],b[M+2];
   double exks8;
};
struct dispersion{
   double sigR0[Nag], sigR1[Nag], sigR2[Nag];  
   double sigT0[Nag], sigT1[Nag], sigT2[Nag];
   double sigZ0[Nag], sigZ1[Nag], sigZ2[Nag];
   double age[Nag];
};
///====================== TEMPELATE ==========================
int Extinction(extinc & ex,source & s);
void read_cmd(CMD & cm);
void func_source(source & s,CMD & cm, extinc & ex);
void func_lens(lens & l, source & s, dispersion & di);
void vrel(source & s,lens & l, dispersion & di);
void Disk_model(source & s);
void optical_depth(source & s);
void lensing(source & s, lens & l,detection & d,int);
double Interpol(double ds, extinc & ex);
double RandN(double sigma,int);
///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//


// SEED Generation By signlesskid - Part 0
	time_t _timeNow;
  	unsigned int _randVal;
	unsigned int _dummyVal;
    FILE * _randStream;
///=========================================

int main(){
///================================================
    // SEED Generation By signlesskid - Part 0 //
    // SRC: http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
   // gettimeofday(&newTime,NULL);
   // ftime(&tp);
	time(&_timeNow);
	_randStream = fopen("/dev/urandom", "r");
	_dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
///================================================
//  SEED Generation By signlesskid //
    _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
    srand(_randVal);
//=================================================

   time(&_timeNow);
   printf("START time:   %s",ctime(&_timeNow));

     source s;
     lens l;
     detection d;
     CMD cm;
     extinc ex;
     dispersion di; 
     read_cmd(cm);

     char filename1[40], filename2[40];
     FILE* result1; FILE* result2;  FILE* result3; 
     FILE* teefi; FILE* TEE;   FILE* fil1;  //FILE* fil2;
     result1=fopen("./files/res120102A.txt", "a+");
     result2=fopen("./files/res120102B.txt", "a+");
     result3=fopen("./files/res120102C.txt", "a+");
     fclose(result1); 
     fclose(result2);
     fclose(result3);
     
     
     
     
     for(int i=0; i<M; ++i){
     if(texp!=30.0) d.Dcm[i]=Dci[i]-1.25*log10(1.0+(pow(10.0,0.8*Dci[i])-1.0)/(texp/30.0));
     else d.Dcm[i]=0.0;}
     
     double nbb,effe,maga[M],lte,Fax,Fbx,sumd=0.0; 
     int vf,ffd[M],flag_det, flagt;
     double Nmicro=0.0,N_bla=0.0,lonn,countt,countd, we;


     FILE* Mont;
     Mont=fopen("./files/Monte2.txt","r");
     for(int i=0; i<NM; ++i){
     fscanf(Mont,"%lf   %d\n",&d.m5c[i],&d.filter[i]);  
     if(d.m5c[i]<20.0 or d.m5c[i]>27.0 or d.filter[i]>=6 or d.filter[i]<0){ 
     cout<<"ERROR airmass: "<<d.m5c[i]<<"\t filter: "<<d.filter[i]<<endl;  int yye; cin>>yye; }} 
     fclose(Mont);

     double err1, err2, err3;
 
     int counter; 
     FILE* Mont1;
     Mont1=fopen("./files/Fitt_coef.txt","r");
     for(int i=0; i<Nag; ++i){
     fscanf(Mont1,"%lf %d %lf %lf %lf %lf %lf %lf\n",&di.age[i],&counter,&di.sigR0[i],&err1, &di.sigR1[i],&err2, &di.sigR2[i],&err3);  
     fscanf(Mont1,"%lf %d %lf %lf %lf %lf %lf %lf\n",&di.age[i],&counter,&di.sigT0[i],&err1, &di.sigT1[i],&err2, &di.sigT2[i],&err3);  
     fscanf(Mont1,"%lf %d %lf %lf %lf %lf %lf %lf\n",&di.age[i],&counter,&di.sigZ0[i],&err1, &di.sigZ1[i],&err2, &di.sigZ2[i],&err3);
     
     di.sigR1[i]=di.sigR1[i]*0.01;  
     di.sigT1[i]=di.sigT1[i]*0.01;  
     di.sigZ1[i]=di.sigZ1[i]*0.01;  
     di.sigR2[i]=di.sigR2[i]*0.0001;  
     di.sigT2[i]=di.sigT2[i]*0.0001;  
     di.sigZ2[i]=di.sigZ2[i]*0.0001;  
     
     cout<<"sigR0:  "<<di.sigR0[i]<<"\t sigR1:  "<<di.sigR1[i]<<"\t sigR2:  "<<di.sigR2[i]<<"\t i:  "<<i<<endl;
     cout<<"sigT0:  "<<di.sigT0[i]<<"\t sigT1:  "<<di.sigT1[i]<<"\t sigT2:  "<<di.sigT2[i]<<"\t i:  "<<i<<endl;
     cout<<"sigZ0:  "<<di.sigZ0[i]<<"\t sigZ1:  "<<di.sigZ1[i]<<"\t sigZ2:  "<<di.sigZ2[i]<<"\t i:  "<<i<<endl;
     cout<<"*********************************************"<<endl;
     } 
     fclose(Mont1);
     cout<<"File of coefficients of polynomials are read !!!"<<endl;




    //// checked
    for (int i=0; i<(M+2); ++i){
    ex.x[i]=double(1.0/wav[i]); 
    if(ex.x[i]>=0.3 and ex.x[i]<=1.1){
    ex.a[i]=+0.574*pow(ex.x[i],1.61);
    ex.b[i]=-0.527*pow(ex.x[i],1.61);}
    else if(ex.x[i]>1.1 and ex.x[i]<=3.3){
    double y=ex.x[i]-1.82;
    ex.a[i]=1.0+0.17699*y-0.50447*y*y-0.02427*y*y*y+0.72085*y*y*y*y+
          0.01979*pow(y,5.0)-0.7753*pow(y,6.0)+0.32999*pow(y,7.0);
    ex.b[i]=1.41338*y+2.28305*y*y+1.07233*y*y*y-5.38434*y*y*y*y-0.62251*pow(y,5.0)+
          5.3026*pow(y,6.0)-2.09002*pow(y,7.0);}
    else{
    cout<<"ERROR(x) is larger than 3.3 "<<ex.x[i]<<endl;
    int yye;   cin>>yye;
    if(ex.x[i]<5.9){Fax=0.0; Fbx=0.0; }
    else{
    Fax=-0.04473*pow(ex.x[i]-5.9,2.0)-0.009779*pow(ex.x[i]-5.9,3.0);
    Fbx= 0.21300*pow(ex.x[i]-5.9,2.0)+0.1207*pow(ex.x[i]-5.9,3.0); }
    ex.a[i]=1.752-0.316*ex.x[i]-0.104/(pow(ex.x[i]-4.67,2.0)+0.341)+Fax;
    ex.b[i]=-3.090+1.825*ex.x[i]+1.206/(pow(ex.x[i]-4.62,2.0)+0.263)+Fbx;}
    cout<<"wave-length(micrometer):  "<<wav[i]<<"\t x[i]: "<<ex.x[i]<<endl;  
    cout<<"a[i]: "<<ex.a[i]<<"\t b[i]: "<<ex.b[i]<<endl; 
    cout<<"***************************************"<<endl;}
  

    //FILE* TEE;  
    //TEE=fopen("./files/tef_last.txt","w"); 
    //fclose(TEE);
    for(int i=0; i<(GG+1); ++i){d.effs_com[i]=0.0;   d.effd_com[i]=0.0;}
    double Nogle[Num]={0.0}; 

///===================== Monte Carlo Simulation ================================
    for(s.lat=-3;  s.lat<=3.0;   s.lat+=ddeg){///degree
    for(lonn=-9.0;  lonn<=-1.0;    lonn +=ddeg){//degree      
    if(lonn<=0.0)   s.lon=360.0+lonn;
    else            s.lon=lonn;
    s.TET=double(360.0-s.lon)/RA;///radian
    s.FI= double(s.lat/RA);///radian
    
    
    cout<<">>>>>>>>>>>>>>>> NEW STEP <<<<<<<<<<<<<<<<"<<endl;
    cout<<"latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
    
    
    d.nstar=0.00001745693;  nbb=0.0; 
    d.ave_Ds=d.ave_Dl=d.ave_vt=d.ave_ml=d.ave_u0=d.ave_tei=0.0;
    d.ave_re=d.ave_vl=d.ave_vs=d.ave_te=d.nevent=d.ave_col=0.0;
    d.ave_opt=d.ave_opt2=d.aveI=s.nstart=s.nstarti=d.ratio=d.nfin=0.0; d.ave_fred=0.0;
    for(int j=0; j<M; ++j){ 
    d.ave_apb[j]=d.ave_apbd[j]=d.ave_dmag[j]=d.ave_npsf[j]=0.0; 
    d.ave_bl[j]=d.ave_aps[j]=d.ave_ex[j]=0.0; d.nsdet[j]=d.neve[j]=0.00000645974;}
    for(int i=0;  i<Num;  ++i){s.nsdis[i]=0.0000145776534762; s.nddis[i]=Nogle[i]=0.0;}
    for(int i=0; i<(GG+1);++i){d.effd[i]=0.0;  d.effs[i]=0.0000000754765;  d.te[i]=double(i*stept);}
    d.N_sim=d.N_det=0; 
    flagt=0; 

    Disk_model(s);///checked
    if(Extinction(ex,s)==1){
    sprintf(filename1,"./files/distribution/%c%.2lf%c%.2lf.dat",'w',s.lat,'_',s.lon);
    fil1=fopen(filename1,"w");
    if(!fil1){cout<<"cannot make file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl; exit(0);}

    //sprintf(filename2,"./files/distribution/%c%.2lf%c%.2lf.dat",'Q',s.lat,'_',s.lon);
    //fil2=fopen(filename2,"w");
    //if(!fil2){cout<<"cannot make file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl; exit(0);}

    teefi=fopen("./files/tefzzz.txt","w"); 
    fclose(teefi);
   
    _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
    srand(_randVal);


   do{
   d.nstar+=1.0;///number of simulated stars
   int uue=0; 
   do{
   ++uue;  
   func_source(s,cm,ex);///checked
   optical_depth(s);///checked
   func_lens(l,s, di);///checked
   if(int(uue)%10000==0){
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);  
   cout<<"uue: "<<uue<<"\t l.Vt: "<<l.Vt<<"\t tE: "<<l.tE<<"\t Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl; }
   }while(l.Vt<10.0 or  l.Vt>250.0 or  l.tE<=0.5 or l.tE>tE_max); 
  
  
   s.nsdis[s.nums]+=1.0;///No.of simulated stars at Ds
   d.ave_w=0.0;  sumd=0.0; N_bla=0.0;    
    
   for(int i=0; i<M; ++i){
   ffd[i]=0; 
   d.weight[i]=0.0; 
   maga[i]=-2.5*log10(fabs(pow(10.0,-0.4*s.Map[i])*(s.magni-1.0)+s.Fluxb[i]));
   if(s.magb[i]>=satu[i] and maga[i]<=detect[i]){
   d.weight[i]=double(1.0/s.nsbl[i]); 
   
   ffd[i]=1;  
   sumd+= ffd[i]; 
   N_bla+=d.weight[i];  
   
   d.nsdet[i]   +=fabs(1.0000*   d.weight[i]); 
   d.ave_aps[i] +=fabs(s.Map[i]* d.weight[i]);
   d.ave_apb[i] +=fabs(s.magb[i]*d.weight[i]);
   d.ave_ex[i]  +=fabs(ex.Ai[i]* d.weight[i]);
   d.ave_npsf[i]+=fabs(s.nsbl[i]*d.weight[i]);} }

   if(s.Isource<17.5){
   Nogle[s.nums] +=1.0; }



   if(ffd[0]>0 or ffd[1]>0 or ffd[2]>0 or ffd[3]>0 or ffd[4]>0){    
   
   d.ave_w=double(N_bla/sumd);  
   s.nddis[s.nums]+=1.0*d.ave_w; 
   
   nbb +=     fabs(1.00*d.ave_w);
   d.ave_opt+=fabs(s.opd*1.0e6*d.ave_w); 
   d.ave_col+=fabs(s.col*d.ave_w);
   d.aveI +=  fabs(s.Isource*d.ave_w);
   
   d.N_sim+=1;
   lensing(s,l,d,d.N_sim);///checked


   flag_det=0; 
   if(d.det==1 and d.peak==1 and d.dchi>=500.0){///detected lensi
   flag_det= 1;     
   d.N_det +=1; 
   
   for(int i=0 ; i<M; ++i){
   d.neve[i] +=   fabs(1.0000000*  d.weight[i]); 
   d.ave_bl[i]+=  fabs(s.blend[i]* d.weight[i]); 
   d.ave_dmag[i]+=fabs(d.deltam[i]*d.weight[i]);
   d.ave_apbd[i]+=fabs(s.magb[i]*  d.weight[i]);}

   d.nevent+=fabs(1.00*d.ave_w); 
   d.ave_te+=fabs(l.tE*d.ave_w);
   d.ave_Ds+=fabs(s.Ds*d.ave_w);
   d.ave_Dl+=fabs(l.Dl*d.ave_w);
   d.ave_vt+=fabs(l.Vt*d.ave_w);
   d.ave_vl+=fabs(l.vl*d.ave_w);
   d.ave_vs+=fabs(l.vs*d.ave_w);
   d.ave_ml+=fabs(l.Ml*d.ave_w);
   d.ave_u0+=fabs(s.u0*d.ave_w);   
   d.ave_opt2+= s.opd*1.0e6*d.ave_w; 
   d.ave_re+=fabs(l.RE/AU)*d.ave_w;
   if(s.cl<=4) d.ave_fred+=1.0*d.ave_w;}
   

    vf=-1; 
    for(int j=0; j<GG; ++j){
    if((l.tE-d.te[j])*(l.tE-d.te[j+1])<0.0 or l.tE==d.te[j]){vf=j;  j=GG+3;}  }
    if(vf==-1 and (l.tE>tE_max or l.tE==tE_max)) vf=GG; 
    if(vf==-1 and (l.tE<tE_min or l.tE==tE_min)) vf=0;    
    if(vf<0 or l.tE<0.0 or vf>GG){cout<<"ERROR tE: "<<l.tE<<"\t vf: "<<vf<<endl;  int uue; cin>>uue;}
    d.effs[vf]    +=1.0;  
    d.effs_com[vf]+=1.0;
    if(flag_det>0.5){d.effd[vf]+=1.0;  d.effd_com[vf]+=1.0;}
    
     
    teefi=fopen("./files/tefzzz.txt","a+");
    fprintf(teefi,"%.7lf   %e  %d  %d\n",l.tE,d.ave_w,vf,flag_det);
    fclose(teefi);


    if(flag_det==1){
    fprintf(fil1,"%d %d %d %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.6lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %d %d   %.8lf\n",
    d.N_sim,s.struc,l.struc,l.Ml,l.RE/AU, l.Dl, s.Ds, l.Vt, l.tE, l.vl, l.vs,
    s.logl,s.logt,s.Rs,s.ro_star,s.col,ex.Av,ex.AI,s.u0,d.t0,s.opd*1.0e6,s.cl,flag_det,d.ave_w);///23
    //for(int i=0; i<M; ++i)
    //fprintf(fil2,"%d  %.6lf  %.8lf  %.6lf  %.4lf   %.2lf  %.4lf  %.4lf  %d   %.8lf\n",
    //d.N_sim,s.magb[i],ex.Ai[i],d.deltam[i],s.blend[i],s.nsbl[i],s.Map[i],s.Mab[i],flag_det,d.weight[i]);
    }
     
     /*if(int(d.N_det)%100000==0 and  flag_det==1){
     cout<<"**********************************************"<<endl;
     cout<<"detecable: "<<sumd<<"\t flag_det: "<<flag_det<<"\t d.nevent: "<<d.nevent<<endl;
     cout<<"nbb: "<<nbb<<"\t ave_w: "<<d.ave_w<<"\t N_det: "<<d.N_det<<endl;
     cout<<"weight0: "<<d.weight[0]<<"\t weight1: "<<d.weight[1]<<"\t weight2: "<<d.weight[2]<<endl;
     cout<<"weight3: "<<d.weight[3]<<"\t weight4: "<<d.weight[4]<<endl;
     cout<<"latitude: "<<s.lat<<"\t longtitude: "<<s.lon<<endl;
     cout<<"Ml (M_sun): "<<l.Ml<<"\t RE (AU): "<<l.RE/AU<<"\t Dl(Kpc): "<<l.Dl<<endl;
     cout<<"relative velosity (km/s): "<<l.Vt<<"\t Vs (km/s): "<<l.vs<<"\t Vl(km/s): "<<l.vl<<endl;
     cout<<"Ds: "<<s.Ds<<"\t tE (days): "<<l.tE<<"\t ro_star"<<s.ro_star<<endl;
     cout<<"app_mag[r]: "<<s.Map[2]<<"\t app_mag(i): "<<s.Map[3]<<endl;
     cout<<"app_mag_bg[r]: "<<s.magb[2]<<"\t app_mag_bg(i): "<<s.magb[3]<<endl;
     cout<<"belnding_r(factor): "<<s.blend[2]<<"\t belnding_i: "<<s.blend[3]<<endl;
     cout<<"optical_depth(x10^6): "<<s.opd*1.0e6<<"\t u0: "<<s.u0<<endl;
     cout<<"N_blend[r_band]: "<<s.nsbl[2]<<"\t N_blend[i_band]: "<<s.nsbl[3]<<endl;
     cout<<"det: "<<d.det<<"\t peak : "<<d.peak<<"\t dci2: "<<d.dchi<<endl;
     cout<<"maga[0]:"<<maga[0]<<"\t maga[1]: "<<maga[1]<<"\t maga[2]: "<<maga[2]<<endl;
     cout<<"maga[3]: "<<maga[3]<<"\t maga[4]: "<<maga[4]<<endl;
     cout<<"*********************************************************************"<<endl; }*/
     if(sumd>M or d.ave_w>1.000001 or d.ave_w<=0.0 or s.Ds>20.0 or l.tE<0.0 or l.tE>tE_max or 
     s.blend[0]>1.0001 or s.blend[1]>1.0001 or s.blend[2]>1.0001 or s.blend[3]>1.0001 or s.blend[4]>1.0001){
     cout<<"ERROR sumd: "<<sumd<<"\t N_bla: "<<N_bla<<"\t ave_w: "<<d.ave_w<<"\t Ds: "<<s.Ds<<endl;  
     cout<<"b0: "<<s.blend[0]<<"\t b1: "<<s.blend[1]<<"\t b2: "<<s.blend[2]<<"\t b3: "<<s.blend[3]<<"\t b4:"<<s.blend[4]<<endl;
     int uue;  cin>>uue; }
     }
     
     //cout<<"N_det:   "<<d.N_det<<"\t neve_0:  "<<d.neve[0]<<"\t neve_1:  "<<d.neve[1]<<endl;
     //cout<<"neve_2:  "<<d.neve[2]<<"\t neve_3:  "<<d.neve[3]<<"\t neve_4:  "<<d.neve[4]<<"\t neve_5:  "<<d.neve[5]<<endl;
     //if(d.neve[0]>0 and d.neve[1]>0 and d.neve[2]>0 and d.neve[3]>0  and d.neve[4]>0 and d.neve[5]>0) flagt=1;
     
     }while(d.N_det<2500);//// شرط ما باید تعیین کند که در هر بازه ی زمانی یک رویداد شبیه سازی شده باشد.  
     fclose(fil1); 
     //fclose(fil2);


     TEE=fopen("./files/tef_last_finalb.txt","w"); 
     for(int j=0;j<(GG+1); ++j){
     d.effd[j]=double(d.effd[j]/d.effs[j]);
     fprintf(TEE,"%d   %.3lf     %.2lf    %.2lf\n",j,d.te[j],d.effs_com[j],d.effd_com[j]);}     
     fclose(TEE);
     


     d.ave_tei=0.0; countd=0.0; //countt=0.0; 
     teefi=fopen("./files/tefzzz.txt","r"); 
     if(!teefi){cout<<"cannot make file teefi:  "<<endl; exit(0);}
     while(!feof(teefi)){
     fscanf(teefi,"%lf  %lf %d  %d\n",&l.tE,&we,&vf,&flag_det);
     //++countt; 
     if(flag_det>0.5 and l.tE>0.0){
     countd += we;  
     d.ave_tei+=fabs(d.effd[vf]/l.tE)*we;}}
     d.ave_tei=double(d.ave_tei/countd); 
     fclose(teefi);
   //if( int(countd)!=d.nevent){
    //cout<<"Error countt: "<<countt<<"\t countd: "<<countd<<"\t nbb: "<<nbb<<endl;   int yyr; cin>>yyr; }


    ///HHHHHHHHHHHHHHHHHHHHHHHHHH SAVING IN FILE HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    result1=fopen("./files/res120102A.txt", "a+");
    result2=fopen("./files/res120102B.txt", "a+"); 
    result3=fopen("./files/res120102C.txt", "a+"); 
    
    d.ave_te= d.ave_te/d.nevent;
    d.ave_re= d.ave_re/d.nevent;
    d.ave_Ds= d.ave_Ds/d.nevent;
    d.ave_Dl= d.ave_Dl/d.nevent;
    d.ave_vt= d.ave_vt/d.nevent;
    d.ave_vl= d.ave_vl/d.nevent;
    d.ave_vs= d.ave_vs/d.nevent;
    d.ave_ml= d.ave_ml/d.nevent;
    d.ave_u0= d.ave_u0/d.nevent;
    d.ave_opt2=d.ave_opt2/d.nevent;
    d.ave_fred=d.ave_fred*100.0/d.nevent;
    d.ave_opt=double(d.ave_opt/nbb);
    d.ave_col=double(d.ave_col/nbb); 
    d.aveI=   double(d.aveI/nbb); 
    d.ratio= 2.0*d.ave_opt*1.0e-6*d.ave_tei*year/M_PI;///[N_event/year]
    
    for(int i=1; i<Num; ++i){
    effe=double(s.nddis[i]/(s.nsdis[i]+0.0000064646745));
    s.nstart += s.Rostari[i]*(s.Nstart/s.Rostart)*effe*ddeg*ddeg;
    effe=double(Nogle[i]/(s.nsdis[i]+0.0000064646745));
    s.nstarti +=s.Rostari[i]*(s.Nstart/s.Rostart)*effe*(92.0/267.0);}
    
    d.nfin=d.ratio*10.0*s.nstart;////No.event  
    Nmicro+=d.nfin;///No. event_total
    
    
    fprintf(result1,"%.2lf %.2lf %.1lf %.1lf %.1lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %e  %e %.8lf %e %e %e %e  %e  %.5lf  %.3lf   %.5lf   %.5lf  %e\n",
    s.lat,s.lon,d.nstar,nbb,d.nevent,d.ave_te,d.ave_re,d.ave_Ds,d.ave_Dl,d.ave_vt,d.ave_vl,  
    d.ave_vs,d.ave_ml,d.ave_u0,s.Nstart,s.Rostart,d.ave_tei,d.ratio,d.ave_opt,d.nfin,Nmicro,
    s.nstart,d.ave_col,d.ave_fred,d.ave_opt2,d.aveI,s.nstarti);//27
    
 
    for(int j=0; j<M; ++j){
    d.ave_aps[j]= d.ave_aps[j]/d.nsdet[j];
    d.ave_apb[j]= d.ave_apb[j]/d.nsdet[j];
    d.ave_ex[j]=  d.ave_ex[j]/d.nsdet[j];
    d.ave_npsf[j]=d.ave_npsf[j]/d.nsdet[j];
    d.ave_apbd[j]=d.ave_apbd[j]/d.neve[j];
    d.ave_bl[j]=  d.ave_bl[j]/d.neve[j];
    d.ave_dmag[j]=d.ave_dmag[j]/d.neve[j];
    fprintf(result2,"%.2lf %.2lf %d %.6lf  %lf  %.1lf  %.2lf  %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.2lf\n",
  s.lat,s.lon,j,d.nstar,nbb,d.nsdet[j],d.neve[j],d.ave_aps[j],d.ave_apb[j],d.ave_apbd[j],d.ave_ex[j],d.ave_bl[j],d.ave_dmag[j],d.ave_npsf[j]);}//14


///==========================================================================

    sprintf(filename1,"./files/distribution/%c%.2lf%c%.2lf.dat",'w',s.lat,'_',s.lon);
    fil1=fopen(filename1,"r");
    if(!fil1){cout<<"cannot make file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl; exit(0);}
    countt=0.0;   countd=0.0;
    double SD_Dl=0.0, SD_Ds=0.0, SD_ml=0.0,  SD_tE=0.0,  SD_vt=0.0; 
    while(!feof(fil1)){
    countt+=we;  
    if(we>0.0) countd+=1.0;
    fscanf(fil1,"%lf  %d  %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d  %lf\n",
    &nbb,&s.struc,&l.struc,&l.Ml,&l.RE,&l.Dl,&s.Ds,&l.Vt,&l.tE,&l.vl,&l.vs,
    &s.logl,&s.logt,&s.Rs,&s.ro_star,&s.col,&ex.Av,&ex.AI,&s.u0,&d.t0,&s.opd,&s.cl,&flag_det,&we);///24
    SD_Dl+=pow(l.Dl-d.ave_Dl,2.0)*we; 
    SD_Ds+=pow(s.Ds-d.ave_Ds,2.0)*we; 
    SD_ml+=pow(l.Ml-d.ave_ml,2.0)*we;  
    SD_tE+=pow(l.tE-d.ave_te,2.0)*we; 
    SD_vt+=pow(l.Vt-d.ave_vt,2.0)*we;}
    SD_Dl= sqrt(SD_Dl/(countt*(countd-1.0)/countd));
    SD_Ds= sqrt(SD_Ds/(countt*(countd-1.0)/countd));
    SD_ml= sqrt(SD_ml/(countt*(countd-1.0)/countd));
    SD_tE= sqrt(SD_tE/(countt*(countd-1.0)/countd));
    SD_vt= sqrt(SD_vt/(countt*(countd-1.0)/countd));
    fclose(fil1);   
    cout<<"<Dl>: "<<d.ave_Dl<<"\t SD_Dl: "<<SD_Dl<<endl;
    cout<<"<Ds>: "<<d.ave_Ds<<"\t SD_Ds: "<<SD_Ds<<endl;
    cout<<"<ml>: "<<d.ave_ml<<"\t SD_ml: "<<SD_ml<<endl;
    cout<<"<tE>: "<<d.ave_te<<"\t SD_tE: "<<SD_tE<<endl;
    cout<<"<vt>: "<<d.ave_vt<<"\t SD_vt: "<<SD_vt<<endl;
    fprintf(result3,"%.2lf %.2lf   %.5lf  %.5lf  %.5lf  %.5lf  %.5lf\n",s.lat,s.lon,SD_Dl, SD_Ds, SD_ml, SD_tE, SD_vt);//7
    fclose(result1);  fclose(result2);  fclose(result3); 
///============================================================================
    
    cout<<" \n\n******  latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
    for (int  j=0; j<M; ++j){
    cout<<"********** LSST FILTER ***************************: "<<j<<endl; 
    cout<<"nsdet: "<<d.nsdet[j]<<endl;
    cout<<"<M_bg_detected>: "<<d.ave_apb[j]<<"\t <M_bg_lensed>: "<<d.ave_apbd[j]<<endl;
    cout<<"<Map_sources>:"<<d.ave_aps[j]<<"\t <EXtinction>: "<<d.ave_ex[j]<<endl;
    cout<<"effi star-det: "<<d.nsdet[j]/d.nstar<<"\t effi_ML: "<<d.neve[j]/d.nsdet[j]<<endl;
    cout<<"<Delta_mag>[i]: "<<d.ave_dmag[j]<<"\t <blending>: "<<d.ave_bl[j]<<endl;}
    cout<<"============================================================"<<endl;
    cout<<"<Ds>: "<<d.ave_Ds<<"\t <Dl>: "<<d.ave_Dl<<endl;
    cout<<"<vl>: "<<d.ave_vl<<"\t <vs>: "<<d.ave_vs<<endl;
    cout<<"<Ml>: "<<d.ave_ml<<"\t <u0>: "<<d.ave_u0<<endl;
    cout<<"<Vt>: "<<d.ave_vt<<"\t<tE>: "<<d.ave_te<<"\t<color>: "<<d.ave_col<<endl;
    cout<<"blending[2]: "<<d.ave_bl[2]<<"\t apb[2]: "<<d.ave_apb[2]<<endl;
    cout<<"<RE>[AU]: "<<d.ave_re<<"\t effi_Microlensing: "<<d.nevent/nbb<<endl;
    cout<<"simulated stars: "<<d.nstar<<"\t detected events: "<<d.nevent<<endl;
    cout<<"No. lensing event[l,b]: "<<d.nfin<<"\t n_start[Omega_l] "<<s.nstart<<"\t nstart(OGLE): "<<s.nstarti<<endl;
    cout<<"<opt>: "<<d.ave_opt<<"\t <opt2>: "<<d.ave_opt2<<"\t event_ra: "<<d.ratio<<endl;
    cout<<">>>>>>>>>>>>>> No. microlensing total: "<<Nmicro<<endl;
    cout<<"fraction_red_gients[%]: "<<d.ave_fred<<"\t <I-band>: "<<d.aveI<<endl;
    cout<<"cadence: "<<cadence<<"\t nbb: "<<nbb<<endl;
    cout<<"================================================================================"<<endl;
     }//end of if extinction
     else cout<<"Extinction does not exit: lat: "<<s.lat<<"\t lon: "<<s.lon<<endl;
     }}///end of subfield
     fclose(_randStream);
     return(0);
}
///==========================================================================//
///                                                                          //
///                   Lensing                                                //
///                                                                          //
///==========================================================================//
void lensing(source & s, lens & l,detection & d,int ei)
{
    int sign, f;
    char filenam1[40];
    d.peak=0;
    int uue;
    double As,u,Delta1,seei,prob1,prob2,ggh,dela;
    double maga,maga2,x;
    double m5[M]={0.0};
    int ds[M]={0};
    double chi2,chi1;
    d.t0=double(fabs((double)rand()/(double)(RAND_MAX+1.)*Tobs));
    d.tmin=d.t0-3.5*l.tE;
    d.tmax=d.t0+3.5*l.tE;
    if(d.tmin<0.0)   d.tmin=0.0;
    if(d.tmax>=Tobs) d.tmax=Tobs;
    if(((d.tmin/year)-int(d.tmin/year))>double(7.0/12.0)) d.tmin=double(int(d.tmin/year)+1.0)*year;
    
    int Nb=int(1.0+(d.tmax-d.tmin)/cadence);
    if(d.tmin<=d.t0 and d.tmax>=d.t0)  d.peak=1;///it means we do see just one side of the light curve

  // cout<<"No. Data: "<<Nb<<"\t tE: "<<l.tE<<"\t t0: "<<d.t0<<"\t u0: "<<s.u0<<endl;
  // cout<<"tmin: "<<d.tmin<<"\t tmax: "<<d.tmax<<"\t ro_*: "<<s.ro_star<<endl;

    if(Nb<=4){
    d.det=0;
    chi1=chi2=d.dchi=0.0;}
    


    else{
    int flag[Nb];
    int numt=0;
    for(int i=0; i<M; ++i) ds[i]=0;
    d.det=0; chi1=chi2=d.dchi=0.0;
    for(int i=0;i<Nb; ++i) flag[i]=0;
    int randd=int(fabs((double)rand()/(double)(RAND_MAX+1.)*NM));
    if(randd>=(NM-1)) randd=0;
    
    for(d.t=d.tmin; d.t<=d.tmax; d.t=d.t+cadence){
    //cout<<"=================================="<<endl;
    randd+=1;
    if(randd>=(20000-1)) randd=0;
    prob1=double(fabs((double)rand()/(double)(RAND_MAX+1.))*100.0);///probability of bad weather
    prob2=double(fabs((double)rand()/(double)(RAND_MAX+1.))*100.0);///probability of uniform observation
    if(prob1>17.0 and prob2>=10.0){
    
    u=sqrt(s.u0*s.u0+(d.t-d.t0)*(d.t-d.t0)/l.tE/l.tE);
    if(u==0.0) u=1.0e-50;
    As=fabs((u*u+2.0)/(u*sqrt(u*u+4.0)));
    if(As<1.0 or As<0.0){cout<<"ERROR AS: "<<As<<"\t u: "<<u<<"\t R:"<<s.ro_star<<endl; cin>>uue; As=1.0;}
    
    f=d.filter[randd];
    if(f<0 or f>5){cout<<"Big error filter: "<<f<<"\t randd: "<<randd<<endl;  int yye; cin>>yye;}
    if(f<LL){////ugrizy
    maga=-2.5*log10(fabs(pow(10.0,-0.4*s.Map[f])*(As-1.0)+s.Fluxb[f]));
 //cout<<"maga: "<<maga<<"\t Map[f]: "<<s.Map[f]<<"\t Fluxb: "<<s.Fluxb[f]<<endl;
 //cout<<"saturation: "<<satu[f]<<"\t detection: "<<detect[f]<<endl;
    if(maga>=satu[f] and maga<=detect[f]){
    //dels=RandN(0.2,1);///amount of flactuation in seeing
   // dela=RandN(1.0,2);///the amount of flactuation in magnification factor  
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
 /*  if(numt==0 and int(ei)%1000==0){
    FILE *test; 
    sprintf(filenam1,"./files/light/%c%.2lf%c%.2lf%c%d.dat",'l',s.lat,'_',s.lon,'_',ei);
    test=fopen(filenam1,"w");
    fclose(test);}*/
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   // seei=double(seeing[f]+RandN(0.2,1));
   // if(seei<=0.4) seei=0.4; 
    ///m5[f]=Cm[f]+d.Dcm[f]+0.5*(msky[f]-21.0)+2.5*log10(fabs(0.7/seei))+1.25*log10(fabs(texp/30.0))-km[f]*(d.airm[randd]-1.0);
   // cout<<"seeing: "<<seei<<"\t airm: "<<d.airm[randd]<<"\t m5: "<<m5[f]<<"\t msky: "<<msky[f]<<endl;
    x=pow(10.0,+0.4*(maga-d.m5c[randd]));
    Delta1=sqrt(fabs((0.04-gama[f])*x+gama[f]*x*x));///random photometric error
    if(((0.04-gama[f])*x+gama[f]*x*x)<0.0){
    cout<<"BIG ERROR negative: "<<((0.04-gama[f])*x+gama[f]*x*x)<<endl; int yye; cin>>yye; }
    d.deltam[f]=sqrt(Delta2*Delta2+Delta1*Delta1);
    //cout<<"Delta1: "<<Delta1<<"\t x: "<<x<<"\t deltam: "<<d.deltam[f]<<endl;

    dela=fabs(RandN(1.0,1));
    if(numt%2==0) sign=-1;
    else          sign=+1; 
    maga2=maga+sign*dela*d.deltam[f];
    chi1 += pow((maga2-s.magb[f])/d.deltam[f],2.0);
    chi2 += pow((maga2-maga )/d.deltam[f],2.0);///it means that fitting is done by considering blending
    if(d.det==0){
    if(fabs(maga2-s.magb[f])<5.0*d.deltam[f]) flag[numt]=0;
    else{
    if(maga2<s.magb[f])  flag[numt]=+1;//upper 
    else                 flag[numt]=-1;}//downer
    if(numt>=3 and abs(flag[numt]+flag[numt-1]+flag[numt-2]+flag[numt-3])==4) d.det=1;}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
  /*if(int(ei)%1000==0){
    FILE *test; 
    sprintf(filenam1,"./files/light/%c%.2lf%c%.2lf%c%d.dat",'l',s.lat,'_',s.lon,'_',ei);
    test=fopen(filenam1,"a+");
    m5[2]=Cm[2]+d.Dcm[2]+0.5*(msky[2]-21.0)+d.m5c[randd]-Cm[f]-d.Dcm[f]-0.5*(msky[f]-21.0); 
    double magg=-2.5*log10(fabs(pow(10.0,-0.4*s.Map[2])*(As-1.0)+s.Fluxb[2]));
    x=pow(10.0,+0.4*(magg-m5[2]));
    Delta1=sqrt((0.04-gama[2])*x+gama[2]*x*x);///random photometric error
    d.deltam[2]=sqrt(Delta2*Delta2+Delta1*Delta1);
    fprintf(test,"%.3lf    %.3lf   %.4lf  %.4lf   %d\n",d.t/year,(d.t-d.t0)/l.tE,ggh,d.deltam[2],f);  
    fclose(test);}*/
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    ++numt;
    ++ds[f];}}}///end of probability of good weather
    if(((d.t/year)-int(d.t/year))>double(7.0/12.0)){d.t= double(int(d.t/year)+1.0)*year;
    } }///end of for time
   // cout<<"********** END OF time loop ************"<<endl;
    d.dchi= fabs(chi2-chi1);
   // cout<<"dchi: "<<d.dchi<<"\t chi2: "<<chi2<<"\t chi1: "<<chi1<<endl;


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   /* if(int(ei)%1000==0){
    FILE *test; 
    sprintf(filenam1,"./files/light/%c%.2lf%c%.2lf%c%d.dat",'l',s.lat,'_',s.lon,'_',ei);
    test=fopen(filenam1,"a+");
    fprintf(test,"%.8lf   %.5lf  %.5lf  %.5lf %.5lf\n",s.u0,d.t0,l.tE,s.Map[2],s.magb[2]);
    for(int id=0; id<M; ++id)
    fprintf(test,"%d  %lf  %d  %d   %d\n",id,d.dchi,d.det,d.peak,ds[id]);
    fclose(test); }  */
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    }///end of else 
   // cout<<"*********** End of lensig function ***************"<<endl;
  
    for (int i=0; i<M; ++i){
    x=pow(10.0,+0.4*(s.magb[i]-M50[i]));
    Delta1=sqrt(fabs((0.04-gama[i])*x+gama[i]*x*x));///random photometric error
    d.deltam[i]=sqrt(Delta2*Delta2+Delta1*Delta1);}
}
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex)
{
  double F=-1.0;
  if(ds<ex.dis[0])       F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{ 
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] and ds<ex.dis[i+1]){
  F = ex.Extks[i]+(double)(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0||F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; int uut; cin>>uut;}
  return(F);
}
///==================================================================
double RandN(double sigma,int nn){
   double rr,f,frand;
   do{
   rr=double(((double)rand()/(double)(RAND_MAX+1.))*2.0-1.0)*sigma*nn; ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/(sigma*sigma));
   frand=fabs((double)rand()/((double)(RAND_MAX+1.))*1.0);
   }while(frand>f);
   return(rr);
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm){
////mass log10(T) log10(age) log10(L) log10(g) metal U G R I Z Y Bj Vj Rj Ij CL TYP
    int yye;     int h,g,k1,k2,uui; 
    double mass, age, gravity,metal,Bjc,Vjc,Rjc,type; 
    char filename[40];
    FILE *fp2;

    double Age[YZ]={0.0}; double B[YZ]={0.0};  double M[YZ]={0.0};   double mm[YZ]={0.0}; 
    int number[70]={0};   int count[70]={0};   double Metal[70]={0.0}; 


    FILE *meta; 
    meta=fopen("./files/metal.txt","r"); 
    for(int i=0; i<70; ++i) {
    fscanf(meta,"%lf   %d  %d\n",&Metal[i],&count[i],&number[i]);    
    if((Metal[i]<Metal[i-1] and i>0) or Metal[i]<0.0 or number[i]==0 or count[i]<0 or count[i]>YZ) {
    cout<<"ERROR Metal[i]: "<<Metal[i]<<"\t count[i]: "<<count[i]<<"\t number[i]: "<<number[i]<<"\t i: "<<i<<endl; cin>>uui;} }
    fclose(meta); 
    
    FILE *yz; 
    yz=fopen("./files/yzma.txt", "r"); 
    for(int i=0; i<YZ; ++i){
    fscanf(yz,"%lf   %lf   %lf  %lf\n",&Age[i],&mm[i],&B[i],&M[i]); 
    if(Age[i]<0.0 or mm[i]<0.0 or fabs(B[i])>1.5 or M[i]<0.5 or Age[i]>18.0) {   
    cout<<"ERROR Age: "<<Age[i]<<"\t metal: "<<mm[i]<<"\t B[i]"<<B[i]<<"\t M[i]: "<<M[i]<<"\t i: "<<i<<endl; cin>>uui; }}
    fclose(yz); 




////=================================== THIN DISK ============================== 
    int j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c%c.dat",'C','M','D','T','i','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTi.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %d %lf\n",
    &mass,&cm.logt_d[j],&age,&cm.logl_d[j],&gravity,&metal,&cm.Mab_d[0][j],
    &cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.Mab_d[3][j],&cm.Mab_d[4][j],&Bjc,&Vjc,&Rjc,&cm.MI_d[j],&cm.cl_d[j],&type);

///*******************************************************
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0]) h=0; 
    else if(metal>Metal[69] or  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 or metal==Metal[i-1]) { h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui; } 
    //cout<<"metal: "<<metal<<"\t MEtal[h]: "<<Metal[h]<<"\t Metal[h+1]: "<<Metal[h+1]<<"\t h: "<<h<<endl;
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>3578 or mm[k1]!=mm[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else {
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl;
    int ooe; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  or  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>3577 or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_d[5][j]= double(B[g]+M[g]*cm.Mab_d[4][j]);  

    if(fabs(cm.Mab_d[5][j]-cm.Mab_d[4][j])>1.5 or fabs(age-Age[g])>3.0 or cm.Mab_d[5][j]==0.0){
    cout<<"ERROR: Mab_d(y-band): "<<cm.Mab_d[5][j]<<"\t Mab_d(z-band): "<<cm.Mab_d[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B[g]<<"\t M[g]: "<<M[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
///********************************************************
    cm.col_d[j]=Vjc-cm.MI_d[j]; 
    if(mass<0.0||cm.logt_d[j]<0.0||cm.Mab_d[2][j]>20.0||metal>0.12||age>10.0||cm.cl_d[j]>7||type>9.0){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N1){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
   // cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;






////=================================== BULGE ================================== 
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c%d.dat",'C','M','D','B','b',2);
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDB.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %d %lf\n",
    &mass,&cm.logt_b[j],&age,&cm.logl_b[j],&gravity,&metal,&cm.Mab_b[0][j],
    &cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.Mab_b[3][j],&cm.Mab_b[4][j],&Bjc,&Vjc,&Rjc,&cm.MI_b[j],&cm.cl_b[j],&type);
   ///cm.Mab_b[5][j]=0.0; 
///*******************************************
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0]) h=0; 
    else if(metal>Metal[69] or metal==Metal[69]) h=69;
    else { 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 or metal==Metal[i-1]) { h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui; } 
    //cout<<"metal: "<<metal<<"\t MEtal[h]: "<<Metal[h]<<"\t Metal[h+1]: "<<Metal[h+1]<<"\t h: "<<h<<endl;
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>3578 or mm[k1]!=mm[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else {
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl;
    int ooe; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  or  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>3577 or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_b[5][j]= double(B[g]+M[g]*cm.Mab_b[4][j]);  

    if(fabs(cm.Mab_b[5][j]-cm.Mab_b[4][j])>1.5 or fabs(age-Age[g])>3.0  or cm.Mab_b[5][j]==0.0){
    cout<<"ERROR:   Mab_b(y-band): "<<cm.Mab_b[5][j]<<"\t Mab_b(z-band): "<<cm.Mab_b[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B[g]<<"\t M[g]: "<<M[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
///***********************************************
    cm.col_b[j]=Vjc-cm.MI_b[j]; 
    if(mass<0.0||cm.logt_b[j]<0.0||Vjc>18.0||age>10.0||metal>0.9||cm.cl_b[j]>7||type>9.0){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
   // cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;









////=================================== THICK DISK =============================
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c%c.dat",'C','M','D','T','k','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTk.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %d %lf \n",
    &mass,&cm.logt_t[j],&age,&cm.logl_t[j],&gravity,&metal,&cm.Mab_t[0][j],
    &cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.Mab_t[3][j],&cm.Mab_t[4][j],&Bjc,&Vjc,&Rjc,&cm.MI_t[j],&cm.cl_t[j],&type);
///************************************************
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0]) h=0; 
    else if(metal>Metal[69] or  metal==Metal[69]) h=69;
    else { 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 or metal==Metal[i-1]) { h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui; } 
    //cout<<"metal: "<<metal<<"\t MEtal[h]: "<<Metal[h]<<"\t Metal[h+1]: "<<Metal[h+1]<<"\t h: "<<h<<endl;
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>3578 or mm[k1]!=mm[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else {
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl;
    int ooe; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  or  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>3577 or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_t[5][j]=double(B[g]+M[g]*cm.Mab_t[4][j]);  

    if(fabs(cm.Mab_t[5][j]-cm.Mab_t[4][j])>1.5 or fabs(age-Age[g])>3.0 or cm.Mab_t[5][j]==0.0){
    cout<<"ERROR:   Mab_t(y-band): "<<cm.Mab_t[5][j]<<"\t Mab_t(z-band): "<<cm.Mab_t[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B[g]<<"\t M[g]: "<<M[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
///*************************************************
    cm.col_t[j]=Vjc-cm.MI_t[j]; 
    if(mass<0.0||cm.logt_t[j]<0.0||cm.Mab_t[2][j]>20.0||metal>0.025||cm.cl_t[j]>7|| type>9.0){
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
   // cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;






////=================================== STELLAR HALO =========================== 
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c.dat",'C','M','D','H','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDH.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &mass,&cm.logt_h[j],&age,&cm.logl_h[j],&gravity,&metal,&cm.Mab_h[0][j],
    &cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.Mab_h[3][j],&cm.Mab_h[4][j],&Bjc,&Vjc,&Rjc,&cm.MI_h[j],&cm.cl_h[j],&type);
///****************************************************************************
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0]) h=0; 
    else if(metal>Metal[69] or  metal==Metal[69]) h=69;
    else { 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 or metal==Metal[i-1]) { h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui; } 
    //cout<<"metal: "<<metal<<"\t MEtal[h]: "<<Metal[h]<<"\t Metal[h+1]: "<<Metal[h+1]<<"\t h: "<<h<<endl;
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>3578 or mm[k1]!=mm[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else {
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl;
    int ooe; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  or  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>3577 or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_h[5][j]= double(B[g]+M[g]*cm.Mab_h[4][j]);  

    if(fabs(cm.Mab_h[5][j]-cm.Mab_h[4][j])>1.5 or fabs(age-Age[g])>3.0 or cm.Mab_h[5][j]==0.0){
    cout<<"ERROR:   Mab_h(y-band): "<<cm.Mab_h[5][j]<<"\t Mab_h(z-band): "<<cm.Mab_h[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B[g]<<"\t M[g]: "<<M[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
///*****************************************************************************
    cm.col_h[j]=Vjc-cm.MI_h[j]; 
    if(mass<0.0 or cm.logt_h[j]<0.0 or cm.Mab_h[2][j]>20.0 or metal>0.01 or cm.cl_h[j]>7|| type>9.0){
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
  //  cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;


}
///==============================================================//
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*M_sun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opd=0.0;
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;}
    s.opd= s.od_disk+s.od_ThD+s.od_bulge+s.od_halo;///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD & cm, extinc & ex)
{
    int num,struc,nums;
    double rho,rf,ratio,N_sblend[M+1];
    double Akv,Avk,Avks,Alv;
    double Ds,Ai[M],Av;
    double Map[M];      
    double extI; 
    double MIsource, Isource;   


    double maxnb=0.0;
    for(int i=0; i<(M+1); ++i){
    s.Fluxb[i]=0.0;
    N_sblend[i]=s.Nstart*pow(FWHM[i]*0.4,2)*M_PI/(3600.0*3600.0);
    N_sblend[i]=N_sblend[i]+RandN(sqrt(N_sblend[i]),1);
    if(N_sblend[i]<=1.0) N_sblend[i]=1.0;
    if(N_sblend[i]<1.0){cout<<"ERROR N_sblend+ deltaN: "<<N_sblend[i]<<endl; int yye; cin>>yye;}
    if(N_sblend[i]>maxnb) maxnb=N_sblend[i];}
    s.FIsource=0.0,s.Isource=0.0;
    double Mab[M]={0.0};   



    for(int k=1; k<=int(maxnb); ++k){
    
    do{
    num=int(fabs((double)rand()/(double)(RAND_MAX+1.)*Num*1.0));
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    }while(rho>s.Rostari[num] or num<5);///distance larger than 50.0 
    Ds=(double)(num*step);///in kpc
    //if(s.Ds>MaxD or s.Ds<0.1){cout<<"ERROR (1): Ds: "<<s.Ds<<"\t MaxD: "<<MaxD<<endl;  int yye; cin>>yye;}
    nums=num;
    //cout<<"Ds: "<<s.Ds<<"\t nums: "<<s.nums<<endl; 



    rf=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[nums];
         if (rf<= s.rho_disk[nums]) struc=0;///thin disk
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
    else if (rf<=s.Rostar0[nums]) struc=3;///halo
   // cout<<"struc: "<<s.struc<<endl;


    if(struc==0){///thin disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    if(k==1){ 
    s.logt=cm.logt_d[num]; 
    s.col=  cm.col_d[num]; 
    s.logl=cm.logl_d[num];
    s.cl= cm.cl_d[num]; }
    MIsource=cm.MI_d[num];}


    if(struc==1){///bulge
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num];}//cout<<"Mab[i]: "<<s.Mab[i]<<endl; }
    if(k==1){ 
    s.logt=cm.logt_b[num];
    s.col=  cm.col_b[num]; 
    s.logl=cm.logl_b[num];
    s.cl= cm.cl_b[num];}
    MIsource=cm.MI_b[num];}


    if(struc==2){///thick disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num]; }//cout<<"Mab[i]: "<<s.Mab[i]<<endl;}
    if(k==1){ 
    s.logt=cm.logt_t[num];
    s.col=  cm.col_t[num]; 
    s.logl=cm.logl_t[num];
    s.cl= cm.cl_t[num];}
    MIsource=cm.MI_t[num];}


    if(struc==3){///stellar halo
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-1.0));
    for(int i=0; i<M; ++i) { Mab[i]=cm.Mab_h[i][num];}//cout<<"Mab[i]: "<<s.Mab[i]<<endl;}
    if(k==1){ 
    s.logt=cm.logt_h[num];
    s.col=  cm.col_h[num]; 
    s.logl=cm.logl_h[num];
    s.cl= cm.cl_h[num]; }
    MIsource=cm.MI_h[num];}


    ex.Aks=Interpol(Ds,ex)+0.00033435136;///extinction in Ks-band
    Akv=fabs(ex.a[M]+ex.b[M]/Rv[struc] ) ;///Cardelli et al. 1989 A(k-band)/A(v-band)
    Avk=1.0/Akv; ////A(v)/A(k)
    Avks=Avk*0.95;///Marshal 2006 A(v)/A(ks)= Av/Ak Ak/Aks  [8, 9]
    Av=double(ex.Aks*Avks);
    if(Avks>10.0 or Avks<7.5 or ex.Aks<0.0 or ex.Aks==0.0 or Av<0.0 or Av==0.0 ){
    cout<<"Error(1) Avks:  "<<Avks<<"\t Av:  "<<Av<<"\t Ds:  "<<Ds<<"\t Aks:  "<<ex.Aks<<endl;  int iis;   cin>>iis; }
   
    
    
    extI=-3.0; 
    Alv=fabs(ex.a[M+1]+ex.b[M+1]/Rv[struc]);///AIc/AV
    extI=double(Av*Alv+RandN(sigma[M+1],1) ); //extinction in Ic-band
    if(extI<0.0) extI=0.0;
    if(Av<0.0 or Av==0.0 or extI<0.0 or fabs(Alv-0.6)>0.1){
    cout<<"Error(2) Av:  "<<Av<<"\t extI:   "<<Av*Alv<<"\t Alv:  "<<Alv<<"\t Ds:   "<<Ds<<endl;  int uue;  cin>>uue;}



 
    for(int i=0; i<M; ++i){
    Alv=fabs(ex.a[i]+ex.b[i]/Rv[struc]);///Alambda/AV
    Ai[i]= Av*Alv+RandN(sigma[i],1); //extinction in other bands
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i]+5.0*log10(Ds*100.0)+Ai[i];
    Isource= MIsource+5.0*log10(Ds*100.0)+extI; 
    if(N_sblend[i]>=k) s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i] ));
    
    if(Alv<=0.0 or Av<0.0 or Ai[i]<0.0 or Map[i]<Mab[i] or Isource<MIsource or extI<0.0 or s.Fluxb[i]<0.0 or N_sblend[i]<1.0){
    cout<<"ERROR(3) Ds: "<<Ds<<"\t Alv: "<<Alv<<"\t Av: "<<Av<<"\t i:  "<<i<<endl;
    cout<<"A[i]: "<<Ai[i]<<"\t extI: "<<extI<<"\t Nblend: "<<N_sblend[i]<<endl; 
    cout<<"Map[i]: "<<Map[i]<<"\t Mab[i]: "<<Mab[i]<<"\t Fluxb:  "<<s.Fluxb[i]<<endl;
    cout<<"MI_abs: "<<MIsource<<"\t MI_app:  "<<Isource<<endl;
    int yue; cin>>yue;} } 
    if(N_sblend[M]>=k) s.FIsource+=fabs(pow(10.0,-0.4*Isource));/// I-band Johnson standard filters


    if(k==1){
    s.struc=struc;  
    s.Ds=Ds;  
    s.nums=nums;
    ex.AI=extI;
    ex.Av=Av; 
    s.MIsource= MIsource; 
    s.Isource= Isource; 
    s.col= s.col+ex.Av-ex.AI;
    for(int i=0; i<M; ++i){ex.Ai[i]=Ai[i];   s.Map[i]=Map[i];   s.Mab[i]=Mab[i];}}
    }///loop over the stars 


    s.blendI=double(pow(10.0,-0.4*s.Isource)/s.FIsource);
    for(int i=0; i<M; ++i){
    if(s.Fluxb[i]<=0.0){cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<endl; int yye; cin>>yye; }
    s.magb[i]=-2.5*log10(fabs(s.Fluxb[i]));
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    s.nsbl[i]=double(N_sblend[i]);
    if(int(s.nsbl[i])<0.0 or s.nsbl[i]==0.0 or s.blend[i]>1.00001 or s.blend[i]<=0.00 or 
    (s.nsbl[i]==1.0 and s.blend[i]<1.0) or s.blendI>1.0 or s.blendI<0.0 or 
    fabs(s.Map[0]-s.Mab[0]-ex.Ai[0]-s.Map[1]+s.Mab[1]+ex.Ai[1])>0.1){ 
    cout<<"N_sblend: "<<N_sblend[M]<<endl; 
    cout<<"blendI: "<<s.blendI<<"\t FIsource: "<<s.FIsource<<"\t Isource: "<<s.Isource<<endl;
    cout<<"BIGG ERRROR nsbl: "<<s.nsbl[i]<<"\t N_sblend: "<<N_sblend[i]<<"\t blend: "<<s.blend[i]<<endl; int uue; cin>>uue;}}


    s.u0=fabs((double)rand()/(1.0 + RAND_MAX));
    //double signs=fabs((double)rand()/(1.0 + RAND_MAX));
    //if(signs<0.5) s.u0=-1.0*s.u0;
    if(s.u0==0.0) s.u0=1.0e-50; 
    s.magni=fabs((s.u0*s.u0+2.0)/(s.u0*sqrt(s.u0*s.u0+4.0)));
    s.Rs=s.logl-4.0*(s.logt-logT_sun);
    s.Rs=pow(10.0,s.Rs/2.0);/// in Sun radius 
}
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 or delt<0.0) {cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     if(s.lon==0.0)     Lon=360.00;

     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 or delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     if(Lon>360.000 or Lon<0.25 or fabs(Lat)>10.0 or (Lon>100 and Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}


     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");
     cout<<"Extnction file name:  Lat:  "<<Lat<<"\t Lon:    "<<Lon<<endl;


     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<64881; ++i) {
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 and fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     //cout<<"distance: "<<ex.dis[i]<<"\t Extks: "<<ex.Extks[i]<<endl;
     if(ex.dis[i]<0.2  or ex.dis[i]>50.0 or ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; }
     if(ex.dis[i]==8.25) ex.exks8= ex.Extks[i]; 
     }}
     //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<endl;
     fclose(fpd);
     return(flag);
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s, dispersion & di)
{
    double f,test;
    double rholens[s.nums+2];
    l.rhomaxl=0.0;

    for(int k=1;k<=s.nums;++k){
    rholens[k]=0.0;
    l.Dl=k*step;
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds) {cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  int yye; cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}

    do{
    l.numl = (int)((double)rand()*1000./((double)(RAND_MAX+1.)*1000.)*(s.nums-5.0)+1.0);
    test = ((double)rand()/(double)(RAND_MAX+1.)*l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;
     int ue; cin>>ue;}
    }while(test>rholens[l.numl]);
   // cout<<"numl: "<<l.numl<<"\t rhomaxl: "<<l.rhomaxl<<endl; 


   double  randflag=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[l.numl];
       if (randflag<=s.rho_disk[l.numl]) l.struc=0;///thin disk
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1; // bulge structure
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2; //thick disk
  else if (randflag<= s.Rostar0[l.numl]) l.struc=3;//halo
  else {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}
 //  cout<<"struc(lens): "<<l.struc<<endl;


  if(l.struc==0){///thin disk
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(4.5-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*57.0);
  if(l.Ml<=1.0) f=pow(l.Ml,-1.6);
  if(l.Ml>1.0) f=pow(l.Ml,-3.0);
  }while(test>f);}


  if(l.struc==1){///Galactic bulge
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(1.4-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*378.0)+0.3;
  f=pow(l.Ml,-2.35);
  }while(test>f);}


  if(l.struc==2){///thick disk
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(1.4-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*3.0)+0.8;
  f=pow(l.Ml,-0.5);
  }while(test>f);}


  if(l.struc==3){///stellar halo
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(0.8-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*1.0)+0.8;
  f=pow(l.Ml,-0.5);
  }while(test>f);}
 // cout<<"lens_ mass: "<<l.Ml<<endl;


    l.Dl=l.numl*step;///kpc
    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*M_sun*s.Ds*KP)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s,l, di);
    l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day


    if(l.tE<0.0 or  l.tE==0.0 or l.RE<0.0 or  l.Dl>s.Ds  or l.xls>=1.0 or l.Dl<0.0 or s.Ds<0.0 or s.Ds>20.0 or l.Ml<0.0 or l.Ml>5.0){ 
    cout<<"BIG ERROR te: "<<l.tE<<"\t RE: "<<l.RE<<"\t V_rel: "<<l.Vt<<endl;l.tE=0.1; int iie; cin>>iie;} 
    s.ro_star=s.Rs*Rsun*l.xls/l.RE; 
   // cout<<"Dl: "<<l.Dl<<"\t xls: "<<l.xls<<"\t rho_star: "<<s.ro_star<<endl;
  //  cout<<"RE: "<<l.RE/AU<<"\t tE: "<<l.tE<<"\t Vt: "<<l.Vt<<endl;
  //  cout<<"************** End of lens function  ***************"<<endl;
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l, dispersion & di)
{
  
  if (l.Dl==0.0) l.Dl=0.00034735;
  double pi=M_PI;
  double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*l.Dl*cos(s.TET)*cos(s.FI));
  double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*s.Ds*cos(s.TET)*cos(s.FI));
  if(Rlc==0.0) Rlc=0.0000000000034346123;
  if(Rsc==0.0) Rsc=0.000000000004762654134; 
 
  double LVx, SVx;
  double SVT, SVR, SVZ, LVT, LVR, LVZ;
  double fv, testfv, test, age;
  double  VSunx, vls2, vls1;
  double betal, betas, deltal, deltas, tetd ;
  double zl, zs;  
  int gd; 


  double NN=2.5;
  double sigma_R_Disk, sigma_T_Disk, sigma_Z_Disk;
  double sigma_R_DiskL,  sigma_T_DiskL, sigma_Z_DiskL;
  double sigma_R_DiskS,  sigma_T_DiskS, sigma_Z_DiskS;
  double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
  double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
  double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
  double Rho[8]={00.0}; double maxr=0.0;
  
  
  if(s.struc==0 or l.struc==0){
  sigma_R_DiskL=0.0;   sigma_T_DiskL=0.0;  sigma_Z_DiskL=0.0;
  sigma_R_DiskS=0.0;   sigma_T_DiskS=0.0;  sigma_Z_DiskS=0.0;
  
  for(int i=0; i<8; ++i){Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}
  for(int i=0;i<2; ++i){
  test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])                       {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
else if(test<=(Rho[0]+Rho[1]))              {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]))       {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3])){sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4])){sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5])){sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
else if(test<=maxr)                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}
    if(i==0){ sigma_R_DiskS= sigma_R_Disk; sigma_T_DiskS= sigma_T_Disk; sigma_Z_DiskS= sigma_Z_Disk;}
    if(i==1){ sigma_R_DiskL= sigma_R_Disk; sigma_T_DiskL= sigma_T_Disk; sigma_Z_DiskL= sigma_Z_Disk;}
///===========================================================
   gd=-1; 
   if(age<di.age[0] or age==di.age[0])  gd=0;
   else if(age>di.age[Nag-1])           gd=Nag-1;  
   else{
   for(int j=1; j<Nag; ++j){
   if(double((age-di.age[j-1])*(age-di.age[j]))<0.0 or age==di.age[j]){gd=j;  break; }}} 
   if(gd<0 or fabs(age-di.age[gd])>1.5){
   cout<<"Error gd<0:  "<<gd<<"\t age:  "<<age<<endl; int uue;  cin>>uue;
   cout<<"age:  "<<age<<"\t di.age[gd]:  "<<di.age[gd]<<"\t gd:  "<<gd<<endl;}
///=============================================================
    zs= fabs(s.Ds*sin(s.FI)*1000.0);//parsec
    if(zs>200.0) zs=100.0;  
    zl= fabs(l.Dl*sin(s.FI)*1000.0);//parsec
    if(zl>200.0) zl=100.0; 
    if(i==0 and zs<200.0){
    sigma_R_DiskS= double(di.sigR0[gd] +  di.sigR1[gd] * zs+ di.sigR2[gd] *zs*zs);
    sigma_T_DiskS= double(di.sigT0[gd] +  di.sigT1[gd] * zs+ di.sigT2[gd] *zs*zs);
    sigma_Z_DiskS= double(di.sigZ0[gd] +  di.sigZ1[gd] * zs+ di.sigZ2[gd] *zs*zs);}
    if(i==1 and zl<200.0){
    sigma_R_DiskL= double(di.sigR0[gd] +  di.sigR1[gd] * zl+ di.sigR2[gd] *zl*zl);
    sigma_T_DiskL= double(di.sigT0[gd] +  di.sigT1[gd] * zl+ di.sigT2[gd] *zl*zl);
    sigma_Z_DiskL= double(di.sigZ0[gd] +  di.sigZ1[gd] * zl+ di.sigZ2[gd] *zl*zl);}}
    
    if(sigma_R_DiskS<0.0 or sigma_T_DiskS<0.0 or sigma_Z_DiskS<0.0  or sigma_R_DiskL<0.0 or sigma_T_DiskL<0.0 or sigma_Z_DiskL<0.0 ){
    cout<<"Error(sigma_s<0.0)\t  zs: "<<zs<<"\t zl:  "<<zl<<"\t Ds:  "<<s.Ds<<"\t Dl:  "<<l.Dl<<endl;
    cout<<"sigma_R_DsiskS:  "<<sigma_R_DiskS<<"\t sigma_T_DiskS:  "<<sigma_T_DiskS<<"\t sigma_Z_DiskS:  "<<sigma_Z_DiskS<<endl;
    cout<<"sigma_R_DsiskL:  "<<sigma_R_DiskL<<"\t sigma_T_DiskL:  "<<sigma_T_DiskL<<"\t sigma_Z_DiskL:  "<<sigma_Z_DiskL<<endl;
    int uue; cin>>uue; }}
        
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_DiskS, NN);
    SVT= RandN(sigma_T_DiskS, NN);
    SVZ= RandN(sigma_Z_DiskS, NN);}

    else if(s.struc==1){///Galactic bulge
    SVR= RandN(sigma_R_Bulge, NN);
    SVT= RandN(sigma_T_Bulge, NN);
    SVZ= RandN(sigma_Z_Bulge, NN);}

    else if(s.struc==2){///thick disk
    SVR= RandN(sigma_R_TDisk, NN);
    SVT= RandN(sigma_T_TDisk, NN);
    SVZ= RandN(sigma_Z_TDisk, NN);}

    else if(s.struc==3){///stellar halo
    SVR= RandN(sigma_R_halo, NN);
    SVT= RandN(sigma_T_halo, NN);
    SVZ= RandN(sigma_Z_halo, NN);}
    if(s.struc==0 or s.struc==2)  SVT =SVT+ vro_sun*(1.00762*pow(Rsc/Dsun,0.0394) + 0.00712);
    l.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
///======================================================================================
 if(l.struc==0){///Galactic disk
   LVR= RandN(sigma_R_DiskL, NN) ;
   LVT= RandN(sigma_T_DiskL, NN) ;
   LVZ= RandN(sigma_Z_DiskL, NN) ;}

   else if(l.struc==1){///Galactic bulge
   LVR= RandN(sigma_R_Bulge, NN) ;
   LVT= RandN(sigma_T_Bulge, NN) ;
   LVZ= RandN(sigma_Z_Bulge, NN) ;}

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN) ;
   LVT= RandN(sigma_T_TDisk, NN) ;
   LVZ= RandN(sigma_Z_TDisk, NN) ;}

   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN);
   LVT= RandN(sigma_T_halo, NN);
   LVZ= RandN(sigma_Z_halo, NN);}
   if(l.struc==0 or l.struc==2)  LVT = LVT+vro_sun *(1.00762*pow(Rlc/Dsun,0.0394) + 0.00712);
   l.vl=sqrt( LVT*LVT + LVZ*LVZ + LVR*LVR );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    betal=betas=0.0;
    tetd= s.TET;
    test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
    if(fabs(test-1.0)<0.01)       betal= pi/2.0;
    else if(fabs(test+1.0)<0.01)  betal=-pi/2.0;
    else                          betal=asin(test);
    
    test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
    if( fabs(test-1.0)<0.01)     betas=pi/2.0;
    else if(fabs(test+1.0)<0.01) betas=-pi/2.0;
    else                         betas=asin(test);
    
    if(Dsun < fabs(l.Dl*cos(s.FI)*cos(tetd))) betal= pi-betal; 
    if(Dsun < fabs(s.Ds*cos(s.FI)*cos(tetd))) betas= pi-betas; 

   if(fabs(l.Dl*cos(s.FI)*sin(tetd))>double(Rlc*1.0003256463) or fabs(test)>1.0063645){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<betal<<endl;
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<betas<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;
   int ew; cin>>ew;}

///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    deltal= pi - fabs(tetd) -fabs(betal);
    deltas= pi - fabs(tetd) -fabs(betas);  
    if(betal<0.0)  deltal= -1.0*deltal;
    if(betas<0.0)  deltas= -1.0*deltas;
    s.deltao= pi-fabs(tetd);
    if(tetd<0.0)  s.deltao=-1.0*s.deltao; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    s.SV_n1 =+SVR * sin(deltas)- SVT * cos(deltas);
    s.LV_n1 =+LVR * sin(deltal)- LVT * cos(deltal);
    s.VSun_n1=+VSunR*sin(s.deltao)-VSunT*cos(s.deltao);
    
    SVx= -SVR*cos(deltas)- SVT*sin(deltas);
    LVx= -LVR*cos(deltal)- LVT*sin(deltal);
    VSunx= -VSunR*cos(s.deltao) -VSunT*sin(s.deltao);
    
    s.SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    s.LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    s.VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
 
    
    vls1= l.xls*s.SV_n1 - s.LV_n1 +(1.0-l.xls)*s.VSun_n1;  ///Source - lens 
    vls2= l.xls*s.SV_n2 - s.LV_n2 +(1.0-l.xls)*s.VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
   

    if (l.Vt<0.0 or l.Vt>1.0e6 or l.Vt==0.0){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;   int yee; cin>>yee;}
//cout<<"Vt: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;
}
///==================================================================
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nn=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   double fd=1.0; ///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;///No limitation 
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;
/*
I assume that the up-limit of the mass is indicated by the simulation. Because it depends on the structre, .... all details can be found in mass_averaged.cpp  code. */

   
  /*char filename[40];
    FILE *fill;
    sprintf(filename,"./files/density/%c%.2lf%c%.2lf.dat",'D',s.lat,'_',s.lon);
    fill=fopen(filename,"w");
    if(!fill){cout<<"cannot open file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl;  exit(0);}*/



for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;

   x=i*step;
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = Dsun-x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///M_sun/pc^3
///=================================================


///========== Galactic Thick Disk =====================
  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nn)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nn)*exp(-fabs(zb)/0.8)/(1.0+0.5*nn);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44);
   else            s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44);///M_sun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));
   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///M_sun/pc^3
///=================================================
///در اینجا اینکه تعداد ستاره ه
///ا به این بستگی دارد که ما  نمودار قدر رنگ مطلق ستاره ها را چگونه درست کرده باشیم.اگر هیچ گونه محدودیتی برای
///درست کردن آن  در نظر نگرفته ایم،  پس تعداد کل ستاره ها را نظر  میگیریم.
/// ولی بهتر است که ما رابطه بین قدر مطلق و جرم را تعیین کنیم. در این صورت می توانیم  ورودی قدر رنگ
///وارد شده به کد را خودمان محدود به ستاره های روشن کنیم تا سرعت اجرای برنامه بالارود.
///averaged mass are the same as the previous work!!! because we did not change the besancon model


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]);///[M_sun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[M_sun/deg^2]
s.Nstari[i]=binary_fraction*(s.rho_disk[i]*fd/0.403445+s.rho_ThD[i]*fh/0.4542+s.rho_halo[i]*fh/0.4542+s.rho_bulge[i]*fb/0.308571);////[Nt/pc^3] 

s.Nstari[i]= s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]

s.Nstart  +=s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
//fprintf(fill,"%e   %e   %e   %e   %e  %e   %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]);
   }
 // fclose(fill);

 // cout<<"Nstart [Nt/deg^2]: "<<s.Nstart<<"\t Ro_star [Mass/deg^2]: "<<s.Rostart<<endl;
 //cout<<">>>>>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<"<<endl;
}
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
