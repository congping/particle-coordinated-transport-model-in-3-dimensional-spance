#include <math.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>

#include "parameter.h"
double ran2(long *);
double rand_norm_1(long * p);
void rand_disk(long *p, double a[2],double R);
void directMoveRec(std::vector<double> & X, double dt, double theta, std::vector<double> Dm, std::vector<double> Dp, double ca,long * p);
void twoStepsMoveRec(std::vector<double> & X, double dt, std::vector<double> Dm, std::vector<double> Dp, double ca, long *p);
void calculateMSD(std::vector<double> & traj, std::vector<double> & MSD, double dt,double tmax);
void moveDrift(std::vector<double> &X, int ind, bool flag, double dt,long*p);
extern double step;
extern FILE * fdetach;

void moveDirectTransport(std::vector<double> &X, int ind, int & mode,long * p){
        if (X.size()>1){
             X[ind]=X[ind]+vdirect*dt*mode;
        }
        double r=ran2(p);
        if (r<omegat) mode=-mode;
        else if (r<omegat+deltat){
              mode=0;
            // printf("detach %f \n",step);
            //  fprintf(fdetach,"%f \n",step*dt);
        }
}

void moveActiveRec(std::vector<double> &X, int & mode, double theta, long *p){
    if (mode==0){
         if (direct) directMoveRec(X,dt,theta,Dm,Dp,ca,p);
         else twoStepsMoveRec(X,dt,Dm,Dp,ca,p);
         //printf("X=%f %f\n",X[0],X[1]);
         //moveDrift(X,dim-1,mode,p);
         if (X[0]<0){
             if (ran2(p)<alphat){
                 if (ran2(p)<0.5) mode=1;
                 else mode=-1;
             }
         }
    }
    else{
        moveDirectTransport(X,dim-1,mode,p);
    }
}


// for disk or cylinder regions
void moveDisk(std::vector<double> & X,double theta,double dt, long *p){
   double epsL=R0-sqrt(Dm[0]*dt)*ca;
   double epsR=R0+sqrt(Dp[0]*dt)*ca;
   double rx=sqrt(X[0]*X[0]+X[1]*X[1]);
     if (rx<epsL){
        for (int i=0;i<X.size();i++){
            double r=rand_norm_1(p)*sqrt(dt);
            X[i]=X[i]+sqrt(Dm[i])*r;  
        }        
        moveDrift(X,dim-1,true,dt,p);
     }
     else if (rx>epsR){
        for (int i=0;i<X.size();i++){
            double r=rand_norm_1(p)*sqrt(dt);
            X[i]=X[i]+sqrt(Dp[i])*r;  
        }
        moveDrift(X,dim-1,false,dt,p);
     }
     else{
        if (fabs(rx-R0)<dplanar){
             double angle;
             if (X[0]>=0) angle=asin(X[1]/rx);
             else angle=-asin(X[1]/rx)+PI;  
             
             //printf("X before rotation %f %f with angle %f\n",X[0], X[1],angle);
             //double x0=cos(angle)*X[0]+sin(angle)*X[1];
             //double x1=-sin(angle)*X[0]+cos(angle)*X[1];
            // printf("X after rotation %f %f \n",x0,x1);
             
             X[0]=rx-R0;//x0;
             X[1]=0;//x1;   
             if (direct) directMoveRec(X,dt,theta,Dm,Dp,ca,p);
             else twoStepsMoveRec(X,dt,Dm,Dp,ca,p);
             angle=-angle;
             X[0]=X[0]+R0;
             double x0=cos(angle)*X[0]+sin(angle)*X[1];
             double x1=-sin(angle)*X[0]+cos(angle)*X[1]; 
             X[0]=x0;
             X[1]=x1;          
        }
        else{
             printf("reduced dt is wrong");
        }
     }
        
   rx=sqrt(X[0]*X[0]+X[1]*X[1]);
   if (rx>R){
      for (int i=0;i<2;i++) X[i]=X[i]/rx*R*2-X[i];       
   }
}


void moveDiskFull(std::vector<double> & X, int & mode, double theta, double dt, long *p){
   double epsL=R0-sqrt(Dm[0]*dt)*ca;
   double epsR=R0+sqrt(Dp[0]*dt)*ca;
   double rx=sqrt(X[0]*X[0]+X[1]*X[1]);
   if (mode==0){
      if ((rx<epsL || rx>epsR ) && R-rx>epsR-R0) moveDisk(X,theta,dt,p);
      else{
        double dt1=dt/double(rat);
        for (int i=0;i<rat;i++)  moveDisk(X,theta,dt1,p);
     }
     if (rx<R0) {
         double r=ran2(p);
         if (r<alphat/2) mode=1;
         else if (r<alphat) mode=-1;
         if (mode!=0) {
                     // printf("attach %f r=% \n",step*dt);//,r);
                 if (flagbind)     fprintf(fdetach,"%f ",step*dt);
         }
     }   
   }
   else moveDirectTransport(X,dim-1,mode,p);
   if (X[2]<0) {
                             X[2]=-X[2];
                             mode=-mode;
   }
   if (X[2]>L) {
                             X[2]=2*L-X[2]; 
                             mode=-mode;
   }
}

void moveMT(std::vector<double> & X,double theta,double dt, std::vector<double> & center, long *p){
   double epsL=R0-sqrt(Dm[0]*dt)*ca;
   double epsR=R0+sqrt(Dp[0]*dt)*ca;
   double rx=sqrt((X[0]-center[0])*(X[0]-center[0])+(X[1]-center[1])*(X[1]-center[1]));
     if (rx<epsL){
        for (int i=0;i<X.size();i++){
            double r=rand_norm_1(p)*sqrt(dt);
            X[i]=X[i]+sqrt(Dm[i])*r;  
        }        
        moveDrift(X,dim-1,true,dt,p);
     }
     else if (rx>epsR){
        for (int i=0;i<X.size();i++){
            double r=rand_norm_1(p)*sqrt(dt);
            X[i]=X[i]+sqrt(Dp[i])*r;  
        }
        moveDrift(X,dim-1,false,dt,p);
     }
     else{
        if (fabs(rx-R0)<dplanar){
             double angle;
             if ((X[0]-center[0])>=0) angle=asin((X[1]-center[1])/rx);
             else angle=-asin((X[1]-center[1])/rx)+PI;  
                          
             /*printf("X before rotation %f %f with angle %f\n",X[0], X[1],angle);
             X[0]=X[0]-center[0];
             X[1]=X[1]-center[1];
             double x0=cos(angle)*X[0]+sin(angle)*X[1];
             double x1=-sin(angle)*X[0]+cos(angle)*X[1];
             printf("X after rotation %f %f \n",x0,x1);
             */
             X[0]=rx-R0;//x0;
             X[1]=0;//x1;   
             if (direct) directMoveRec(X,dt,theta,Dm,Dp,ca,p);
             else twoStepsMoveRec(X,dt,Dm,Dp,ca,p);
             angle=-angle;
             X[0]=X[0]+R0;
             double x0=cos(angle)*X[0]+sin(angle)*X[1];
             double x1=-sin(angle)*X[0]+cos(angle)*X[1]; 
             X[0]=x0+center[0];
             X[1]=x1+center[1];          
        }
        else{
             printf("reduced dt is wrong");
        }
     }
        
   rx=sqrt(X[0]*X[0]+X[1]*X[1]);
   if (rx>R){
      for (int i=0;i<2;i++) X[i]=X[i]/rx*R*2-X[i];       
   }
}


void moveMTFull(std::vector<double> & X, int & mode, double theta, double dt, std::vector<double> & center, long *p){
     // move MT
     for (int i=0;i<2;i++){
        double r=rand_norm_1(p)*sqrt(dt);
        center[i]=center[i]+sqrt(Dmt[i])*r;  
     }
     double rc=sqrt(center[0]*center[0]+center[1]*center[1]);
     double Rc=R-R0;
     if (rc>Rc){
         for (int i=0;i<2;i++) center[i]=center[i]/rc*Rc*2-center[i];
     }
    // printf("center %f %f\n",center[0],center[1]);
   double epsL=R0-sqrt(Dm[0]*dt)*ca;
   double epsR=R0+sqrt(Dp[0]*dt)*ca;
   double rx=sqrt((X[0]-center[0])*(X[0]-center[0])+(X[1]-center[1])*(X[1]-center[1]));
   if (mode==0){
      if ((rx<epsL || rx>epsR ))// && R-rx>epsR-R0)
            moveMT(X,theta,dt,center,p);
      else{
        double dt1=dt/double(rat);
      //  printf(" within layer X %f %f %f \n", X[0],X[1], rx);
        for (int i=0;i<rat;i++)  moveMT(X,theta,dt1,center,p);
        
     }
     if (rx<R0) {
         double r=ran2(p);
         if (r<alphat/2) mode=1;
         else if (r<alphat) mode=-1;
         if (mode!=0) {
                     // printf("attach %f r=% \n",step*dt);//,r);
                 if (flagbind)     fprintf(fdetach,"%f ",step*dt);
         }
     }   
   }
   else moveDirectTransport(X,dim-1,mode,p);
   if (X[2]<0) {
                             X[2]=-X[2];
                             mode=-mode;
   }
   if (X[2]>L) {
                             X[2]=2*L-X[2]; 
                             mode=-mode;
   }
}
void getDensityRec(long *p){
   if (dim>=2){
      //Dm[1]=Dp[1];
      vdrift[1]=-vdrift[0];  
   }
   std::vector<double> X(dim,0.005); 
   int mode=0;  
   int Nr=int(floor((L2-L1)/hr));
   int Nz=int(floor(L/hr));
   double theta=(sqrt(Dp[0])-sqrt(Dm[0]))/(sqrt(Dm[0])+sqrt(Dp[0]));//discontinuity in the first dimension
   printf("Nr=%d Nz=%d, theta=%f, vdirft %f %f\n",Nr,Nz,theta, vdrift[0], vdrift[1]);
   printf("discountiy at 0, interface layer: left end %f right end %f\n",-sqrt(Dm[0]*dt)*ca, sqrt(Dp[0]*dt)*ca);

   std::vector<int> den(Nr,0);
   std::vector<std::vector<int> > den2(Nz,den);

   std::vector<std::vector<int> > dendirect(Nz,den);
   FILE * fdiff;
   FILE * fdirect;
   char filename[256];
   sprintf(filename,"diff_L_%2.2f_hr%2.4f_T%d_dt%f_Iter%.1f.txt",L2,hr,int(T),dt,Iter);
   fdiff=fopen(filename,"w");
   sprintf(filename,"pro_L_%2.2f_hr%2.4f_T%d_dt%f_Iter%.1f.txt",L2,hr,int(T),dt,Iter);
   fdirect=fopen(filename,"w");
   
   for (double i=0;i<Iter;i++){
          moveActiveRec(X,mode,theta,p);
          if (X[0]<L1) X[0]=2*L1-X[0];
          if (X[0]>L2) X[0]=2*L2-X[0];
          int a=int(floor((X[0]-L1)/hr));
          int b;
          if (X.size()==2){
                 if (X[1]<0) {
                             X[1]=-X[1];
                             mode=-mode;
                 }
                 if (X[1]>L) {
                                 X[1]=2*L-X[1]; 
                                 mode=-mode;
                 }
                 if (mode==0){
                     b=int(floor(X[1]/hr)); 
                     den2[b][a]=den2[b][a]+1;  
                 }
                 else{
                      b=int(floor(X[1]/hr)); 
                      dendirect[b][a]=dendirect[b][a]+1; 
                 } 
          }
          else if(X.size()==1) den[a]=den[a]+1;          
   }
   if (X.size()==1){
      for (int i=0;i<den.size();i++) fprintf(fdiff,"%d ",den[i]);
   }
   else if (X.size()==2){
     for (int i=0;i<den2.size();i++){
       for (int j=0;j<den2[0].size();j++){
           fprintf(fdiff,"%d ",den2[i][j]);
           fprintf(fdirect,"%d ",dendirect[i][j]);
       }
       fprintf(fdiff,"\n");
       fprintf(fdirect,"\n");
     }
   }
   fclose(fdiff);
   fclose(fdirect);
}

void getDensityDisk(long *p){  
   int mode=0;
   std::vector<double> X(dim,0.2);
   X[dim-1]=0.1;
   if (dim==2) X[1]=-0.03;
   if (dim==3) {
          Dm[2]=2*0.015;
          Dp[2]=2*0.003;//pow(10,-4.8);//0.0003;//0.015;
   }
   int Nr=int(ceil(2*R/hr));
   int Nl=int(ceil(L/hr));
   double theta=(sqrt(Dp[0])-sqrt(Dm[0]))/(sqrt(Dm[0])+sqrt(Dp[0]));
   printf("dimension=%d Nr=%d theta=%f drift %f %f Dp %f %f \n",dim,Nr,theta,vdrift[0], vdrift[1], Dp[0], Dp[2]);
   printf("discountiy at 0, interface layer: left end %f right end %f\n",-sqrt(Dm[0]*dt)*ca, sqrt(Dp[0]*dt)*ca);
   std::vector<int> den(Nr,0);
   std::vector<int> denx(Nl,0);
   std::vector<std::vector<int> > den2;//(Nr,den);
   std::vector<std::vector<int> > denR(Nl,den);
   std::vector<std::vector<int> > denRdir(Nl,den);
   std::vector<std::vector<std::vector<int> > > den3;//(4,den2);
   std::vector<std::vector<std::vector<int> > > dendirect;//(4,den2);
   FILE * ftraj;
   FILE * fdiff;
   FILE * fdr;// rectangle r-x
   FILE * fdirect;
   FILE * fRdirect;
   FILE * fdx;
   FILE * fr;
   FILE * fMT;   

   char filename[256];
   if (flagden){
       sprintf(filename,"D.txt");//%d_R0%.1f_dpl%1.5f_ca%d_L_%2.2f_hr%2.4f_drift%1.4f_Dm%.4f_Dp%.4f_T%d_dt%1.4f_vdirect%.2f.txt",dim,R0,dplanar,int(ca),L,hr,int(T),vdrift[0],Dm[0],Dp[0],dt,vdirect);
       fdiff=fopen(filename,"w");
   }
   if (flagdenAL){
       sprintf(filename,"rhoRL_R0%.2f_L_%2.2f_hr%2.4f_drift%1.4f_Dmt%.4f_Dm[0]%.8f_Dm[2]%.5f_Dp[0]%.8f_Dp[2]%.5f_T%d_dt%1.4f_alpha%.6f_wd%.6f_w%.6f.txt",R0,L,hr,vdrift[0],Dmt[0],Dm[0],Dm[2],Dp[0],Dp[2],int(T),dt,alpha,delta,omega);
       fdr=fopen(filename,"w");
       sprintf(filename,"rhoRL_direct_R0%.2f_L_%2.2f_hr%2.4f_drift%1.4f_Dmt%.4f_Dm[0]%.8f_Dm[2]%.5f_Dp[0]%.8f_Dp[2]%.5f_T%d_dt%1.4f_alpha%.6f_wd%.6f_w%.6f.txt",R0,L,hr,vdrift[0],Dmt[0],Dm[0],Dm[2],Dp[0],Dp[2],int(T),dt,alpha,delta,omega);
       fRdirect=fopen(filename,"w");
   }
   if (flagdenLateral){
       sprintf(filename,"rhoR_R0%.2f_L_%2.2f_hr%2.4f_drift%1.4f_Dmt%.4f_Dm[0]%.8f_Dm[2]%.5f_Dp[0]%.8f_Dp[2]%.5f_T%d_dt%1.4f_alpha%.6f_wd%.6f_w%.6f.txt",R0,L,hr,vdrift[0],Dmt[0],Dm[0],Dm[2],Dp[0],Dp[2],int(T),dt,alpha,delta,omega);
       fr=fopen(filename,"w");
   }
   if (flagtraj){
       sprintf(filename,"testTraj.txt");
       ftraj=fopen(filename,"w");
   }
   if (flagMT){
       sprintf(filename,"testMT.txt");
       fMT=fopen(filename,"w");
   }
   if (flagdenDirect){
       sprintf(filename,"rhodirect_x_R0%.2f_L_%2.2f_hr%2.4f_drift%1.4f_Dmt%.4f_Dm[0]%.8f_Dm[2]%.5f_Dp[0]%.8f_Dp[2]%.5f_T%d_dt%1.4f_alpha%.6f_wd%.6f_w%.6f.txt",R0,L,hr,vdrift[0],Dmt[0],Dm[0],Dm[2],Dp[0],Dp[2],int(T),dt,alpha,delta,omega);
       fdirect=fopen(filename,"w");
   }
   if (flagdenAxis){
       sprintf(filename,"rhox_R0%.2f_L_%2.2f_hr%2.4f_drift%1.4f_Dmt%.4f_Dm[0]%.8f_Dm[2]%.5f_Dp[0]%.8f_Dp[2]%.5f_T%d_dt%1.4f_alpha%.6f_wd%.6f_w%.6f.txt",R0,L,hr,vdrift[0],Dmt[0],Dm[0],Dm[2],Dp[0],Dp[2],int(T),dt,alpha,delta,omega);
       fdx=fopen(filename,"w");
   }
   std::vector<double> center(2,0);
   for (step=0;step<Iter;step++){
          
   //       moveDiskFull(X,mode,theta,dt,p);
          moveMTFull(X,mode,theta,dt,center,p);
          //printf("X=%f %f %f\n",X[0],X[1],X[2]);
          double rx=sqrt(X[0]*X[0]+X[1]*X[1]);
          int a=int(floor((X[0]+R)/hr));
          int b=int(floor((X[1]+R)/hr));
          int c,c3;
          
          if (X.size()==3){
               /*  if (X[2]<0) {
                             X[2]=-X[2];
                             mode=-mode;
                 }
                 if (X[2]>L) {
                             X[2]=2*L-X[2]; 
                             mode=-mode;
                 }
                 */
                 c=int(floor(X[2]/hr));//*3/Nl); 
                 //printf("c=%d \n",c);
                 c3=int(floor(X[2]/hr*den3.size()/Nl));
                 int rr=int(floor(rx/hr));
                 if (flagden) den3[c3][b][a]=den3[c3][b][a]+1;
                 if (flagtraj) fprintf(ftraj,"%f, %f %f\n",X[0],X[1],X[2]);  
                 if (flagMT) fprintf(fMT, "%f %f\n", center[0],center[1]); 
                 
                 if (flagdenAL & mode==0){                    
                    denR[c][rr]=denR[c][rr]+1; 
                 }   
                 if (mode!=0 && flagdenAL ){//flagdenDirect){
                     denRdir[c][rr]=denRdir[c][rr]+1;
                     // dendirect[c3][b][a]=dendirect[c3][b][a]+1;
                 }
                 if (flagdenAxis){ 
                   int xx=int(floor(X[2]/hr));
                   denx[xx]=denx[xx]+1;
                 }
                 if (flagdenLateral & mode==0) den[rr]=den[rr]+1;
          }
          else if(X.size()==2) {
               den2[b][a]=den2[b][a]+1;
               //fprintf(ftraj,"%f, %f \n",X[0],X[1]);
               int rr=int(floor(rx/hr));
               den[rr]=den[rr]+1;
          }
          
          //fprintf(fden,"%f \n",rx);
   }
   if (dim==2){
               for (int i=0;i<Nr/2;i++) fprintf(fdr,"%f \n",double(den[i])/double(i+1));
   }
   if (dim==3){
                for (int i=0;i<Nl;i++) {
                    if (flagdenAL){
                        for (int j=0;j<Nr/2;j++){
                             fprintf(fdr,"%f ",double(denR[i][j])/double(j+0.5));
                             fprintf(fRdirect,"%f ",double(denRdir[i][j])/double(j+0.5));
                        }  
                        fprintf(fdr,"\n");
                        fprintf(fRdirect,"\n");
                    }
                    if (flagdenAxis)   fprintf(fdx, "%d ",denx[i]);
                }
                if (flagdenLateral) {
                     for (int i=0;i<den.size();i++) fprintf(fr,"%f ",den[i]/double(i+0.5));
                }
   }
   
   
  
   if (X.size()==2){
     for (int i=0;i<den2.size();i++){
       for (int j=0;j<den2[0].size();j++){
           fprintf(fdiff,"%d ",den2[i][j]);
       }
       fprintf(fdiff,"\n");
     }
   }
   else if (X.size()==3 && flagden){
       for (int ic=0;ic<den3.size();ic++){
            for (int ib=0;ib<Nr;ib++){
                for (int ia=0;ia<Nr;ia++){
                    fprintf(fdiff,"%d ",den3[ic][ib][ia]);
                    if (vdirect>0 && flagdenDirect) fprintf(fdirect,"%d ",dendirect[ic][ib][ia]);
                }
            }
             fprintf(fdiff,"\n");
            if (flagdenDirect) fprintf(fdirect,"\n");
        }
   }
   if (flagdenLateral) fclose(fr);
   if (flagden) fclose(fdiff);
   if (flagtraj) fclose(ftraj);
   if (flagMT) fclose(fMT);
   if (flagdenDirect) fclose(fdirect);
   if (flagdenAxis) fclose(fdx);
   if (flagdenAL) {
      fclose(fdr);
      fclose(fRdirect);
   }
}


void getHitTime(long *p){
          Dm[2]=2*0.0152;
          Dp[2]=2*0.003;          
   double theta=(sqrt(Dp[0])-sqrt(Dm[0]))/(sqrt(Dm[0])+sqrt(Dp[0]));
   FILE * fhit;
   char filename[256];
   if (flaghit){
   sprintf(filename,"hittime_R0%.1f_dpl%1.5f_L_%2.2f_drift%1.4f_Dmt%.4f_Dm[0]%.4f_Dp[0]%.4f_T%d_dt%1.4f_valpha%.5f.txt",R0,dplanar,L,vdrift[0],Dmt[0],Dm[0],Dp[0],int(T),dt,alpha);
   fhit=fopen(filename,"w");
   }
   int runs=1000;
   for (int run=0;run<runs;run++){
        int Nhit=0;
        double a[2]={0,0};
        double rx=sqrt(a[0]*a[0]+a[1]*a[1]);
        while(rx<0){
          rand_disk(p,a,R);
         rx=sqrt(a[0]*a[0]+a[1]*a[1]);
        }
        
        std::vector<double> X(dim,0.2);
        X[0]=a[0];
        X[1]=a[1];
        X[2]=1;
        std::vector<double> center(2,0);
        int mode=0;
        for (step=0;step<Iter;step++){      
            moveMTFull(X,mode,theta,dt,center,p);
            double rc=sqrt((X[0]-center[0])*(X[0]-center[0])+(X[1]-center[1])*(X[1]-center[1]));
            if (rc<R0 && flaghit) Nhit=Nhit+1;
        }
        fprintf(fhit,"%f ",double(Nhit)/double(Iter));
   }
   if (flaghit) fclose(fhit);
}


void calculateBindProb(long *p){
   extern bool flagexit;
   
   Dm[2]=2*0.0152;
   Dp[2]=2*0.003;
          
   int Nexit=0;
   int Nbind=0; 
   std::vector<double> Tbind;
   double theta=(sqrt(Dp[0])-sqrt(Dm[0]))/(sqrt(Dm[0])+sqrt(Dp[0]));
   FILE * fbinding;
   
   char filename[256];
   int pxt=5000;
   sprintf(filename,"binding3_px%d_MT%1.4f_R0%.1f_L%2.2f_drift%1.4f_Dm[0]%.4f_Dp[0]%.4f_dt%1.4f_alpha%.4f.txt",pxt,Dmt[0],R0,L,vdrift[0],Dm[0],Dp[0],dt,alpha);
   fbinding=fopen(filename,"w");
  // extern double alpha;
   //printf("alpha=%f",alpha);
   int px;
   for (px=0;px<pxt;px++){
      int mode=0;
      flagexit=false;
      double a[2]={0,0};
      double rx=sqrt(a[0]*a[0]+a[1]*a[1]);
      //while(rx<0){
      rand_disk(p,a,R);
      rx=sqrt(a[0]*a[0]+a[1]*a[1]);
      //}
      step=0;
      std::vector<double> X(dim,0.2);
      X[0]=a[0];
      X[1]=a[1];
      X[2]=1;
     // printf("X=%f %f %f rx=%f \n",X[0],X[1],X[2],rx);
      std::vector<double> center(2,0);       
      while (mode==0){//!flagexit  && mode==0 && step<3000){           
          //moveDiskFull(X,mode,theta,dt,p);
          moveMTFull(X,mode,theta,dt,center,p); 
          step=step+1;         
      }  
    //  printf("mode=%d \n",mode);
      if (mode==1 || mode==-1) {
                   //Tbind.insert(Tbind.end(),step*dt);
                   Nbind=Nbind+1;
                   fprintf(fbinding,"%f ",step*dt);
                  //  printf("bind time %f \n",step*dt);
      }
      else if(flagexit) {
            Nexit=Nexit+1;
           // printf("exit time %f \n",step*dt);
            double rx=sqrt(X[0]*X[0]+X[1]*X[1]);
           // printf("rx=%f \n",rx);
      }
   }    
   printf("total=%d bind=%d exit=%d\n",px,Nbind,Nexit);
          
}



void getMSD(long *p){
   Dm[2]=2*0.015;
   Dp[2]=2*0.003;
          
   double theta=(sqrt(Dp[0])-sqrt(Dm[0]))/(sqrt(Dm[0])+sqrt(Dp[0]));
   FILE * fMSD;
   
   char filename[256];
   int pxt=5000;
   sprintf(filename,"MSD_dim3_px%d_R0%.1f_L%2.2f_drift%1.4f_Dmt%.4f_Dm[0]%.4f_Dm[2]%.4f_Dp[0]%.4f_Dp[2]%.4f_dt%1.4f_alpha%.4f.txt",pxt,R0,L,vdrift[0],Dmt[0],Dm[0],Dm[2],Dp[0],Dp[2],dt,alpha);
   fMSD=fopen(filename,"w");
   extern double alpha;
   alpha=0;
   int px;
   double tmax=50;//10;
   for (double i=1;i<tmax/dt;i++) fprintf(fMSD,"%f ",dt*i);
   fprintf(fMSD,"\n");
   
   for (px=0;px<pxt;px++){
      int mode=0;
      std::vector<double> traj;
      std::vector<double> MSD;
      double a[2]={0,0};
      rand_disk(p,a,R);
      std::vector<double> X(dim,0.2);
      X[0]=a[0];
      X[1]=a[1];
      X[2]=1;
      std::vector<double> center(2,0);
      for (double step=0;step<2.0*tmax/dt;step++){
          //moveDiskFull(X,mode,theta,dt,p);
          moveMTFull(X,mode,theta,dt,center,p);
          traj.insert(traj.end(),X[2]);
      }
      calculateMSD(traj,MSD,dt,tmax);
      for (int i=0;i<MSD.size();i++) fprintf(fMSD,"%f ",MSD[i]);
      fprintf(fMSD,"\n");
   }
   fclose(fMSD);
}
void testpdf(long *p){
     //std::vector<double> X(dim,0.5); 
     int mode=0;  
     int N=int(floor(1.0/dt));
     printf("N=%d", N); 
     FILE * ftest;
     char filename[256];
     sprintf(filename,"pdf.txt");
     ftest=fopen(filename,"w");
     double theta=(sqrt(Dp[0])-sqrt(Dm[0]))/(sqrt(Dm[0])+sqrt(Dp[0]));
     for (double i=0;i<1;i++){
       printf("i=%f \n",i);
       std::vector<double> X(1,0.5);
       for (int j=0;j<N;j++) {
          //twoStepsMoveRec(X,dt,Dm,Dp,ca,p);
          moveActiveRec(X,mode,theta,p);
          if (X[0]>1) X[0]=2-X[0];
          else if (X[0]<-1) X[0]=-2-X[0];
          printf("mode=%d X=%f \n",mode,X[0]);
          fprintf(ftest,"%f ",X[0]);    
       }
       //fprintf(ftest,"%f ",X[0]);       
     }
     fclose(ftest);
}
