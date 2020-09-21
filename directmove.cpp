#include <math.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>



double ran2(long *);

double rand_norm_1(long * p);
void exitTime(double & x, double & dt,long * p);

extern bool flagexit;
extern std::vector<double> vdrift;
extern int dim;

double Zplus(double z,long *p){
   if (z<0){
       double r=rand_norm_1(p);
       while (r<-z)
           r=rand_norm_1(p);
       return z+r;             
   }
   else{
        printf("positive z");
        system("pause");
   }
}


double Zminus(double z, double theta, long *p){
   if (z<=0 && theta>0){
        bool flag=true;
        double H;
        while(flag){
            double r=rand_norm_1(p);
            while(r<z) r=rand_norm_1(p);
            H=z-r;
            double r1=ran2(p);
            if (r1<theta*exp(-2*H*z)) flag=true;
            else flag=false;        
        }
        return H;
   }
   else if (z<=0 && theta<=0){
        double y;
        double r1=ran2(p);
        double phi=(1+erf(-z/sqrt(2)))/2;
        double u=1-phi;
        double alpha=(1-u)/(1-(1+theta)*u);
        if (r1<alpha){
           double r=rand_norm_1(p);
           while(r<z) r=rand_norm_1(p);
           y=z-r;
        }
        else{
           double r=rand_norm_1(p);
           while(r<-z) r=rand_norm_1(p);
           y=-z-r;  
        }
        return y;
   }
   else{
        printf("positive z");
        system("pause");
   }
}


void moveInterfaceDirect(double & x, double dt, double theta, long *p){
    double H;
    double x1=x;
    if (x>0){
             x1=-x;
             theta=-theta;
    }
    double u=ran2(p);
    double z=x1/sqrt(dt);
    double phi=(1+erf(-z/sqrt(2)))/2;
    phi=(1+theta)*(1-phi);
    if (u<=phi) H=sqrt(dt)*Zplus(z,p);
    else H=sqrt(dt)*Zminus(z,theta,p);
    
    if (x>0) H=-H;
    x=H;
}

void moveDrift(std::vector<double> &X, int ind, bool flag, double dt,long*p){
         if (flag && X.size()>1){
              X[ind]=X[ind]+vdrift[0]*dt;   
         }
         else if (!flag && X.size()>1) {
             // printf("before drift X=%f %f\n",X[0],X[1]);
              X[ind]=X[ind]+vdrift[1]*dt;
             // printf("after drift X=%f %f\n",X[0],X[1]);
         }
}

void directMoveRec(std::vector<double> & X, double dt, double theta, std::vector<double> Dm, std::vector<double> Dp, double ca,long * p){     
    double x1=X[0];      
    if (X[0]<-sqrt(Dm[0]*dt)*ca){
         double r=rand_norm_1(p)*sqrt(dt);                           
         X[0]=X[0]+sqrt(Dm[0])*r;
         //printf("left X=%f \n",X[0]);
    }
    else if (X[0]>sqrt(Dp[0]*dt)*ca){
         double r=rand_norm_1(p)*sqrt(dt);   
         X[0]=X[0]+sqrt(Dp[0])*r;
    }
    else {
         if (X[0]>0) X[0]=X[0]/sqrt(Dp[0]);
         else X[0]=X[0]/sqrt(Dm[0]);
         moveInterfaceDirect(X[0],dt,theta,p);//here X[0] is dimensionless
         if (X[0]>0) X[0]=X[0]*sqrt(Dp[0]);
         else X[0]=X[0]*sqrt(Dm[0]);
    }
    //printf("X size %d",X.size());
    if (X.size()>1){  
        //printf("before X=%f %f\n",X[0],X[1],X[2]);
        double Dmix;
        double tt=dt;
        if (x1>=-sqrt(Dm[0]*dt)*ca && x1<=sqrt(Dp[0]*dt)*ca){//cgabge X[0] to x1 --CL 06_10_14
            double x=x1;
            if (x1>0) x=x/sqrt(Dp[0]);
            else x=x/sqrt(Dm[0]); 
            exitTime(x,tt,p);
        }
        if (x1<0) moveDrift(X,dim-1,true,tt,p); 
        else  moveDrift(X,dim-1,false,tt,p);          
        if (tt==dt){
             for (int i=1;i<X.size();i++){  
                  Dmix=Dm[i]*(x1<0)+Dp[i]*(x1>=0); 
                  X[i]=X[i]+sqrt(Dmix*dt)*rand_norm_1(p); 
             }    
            // flagexit=false;                               
        }  
        else{//crossing
            if (x1<0) flagexit=true;  
           // else flagexit=false; 
            double alpha=sqrt(Dm[0])/(sqrt(Dm[0])+sqrt(Dp[0]));
            for (int i=1;i<X.size();i++){                   
                 Dmix=Dm[i]*alpha+Dp[i]*(1-alpha);   
                 if (x1>0) {
                       X[i]=X[i]+sqrt(tt*Dp[i])*rand_norm_1(p);
                       X[i]=X[i]+sqrt((dt-tt)*Dmix)*rand_norm_1(p);                                  
                  }   
                  else{
                       X[i]=X[i]+sqrt(tt*Dm[i])*rand_norm_1(p);
                       X[i]=X[i]+sqrt((dt-tt)*Dmix)*rand_norm_1(p);
                  }                         
            }   
            double v=vdrift[0]*alpha+(1-alpha)*vdrift[1];
            X[dim-1]=X[dim-1]+v*(dt-tt);
        }     
        //printf("after X=%f %f\n",X[0],X[1],X[2]);       
    }
}


void testdirect(long *p){
    int N=100000;
    double theta=0.1;
    double dt=0.01;
    FILE * ftest;
    char filename[256];
    sprintf(filename,"test.txt");
    ftest=fopen(filename,"w");
    for (int i=0;i<N;i++){
         double z=0.5;//-ran2(p);
         //double z=Zplus(r,p);
         //double z=Zminus(r,theta,p);
         for (int j=0;j<100;j++) moveInterfaceDirect(z,dt,theta,p);
         fprintf(ftest,"%f ",z);
    }
    fclose(ftest);
}
