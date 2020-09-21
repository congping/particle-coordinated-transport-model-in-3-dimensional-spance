#include <math.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>

double ran2(long *);
double rand_norm_1(long * p);
double inverseGaussian(double lam,double mu, long * p){
       double v=rand_norm_1(p);
       v=v*v;
       double x1=mu+mu/(2*lam)*(mu*v-sqrt(4*mu*lam*v+v*v*mu*mu));
       double r=ran2(p);
       if (r<mu/(mu+x1)) return x1;
       else return mu*mu/x1;
}


   
void exitTime(double & x, double & dt,long * p){
      double r=rand_norm_1(p);
      double y=x+sqrt(dt)*r;
      
      if (x*y<0){
           r=inverseGaussian(fabs(x/y),x*x/dt,p);
           dt=dt*r/(1+r);
           x=0;
      }
      else{
           r=ran2(p);
           if (r<exp(-2*x*y/dt)){
               double r3=inverseGaussian(fabs(x/y),x*x/dt,p);
               dt=dt*r3/(1+r3);
               x=0;
           }
           else {
                x=y;
           }
      }
}


double moveInterface2steps(double & x, double dt, double theta, long *p){//dimensionless
      double t=dt;
      exitTime(x,t,p);
     // if (t==dt) return x
      //else{
      if (dt!=t){
           double eps=-1;
           if (ran2(p)<(1+theta)/2) eps=1;   
           double r=rand_norm_1(p);
           x=eps*sqrt(dt-t)*fabs(r);              
      }    
      return t;   
}

// dimensionless in x, interface x=0, Dm, Dp is for diffusion in y axis, theta is parameter from diffusion in x-axis
void moveInterface2stepsMD(std::vector<double> & X, double dt, double theta, std::vector<double> Dm, std::vector<double> Dp, long * p){     
     double t=moveInterface2steps(X[0],dt,theta,p);
     if (t==dt){// no cross
       for (int i=1;i<X.size();i++){
          double r=rand_norm_1(p)*sqrt(dt);
          if (X[0]<0) X[i]=X[i]+sqrt(Dm[i])*r;
          else X[i]=X[i]+sqrt(Dp[i])*r; 
       }
     }
     else{
        for (int i=1;i<X.size();i++){
           double Dalpha=sqrt(Dm[i])/(sqrt(Dm[i])+sqrt(Dp[i]));
           double Dmix=Dm[i]*Dalpha+Dp[i]*(1-Dalpha);
        
           double r=rand_norm_1(p)*sqrt(t);
           if (X[0]<0) X[i]=X[i]+sqrt(Dm[i])*r;
           else X[i]=X[i]+sqrt(Dp[i])*r;            
           r=rand_norm_1(p)*sqrt(dt-t);
           X[i]=X[i]+sqrt(Dmix)*r;
        }
     }
     
}

// x=0 interface, Dm, Dp is x<0 and x>0   for rectangle region
void twoStepsMoveRec(std::vector<double> & X, double dt, std::vector<double> Dm, std::vector<double> Dp, double ca, long *p){
   if (X[0]<-sqrt(Dm[0]*dt)*ca){
      for (int i=0;i<X.size();i++){
          double r=rand_norm_1(p)*sqrt(dt);
          X[i]=X[i]+sqrt(Dm[i])*r;                           
      }
   }
   else if (X[0]>sqrt(Dp[0]*dt)*ca){
      for (int i=0;i<X.size();i++){
          double r=rand_norm_1(p)*sqrt(dt);
          X[i]=X[i]+sqrt(Dp[i])*r;                           
      }
   }
   else{      
      if (X[0]>0) X[0]=X[0]/sqrt(Dp[0]);
      else X[0]=X[0]/sqrt(Dm[0]);
      double theta=(sqrt(Dp[0])-sqrt(Dm[0]))/(sqrt(Dp[0])+sqrt(Dm[0]));
      moveInterface2stepsMD(X,dt,theta,Dm,Dp,p);
      if (X[0]>0)  X[0]=X[0]*sqrt(Dp[0]);
      else  X[0]=X[0]*sqrt(Dm[0]);
   }
}



void test(long *p){
   FILE * ftest;
   ftest=fopen("test.txt","w");
   double Dp=0.10;
   double Dm=0.005;
   double dt=0.005;
   
   double theta=(sqrt(Dp)-sqrt(Dm))/(sqrt(Dp)+sqrt(Dm));
   for (int i=0;i<500000000;i++){
       double x=-0.05;       
       double t=moveInterface2steps(x,dt,theta,p);
       //exitTime(x,dt,p);
       //double x=inverseGaussian(1,1,p);
       fprintf(ftest,"%f ",x);
   }
}
