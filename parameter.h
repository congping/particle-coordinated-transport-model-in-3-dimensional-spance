#ifndef PARAMETER_H
#define PARAMETER_H
#include <math.h>
#include <vector>
#define PI 3.141592653589793
bool flagdenAxis=true;//false;//true;
bool flagdenLateral=true;
bool flagdenAL=true;
bool flagden=false;
bool flagdenDirect=true;
bool flagtraj=false;//true;
bool flagbind=false;
bool flagMT=false;//true;
bool flaghit=false;//true;

bool direct=true;
int dim=3;
std::vector<double> Dm(dim,2*0.0000014);
std::vector<double> Dp(dim,2*0.0014);//pow(10,-4.4));//0/.000014)
std::vector<double> Dmt(2,2*0.0);
double ca=5;
double dt=0.01; 
double T=100000000;//000;
double Iter=floor(T/dt);
std::vector<double> vdrift(2,0.0006);//(2,0.00044*2);

double hr=0.01;//0.0645;
double hx=hr;
double L1=-0.4; // for the rectangle domain : discontinuity at 0
double L2=-L1;
double L=30;// [0,L] for the other dimension besides [L1,L2] in rectangle or D(0,R) in disk
double R=0.86;
double R0=0.25;//0.4*R;//*0.4; // for the discontinuity in disk domain
double dplanar=R0*0.02;
double dt1=pow(dplanar/ca,2)/std::max(Dm[0],Dp[0]);//reduced time step in interface layer
int rat=int(dt/dt1)+1;//reduced time step in interface layers
double vdirect=1.9;
    //attach rate
double alpha=0.0034*100;//1.0;
double alphat=alpha*dt;
    //turning rate
double omega = 0.12;
double omegat =omega*dt;
    // detach rate
double delta = 0.17;
double deltat=delta*dt;




#endif // PARAMETER_H
