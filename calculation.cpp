// calculations
#include <vector>
#include <stdio.h>      /* printf */
#include <math.h>       /* floor */
void calculateMSD(std::vector<double> & traj, std::vector<double> & MSD, double dt,double tmax){
     int T=int(floor(tmax/dt));
     MSD.clear();
     for (int i=1;i<T;i++){
        double s=0;
        for (int j=0;j<traj.size()-i;j++){
              s=s+(traj[j+i]-traj[j])*(traj[j+i]-traj[j]);
        }
        s=s/(traj.size()-i);
        //printf("s=%f size=%d",s,traj.size()-i);
        MSD.insert(MSD.end(),s);
     }
}
