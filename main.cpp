#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>

void test(long *p);
void test2steps(long *p);
void getDensityRec(long *p);
void getDensityDisk(long *p);
double ran2(long *p);
void testdirect(long *p);
void testpdf(long *p);
void calculateBindProb(long *p);
void getMSD(long *p);
void getHitTime(long *p);
extern double R0;   
double step;

extern bool flagbind;
FILE * fdetach; 
bool flagexit;
int main(){
    clock_t t;
    t = clock();
    long idum;
    long *P2Idum=&idum;
    idum=-(long)time(NULL);
    printf("Ro=%f",R0);
    char filename[256];
    sprintf(filename,"detach_time_R0%f.txt",R0);
    if (flagbind) fdetach=fopen(filename,"w");
//    getMSD(P2Idum);
    //getHitTime(P2Idum);//
    //calculateBindProb(P2Idum);
    getDensityDisk(P2Idum);
    //getDensityRec(P2Idum);
    //test2steps(P2Idum);
    //test(P2Idum);
    //testpdf(P2Idum);
    //testdirect(P2Idum);
    t=clock()-t;
    if (flagbind) fclose(fdetach);
    printf("time processed %f second\n",float(t)/CLOCKS_PER_SEC);
    system("pause");
    return 0;
}

