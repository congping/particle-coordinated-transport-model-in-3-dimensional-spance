#include<time.h>
#include<iostream>
#include <stdlib.h>
#include <math.h>
//#include <stdio.h>
#define PI 3.14159265359
double ran2(long *);

double rand_norm_1(long * p){
      double x,y;
      x=ran2(p);
      y=ran2(p);
      double s;
      //double PI=3.14159265359;
      s=sqrt(-2*log(x))*cos(2*PI*y);
      return s;
}
double rand_norm_2(long *p){
      double x1, x2, w, y1, y2;    
 
         do {
                 x1 = 2.0 * ran2(p) - 1.0;
                 x2 = 2.0 * ran2(p) - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
         y1 = x1 * w;
         y2 = x2 * w;
         return y1;
}


void rand_disk(long *p, double a[2],double R){// uniform distribution in a disk of radius R
   double t = 2*PI*ran2(p);
   double u = ran2(p)+ran2(p);
   double r;
   if (u>1) r=2-u;
   else  r=u;
   r=r*R;
      
   a[0]=r*cos(t);
   a[1]=r*sin(t);
   //printf("a=%f %f \n",a[0],a[1]);;
   
}
 void test_norm(long * p){   
      /* 
    FILE * rand_gen;
    rand_gen=fopen("norm.txt","w");
    for (int i=0;i<100000;i++) fprintf(rand_gen, "%4f ",rand_norm_1(p));
    std::cout<<"test";
    fclose(rand_gen);
    */
}
      
 
