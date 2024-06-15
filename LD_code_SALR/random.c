//Parallised random number generator, provided by Adam Wysocki
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include "system.h"

/* fix RNOR */

float fix(void)
{const float r=3.442620;
 float x,y;

 for (;;)
    {x=(float) (HZ)*WN[IZ];
     if (IZ==0)
       {do
          {x=-log(RandomNumber)*0.2904764; // 1/r=0.2904764
           y=-log(RandomNumber);
          }
        while(y+y<x*x);
        return (HZ>0)?r+x:-r-x;
       }
      if(FN[IZ]+RandomNumber*(FN[IZ-1]-FN[IZ])<exp(-0.5*x*x)) return x;
      HZ=SHR0;
      IZ=HZ&127;
      if(abs(HZ)<KN[IZ]) return ((float) (HZ)*WN[IZ]);
     }
}

/* initiate RNOR */

void zigset(void)
{int i;
 double q,dn=3.442619855899,tn=3.442619855899;
 const double m1=2147483648.0,vn=9.91256303526217e-3;

 q=vn/exp(-0.5*dn*dn);

 KN[0]=(int) ((dn/q)*m1);
 KN[1]=0;

 WN[0]=(float) (q/m1);
 WN[127]=(float) (dn/m1);

 FN[0]=1.0;
 FN[127]=(float) (exp(-0.5*dn*dn));

 for (i=126;i>=1;i--)
    {dn=sqrt(-2.0*log(vn/dn+exp(-0.5*dn*dn)));
     KN[i+1]=(int) ((dn/tn)*m1);
     tn=dn;
     FN[i]=(float) (exp(-0.5*dn*dn));
     WN[i]=(float) (dn/m1);
    }
}


/* 

int InitializeRandomNumberGenerator(void)
{int thread;

 zigset();
 JSR=time(NULL);
 for (thread=0;thread<N_threads;thread++) {seed[thread]=SHR3;}

return 0;
}
*/

int InitializeRandomNumberGenerator(double seedseed)
{int thread;

 zigset();
// JSR=time(NULL);
 JSR=seedseed;
 for (thread=0;thread<N_threads;thread++) {seed[thread]=SHR3;}

return 0;
}

/* */

// Generates A Random Velocity According To A Boltzmann Distribution
double RandomVelocity(double temperature)
{
  return sqrt(temperature)*RandomGaussianNumber;
}

