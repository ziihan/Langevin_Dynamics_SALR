#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "random.h"
#include "system.h"
#include <omp.h>

//Some fixed parameters, constants. Before initialization 
//the second half is the Main function 

void Preinit(void)
{
dq=2.0*PI/MIN(BoxSize_x,BoxSize_y);

Nx=24;
Ny=24;
Nz=1;
SampleSteps=10;
PrintStep=100;
RecordStep=0;
cCon=1.0;//initial value for bending (cosine(PI-bondangle) value)

Rcelli=3;
Rcello=3;
//set up the locations of walls
Wall0=1.0;
Wallz=BoxSize_z-1.0;
Nobs=1;
//MDt=0.02*Deltat;
NumberOfColloids=Nobs*Nx*Ny*Nz;
NumberChainColloids=0;

Rcol=2.5;//colloid radius
Rcc=2.0*Rcol;
BondLength=8.5;//distance between two consecutive colloids
RodLength=2.0*Rcol;//actually here is the spring length
//how to chose these parameters
ColloidBending=500000.0;
ColloidSpring=50000.0;//rather large

Mcol=654.498;
invMcol=1.0/Mcol;
HeatOn=0.2;

xif=6.0*PI*8.7*Rcol/Mcol;

//CutOff=1.2599210498948732*Rcol;//1.12246204831*Rcol;

CutOffi=1.12246204831*Rcol;//12-6 rLJ
CutOffo=1.121353391970139*Rcol;

CutOffcc=Rcc*pow(2.0,1.0/50.);//1.12246204831*Rcc;
CutOffYK=5.0*Rcc;

CutOffw=1.12246204831*BoxSize_y;

qq=2.0;
DebyeLength=1.794*Rcc;

//CutOff=pow(2.0,0.166666666667)*Rcol;
Eshif=1.0;//4*(pow(0.4,12.0)-pow(0.4,6.0));

Ecut=4.0*(pow((1.0/2.5),12.0)-pow((1.0/2.5),6.0));//standard LJ potetial for colloid-colloid interaction
EcutYK= 4.0*epsilon*(pow((Rcc/CutOffYK),100.0)-pow((Rcc/CutOffYK),50.0))+qq*(DebyeLength/CutOffYK)*exp(-CutOffYK/DebyeLength);
FcutYK= 8.0*50.0*epsilon*(pow((Rcc/CutOffYK),100.0)-0.5*pow((Rcc/CutOffYK),50.0))*(1.0/CutOffYK)+qq*DebyeLength*exp(-CutOffYK/DebyeLength)*(1.0/CutOffYK+1.0/DebyeLength)*(1.0/CutOffYK);
//initialization the forces for the velocity Verlet algorithm
beta=(90.0/180.0)*PI;

c0=(1.0/xif)*(1.0-exp(-xif*MDt));
c1=(1.0/(xif*xif*Mcol))*(xif*MDt-1.0+exp(-xif*MDt));
c2=xif/(exp(xif*MDt)-1.0);
c3=(1.0/(xif*xif*Mcol))*(exp(xif*MDt)-xif*MDt-1.0);
c2c3=c2*c3;

Cx=2.0*xif*MDt-3.0+4.0*exp(-xif*MDt)-exp(-2.0*xif*MDt);
Gx=exp(xif*MDt)-2.0*xif*MDt-exp(-xif*MDt);
Ex=-Cx*(-2.0*xif*MDt-3.0+4.0*exp(xif*MDt)-exp(2.0*xif*MDt))-SQR(Gx);
}


//Main function

//int main(void)
int main(int argc, char *argv[])
{
  char command[200];
  // Multiple Particle Collision Dynamics Program 
    printf("Starting!\n");
  Counter=atoi(argv[1]);
  NUM_THREADS=atoi(argv[2]);
 (void) omp_set_num_threads(NUM_THREADS);//tell the code how many threads will be used
  N_threads=NUM_THREADS;//in the following parallelization sections
  epsilon=atof(argv[3]);//charged particles

  InitializeRandomNumberGenerator(time(0l));
printf("how many threads? %d\n",N_threads);//here we just set number of threads, program starts with #pragma...
  Readdat();//update restart data
  Preinit();//initialize parameters every (re)start program
 printf("Temperature: %f\n",Temperature);

if(RunStep==0)
{
//sprintf(command,"mkdir data_run%d",Counter);
sprintf(command,"mkdir /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d",NumberOfColloids,Counter);
system(command);
Init();
BDloop();
} 
else 
{
initial_pos_vel_restart(Counter,(RunStep));
BDloop();
}
//restart procedure?  
  return (0); 
}
