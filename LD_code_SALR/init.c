#include <stdio.h>
#include <stdlib.h>
//#include <stdbool.h>
#include <math.h> 
#include "system.h"
#include "random.h"

// Generating initial configurations: positions/velocities
void Init(void)
{
double dx,dy,dz;
int j,ithread,kj;
double RowLength=2.0*BondLength+(Nx-1)*BondLength;
double rdbeta,theta;
VECTOR tpPosition; //define the colloid centers
int di,dj,dk;
dx=BoxSize_x/((double)Nx);dy=BoxSize_y/((double)Ny);dz=BoxSize_z/((double)Nz);
for (j=0;j<NumberOfColloids;j++)
{
rdbeta=RandomNumber*2.0*PI;//an arbitary orientation
theta=acos(RandomNumber*2.0-1.0);
//theta=RandomNumber*PI;
di=j%Nx;
ColloidPositions[j].x=(di+0.5)*dx;
//PatchPositions[j].x=ColloidPositions[j].x-Rdipole*sin(theta)*cos(rdbeta);
//PatchPositions[j+NumberOfColloids].x=ColloidPositions[j].x+Rdipole*sin(theta)*cos(rdbeta);

dj=(j-di)/Nx%Ny;
ColloidPositions[j].y=(dj+0.5)*dy;
//PatchPositions[j].y=ColloidPositions[j].y-Rdipole*sin(theta)*sin(rdbeta);
//PatchPositions[j+NumberOfColloids].y=ColloidPositions[j].y+Rdipole*sin(theta)*sin(rdbeta);

dk=((j-di)/Nx-dj)/Ny;
ColloidPositions[j].z =(dk+0.5)*dz;
//PatchPositions[j].z = ColloidPositions[j].z-Rdipole*cos(theta);
//PatchPositions[j+NumberOfColloids].z = ColloidPositions[j].z+Rdipole*cos(theta);
}

for (j=0;j<NumberOfColloids;j++)
{
        CountBox_x[j]=0;
        CountBox_y[j]=0;
        CountBox_z[j]=0;

for(ithread=0;ithread<N_threads;ithread++){
 ColloidForcesPriv[j][ithread].x=0.0;
 ColloidForcesPriv[j][ithread].y=0.0;
 ColloidForcesPriv[j][ithread].z=0.0;
   }
}
//
    
  int i;//3 Dimensional
  double scale;
  VECTOR Momentum;

  Momentum.x=0.0;
  Momentum.y=0.0;
  Momentum.z=0.0;
  //Uold=0.0;
  UKinetic=0.0;
 
  for(i=0;i<NumberOfColloids;i++)
  {
/*
    ColloidVelocities[i].x=sqrt(Temperature/Mcol)*RandomGaussianNumber;
    ColloidVelocities[i].y=sqrt(Temperature/Mcol)*RandomGaussianNumber;
    ColloidVelocities[i].z=sqrt(Temperature/Mcol)*RandomGaussianNumber;
*/
    ColloidVelocities[i].x=RandomGaussianNumber;
    ColloidVelocities[i].y=RandomGaussianNumber;
    ColloidVelocities[i].z=RandomGaussianNumber;

    Momentum.x+=ColloidVelocities[i].x;
    Momentum.y+=ColloidVelocities[i].y;
    Momentum.z+=ColloidVelocities[i].z;
  }
 
  Momentum.x/=(double)NumberOfColloids;//
  Momentum.y/=(double)NumberOfColloids;
  Momentum.z/=(double)NumberOfColloids;

  // calculate the kinetic energy UKinetic
 
  for(i=0;i<NumberOfColloids;i++)
  {
    ColloidVelocities[i].x-=Momentum.x;
    ColloidVelocities[i].y-=Momentum.y;
    ColloidVelocities[i].z-=Momentum.z;
 
    UKinetic+=(SQR(ColloidVelocities[i].x)+SQR(ColloidVelocities[i].y)+SQR(ColloidVelocities[i].z));
  }
  // scale all velocities to the correct temperature
  scale=sqrt(Temperature*(3.0*NumberOfColloids-3.0)/UKinetic);
  for(i=0;i<NumberOfColloids;i++)
  {
    ColloidVelocities[i].x*=scale;
    ColloidVelocities[i].y*=scale;
    ColloidVelocities[i].z*=scale;
  }

printf("initializaiton checked!\n");
}

 
