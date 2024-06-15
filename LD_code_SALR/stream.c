#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
#include "random.h"
#include <omp.h>
 
// Calculate Collision positions
void Stream(void)
{
  int i;
 // int Null;
 // Null=0;
  double WallTimez,WallTimey;
  //VECTOR dr;
#pragma omp parallel private(i,WallTimey,WallTimez)
{
//printf("how many threads?:%d\n",omp_get_num_threads());

   #pragma omp for

for(i=0;i<NumberOfParticles;i++){
//streaming all of the out-of-list solvent particles 
  if(flag[i]==0){
    //Null++;
    Positions[i].x +=Deltat*Velocities[i].x;
    Positions[i].y +=Deltat*Velocities[i].y;
    Positions[i].z +=Deltat*Velocities[i].z;
  // Bounce back rules
///////////////////////////////////////////////////////////////////////////////
  if((Positions[i].z<0.0)&&(Positions[i].y<BoxSize_y)&&(Positions[i].y>0.0))
     {
       WallTimez=Positions[i].z/Velocities[i].z;
       Positions[i].x -=2.0*Velocities[i].x*WallTimez;
       Positions[i].y -=2.0*Velocities[i].y*WallTimez;
       Positions[i].z -=2.0*Velocities[i].z*WallTimez;

       Velocities[i].x *=-1.0;
       Velocities[i].y *=-1.0;
       Velocities[i].z *=-1.0;
     }
   if((Positions[i].z>BoxSize_z)&&(Positions[i].y<BoxSize_y)&&(Positions[i].y>0.0))
     {
       WallTimez=(Positions[i].z-BoxSize_z)/Velocities[i].z;
       Positions[i].x -=2.0*Velocities[i].x*WallTimez;
       Positions[i].y -=2.0*Velocities[i].y*WallTimez;
       Positions[i].z -=2.0*Velocities[i].z*WallTimez;

       Velocities[i].x *=-1.0;
       Velocities[i].y *=-1.0;
       Velocities[i].z *=-1.0;
     }

/////////////////////////////////////////////////////////////////////////////////
  if((Positions[i].y<0.0)&&(Positions[i].z<BoxSize_z)&&(Positions[i].z>0.0))
     {   
       WallTimey=Positions[i].y/Velocities[i].y;
       Positions[i].x -=2.0*Velocities[i].x*WallTimey;
       Positions[i].y -=2.0*Velocities[i].y*WallTimey;
       Positions[i].z -=2.0*Velocities[i].z*WallTimey;

       Velocities[i].x *=-1.0;
       Velocities[i].y *=-1.0;
       Velocities[i].z *=-1.0;
     }   
   if((Positions[i].y>BoxSize_y)&&(Positions[i].z<BoxSize_z)&&(Positions[i].z>0.0))
     {   
       WallTimey=(Positions[i].y-BoxSize_y)/Velocities[i].y;
       Positions[i].x -=2.0*Velocities[i].x*WallTimey;
       Positions[i].y -=2.0*Velocities[i].y*WallTimey;
       Positions[i].z -=2.0*Velocities[i].z*WallTimey;

       Velocities[i].x *=-1.0;
       Velocities[i].y *=-1.0;
       Velocities[i].z *=-1.0;
     }   
////////////////////////////////////////////////////////////////////////

   if((Positions[i].y<0.0)&&(Positions[i].z<0.0))
     {
       WallTimey=(Positions[i].y)/Velocities[i].y;
       WallTimez=(Positions[i].z)/Velocities[i].z;
       Positions[i].x -=2.0*Velocities[i].x*MAX(WallTimey,WallTimez);
       Positions[i].y -=2.0*Velocities[i].y*MAX(WallTimey,WallTimez);
       Positions[i].z -=2.0*Velocities[i].z*MAX(WallTimey,WallTimez);

       Velocities[i].x *=-1.0;
       Velocities[i].y *=-1.0;
       Velocities[i].z *=-1.0;
     }

   if((Positions[i].y<0.0)&&(Positions[i].z>BoxSize_z))
     {   
       WallTimey=(Positions[i].y)/Velocities[i].y;
       WallTimez=(Positions[i].z-BoxSize_z)/Velocities[i].z;
       Positions[i].x -=2.0*Velocities[i].x*MAX(WallTimey,WallTimez);
       Positions[i].y -=2.0*Velocities[i].y*MAX(WallTimey,WallTimez);
       Positions[i].z -=2.0*Velocities[i].z*MAX(WallTimey,WallTimez);

       Velocities[i].x *=-1.0;
       Velocities[i].y *=-1.0;
       Velocities[i].z *=-1.0;
     }   

   if((Positions[i].y>BoxSize_y)&&(Positions[i].z<0.0))
     {   
       WallTimey=(Positions[i].y-BoxSize_y)/Velocities[i].y;
       WallTimez=(Positions[i].z)/Velocities[i].z;
       Positions[i].x -=2.0*Velocities[i].x*MAX(WallTimey,WallTimez);
       Positions[i].y -=2.0*Velocities[i].y*MAX(WallTimey,WallTimez);
       Positions[i].z -=2.0*Velocities[i].z*MAX(WallTimey,WallTimez);

       Velocities[i].x *=-1.0;
       Velocities[i].y *=-1.0;
       Velocities[i].z *=-1.0;
     }   

   if((Positions[i].y>BoxSize_y)&&(Positions[i].z>BoxSize_z))
     {   
       WallTimey=(Positions[i].y-BoxSize_y)/Velocities[i].y;
       WallTimez=(Positions[i].z-BoxSize_z)/Velocities[i].z;
       Positions[i].x -=2.0*Velocities[i].x*MAX(WallTimey,WallTimez);
       Positions[i].y -=2.0*Velocities[i].y*MAX(WallTimey,WallTimez);
       Positions[i].z -=2.0*Velocities[i].z*MAX(WallTimey,WallTimez);

       Velocities[i].x *=-1.0;
       Velocities[i].y *=-1.0;
       Velocities[i].z *=-1.0;
     }   

   //Apply PBC
}
    Positions[i].x -= BoxSize_x*(floor(Positions[i].x/BoxSize_x));if(Positions[i].x==BoxSize_x){Positions[i].x=0.0;}
    Positions[i].y -= BoxSize_y*(floor(Positions[i].y/BoxSize_y));if(Positions[i].y==BoxSize_y){Positions[i].y=0.0;}
    Positions[i].z -= BoxSize_z*(floor(Positions[i].z/BoxSize_z));if(Positions[i].z==BoxSize_z){Positions[i].z=0.0;}
// even though you put a wall at z direction, this expression is still valid
    }
  }
//printf("StreamingThis=%d\tNumberOfParticles=%d\n",Null,NumberOfParticles);
}

void ReScale(void)
{
int i;
VECTOR Momentum;
double scale;
//double TemperatureOfSlab[NumberOfCells_z];
  Momentum.x=0.0;
  Momentum.y=0.0;
  Momentum.z=0.0;
  //Uold=0.0;
  UKinetic=0.0;
  for(i=0;i<NumberOfParticles;i++)
  {
    Momentum.x+=Velocities[i].x;
    Momentum.y+=Velocities[i].y;
    Momentum.z+=Velocities[i].z;
  }
  Momentum.x/=(double)NumberOfParticles;//
  Momentum.y/=(double)NumberOfParticles;
  Momentum.z/=(double)NumberOfParticles;
  // calculate the kinetic energy UKinetic
  for(i=0;i<NumberOfParticles;i++)
  {
    Velocities[i].x-=Momentum.x;
    Velocities[i].y-=Momentum.y;
    Velocities[i].z-=Momentum.z;

    UKinetic+=(SQR(Velocities[i].x)+SQR(Velocities[i].y)+SQR(Velocities[i].z));
  }
  // scale all velocities to the correct temperature
  scale=sqrt(Temperature*(3.0*NumberOfParticles-3.0)/UKinetic);
  for(i=0;i<NumberOfParticles;i++)
  {
    Velocities[i].x*=scale;
    Velocities[i].y*=scale;
    Velocities[i].z*=scale;
  }
}
