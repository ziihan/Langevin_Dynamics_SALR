#include <stdio.h>
#include <stdlib.h>
//#include <stdint.h>
#include <math.h>
#include "system.h"
#include "random.h"
#include <omp.h>
//
//
//Force calculation for the direct interaction of Brownian particles
//Integrating the many-particle Langevin Equation (inertia effects included)
//
// Force calculation
 void Force(void)
{
int thread,ithread;
int i,j,jj,k;
int cod;
// int StartOfNeighbour;
double R2,Rp,Rp2,Fcc,Fs,R2i,R50i,R6i;
double r2,Ff,r6i,r48i,r2i;
double sqrcc;
//double Ff,ri,r3i,r2;
VECTOR dr;
VECTOR dR;

UPotential=0.0;

USpring=0.0;

cod=ceil((double)(NumberOfColloids)/(double)(N_threads));
#pragma omp parallel\
 private(j)
{
#pragma omp for

for (j=0;j<NumberOfColloids;j++)
{
ColloidForces[j].x=0.0;
ColloidForces[j].y=0.0;
ColloidForces[j].z=0.0;
}
}
//colloid-colloid interactions, the first step is brutal force
// private(j,jj,R2,R2i,R6i,Fcc,dR,sqrcc,Fs)
#pragma omp parallel\
 private(j,jj,R2,R2i,R50i,Fcc,dR,sqrcc,Rp)
{
#pragma omp for reduction (+:UPotential)
for (j=0;j<(NumberOfColloids-1);j++) {
   for(jj=j+1;jj<NumberOfColloids;jj++){
      dR.x=ColloidPositions[j].x-ColloidPositions[jj].x;
      dR.y=ColloidPositions[j].y-ColloidPositions[jj].y;
      dR.z=ColloidPositions[j].z-ColloidPositions[jj].z;

      if (dR.x>=0.5*BoxSize_x)dR.x-=BoxSize_x;
      else if(dR.x<-0.5*BoxSize_x)dR.x+=BoxSize_x;
      if (dR.y>=0.5*BoxSize_y)dR.y-=BoxSize_y;
      else if(dR.y<-0.5*BoxSize_y)dR.y+=BoxSize_y;
      if (dR.z>=0.5*BoxSize_z)dR.z-=BoxSize_z;
      else if(dR.z<-0.5*BoxSize_z)dR.z+=BoxSize_z;
      R2=SQR(dR.x)+SQR(dR.y)+SQR(dR.z);
//printf("Vx=%lf\n",ColloidVelocities[j].x);
//if (jj!=(j+NumberChainColloids))
//{

if(R2<SQR(CutOffYK)){
Rp=sqrt(R2);

        R2i=SQR(Rcc)/R2;
/*
        R6i=pow(R2i,3.0);
        UPotential +=4.0*R6i*(R6i-1.0)+Eshif;
        Fcc=8.0*6.0*R6i*(R6i-0.5);
*/
        R50i=pow(R2i,25.0);
        UPotential +=(epsilon*4.0*R50i*(R50i-1.0)+(qq*DebyeLength/Rp)*exp(-Rp/DebyeLength)-EcutYK);
        
        Fcc=(epsilon*400.0*R50i*(R50i-0.5)+qq*DebyeLength*exp(-Rp/DebyeLength)*(1.0/Rp+1.0/DebyeLength))*(1.0/Rp);
        Fcc=Fcc/Rp;//pair wise force action on colloid (net force)

        ColloidForces[j].x +=Fcc*dR.x;
        ColloidForces[j].y +=Fcc*dR.y;
        ColloidForces[j].z +=Fcc*dR.z;
        //reaction on mpc solvent particles
        ColloidForces[jj].x -=Fcc*dR.x;
        ColloidForces[jj].y -=Fcc*dR.y;
        ColloidForces[jj].z -=Fcc*dR.z;
                    }

//}
/*
else if (jj==(j+NumberChainColloids))
{
    sqrcc=sqrt(R2);
    Fs=-ColloidSpring*(sqrcc-(RodLength));//initial length RodLength
    Fs=Fs/sqrcc;

    ColloidForces[j].x += Fs*dR.x;
    ColloidForces[j].y += Fs*dR.y;
    ColloidForces[j].z += Fs*dR.z;

    ColloidForces[jj].x -= Fs*dR.x;
    ColloidForces[jj].y -= Fs*dR.y;
    ColloidForces[jj].z -= Fs*dR.z;

    USpring +=0.5*ColloidSpring*SQR(sqrcc-(RodLength));
}
*/
}
}
//printf("USpring=%f\n",USpring);
}

//implementation of walls pay more attention on the sign of the forces
double Fwy,yi,y6i;
double Fwz,zi,z6i;
UWalls=0.0;
/*
#pragma omp parallel\
 private(j,yi,y6i,Fwy,zi,z6i,Fwz)
{
#pragma omp for reduction (+:UWalls)

for(j=0;j<NumberOfColloids;j++)
   {
   if((0.0<(ColloidPositions[j].y))&&((ColloidPositions[j].y)<CutOffw)){
       yi=0.5*BoxSize_y/(ColloidPositions[j].y);
       y6i=CUBE(yi)*CUBE(yi);
       UWalls +=4.0*y6i*(y6i-1.0)+Eshif;
     Fwy=48.0*y6i*(y6i-0.5);
     Fwy=Fwy/(ColloidPositions[j].y);
     ColloidForces[j].y +=Fwy;
    }

   if((0.0<(BoxSize_y-ColloidPositions[j].y))&&((BoxSize_y-ColloidPositions[j].y)<CutOffw)){
       yi=0.5*BoxSize_y/(BoxSize_y-ColloidPositions[j].y);
       y6i=CUBE(yi)*CUBE(yi);
       UWalls +=4.0*y6i*(y6i-1.0)+Eshif;
     Fwy=-48.0*y6i*(y6i-0.5);
     Fwy=Fwy/(BoxSize_y-ColloidPositions[j].y);
     ColloidForces[j].y +=Fwy;
    }

   if((0.0<(ColloidPositions[j].z-Wall0))&&((ColloidPositions[j].z-Wall0)<CutOffcc)){
       zi=1.5*Rcc/(ColloidPositions[j].z-Wall0);
       z6i=CUBE(zi)*CUBE(zi);
       UWalls +=4.0*z6i*(z6i-1.0)+Eshif;
     Fwz=48.0*z6i*(z6i-0.5);
     Fwz=Fwz/(ColloidPositions[j].z-Wall0);
     ColloidForces[j].z +=Fwz;
    }

   if((0.0<(Wallz-ColloidPositions[j].z))&&((Wallz-ColloidPositions[j].z)<CutOffcc)){
       zi=1.5*Rcc/(Wallz-ColloidPositions[j].z);
       z6i=CUBE(zi)*CUBE(zi);
       UWalls +=4.0*z6i*(z6i-1.0)+Eshif;
     Fwz=-48.0*z6i*(z6i-0.5);
     Fwz=Fwz/(Wallz-ColloidPositions[j].z);
     ColloidForces[j].z +=Fwz;
    }

    }
   //  printf("UPotential2=%f\n",UWalls);
}
*/
//Applying colloids harmonic potential
//   USpring=0.0;
/*
for(j=0;j<(NumberOfColloids-1);j++)
  {
    VECTOR drcc;
    double rcc,Fc,sqrcc;
    drcc.x=ColloidPositions[j+1].x-ColloidPositions[j].x;//-2.0*Rcol;
    drcc.y=ColloidPositions[j+1].y-ColloidPositions[j].y;//-2.0*Rcol;
    drcc.z=ColloidPositions[j+1].z-ColloidPositions[j].z;//-2.0*Rcol;
    
      if (drcc.x>0.5*BoxSize_x)drcc.x-=BoxSize_x;
      else if(drcc.x<-0.5*BoxSize_x)drcc.x+=BoxSize_x;
      if (drcc.y>0.5*BoxSize_y)drcc.y-=BoxSize_y;
      else if(drcc.y<-0.5*BoxSize_y)drcc.y+=BoxSize_y;
      if (drcc.z>0.5*BoxSize_z)drcc.z-=BoxSize_z;
      else if(drcc.z<-0.5*BoxSize_z)drcc.z+=BoxSize_z;

    rcc=SQR(drcc.x)+SQR(drcc.y)+SQR(drcc.z);
    sqrcc=sqrt(rcc);
    Fc=-ColloidSpring*(sqrcc-(BondLength));//initial length
    Fc=Fc/sqrcc;

    ColloidForces[j].x -= Fc*drcc.x;
    ColloidForces[j].y -= Fc*drcc.y;
    ColloidForces[j].z -= Fc*drcc.z;

    ColloidForces[j+1].x += Fc*drcc.x;
    ColloidForces[j+1].y += Fc*drcc.y;
    ColloidForces[j+1].z += Fc*drcc.z;

    USpring +=0.5*ColloidSpring*SQR(sqrcc-(BondLength));
}
    // printf("UPotential3=%f\n",USpring);
//Applying colloids bending potential
*/
//Bond angle potential coding is from "The Art of Molecular Dynamics Simulation"

   UBending=0.0;
/*
for(j=0;j<(NumberOfColloids-2);j++)
  {
    VECTOR drbc1,drbc2,w1,w2;
    double Fb,c11,c12,c22,c,cd;

    drbc2.x=(ColloidPositions[j+2].x-ColloidPositions[j+1].x);
    drbc1.x=(ColloidPositions[j+1].x-ColloidPositions[j].x);
    drbc2.y=(ColloidPositions[j+2].y-ColloidPositions[j+1].y);
    drbc1.y=(ColloidPositions[j+1].y-ColloidPositions[j].y);
    drbc2.z=(ColloidPositions[j+2].z-ColloidPositions[j+1].z);
    drbc1.z=(ColloidPositions[j+1].z-ColloidPositions[j].z);

      if (drbc1.x>0.5*BoxSize_x)drbc1.x-=BoxSize_x;
      else if(drbc1.x<-0.5*BoxSize_x)drbc1.x+=BoxSize_x;
      if (drbc1.y>0.5*BoxSize_y)drbc1.y-=BoxSize_y;
      else if(drbc1.y<-0.5*BoxSize_y)drbc1.y+=BoxSize_y;
     if (drbc1.z>0.5*BoxSize_z)drbc1.z-=BoxSize_z;
      else if(drbc1.z<-0.5*BoxSize_z)drbc1.z+=BoxSize_z;

      if (drbc2.x>0.5*BoxSize_x)drbc2.x-=BoxSize_x;
      else if(drbc2.x<-0.5*BoxSize_x)drbc2.x+=BoxSize_x;
      if (drbc2.y>0.5*BoxSize_y)drbc2.y-=BoxSize_y;
      else if(drbc2.y<-0.5*BoxSize_y)drbc2.y+=BoxSize_y;
     if (drbc2.z>0.5*BoxSize_z)drbc2.z-=BoxSize_z;
      else if(drbc2.z<-0.5*BoxSize_z)drbc2.z+=BoxSize_z;
 
    c11=SQR(drbc1.x)+SQR(drbc1.y)+SQR(drbc1.z);
 //   sqrbc1=sqrt(c11);

    c12=drbc1.x*drbc2.x+drbc1.y*drbc2.y+drbc1.z*drbc2.z;
 
    c22=SQR(drbc2.x)+SQR(drbc2.y)+SQR(drbc2.z);
 //   sqrbc2=sqrt(c22);

    cd=sqrt(c11*c22);

    c=c12/cd;
//printf("Cosine=%f\n",c);
    Fb=-ColloidBending*(c-cCon);//bending force component

    w1.x=(c12/c11)*drbc1.x-1.0*drbc2.x;
    w1.y=(c12/c11)*drbc1.y-1.0*drbc2.y;
    w1.z=(c12/c11)*drbc1.z-1.0*drbc2.z;

    w1.x *=(Fb/cd);
    w1.y *=(Fb/cd);
    w1.z *=(Fb/cd);

    w2.x=1.0*drbc1.x-(c12/c22)*drbc2.x;
    w2.y=1.0*drbc1.y-(c12/c22)*drbc2.y;
    w2.z=1.0*drbc1.z-(c12/c22)*drbc2.z;

    w2.x *=(Fb/cd);
    w2.y *=(Fb/cd);
    w2.z *=(Fb/cd);

    ColloidForces[j].x +=w1.x;
    ColloidForces[j].y +=w1.y;
    ColloidForces[j].z +=w1.z;

    ColloidForces[j+1].x -=w1.x;
    ColloidForces[j+1].y -=w1.y;
    ColloidForces[j+1].z -=w1.z;

    ColloidForces[j+1].x -=w2.x;
    ColloidForces[j+1].y -=w2.y;
    ColloidForces[j+1].z -=w2.z;

    ColloidForces[j+2].x +=w2.x;
    ColloidForces[j+2].y +=w2.y;
    ColloidForces[j+2].z +=w2.z;

    UBending +=0.5*ColloidBending*SQR(c-cCon);
}
*/
 //    printf("UPotential4=%f\n",UBending);
UPotential +=UWalls+USpring+UBending;
}

////////////////////////////////////////////////////
//////BROWNIAN DYNAMICS INTEGRATION/////////////////////////////
///////////////////////////////////////////////////
//Brownian Dynamics require effective solvent (friction), hydrodynamics and memory effects
void Integrate(int Switch)
{
int k,j,i;
int Ing=1;
VECTOR g1;
VECTOR g2;
VECTOR randomF;
VECTOR randomV;
VECTOR randomV1;
VECTOR randomV2;
VECTOR delta_r;
if(Switch==0){Ing=0;}//do nothing
else if(Switch==1){
  for(j=0;j<NumberOfColloids;j++)//loop twin colloids
  {
    g1.x=RandomGaussianNumber;
    g1.y=RandomGaussianNumber;
    g1.z=RandomGaussianNumber;

    g2.x=RandomGaussianNumber;
    g2.y=RandomGaussianNumber;
    g2.z=RandomGaussianNumber;

//random force calculations   randomF and randomV
   randomF.x=sqrt((Temperature/(Mcol*xif*xif))*Cx)*g1.x;
   randomF.y=sqrt((Temperature/(Mcol*xif*xif))*Cx)*g1.y;
   randomF.z=sqrt((Temperature/(Mcol*xif*xif))*Cx)*g1.z;

   randomV1.x=sqrt((Temperature/(Mcol*xif*xif))*(Ex/Cx))*g2.x;
   randomV1.y=sqrt((Temperature/(Mcol*xif*xif))*(Ex/Cx))*g2.y;
   randomV1.z=sqrt((Temperature/(Mcol*xif*xif))*(Ex/Cx))*g2.z;

   randomV2.x=(Gx/Cx)*randomF.x;
   randomV2.y=(Gx/Cx)*randomF.y;
   randomV2.z=(Gx/Cx)*randomF.z;

   VAdd(randomV,randomV1,randomV2);

//updating the positions of colloids, akin to Velocity-Verlet algorithm
   delta_r.x =c0*ColloidVelocities[j].x + c1*ColloidForces[j].x+randomF.x;
   delta_r.y =c0*ColloidVelocities[j].y + c1*ColloidForces[j].y+randomF.y;
//   delta_r.z =c0*ColloidVelocities[j].z + c1*ColloidForces[j].z+randomF.z;

   ColloidPositions[j].x +=delta_r.x;
   ColloidPositions[j].y +=delta_r.y;
//   ColloidPositions[j].z +=delta_r.z;


    if(ColloidPositions[j].x>BoxSize_x)
   {
    CountBox_x[j]++;
   }
   if(ColloidPositions[j].x<0.0)
   {
   CountBox_x[j]--;
   }
   if(ColloidPositions[j].y>BoxSize_y)
   {
    CountBox_y[j]++;
   }
   if(ColloidPositions[j].y<0.0)
   {
   CountBox_y[j]--;
   }
/*
    if(ColloidPositions[j].z>BoxSize_z)
   {
    CountBox_z[j]++;
   }
   if(ColloidPositions[j].z<0.0)
   {
   CountBox_z[j]--;
   }
*/
//Apply PBC
   ColloidPositions[j].x -= BoxSize_x*(floor(ColloidPositions[j].x/BoxSize_x));if(ColloidPositions[j].x==BoxSize_x){ColloidPositions[j].x=0.0;}
   ColloidPositions[j].y -= BoxSize_y*(floor(ColloidPositions[j].y/BoxSize_y));if(ColloidPositions[j].y==BoxSize_y){ColloidPositions[j].y=0.0;}
//   ColloidPositions[j].z -= BoxSize_z*(floor(ColloidPositions[j].z/BoxSize_z));if(ColloidPositions[j].z==BoxSize_z){ColloidPositions[j].z=0.0;}

//updating the first half of the velocity
   ColloidVelocities[j].x = c2*(delta_r.x+randomV.x);
   ColloidVelocities[j].y = c2*(delta_r.y+randomV.y);
//   ColloidVelocities[j].z = c2*(delta_r.z+randomV.z);
 
}
}
//calculate the NEW force
Force();//the force at t+dt
if(Switch==0){Ing=0;}
else if(Switch==1){
//update the second-half-time Velocities
for (j=0;j<NumberOfColloids;j++){
   ColloidVelocities[j].x +=c2c3*ColloidForces[j].x;
   ColloidVelocities[j].y +=c2c3*ColloidForces[j].y;
//   ColloidVelocities[j].z +=c2c3*ColloidForces[j].z;
   }
}
}
