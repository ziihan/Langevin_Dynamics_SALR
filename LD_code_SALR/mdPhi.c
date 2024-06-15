#include <stdio.h>
#include <stdlib.h>
//#include <stdint.h>
#include <math.h>
#include "system.h"
#include "random.h"
#include <omp.h>

void Neighbours(void)
{
int i,j,tempX,tempY,tempZ;
int nNeighbour,NiCell;
int StartNeighbour;
double ro2;
int p1,p2,p3;
//int label1=0;
int p;

VECTOR dri;
GridCells(0);//sorting particles

for (p=0;p<NumberOfParticles;p++)
{
flag[p]=0;
}

for (j=0;j<NumberOfColloids;j++){
 StartOfNeighbour[j]=0;//the sequence of the Starting Neighbour Index
}

iNeighbour=0;
StartNeighbour=0;//the sequence of the Starting Neighbour Index

for (j=0;j<NumberOfColloids;j++)
{

StartOfNeighbour[j]=StartNeighbour;
//From StartOfNeighbour[j]+1 to StartOfNeighbour[j+1]-1
if(j<NumberPhiColloids){
// printf("Good!\n");
// StartOfNeighbour[j]=StartNeighbour;
  nNeighbour=0;//neighbour counter
//printf("Cx=%d\tCy=%d\tCy=%d\n",Nix[j],Niy[j],Niz[j]);

    for(p3=-Rcelli;p3<=Rcelli;p3++){
     for(p2=-Rcelli;p2<=Rcelli;p2++){
      for(p1=-Rcelli;p1<=Rcelli;p1++){
      tempX=Nix[j]+p1;
      tempY=Niy[j]+p2;
      tempZ=Niz[j]+p3;
//printf("Cx=%d\tCy=%d\tCy=%d\n",Nix[j],Niy[j],Niz[j]);
      if(tempX<0) tempX +=NumberOfCells_x;
      else if(tempX>=NumberOfCells_x) tempX -=NumberOfCells_x;
      if(tempY<0) tempY +=NumberOfCells_y;
      else if(tempY>=NumberOfCells_y) tempY -=NumberOfCells_y;
      if(tempZ<0) tempZ +=NumberOfCells_z;
      else if(tempZ>=NumberOfCells_z) tempZ -=NumberOfCells_z;
    NiCell=tempX+tempY*NumberOfCells_x+tempZ*NumberOfCells_x*NumberOfCells_y;
//printf("Good,j=%d\tNiX=%d\tNiY=%d\tNiZ=%d\n",j,Nix[j],Niy[j],Niz[j]);
    i=hoc[NiCell];
//printf("Good!!\n");
//distinguishes the colloids belongs to which cell  
  while(i!=-1)         //i==-1 ?         /*  particles in cell */
        {
         if (i < NumberOfParticles) {
 dri.x=Positions[i].x-ColloidPositions[j].x;
 dri.y=Positions[i].y-ColloidPositions[j].y;
 dri.z=Positions[i].z-ColloidPositions[j].z;
//Apply PBCs no overlapping
 if(dri.x>0.5*BoxSize_x)dri.x -=BoxSize_x;
 else if(dri.x<-0.5*BoxSize_x)dri.x +=BoxSize_x;
 if (dri.y>0.5*BoxSize_y)dri.y -=BoxSize_y;
 else if(dri.y<-0.5*BoxSize_y)dri.y +=BoxSize_y;
 if (dri.z>0.5*BoxSize_z)dri.z -=BoxSize_z;
 else if(dri.z<-0.5*BoxSize_z)dri.z +=BoxSize_z;

 ro2=SQR(dri.x)+SQR(dri.y)+SQR(dri.z);
 if (ro2<SQR(OutOffi)){
  Neighbour[StartNeighbour+nNeighbour]=i;//store the neighbour index 
  nNeighbour++;
  flag[i]++;//particles in the list-if flag larger than 1 - count once in integrate array
  if(flag[i]==1){
  IntegrateNeighbour[iNeighbour]=i;
  iNeighbour++;
    }
//  label1++;
      }
       }
i=ll[i];
        }
}
}
}}
if(j>=NumberPhiColloids){
//StartOfNeighbour[j]=StartNeighbour;
  nNeighbour=0;//neighbour counter
    for(p3=-Rcello;p3<=Rcello;p3++){
     for(p2=-Rcello;p2<=Rcello;p2++){
      for(p1=-Rcello;p1<=Rcello;p1++){
      tempX=Nix[j]+p1;
      tempY=Niy[j]+p2;
      tempZ=Niz[j]+p3;
//printf("Cx=%d\tCy=%d\tCy=%d\n",Nix[j],Niy[j],Niz[j]);
      if(tempX<0) tempX +=NumberOfCells_x;
      else if(tempX>=NumberOfCells_x) tempX -=NumberOfCells_x;
      if(tempY<0) tempY +=NumberOfCells_y;
      else if(tempY>=NumberOfCells_y) tempY -=NumberOfCells_y;
      if(tempZ<0) tempZ +=NumberOfCells_z;
      else if(tempZ>=NumberOfCells_z) tempZ -=NumberOfCells_z;
    NiCell=tempX+tempY*NumberOfCells_x+tempZ*NumberOfCells_x*NumberOfCells_y;
    i=hoc[NiCell];
//distinguishes the colloids belongs to which cell
  while(i!=-1)         //i==-1 ?         /*  particles in cell */
        {
         if (i < NumberOfParticles) {
 dri.x=Positions[i].x-ColloidPositions[j].x;
 dri.y=Positions[i].y-ColloidPositions[j].y;
 dri.z=Positions[i].z-ColloidPositions[j].z;
//Apply PBCs no overlapping
 if(dri.x>0.5*BoxSize_x)dri.x -=BoxSize_x;
 else if(dri.x<-0.5*BoxSize_x)dri.x +=BoxSize_x;
 if (dri.y>0.5*BoxSize_y)dri.y -=BoxSize_y;
 else if(dri.y<-0.5*BoxSize_y)dri.y +=BoxSize_y;
 if (dri.z>0.5*BoxSize_z)dri.z -=BoxSize_z;
 else if(dri.z<-0.5*BoxSize_z)dri.z +=BoxSize_z;

 ro2=SQR(dri.x)+SQR(dri.y)+SQR(dri.z);
 if (ro2<SQR(OutOffo)){
  Neighbour[StartNeighbour+nNeighbour]=i;//store the neighbour index
  nNeighbour++;
  flag[i]++;//particles in the list-if flag larger than 1 - count once in integrate array
  if(flag[i]==1){
  IntegrateNeighbour[iNeighbour]=i;
  iNeighbour++;
    }
//  label1++;
      }
       }
i=ll[i];
        }
}
}
}}

NumberOfNeighbour[j]=nNeighbour;//Number of neighbours of jth particles(colloids)
StartNeighbour += nNeighbour;//the Starting Neighbour Index of the next particle(colloid)
}
//printf("NeighbourCount=%d\n",label1);
}



/*
void Neighbours(void)
{
int i,j,tempX,tempY,tempZ;
int nNeighbour,NiCell;
int StartNeighbour;
double ro2;
int p1,p2,p3;
//int label1=0;
int p;
int cod;

VECTOR dri;
GridCells(0);//sorting particles
#pragma omp parallel private(p) 
{
   #pragma omp for
 for (p=0;p<NumberOfParticles;p++)
  {
   flag[p]=0;
  }
}

for (j=0;j<NumberOfColloids;j++){
StartOfNeighbour[j]=0;//the sequence of the Starting Neighbour Index
}

iNeighbour=0;
StartNeighbour=0;//the sequence of the Starting Neighbour Index

for (j=0;j<NumberOfColloids;j++)
{
//From StartOfNeighbour[j]+1 to StartOfNeighbour[j+1]-1
  StartOfNeighbour[j]=StartNeighbour;
  nNeighbour=0;//neighbour counter
    for(p3=-Rcelli;p3<=Rcelli;p3++){
     for(p2=-Rcelli;p2<=Rcelli;p2++){
      for(p1=-Rcelli;p1<=Rcelli;p1++){
      tempX=Nix[j]+p1;
      tempY=Niy[j]+p2;
      tempZ=Niz[j]+p3;
//printf("Cx=%d\tCy=%d\tCy=%d\n",Nix[j],Niy[j],Niz[j]);
      if(tempX<0) tempX +=NumberOfCells_x;
      else if(tempX>=NumberOfCells_x) tempX -=NumberOfCells_x;
      if(tempY<0) tempY +=NumberOfCells_y;
      else if(tempY>=NumberOfCells_y) tempY -=NumberOfCells_y;
      if(tempZ<0) tempZ +=NumberOfCells_z;
      else if(tempZ>=NumberOfCells_z) tempZ -=NumberOfCells_z;
    NiCell=tempX+tempY*NumberOfCells_x+tempZ*NumberOfCells_x*NumberOfCells_y;
    i=hoc[NiCell];
//distinguishes the colloids belongs to which cell  
  while(i!=-1)         //i==-1 ?    
        {
         if (i < NumberOfParticles) {
 dri.x=Positions[i].x-ColloidPositions[j].x;
 dri.y=Positions[i].y-ColloidPositions[j].y;
 dri.z=Positions[i].z-ColloidPositions[j].z;
//Apply PBCs no overlapping
 if(dri.x>0.5*BoxSize_x)dri.x -=BoxSize_x;
 else if(dri.x<-0.5*BoxSize_x)dri.x +=BoxSize_x;
 if (dri.y>0.5*BoxSize_y)dri.y -=BoxSize_y;
 else if(dri.y<-0.5*BoxSize_y)dri.y +=BoxSize_y;
 if (dri.z>0.5*BoxSize_z)dri.z -=BoxSize_z;
 else if(dri.z<-0.5*BoxSize_z)dri.z +=BoxSize_z;

 ro2=SQR(dri.x)+SQR(dri.y)+SQR(dri.z);
 if (ro2<SQR(OutOff)){
  Neighbour[StartNeighbour+nNeighbour]=i;//store the neighbour index 
  nNeighbour++;
  flag[i]++;//particles in the list-if flag larger than 1 - count once in integrate array
  if(flag[i]==1){
  IntegrateNeighbour[iNeighbour]=i;
  iNeighbour++;
    }
//  label1++;
      }
       }
i=ll[i];
        }
}
}
}
NumberOfNeighbour[j]=nNeighbour;//Number of neighbours of jth particles(colloids)
//#pragma omp barrier
StartNeighbour += nNeighbour;//the Starting Neighbour Index of the next particle(colloid)
}
//printf("NeighbourCount=%d\n",label1);
}
*/


// Force calculation
 void Force(void)
{
int thread,ithread;
int i,j,jj,k;
int cod;
// int StartOfNeighbour;
double R2,Fcc,R2i,R6i;
double r2,Ff,r6i,r48i,r2i;
//double Ff,ri,r3i,r2;
VECTOR dr;
VECTOR dR;

UPotential=0.0;
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
#pragma omp parallel\
 private(j,jj,R2,R2i,R6i,Fcc,dR)
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
if(R2<SQR(CutOffcc)){
        R2i=SQR(Rcc)/R2;
        R6i=CUBE(R2i);
        UPotential +=4.0*R6i*(R6i-1.0)+Eshif;
        Fcc=48.0*R6i*(R6i-0.5);
        Fcc=Fcc/R2;//pair wise force
//action on colloid (net force)
        ColloidForces[j].x +=Fcc*dR.x;
        ColloidForces[j].y +=Fcc*dR.y;
        ColloidForces[j].z +=Fcc*dR.z;
//reaction on mpc solvent particles
        ColloidForces[jj].x -=Fcc*dR.x;
        ColloidForces[jj].y -=Fcc*dR.y;
        ColloidForces[jj].z -=Fcc*dR.z;
}
}
}

}
#pragma omp parallel\
 private(i,k)
{
#pragma omp for
for(k=0;k<(iNeighbour);k++){
    i=IntegrateNeighbour[k];
     Forces[i].x = 0.0;
     Forces[i].y = 0.0;
     Forces[i].z = 0.0;
}
}
//#pragma omp barrier
#pragma omp parallel\
 private(thread,i,j,k,dr,r2,r6i,r48i,Ff)
{
//printf("n_threads=%d\n",omp_get_num_threads());
 thread=omp_get_thread_num();
#pragma omp for schedule(static, cod)\
 reduction(+:UPotential)
//sorting particles into cells
//  StartOfNeighbour=0;
  for(j=0;j<NumberOfColloids;j++)//loop twin colloids
  {
   for(k=StartOfNeighbour[j];k<(StartOfNeighbour[j]+NumberOfNeighbour[j]);k++){
         i=Neighbour[k];
  //Only interactions of Colloids-SolventParticles are considered
      dr.x=ColloidPositions[j].x-Positions[i].x;
      dr.y=ColloidPositions[j].y-Positions[i].y;
      dr.z=ColloidPositions[j].z-Positions[i].z;
//notice the directions
      // apply boundary conditions--Minimum image convension for the force calculation
      if (dr.x>=0.5*BoxSize_x)dr.x-=BoxSize_x;
      else if(dr.x<-0.5*BoxSize_x)dr.x+=BoxSize_x;
      if (dr.y>=0.5*BoxSize_y)dr.y-=BoxSize_y;
      else if(dr.y<-0.5*BoxSize_y)dr.y+=BoxSize_y;
      if (dr.z>=0.5*BoxSize_z)dr.z-=BoxSize_z;
      else if(dr.z<-0.5*BoxSize_z)dr.z+=BoxSize_z;
      r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      // check if the distance is within the cutoff radius
      if((r2<SQR(CutOffi))&&(j<NumberPhiColloids))
      {

        r2i=SQR(Rcol)/r2;
        r6i=CUBE(r2i);
        UPotential +=4.0*r6i*(r6i-1.0)+Eshif;
        Ff=48.0*r6i*(r6i-0.5);
        Ff=Ff/r2;//pair wise force

/*
        r2i=SQR(Rcol)/r2;//(sigma/R)-square
        r48i=pow(r2i,24.0);//

        UPotential +=4.0*r48i*(r48i-1.0)-Ecut;
        Ff=8.0*48.0*r48i*(r48i-0.5);
        Ff=Ff/r2;//pair wise force
*/
/*
        ri=Rcol/sqrt(r2);
        r3i=CUBE(ri);
        UPotential +=4.0*r3i*(r3i-1.0)+Eshif;
    // printf("UPotential=%15le\n",UPotential);
        Ff=24.0*r3i*(r3i-0.5);
        Ff=Ff/r2;//pair wise force
*/
/*
//action on colloid (net force)
        ColloidForces[j].x +=Ff*dr.x;
        ColloidForces[j].y +=Ff*dr.y;
        ColloidForces[j].z +=Ff*dr.z;

//reaction on mpc solvent particles
        Forces[i].x -=Ff*dr.x;
        Forces[i].y -=Ff*dr.y;
        Forces[i].z -=Ff*dr.z;
*/
//action on colloid (net force)
        ColloidForcesPriv[j][thread].x +=Ff*dr.x;
        ColloidForcesPriv[j][thread].y +=Ff*dr.y;
        ColloidForcesPriv[j][thread].z +=Ff*dr.z;
//reaction on mpc solvent particles
        ForcesPriv[i][thread].x -=Ff*dr.x;
        ForcesPriv[i][thread].y -=Ff*dr.y;
        ForcesPriv[i][thread].z -=Ff*dr.z;
//printf("ColloidForces[i].x:%lf\n,Forces[i].x:%lf\n",ColloidForces[i].x,Forces[i].x);
       }

   if((r2<SQR(CutOffo))&&(j>=NumberPhiColloids))
      {
        r2i=SQR(Rcol)/r2;//(sigma/R)-square
        r48i=pow(r2i,24.0);//

        UPotential +=4.0*r48i*(r48i-1.0)-Ecut;
        Ff=8.0*48.0*r48i*(r48i-0.5);
        Ff=Ff/r2;//pair wise force

/*
        r2i=SQR(Rcol)/r2;//(sigma/R)-square
        r12i=pow(r2i,6.0);//
        //48-96-aLJ Potential
        UPotential +=4.0*r12i*(r12i-1.0)+Eshif;
        Ff=8.0*12.0*r12i*(r12i-0.5);
        Ff=Ff/r2;//pair wise force
*/
/*
        r2i=SQR(Rcol)/r2;
        r6i=CUBE(r2i);
        r24i=r6i*r6i*r6i*r6i;

        UPotential +=4.0*r24i*(r24i-1.0)-Ecut;
        Ff=192.0*r24i*(r24i-0.5);
        Ff=Ff/r2;//pair wise force
*/
//action on colloid (net force)
        ColloidForcesPriv[j][thread].x +=Ff*dr.x;
        ColloidForcesPriv[j][thread].y +=Ff*dr.y;
        ColloidForcesPriv[j][thread].z +=Ff*dr.z;
//reaction on mpc solvent particles
        ForcesPriv[i][thread].x -=Ff*dr.x;
        ForcesPriv[i][thread].y -=Ff*dr.y;
        ForcesPriv[i][thread].z -=Ff*dr.z;
/*
//action on colloid (net force)
        ColloidForces[j].x +=Ff*dr.x;
        ColloidForces[j].y +=Ff*dr.y;
        ColloidForces[j].z +=Ff*dr.z;
//reaction on mpc solvent particles
        Forces[i].x -=Ff*dr.x;
        Forces[i].y -=Ff*dr.y;
        Forces[i].z -=Ff*dr.z;
*/
       }

      }
//StartOfNeighbour +=NumberOfNeighbour[j];
     }
}

for(j=0;j<NumberOfColloids;j++)
   {for (ithread=0;ithread<N_threads;ithread++)
     {
       ColloidForces[j].x +=ColloidForcesPriv[j][ithread].x;
       ColloidForces[j].y +=ColloidForcesPriv[j][ithread].y;
       ColloidForces[j].z +=ColloidForcesPriv[j][ithread].z;

       ColloidForcesPriv[j][ithread].x=0.0;
       ColloidForcesPriv[j][ithread].y=0.0;
       ColloidForcesPriv[j][ithread].z=0.0;
     }
   }
#pragma omp parallel private(k,i,ithread) 
{
   #pragma omp for
//think more!!!
   for(k=0;k<(iNeighbour);k++){
         i=IntegrateNeighbour[k];
   {for (ithread=0;ithread<N_threads;ithread++)
      {
       Forces[i].x +=ForcesPriv[i][ithread].x;
       Forces[i].y +=ForcesPriv[i][ithread].y;
       Forces[i].z +=ForcesPriv[i][ithread].z;

       ForcesPriv[i][ithread].x=0.0;
       ForcesPriv[i][ithread].y=0.0;
       ForcesPriv[i][ithread].z=0.0;
      }
    }
 }
}

  //   printf("UPotential1=%f\n",UPotential);
//printf("ColloidIndex%d\tCheckNumber=%d\n",StartOfNeighbour,Check);

//implementation of walls pay more attention on the sign of the forces
double Fwy,yi,y6i;
double Fwz,zi,z6i;
UWalls=0.0;

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
//Applying colloids harmonic potential
   USpring=0.0;
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

void Integrate(int Switch)
{
int k,j,i;
int Ing=1;
if(Switch==0){Ing=0;}//do nothing
else if(Switch==1){
  for(j=0;j<NumberOfColloids;j++)//loop twin colloids
  {
//updating the positions of colloids
   ColloidPositions[j].x +=MDt*ColloidVelocities[j].x + (SQR(MDt)*0.5*invMcol)*ColloidForces[j].x;
   ColloidPositions[j].y +=MDt*ColloidVelocities[j].y + (SQR(MDt)*0.5*invMcol)*ColloidForces[j].y;
   ColloidPositions[j].z +=MDt*ColloidVelocities[j].z + (SQR(MDt)*0.5*invMcol)*ColloidForces[j].z;
//Apply PBC
   ColloidPositions[j].x -= BoxSize_x*(floor(ColloidPositions[j].x/BoxSize_x));if(ColloidPositions[j].x==BoxSize_x){ColloidPositions[j].x=0.0;}
   ColloidPositions[j].y -= BoxSize_y*(floor(ColloidPositions[j].y/BoxSize_y));if(ColloidPositions[j].y==BoxSize_y){ColloidPositions[j].y=0.0;}
   ColloidPositions[j].z -= BoxSize_z*(floor(ColloidPositions[j].z/BoxSize_z));if(ColloidPositions[j].z==BoxSize_z){ColloidPositions[j].z=0.0;}

   ColloidVelocities[j].x +=(MDt*0.5*invMcol)*ColloidForces[j].x;
   ColloidVelocities[j].y +=(MDt*0.5*invMcol)*ColloidForces[j].y;
   ColloidVelocities[j].z +=(MDt*0.5*invMcol)*ColloidForces[j].z;
  }
}

#pragma omp parallel private(k,i) 
{
   #pragma omp for

   for(k=0;k<(iNeighbour);k++){
         i=IntegrateNeighbour[k];
   Positions[i].x +=MDt*Velocities[i].x + 0.5*SQR(MDt)*Forces[i].x;
   Positions[i].y +=MDt*Velocities[i].y + 0.5*SQR(MDt)*Forces[i].y;
   Positions[i].z +=MDt*Velocities[i].z + 0.5*SQR(MDt)*Forces[i].z;

   Positions[i].x -= BoxSize_x*(floor(Positions[i].x/BoxSize_x));if(Positions[i].x==BoxSize_x){Positions[i].x=0.0;}
   Positions[i].y -= BoxSize_y*(floor(Positions[i].y/BoxSize_y));if(Positions[i].y==BoxSize_y){Positions[i].y=0.0;}
   Positions[i].z -= BoxSize_z*(floor(Positions[i].z/BoxSize_z));if(Positions[i].z==BoxSize_z){Positions[i].z=0.0;}

   Velocities[i].x +=0.5*MDt*Forces[i].x;
   Velocities[i].y +=0.5*MDt*Forces[i].y;
   Velocities[i].z +=0.5*MDt*Forces[i].z;
       }
}
//calculate the NEW force
Force();//the force at t+dt
if(Switch==0){Ing=0;}
else if(Switch==1){
//update the second-half-time Velocities
for (j=0;j<NumberOfColloids;j++){
   ColloidVelocities[j].x +=(0.5*MDt*invMcol)*ColloidForces[j].x;
   ColloidVelocities[j].y +=(0.5*MDt*invMcol)*ColloidForces[j].y;
   ColloidVelocities[j].z +=(0.5*MDt*invMcol)*ColloidForces[j].z;
   }
}
//#pragma omp barrier

//loop all 27 neighbour&center cells, initialize their forces
#pragma omp parallel private(k,i)
{
   #pragma omp for
   for(k=0;k<(iNeighbour);k++){
         i=IntegrateNeighbour[k];
   Velocities[i].x +=(0.5*MDt)*Forces[i].x;
   Velocities[i].y +=(0.5*MDt)*Forces[i].y;
   Velocities[i].z +=(0.5*MDt)*Forces[i].z;
   }
//printf("that=%d\n",that);
}
}


void Drdf(void)
{
int i,j,kg;
VECTOR DR;
double R2;
for(j=0;j<NumberOfColloids;j++){
for(i=0;i<NumberOfParticles;i++){
DR.x=ColloidPositions[j].x-Positions[i].x;
DR.y=ColloidPositions[j].y-Positions[i].y;
DR.z=ColloidPositions[j].z-Positions[i].z;

      if (DR.x>0.5*BoxSize_x)DR.x-=BoxSize_x;
      else if(DR.x<-0.5*BoxSize_x)DR.x+=BoxSize_x;
      if (DR.y>0.5*BoxSize_y)DR.y-=BoxSize_y;
      else if(DR.y<-0.5*BoxSize_y)DR.y+=BoxSize_y;
      if (DR.z>0.5*BoxSize_z)DR.z-=BoxSize_z;
      else if(DR.z<-0.5*BoxSize_z)DR.z+=BoxSize_z;
      R2=SQR(DR.x)+SQR(DR.y)+SQR(DR.z);
      if((R2<=SQR(0.5*BoxSize_x))&&(R2<=SQR(3.0*Rcol))){
     kg= (int)(sqrt(R2)/BinSize);
     dg[kg]++;
       }
}
}
DensCount++;
}

void PrintDrdf(void)
{
double r,vb,nid;
FILE* BinDens;
char binname[100];
    sprintf(binname, "flow45field%d/SomeRDF%d-%d.dat", NumberOfColloids,Counter,Needle);
    BinDens=fopen(binname,"wt");
int k;
for (k=0;k<NumberOfBins;k++){
r=BinSize*(k+0.5);
vb=((double)(CUBE(k+1))-(double)(CUBE(k)))*CUBE(BinSize);
nid=(4.0/3.0)*PI*vb*Rho;
//dg[k]=(double)(dg[k])/((double)(DensCount*NumberOfColloids)*nid);
fprintf(BinDens,"%lf\t%lf\n",r,((double)dg[k])/(((double)(DensCount*NumberOfColloids))*nid));
}
fclose(BinDens);
} 
