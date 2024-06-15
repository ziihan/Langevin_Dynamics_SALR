#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "random.h"
#include "system.h"
// Langevin Dynamics loop and print output files
void BDloop(void)
{ 
   int j,k;
   double UKinetic_s;
   FILE *f;
   char filename[200];

   FILE *f0;
   char filename0[200];

   FILE *FilePtr;
   char moviename[200];
   sprintf(moviename,"Particlediffusion%d/BDtraj%d-StartStep-%d.xyz",NumberOfColloids,Counter,RunStep);
   FilePtr=fopen(moviename,"wt");

  sprintf(filename,"Particlediffusion%d/BD-%d-%d.dat",NumberOfColloids,Counter,RunStep);
  f=fopen(filename,"wt");

 //Initialize Diffusion Coefficient
   SampleDiff(INITIALIZE);
SampleRDF(INITIALIZE);   
//   Neighbours();
   Force();
for(Step=(RunStep+1);Step<=NumberOfSteps;Step++)
  {
  Integrate(1);
 if((Step%500==0)&&(Step>=500))
  {
 int i;
 UKinetic=0.0;
 UKinetic_s=0.0;
// calculate the kinetic energy UKinetic
 for(i=0;i<NumberOfColloids;i++)
      {
 UKinetic_s +=0.5*Mcol*(SQR(ColloidVelocities[i].x)+SQR(ColloidVelocities[i].y));
      }
 //UKinetic +=UKinetic_s;
 UTotal = UKinetic+UKinetic_s+UPotential;
 UTotal = 2.0*UTotal/(3.0*(double)(NumberOfColloids));
 UKinetic_s = 2.0*UKinetic_s/(3.0*(double)(NumberOfColloids));
 printf("%.1f\t%lf\t%lf\t%lf\n",(Step*MDt),UPotential,UTotal,UKinetic_s);
                                           //if(Step>RecordStep){
 fprintf(f,"%.1f\t%lf\t%lf\n",(Step*MDt),UTotal,UKinetic_s);
}
                                              //                   }
//                                                      }//print out
if(Step%10000==0){
 WritePdb(FilePtr);
}
if((Step%SampleSteps==0)&&(Step>=SampleSteps)){
//  WritePdb(FilePtr);
  SampleDiff(SAMPLE);
  SampleRDF(SAMPLE);
  }
//**********************Output data for restart****************//
//repeat the precedure like before  md force integrate, stream collision
if((Step>=(output_restart_steps))&&((Step%output_restart_steps)==0))
{
//SampleDiff(WRITE_RESULTS);
SampleRDF(WRITE_RESULTS);
output_pos_vel_restart(Counter,Step);
// fclose(f);
SampleDiff(WRITE_RESULTS);
}
}
 // SampleDiff(WRITE_RESULTS);
  fclose(f);
  fclose(FilePtr);//write a movie    
   printf("Simulation finished!\n");
}
