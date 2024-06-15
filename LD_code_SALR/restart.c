#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

//Function one: writing restart binary files, save for the restart.
//Function two: reading restart files 
char filename[400],dirname[400],command[400];
FILE *file;

//writing binary files//
void output_pos_vel_restart(int run,int t)
{int i,j;
 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_coll_run%d.tmp",NumberOfColloids,run,t);
 file=fopen(filename,"w");

 for (j=0;j<NumberOfColloids;j++)
   {
fprintf(file,"%lf %lf %lf %lf %lf %lf %d %d %d\n",ColloidPositions[j].x,ColloidPositions[j].y,ColloidPositions[j].z,ColloidVelocities[j].x,ColloidVelocities[j].y,ColloidVelocities[j].z,CountBox_x[j],CountBox_y[j],CountBox_z[j]);
//fprintf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d\n",ColloidPositions[j].x,ColloidPositions[j].y,ColloidPositions[j].z,ColloidVelocities[j].x,ColloidVelocities[j].y,ColloidVelocities[j].z,Omega[j].x,Omega[j].y,Omega[j].z,PatchPositions[j].x,PatchPositions[j+NumberOfColloids].x,PatchPositions[j].y,PatchPositions[j+NumberOfColloids].y,PatchPositions[j].z,PatchPositions[j+NumberOfColloids].z,CountBox_x[j],CountBox_y[j],CountBox_z[j]);
   }
 fclose(file);

 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_msd_run%d.tmp",NumberOfColloids,run,t);
 file=fopen(filename,"wb");
 fwrite(R2,sizeof(R2[0]),MAXT,file);//note that cells are more than flow field grids
 fclose(file);

 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_vacf_run%d.tmp",NumberOfColloids,run,t);
 file=fopen(filename,"wb");
 fwrite(Vacf,sizeof(Vacf[0]),MAXT,file);//note that cells are more than flow field grids
 fclose(file);

/*
 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_oacf_run%d.tmp",NumberOfColloids,run,t);
 file=fopen(filename,"wb");
 fwrite(U2,sizeof(U2[0]),MAXT,file);//note that cells are more than flow field grids
 fclose(file);
 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_wacf_run%d.tmp",NumberOfColloids,run,t);
 file=fopen(filename,"wb");
 fwrite(Wacf,sizeof(Wacf[0]),MAXT,file);//note that cells are more than flow field grids
 fclose(file);
*/

 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_samplecounter_run%d.tmp",NumberOfColloids,run,t);
 file=fopen(filename,"wb");
 fwrite(SampleCounter,sizeof(SampleCounter[0]),MAXT,file);//note that cells are more than flow field grids
 fclose(file);

////writing the restart files of the SQT
/*
 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_sqt_run%d.tmp",NumberOfColloids,run,t);
 file=fopen(filename,"wb");
 fwrite(S_qt,sizeof(S_qt[0][0]),MAX_q*MAXT,file);//note that cells are more than flow field grids
 fclose(file);
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 sprintf(command,"mv /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_coll_run%d.tmp /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_coll_run%d.dat",NumberOfColloids,run,t,NumberOfColloids,run,t);
 system(command);

 sprintf(command,"mv /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_msd_run%d.tmp /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_msd_run%d.dat",NumberOfColloids,run,t,NumberOfColloids,run,t);
 system(command);

 sprintf(command,"mv /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_vacf_run%d.tmp /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_vacf_run%d.dat",NumberOfColloids,run,t,NumberOfColloids,run,t);
 system(command);
/*
 sprintf(command,"mv /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_oacf_run%d.tmp /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_oacf_run%d.dat",NumberOfColloids,run,t,NumberOfColloids,run,t);
 system(command);

 sprintf(command,"mv /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_wacf_run%d.tmp /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_wacf_run%d.dat",NumberOfColloids,run,t,NumberOfColloids,run,t);
 system(command);
*/
 sprintf(command,"mv /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_samplecounter_run%d.tmp /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_samplecounter_run%d.dat",NumberOfColloids,run,t,NumberOfColloids,run,t);
 system(command);

//writing the restart files of the SQT
/*
 sprintf(command,"mv /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_sqt_run%d.tmp /p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_sqt_run%d.dat",NumberOfColloids,run,t,NumberOfColloids,run,t);
 system(command);
*/
//////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 sprintf(filename,"param/run%d.dat",run);
 file=fopen(filename,"w");
 fprintf(file,"%lf\n",BoxSize_x);
 fprintf(file,"%lf\n",BoxSize_y);
 fprintf(file,"%lf\n",BoxSize_z);
 fprintf(file,"%lf\n",Temperature);
 fprintf(file,"%lf\n",MDt);
 fprintf(file,"%d\n",NumberOfSteps);
 fprintf(file,"%d\n",output_restart_steps);
 fprintf(file,"%d\n",t);//restart time step
 fclose(file);
}


void initial_pos_vel_restart(int run, int t)
{int i,j;
double invOsmod;

 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_coll_run%d.dat",NumberOfColloids,run,t);
 file=fopen(filename,"r");
 for (j=0;j<NumberOfColloids;j++)
    {
fscanf(file,"%lf %lf %lf %lf %lf %lf %d %d %d\n",&ColloidPositions[j].x,&ColloidPositions[j].y,&ColloidPositions[j].z,&ColloidVelocities[j].x,&ColloidVelocities[j].y,&ColloidVelocities[j].z,&CountBox_x[j],&CountBox_y[j],&CountBox_z[j]);

//fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d\n",&ColloidPositions[j].x,&ColloidPositions[j].y,&ColloidPositions[j].z,&ColloidVelocities[j].x,&ColloidVelocities[j].y,&ColloidVelocities[j].z,&Omega[j].x,&Omega[j].y,&Omega[j].z,&PatchPositions[j].x,&PatchPositions[j+NumberOfColloids].x,&PatchPositions[j].y,&PatchPositions[j+NumberOfColloids].y,&PatchPositions[j].z,&PatchPositions[j+NumberOfColloids].z,&CountBox_x[j],&CountBox_y[j],&CountBox_z[j]);
//renormalize the unit vector to avoid imaginary part comes from numerical errors
/*
Os[j].x=(PatchPositions[NumberOfColloids+j].x-PatchPositions[j].x);
Os[j].y=(PatchPositions[NumberOfColloids+j].y-PatchPositions[j].y);
Os[j].z=(PatchPositions[NumberOfColloids+j].z-PatchPositions[j].z);
invOsmod=1.0/(sqrt(SQR(Os[j].x)+SQR(Os[j].y)+SQR(Os[j].z)));
Os[j].x *=invOsmod;
Os[j].y *=invOsmod;
Os[j].z *=invOsmod;

VCross(OsDot[j],Omega[j],Os[j]);

Arms[j]=-Rdipole;
Arms[j+NumberOfColloids]=Rdipole;
*/ 
   }
 fclose(file);

 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_msd_run%d.dat",NumberOfColloids,run,t);
 file=fopen(filename,"rb");
 fread(R2,sizeof(R2[0]),MAXT,file);
 fclose(file);

 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_vacf_run%d.dat",NumberOfColloids,run,t);
 file=fopen(filename,"rb");
 fread(Vacf,sizeof(Vacf[0]),MAXT,file);
 fclose(file);
/*
 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_oacf_run%d.dat",NumberOfColloids,run,t);
 file=fopen(filename,"rb");
 fread(U2,sizeof(U2[0]),MAXT,file);
 fclose(file);

 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_wacf_run%d.dat",NumberOfColloids,run,t);
 file=fopen(filename,"rb");
 fread(Wacf,sizeof(Wacf[0]),MAXT,file);
 fclose(file);
*/

 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_samplecounter_run%d.dat",NumberOfColloids,run,t);
 file=fopen(filename,"rb");
 fread(SampleCounter,sizeof(SampleCounter[0]),MAXT,file);
 fclose(file);

//read the Sqt binary file
/*
 sprintf(filename,"/p/scratch/jics30/BDfiles/Particlediffusion%d/salr_restart_files/data_run%d/restart_sqt_run%d.dat",NumberOfColloids,run,t);
 file=fopen(filename,"rb");
 fread(S_qt,sizeof(S_qt[0][0]),MAX_q*MAXT,file);
 fclose(file);
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
}

