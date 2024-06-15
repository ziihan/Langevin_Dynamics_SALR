#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include "system.h"

// read parameters: simulation box sizes, temperature, time steps, check point for restart files, etc..
void Readdat(void)
{ 
  FILE *FilePtr; 
  char filename[200];
  sprintf(filename,"param/run%d.dat",Counter);
  FilePtr=fopen(filename,"r");

 fscanf(FilePtr,"%lf\n",&BoxSize_x);
 fscanf(FilePtr,"%lf\n",&BoxSize_y);
 fscanf(FilePtr,"%lf\n",&BoxSize_z);
 fscanf(FilePtr,"%lf\n",&Temperature);
 //fscanf(FilePtr,"%lf\n",&Deltat);
 fscanf(FilePtr,"%lf\n",&MDt);
 fscanf(FilePtr,"%d\n",&NumberOfSteps);
 fscanf(FilePtr,"%d\n",&output_restart_steps);//
 fscanf(FilePtr,"%d\n",&RunStep);//running step-for restart file-can change it by hands
 fclose(FilePtr);

  // print information to the screen
 
  printf("Langevein Dynamics Program:\n");
  printf("\n");
//  printf("Number of Cells       : %d\n",NumberOfCells_x*NumberOfCells_y*NumberOfCells_z);
  printf("Boxlength             : %f\n",BoxSize_x);
  printf("Temperature           : %f\n",Temperature);
  printf("Density                : %f\n", Rho);
  printf("Number of steps       : %d\n",NumberOfSteps);
  printf("Timestep              : %f\n",MDt);
}
