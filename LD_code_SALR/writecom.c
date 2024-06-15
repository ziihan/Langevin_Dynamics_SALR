#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"
#include <omp.h>
// make a movie file of the simulation box
// use vmd to view it..
void WritePdb(FILE *FilePtr)
{
  int j;
  static int Countstep=0,CountCol=0;

fprintf(FilePtr,"%9d\n",Countstep*1000);
    Countstep++;
//  fprintf(FilePtr,"%9d\n",3*NumberOfColloids);
//  fprintf(FilePtr,"#\tX\tY\tZ\n");
#pragma omp parallel private(j)
{
#pragma omp for
  for(j=0;j<NumberOfColloids;j++)
  {
    //CountCol++;
    fprintf(FilePtr,"%8.3lf\t%8.3lf\t%8.3lf\n",ColloidPositions[j].x*1.0,ColloidPositions[j].y*1.0,ColloidPositions[j].z*1.0);

//    fprintf(FilePtr,"%s\t%8.3lf\t%8.3lf\t%8.3lf\n%s\t%8.3lf\t%8.3lf\t%8.3lf\n%s\t%8.3lf\t%8.3lf\t%8.3lf\n",
//      "C",ColloidPositions[j].x*1.0,ColloidPositions[j].y*1.0,ColloidPositions[j].z*1.0,"O",PatchPositions[j].x*1.0,PatchPositions[j].y*1.0,PatchPositions[j].z*1.0,
//      "N",PatchPositions[j+NumberOfColloids].x*1.0,PatchPositions[j+NumberOfColloids].y*1.0,PatchPositions[j+NumberOfColloids].z*1.0);
/*
    fprintf(FilePtr,"%s\t%8.3lf\t%8.3lf\t%8.3lf\n",
      "N",PatchPositions[j].x*1.0,PatchPositions[j].y*1.0,PatchPositions[j].z*1.0);

    fprintf(FilePtr,"%s\t%8.3lf\t%8.3lf\t%8.3lf\n",
      "O",PatchPositions[j+NumberOfColloids].x*1.0,PatchPositions[j+NumberOfColloids].y*1.0,PatchPositions[j+NumberOfColloids].z*1.0);
*/
//firintf(FilePtr,"%s\n","ENDMDL");
  }
}
}
