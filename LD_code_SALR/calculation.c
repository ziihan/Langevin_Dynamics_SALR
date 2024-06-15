/*This script calculates the local density, local orientation of concentrated rods from vmd xyz trajectory files */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#define MAXFRAME 1000000     //preset how many frames in the xyz files
#define MAXCOLL 500          //how many atoms or beads in each frame
#define MAXSLAB 100
#define SQR(x) ((x)*(x))    //define the square
#define CUBE(x) ((x)*(x)*(x)) //define cubic

//#define MIN(x,y) ((x)<(y)?(x):(y))
//#define MAX(x,y) ((x)>(y)?(x):(y))

typedef struct  //define x.y.z components
{
  double x, y, z;
} VECTOR;

VECTOR ColloidPositions[MAXFRAME][MAXCOLL];
//double LocalTimeDensity_Z[MAXFRAME][MAXSLAB];
//VECTOR LocalTimeAlignment[MAXFRAME][MAXSLAB];

extern VECTOR ColloidPositions[MAXFRAME][MAXCOLL];

//double LocalTimeDensity_Z[MAXFRAME][MAXSLAB];
VECTOR LocalTimeAlignment[MAXFRAME][MAXSLAB];
//
//extern double LocalTimeDensity_Z[MAXFRAME][MAXSLAB];
extern VECTOR LocalTimeAlignment[MAXFRAME][MAXSLAB];


//calculate how many lines of the file
int numLines(char *fileName) {
    FILE *f;
    char c;
    int lines = 0;

    f = fopen(fileName, "r");

    if(f == NULL)
        return 0;

    while((c = fgetc(f)) != EOF)
        if(c == '\n')
            lines++;

    fclose(f);
/*
    if(c != '\n')
        lines++;
*/
    return lines;
}




int main()
{ 
FILE *file;
char filename[400];//,dirname[400],command[400];
unsigned int NumberOfColloids;   //number of atoms
unsigned int NumberOfRods;  //number of rods
unsigned int i,j;
unsigned int filelines,Nframe;
unsigned int BoxSize_z=150;
unsigned int slabZ=5;
unsigned int NUMslab;


NumberOfColloids=300;////specify how many atoms 
NumberOfRods=150;
NUMslab=BoxSize_z/slabZ;


sprintf(filename,"/usr/users/iff_th2/ztan/RETOOL/MyMPC/DimerColl/fixednew.dat");
file=fopen(filename,"r");
printf("what!\n");

filelines=numLines(filename);
printf("lines of file=%d\n",filelines);
Nframe=filelines/NumberOfColloids;
printf("lines of Frame=%d\n",Nframe);


/*Here we initialize the local density and alignment at different frame*/
//VECTOR ColloidPositions[Nframe][NumberOfColloids];

//double LocalTimeDensity_Z[Nframe][NUMslab];
//VECTOR LocalTimeAlignment[Nframe][NUMslab];

for (i=0;i<Nframe;i++){
for (j=0;j<NumberOfColloids;j++)
{
fscanf(file,"%lf\t%lf\t%lf\n",&ColloidPositions[i][j].x,&ColloidPositions[i][j].y,&ColloidPositions[i][j].z);
}
}




printf("caculation finished!\n");
fclose(file);
}
