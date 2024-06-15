#include <stdio.h>
#include <stdlib.h>
#include <math.h>





//#define MAXLINE 1000000000
#define MAXFRAME 200
#define MAXCOLL 500
#define SQR(x) ((x)*(x))    //define the square
#define CUBE(x) ((x)*(x)*(x)) //define cubic
#define MIN(x,y) ((x)<(y)?(x):(y)) //define 
#define MAX(x,y) ((x)>(y)?(x):(y))


typedef struct  //define x.y.z components
{
  double x, y, z;
} VECTOR;

VECTOR ColloidPositions[MAXFRAME][MAXCOLL];

extern VECTOR ColloidPositions[MAXFRAME][MAXCOLL];



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

    if(c != '\n')
        lines++;

    return lines;
}




int main()
{ 
FILE *file;
char filename[400];//,dirname[400],command[400];
int NumberOfColloids;
int NumberPhiColloids;
int i,j;
int filelines,Nframe;

NumberOfColloids=300;
NumberPhiColloids=150;

sprintf(filename,"/usr/users/iff_th2/ztan/RETOOL/MyMPC/DimerColl/fixed.dat");
filelines=numLines(filename)-1;
printf("lines of file=%d\n",filelines);
Nframe=filelines/NumberOfColloids;
printf("lines of frame=%d\n",Nframe);
file=fopen(filename,"r");


for (i=0;i<Nframe;i++){
for (j=0;j<NumberOfColloids;j++)
{
fscanf(file,"%lf\t%lf\t%lf\n",&ColloidPositions[i][j].x,&ColloidPositions[i][j].y,&ColloidPositions[i][j].z);
}

//if(i==299){
printf("Px=%lf\tPy=%lf\tPz=%lf\n",ColloidPositions[i][299].x,ColloidPositions[i][299].y,ColloidPositions[i][299].z);

//}
}


fclose(file);
}
