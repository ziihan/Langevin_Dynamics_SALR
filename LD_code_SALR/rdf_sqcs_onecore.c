#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

//#define Maxx 200
//#define MAX_q 3600 //wave vector two dimensions
//#define q_max 60
//#define Maxo 1
// samples the radial distribution function
void SampleRDF(int Ich)
{
  int i,j,k;
  double r2,o2,f2;//calculate the relative distance and relative orientation
 
  static double DeltaR,Gg[Maxx];//DeltaOs,DeltaOf;
  static double S_qc[MAX_q],S_qs[MAX_q];

  static int Bla;//sample counter
  VECTOR dr,dw;
  unsigned int qx,qy,qz;
  double sqc,sqs;
//  static double Unz,Sn;
  //Sn, measuring the orientational order parameter

  char filename[200];
  FILE *File;

  char filenamegr[200];
  FILE *Filegr;

  char filename0[200];
  static FILE *FileNeO; 
  switch(Ich)
  {
    case INITIALIZE:
      for(i=0;i<(Maxx);i++)
        Gg[i]=0.0;

    for (k=0;k<MAX_q;k++)
      {
        S_qc[k]=0.0;// + 0.0*I;          
        S_qs[k]=0.0; 
      //rho_qt[k][i]=0.0 + 0.0*I;
      }
                        
     //   Unz=0.0;
      Bla=0;//sample counter initialization
      DeltaR=(10.0*Rcol)/((double)Maxx);
      //DeltaOs=DeltaOf=PI/((double)Maxx);

//      sprintf(filename0,"Particlediffusion%d/Protein-NematicOrder-%d.dat",NumberOfColloids,Counter);
//      FileNeO=fopen(filename0,"wt");

      break;
    case SAMPLE:
      Bla++;
// printf("Bla=%d\n",Bla);     
 // loop over all pairs


         for (k=0;k<MAX_q;k++)
{
       qx = k % q_max;
       qy= (k-qx)/q_max % q_max;

   sqc=0.0;
   sqs=0.0;
      for(i=0;i<NumberOfColloids;i++)
      {
       sqc +=cos(dq*((qx+1.0)*ColloidPositions[i].x+(qy+1.0)*ColloidPositions[i].y));
       sqs +=sin(dq*((qx+1.0)*ColloidPositions[i].x+(qy+1.0)*ColloidPositions[i].y));
      }
   S_qc[k] +=SQR(sqc);
   S_qs[k] +=SQR(sqs);
  //printf("Sqc=%f\tSqs=%f\n",S_qc[k],S_qs[k]);
}
 // loop over all pairs
      for(i=0;i<NumberOfColloids-1;i++)
      {
        for(j=i+1;j<NumberOfColloids;j++)
        {
          dr.x=ColloidPositions[i].x-ColloidPositions[j].x;
          dr.y=ColloidPositions[i].y-ColloidPositions[j].y;
      //    dr.z=ColloidPositions[i].z-ColloidPositions[j].z;

          // apply boundary conditions
      if (dr.x>0.5*BoxSize_x)dr.x-=BoxSize_x;
      else if(dr.x<-0.5*BoxSize_x)dr.x+=BoxSize_x;
      if (dr.y>0.5*BoxSize_y)dr.y-=BoxSize_y;
      else if(dr.y<-0.5*BoxSize_y)dr.y+=BoxSize_y;

          r2=sqrt(SQR(dr.x)+SQR(dr.y));//+SQR(dr.z));

          // calculate in which bin this interaction is in
//          if(r2<0.5*MIN(BoxSize_x,BoxSize_y))
           if(r2<8.0*Rcol)   
         Gg[(int)(r2/DeltaR)]+=2.0;

        }
      }
/*
     for(i=0;i<NumberOfColloids;i++)
        {
//        Un.x +=Os[i].x;
//        Un.y +=Os[i].y;
        Unz +=SQR(Os[i].z);
        }
        o2=(double)(NumberOfColloids*Bla);
        Sn = 0.5*(3.0*(Unz/o2)-1.0);//       
//printf("o2=%f\tSn=%f\n",o2,Sn);//
      fprintf(FileNeO,"%f %f\n",((double)Bla)*SampleSteps*MDt,Sn);
*/

      break;
    case WRITE_RESULTS:
      // Write Results To Disk
      sprintf(filename,"Particlediffusion%d/Protein-Sq-%d-Step%d.dat",NumberOfColloids,Counter,Step);
      File=fopen(filename,"wt");
        for (k=0;k<MAX_q;k++)
             {
         // S_qt[k][i] /=(double)(NumberOfColloids*SampleCounter[i]);
            qx = k % q_max;
            qy= (k-qx)/q_max % q_max;
              fprintf(File,"%lf %lf\n",Rcc*dq*sqrt(SQR(qx+1.0)+SQR(qy+1.0)),(S_qc[k]+S_qs[k])/((double)(NumberOfColloids*Bla)));
//            fprintf(FilePtrSqt,"%lf %lf %lf %lf\n",dq*(k+1.0),SampleSteps*i*MDt,creal(S_qt[k][i]),cimag(S_qt[k][i]));
             }

      sprintf(filenamegr,"Particlediffusion%d/Protein-rodf-%d-Step%d.dat",NumberOfColloids,Counter,Step);
      Filegr=fopen(filenamegr,"wt");
      for(i=0;i<Maxx-1;i++)
      {
        r2=PI*(NumberOfColloids/(BoxSize_x*BoxSize_y))*SQR(DeltaR)*(SQR(i+1)-SQR(i));
        fprintf(Filegr,"%f %f\n",(i+0.5)*DeltaR/Rcc,Gg[i]/((double)(Bla*NumberOfColloids)*r2));
      }
      fclose(File);
      fclose(Filegr);
//      fclose(FileNeO);
      break;
  }
}
