#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
//#include <complex.h>

//#define MAX_q 10
//
//Computing mean-squared displacement, velocity autocorrelation functions. cf. the algorithm from the book "Unerstanding Molecular Simuations: from alogrithms to applications" by Daan Frenkel and Berend Smit.
//
void SampleDiff(int Switch)
{
  int t0index,i,j,jj,k,CorrelTime;
  double CumIntegrationV,CumIntegrationW;
  static int time,t0Counter,t0time[MAXT0];
  double Rsqr,Usqr,VCF,WCF;
  double dcc;
//temperary varables
  static double Vxt0[MAXCOLL][MAXT0],Vyt0[MAXCOLL][MAXT0],Vzt0[MAXCOLL][MAXT0];
  static double Rx0[MAXCOLL][MAXT0],Ry0[MAXCOLL][MAXT0],Rz0[MAXCOLL][MAXT0];
  static double Ux0[MAXCOLL][MAXT0],Uy0[MAXCOLL][MAXT0],Uz0[MAXCOLL][MAXT0];
  static double Wxt0[MAXCOLL][MAXT0],Wyt0[MAXCOLL][MAXT0],Wzt0[MAXCOLL][MAXT0];
      char filename1[200];
      char filename2[200];
      char filename3[200];
      char filename4[200];

  FILE *FilePtrMsd;
  FILE *FilePtrVacf;
  FILE *FilePtrUacf;//orientation auto-correlation function
  FILE *FilePtrWacf;//angular velocity auto-correlation function

  switch(Switch)
  {
    // initialize everything
    case INITIALIZE:
      time=0;
      t0Counter=0;
      
      //dq=2.0*PI/MIN(BoxSize_x,BoxSize_y);
    if (RunStep<(output_restart_steps)) 
    {  for(i=0;i<MAXT;i++)
      {
        R2[i]=0.0;

//	U2[i]=0.0;
        SampleCounter[i]=0;
        Vacf[i]=0.0;
//	Wacf[i]=0.0;
      }
    }
      break;
      case SAMPLE:
      time++;
      if((time%FREQT0)==0)//decide to take a new t = 0
      {
        t0Counter++;
        t0index=(t0Counter-1)%MAXT0;
        t0time[t0index]=time;
//because rods are rigid, two positions at ends determine the center of mass and orientation vector 
       for(j=0;j<NumberOfColloids;j++){
       Rx0[j][t0index]=ColloidPositions[j].x+BoxSize_x*CountBox_x[j];
       Ry0[j][t0index]=ColloidPositions[j].y+BoxSize_y*CountBox_y[j];
       Rz0[j][t0index]=ColloidPositions[j].z+BoxSize_z*CountBox_z[j];
       Vxt0[j][t0index]=ColloidVelocities[j].x;
       Vyt0[j][t0index]=ColloidVelocities[j].y;
       Vzt0[j][t0index]=ColloidVelocities[j].z;

/*
       Ux0[j][t0index]=Os[j].x;
       Uy0[j][t0index]=Os[j].y;
       Uz0[j][t0index]=Os[j].z;
       
       Wxt0[j][t0index]=Omega[j].x;
       Wyt0[j][t0index]=Omega[j].y;
       Wzt0[j][t0index]=Omega[j].z;
*/
       }
      }
      // loop over all time origins that have been stored
      // MAXT0 is the number of time origins

  for(t0index=0;t0index<MIN(t0Counter,MAXT0);t0index++){
         CorrelTime=time-t0time[t0index];//?????? the index of sample in time_origin----tmax
    if(CorrelTime<MAXT){
          SampleCounter[CorrelTime]++;

       for(j=0;j<NumberOfColloids;j++){

Vacf[CorrelTime]+=(ColloidVelocities[j].x*Vxt0[j][t0index]+ColloidVelocities[j].y*Vyt0[j][t0index]+ColloidVelocities[j].z*Vzt0[j][t0index]);
//     Wacf[CorrelTime]+=(Omega[j].x*Wxt0[j][t0index]+Omega[j].y*Wyt0[j][t0index]+Omega[j].z*Wzt0[j][t0index]);
 
       ColloidNewPositions[j].x = ColloidPositions[j].x+BoxSize_x*CountBox_x[j];
       ColloidNewPositions[j].y = ColloidPositions[j].y+BoxSize_y*CountBox_y[j];
       ColloidNewPositions[j].z = ColloidPositions[j].z+BoxSize_z*CountBox_z[j]; 

//       U2[CorrelTime] +=Os[j].x*Ux0[j][t0index]+Os[j].y*Uy0[j][t0index]+Os[j].z*Uz0[j][t0index];

       R2[CorrelTime] +=(SQR(ColloidNewPositions[j].x-Rx0[j][t0index])+SQR(ColloidNewPositions[j].y-Ry0[j][t0index])+SQR(ColloidNewPositions[j].z-Rz0[j][t0index]));
                        }
          }
      }
      // end modification
    break;
    case WRITE_RESULTS:
      CumIntegrationV=0.0;
 //     CumIntegrationW=0.0;
sprintf(filename2,"Particlediffusion%d/Protein-vacf-%d-Step%d.dat",NumberOfColloids,Counter,Step);
//sprintf(filename4,"Particlediffusion%d/Protein-wacf-%d-Step%d.dat",NumberOfColloids,Counter,Step);
//sprintf(filename3,"Particlediffusion%d/Protein-oacf-%d-Step%d.dat",NumberOfColloids,Counter,Step);
sprintf(filename1,"Particlediffusion%d/Protein-msd-%d-Step%d.dat",NumberOfColloids,Counter,Step);

      FilePtrVacf=fopen(filename2,"wt");
  //    FilePtrWacf=fopen(filename4,"wt");
 //     FilePtrUacf=fopen(filename3,"wt");
      FilePtrMsd=fopen(filename1,"wt");
      for(i=0;i<(MAXT-0);i++)
      {
        if(SampleCounter[i]>0)
        {
          VCF=Vacf[i]/(double)(NumberOfColloids*SampleCounter[i]);
  //        WCF=Wacf[i]/(double)(NumberOfColloids*SampleCounter[i]);
   //       Usqr=U2[i]/(double)(NumberOfColloids*SampleCounter[i]);
          Rsqr=R2[i]/(double)(NumberOfColloids*SampleCounter[i]);
        }
        else
        {
          VCF=0.0;
     //     WCF=0.0;
   //       Usqr=0.0;
          Rsqr=0.0; 
        }

        CumIntegrationV +=SampleSteps*MDt*VCF/3.0;//the prefactor indicates how many MD time steps per sampling performence
        fprintf(FilePtrVacf,"%lf %lf %lf\n",SampleSteps*(i)*MDt,VCF,CumIntegrationV);
       // CumIntegrationW +=SampleSteps*MDt*WCF/3.0;
      //  fprintf(FilePtrWacf,"%lf %lf %lf\n",SampleSteps*(i)*MDt,WCF,CumIntegrationW);
      //  fprintf(FilePtrUacf,"%lf %lf %lf\n",SampleSteps*(i)*MDt,Usqr,(-0.5*log(Usqr)/(SampleSteps*(i+1)*MDt)));
        if(i>=1){
            fprintf(FilePtrMsd,"%lf %lf %lf\n",SampleSteps*(i)*MDt,Rsqr,Rsqr/(6.0*(i+1)*MDt*SampleSteps));
          }
        else{
            fprintf(FilePtrMsd,"%lf %lf %lf\n",SampleSteps*(i)*MDt,Rsqr,0.0);
          }
      }
      fclose(FilePtrVacf);
    //  fclose(FilePtrWacf);
    //  fclose(FilePtrUacf);
      fclose(FilePtrMsd);
  }
}

