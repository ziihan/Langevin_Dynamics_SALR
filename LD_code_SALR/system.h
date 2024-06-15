//All of the global variables
#include <stdlib.h>
#include <stdio.h>

#define PI 3.1415926535898
#define MAXCOLL 5000       //maximum number of colloids
#define VMax 10            ////maximum velocity of solvent particles(used in MPC code)

//Macros for moving window averaging
#define Maxx 200
#define MAXT 20000         //maximum correlation time
#define MAXT0 2000         //maximum number of correlation time origin
#define FREQT0 10          //the frequency for setting the time origin
#define MAX_q 3600         //number of points in q-space
#define q_max 60           //maximum number of wavenumbers

//FOR random number generator
#define SHR0 (JSR^=(JSR<<21),JSR^=(JSR>>35),JSR^=(JSR<<4)) // xorshift RNG
#define SHR3 (JZ=JSR,JSR^=(JSR<<21),JSR^=(JSR>>35),JSR^=(JSR<<4),JZ+JSR) // xorshift RNG
#define RandomNumber (0.5+(signed) SHR0*0.2328306e-09) // uniform RNG
#define RandomGaussianNumber (HZ=SHR0,IZ=(HZ&127),(abs(HZ)<KN[IZ])?(float) (HZ)*WN[IZ]:fix()) // ziggurat normal RNG
#define MAXTHREAD 50
//#define NUM_THREADS 8
extern unsigned int NUM_THREADS;   //number of threads
extern unsigned long int JZ,JSR;   //these three lines are related to random number generator
extern int HZ,IZ,KN[128];          //
extern float WN[128],FN[128];      //

extern unsigned long int seed[MAXTHREAD];    //seed for (parallelized) random number generator
extern unsigned int N_threads;               //the index of thread?

#define SQR(x) ((x)*(x))    //define the square
#define CUBE(x) ((x)*(x)*(x)) //define cubic
#define MIN(x,y) ((x)<(y)?(x):(y)) //define 
#define MAX(x,y) ((x)>(y)?(x):(y))
#define BINAry(x,y) ((x)<(y)?(0):(1))
enum{INITIALIZE,SAMPLE,WRITE_RESULTS};         //for the measurement of g(r),S(q), mean-squared displacement etc..

typedef struct  //define x.y.z components
{
  double x, y, z;
} VECTOR;

//set certain values to vector
#define VSet(v,sx,sy,sz)	\
(v).x = sx,	\
(v).y = sy,     \
(v).z = sz    
//set same values to vector
#define VSetAll(v,s)	\
VSet(v, s, s, s)

#define VScale(v,s)	\
(v).x *= s,	\
(v).y *= s,     \
(v).z *= s

#define VAdd(v1,v2,v3)	\
(v1).x = (v2).x + (v3).x,	\
(v1).y = (v2).y + (v3).y,        \
(v1).z = (v2).z + (v3).z

#define VSub(v1,v2,v3)	\
(v1).x = (v2).x - (v3).x,        \
(v1).y = (v2).y - (v3).y,        \
(v1).z = (v2).z - (v3).z

#define VDot(v1 , v2)	\
((v1).x*(v2).x + (v1).y*(v2).y + (v1).z*(v2).z)

#define VCross(v1, v2, v3)	\
(v1).x = (v2).y * (v3).z - (v2).z * (v3).y,	\
(v1).y = (v2).z * (v3).x - (v2).x * (v3).z,     \
(v1).z = (v2).x * (v3).y - (v2).y * (v3).x

#define VZero(v)	VSetAll(v, 0)	//set vector components to zero
#define VLenSq(v)	VDot(v,v)	//the square of length of vector
#define VLen(v)		sqrt(VDot(v,v))	//the length of vector
#define VVAdd(v1,v2)	VAdd(v1, v1, v2)	//update compoents of vector
#define VVSub(v1,v2)    VSub(v1, v1, v2)	//update compoents of vector


extern double Temperature;
extern int SampleSteps;             //how often the sampling takes place for correlation functions

extern double xif; //for friction, the other 8 lines are 
extern double c0;
extern double c1;
extern double c2;
extern double c3;
extern double c2c3;

extern double Cx;
extern double Gx;
extern double Ex;

extern double dq;
extern unsigned int COUNT;
extern unsigned int PrintStep;
extern unsigned int Step;
extern unsigned int MPCstep;
extern unsigned int RecordStep;

extern int SampleCounter[MAXT];
extern double Vacf[MAXT];
extern double Wacf[MAXT];
extern double R2[MAXT];
extern double U2[MAXT];
extern double S_qt[MAX_q][MAXT];  //turn it on when it is measured
//extern double dq_sigma;

extern unsigned int PrintStep;
extern unsigned int Step;
extern unsigned int RecordStep;

//extern int mark[MAXCOLL];//if Philic Colloid mark 0, if Phobic Colloid mark 1
//variables for linked-list sorting of the collision cells
extern int Counter;
extern int Needle;
extern double beta;

extern int Rcelli;
extern int Rcello;
//linked-list for MD calculation
extern double OutOffi;//the Radius of the outer shell for Verlet-list or Cell-list approaches
extern double OutOffo;
extern double OutOffcc;
extern double CellSize;
//cell-index of colloids
extern int Nix[MAXCOLL];
extern int Niy[MAXCOLL];
extern int Niz[MAXCOLL];

extern int iNeighbour;//number of neighbour without overlapping
extern int IntegrateNeighbour[MAXCOLL*100000];
extern int Neighbour[MAXCOLL*100000];//big enough Rho*(4/3)*PI*(Rcut^3-Rcol^3)
extern int HeatNeighbour[MAXCOLL*100000];
extern int StartOfNeighbour[MAXCOLL];
extern int NumberOfNeighbour[MAXCOLL];

//extern int *ll;
//extern int *hoc;
 
extern int icell;
extern int FlowCount;//counter for storing flowfield
extern int TempCount;//counter for storing temperature profile
extern int DensCount;

//Parameters of Colloids at the MD part 
extern double NaT;

extern double Wall0;
extern double Wallz;
extern double UWalls;
extern double USpring;
extern double UBending;

extern VECTOR CNP;
extern VECTOR COP;

extern double cCon;
extern double ColloidBending;
extern double ColloidSpring;
extern double BondLength;
extern double RodLength;


extern VECTOR ColloidNewPositions[MAXCOLL];
extern VECTOR ColloidOldPositions[MAXCOLL];
extern VECTOR Colloidt0Positions[MAXCOLL];

extern VECTOR ColloidPositions[MAXCOLL];
extern VECTOR ColloidVelocities[MAXCOLL];
extern VECTOR ColloidForces[MAXCOLL];

extern VECTOR ColloidPositionsPriv[MAXCOLL][MAXTHREAD];
extern VECTOR ColloidVelocitiesPriv[MAXCOLL][MAXTHREAD];
extern VECTOR ColloidForcesPriv[MAXCOLL][MAXTHREAD];

extern double Mcol;//colloid mass
extern double invMcol;//colloid inversed mass
extern double CutOffi;        // Cut-Off Radius
extern double Eshif;         // Shifted energy
extern double HeatOn;
extern double CutOffo;
extern double Ecut;
extern double EcutYK;
extern double FcutYK;
extern double CutOffcc;
extern double CutOffYK;
//extern double CutOffwy;
//extern double CutOffwz;
extern double CutOffw;
extern double DebyeLength;
extern double qq;
extern double epsilon;

//extern VECTOR TempPosition;//
extern double Rcol;          //colloid radius
extern double Rcc;//colloid-colloid interaction radius
extern double Inertia;
extern double invInertia;



/*
extern VECTOR COM;
extern VECTOR OmegaDt0;
extern VECTOR Omega[MAXCOLL];
extern VECTOR OmegaDt[MAXCOLL];
//extern VECTOR TorqueAll;
extern double Arms[2*MAXCOLL];//lever arms
//extern VECTOR Torque[MAXCOLL];
extern VECTOR Os[MAXCOLL];//unit vector along bond
extern VECTOR OsDot[MAXCOLL];//d(Os)/dt
extern VECTOR Os2Dot[MAXCOLL];//d2(Os)/dt2
extern VECTOR SignedForce[MAXCOLL];
extern VECTOR OmegaxOs[MAXCOLL];
extern VECTOR Os2Dot1[MAXCOLL];
extern VECTOR Os2Dot2[MAXCOLL];
*/

//variables for Mean-Square-Displacement
extern int CountBox_x[MAXCOLL];
extern int CountBox_y[MAXCOLL];
extern int CountBox_z[MAXCOLL];

//extern double force;

extern double Rho;                //density

extern double MDt;            // MD Timestep

extern double BoxSize_x;           // Boxlengths
extern double BoxSize_y;
extern double BoxSize_z;
//some MD things

//add this to test the conservation of energy
extern double UKinetic;      // Kinetic Energy
extern double UPotential;    // Potential Energy
extern double UTotal;        // Total Energy

extern unsigned int output_restart_steps;
extern unsigned int RunStep;
extern unsigned int NumberOfSteps;   //number of MD steps for MD-MPC coupling
extern unsigned int NumberOfCells;   //number of collision cells
extern unsigned int NumberOfParticles;//
extern unsigned int NumberHeatColloids;
extern unsigned int NumberOfColloids;
extern unsigned int Nx;
extern unsigned int Ny;
extern unsigned int Nz;
extern unsigned int Nobs;
extern unsigned int NumberChainColloids;
//number of collision cells along x-y-z directions

extern int NumberOfBins;
extern double BinSize;
//some subroutines
void WritePdb(FILE *FilePtr);
void Readdat(void);
void Init(void);
void Stream(void);
void Collision(void);
void BDloop(void);
void MSD(void);
void Flow(void);
void GridCells(int);
void StoreFlow(void);
void SlabTemp(void);
void TempK(void);
void PrintTemp(void);
void Exchange(int EXPart);
void Force(void);
//void Integrate(void);
void Integrate(int Switch);
void Neighbours(void);
//void CountNeighbours(void);
void LocalHeat(void);
void MeasureHeat(void);
void Thermostats(void);
//void ReSequence(void);
//void ReScale(void);

void SampleRDF(int Switch);
void SampleDiff(int Switch);
void output_pos_vel_restart(int run, int t);
void initial_pos_vel_restart(int run, int t);
