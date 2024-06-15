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

#define VZero(v)	VSetAll(v, 0)			//set vector components to zero
#define VLenSq(v)	VDot(v,v)			//the square of length of vector
#define VLen(v)		sqrt(VDot(v,v))			//the length of vector
#define VVAdd(v1,v2)	VAdd(v1, v1, v2)		//update compoents of vector
#define VVSub(v1,v2)    VSub(v1, v1, v2)		//update compoents of vector


extern double Temperature;
extern int SampleSteps;             //how often the sampling takes place for correlation functions

extern double xif; //for friction, the other 8 lines are all for Langevin Dynamics 
extern double c0;  //
extern double c1;  //
extern double c2;  //
extern double c3;  //
extern double c2c3;//

extern double Cx;  //
extern double Gx;  //
extern double Ex;  //

extern double dq;  //step size in q-space
extern unsigned int COUNT;    
extern unsigned int PrintStep;       //Every "PrintStep", the program prints out to the screen or write to the output files
extern unsigned int Step;            //simulation step counter
extern unsigned int MPCstep;         ////MPCstep, used for MPC code
extern unsigned int RecordStep;      //wait for RecordStep to equilibrate

extern int SampleCounter[MAXT];      //count number of effective samplings at each correlation time
extern double Vacf[MAXT];            //velocity autocorrelation function
extern double Wacf[MAXT];            //angular velocity autocorrelation function
extern double R2[MAXT];              //mean-squared displacement
extern double U2[MAXT];              //orientaiton mean-squared displacement
extern double S_qt[MAX_q][MAXT];     //DYNAMIC structure factor
//extern double dq_sigma;

//variables for linked-list sorting of the collision cells
extern int Counter;                  //the index of the program
extern int Needle;                   //////used for old code, maybe useless
extern double beta;                  //////old code

extern int Rcelli;                   //
extern int Rcello;                   //
//linked-list for MD calculation
extern double OutOffi;               //the Radius of the outer shell for Verlet-list or Cell-list approaches, useful for MPC
extern double OutOffo;             
extern double OutOffcc;
extern double CellSize;
//cell-index of colloids
extern int Nix[MAXCOLL];             ////integer coordinate_x of colloids
extern int Niy[MAXCOLL];             //// *y*
extern int Niz[MAXCOLL];             //// *z*


////Used for MPC code, variables are used for building neighbour (verlet) list based on cell-list 
extern int iNeighbour;//number of neighbour without overlapping
extern int IntegrateNeighbour[MAXCOLL*100000];
extern int Neighbour[MAXCOLL*100000];//big enough Rho*(4/3)*PI*(Rcut^3-Rcol^3)
extern int HeatNeighbour[MAXCOLL*100000];
extern int StartOfNeighbour[MAXCOLL];
extern int NumberOfNeighbour[MAXCOLL];

//extern int *ll;
//extern int *hoc;
 
extern int icell;
extern int FlowCount;
extern int TempCount;
extern int DensCount;

//Parameters of Colloids at the MD part 
extern double NaT;

extern double Wall0;                  ////seting wall (along-z) positions if needed
extern double Wallz;                  ////
extern double UWalls;                 ////Wall potential
extern double USpring;                ////spring potential
extern double UBending;               ////bending potential
                                     
extern VECTOR CNP;
extern VECTOR COP;

extern double cCon;
extern double ColloidBending;         //force derived from bending
extern double ColloidSpring;          //force derived form chain molecules (colloids)
extern double BondLength;             
extern double RodLength;

extern VECTOR ColloidNewPositions[MAXCOLL];
extern VECTOR ColloidOldPositions[MAXCOLL];
extern VECTOR Colloidt0Positions[MAXCOLL];                      //colloid positions stored at the origin of correlation time

extern VECTOR ColloidPositions[MAXCOLL];                        //self-explained
extern VECTOR ColloidVelocities[MAXCOLL];
extern VECTOR ColloidForces[MAXCOLL];

extern VECTOR ColloidPositionsPriv[MAXCOLL][MAXTHREAD];         //for parallelization
extern VECTOR ColloidVelocitiesPriv[MAXCOLL][MAXTHREAD];
extern VECTOR ColloidForcesPriv[MAXCOLL][MAXTHREAD];

extern double Mcol;                                             //colloid mass
extern double invMcol;                                          //colloid inversed mass
extern double CutOffi;                                          // Cut-Off Radius
extern double Eshif;                                            // WCA potential Shifted energy
extern double HeatOn;
extern double CutOffo;
extern double Ecut;
extern double EcutYK;                                           //Yukawa cutoff energy
extern double FcutYK;                                           //Yukawa cutoff force
extern double CutOffcc;                                         ////cutoff for colloid-colloid 
extern double CutOffYK;                                         //

extern double CutOffw;                                          //cutoff for wall potential
extern double DebyeLength;           
extern double qq;                                               //related to effective charge
extern double epsilon;                                          //attraction strength

//extern VECTOR TempPosition;//
extern double Rcol;                                             //colloid radius
extern double Rcc;                                              //colloid-colloid interaction radius
extern double Inertia;                                          ////moment of inertia
extern double invInertia;                                       ////the inverse of the moment of inertia

//variables for Mean-Square-Displacement
extern int CountBox_x[MAXCOLL];
extern int CountBox_y[MAXCOLL];
extern int CountBox_z[MAXCOLL];

//extern double force;

extern double Rho;                 //density
extern double MDt;                 // MD Timestep

extern double BoxSize_x;           // Boxlengths
extern double BoxSize_y;
extern double BoxSize_z;
//some MD things

//add this to test the conservation of energy
extern double UKinetic;      // Kinetic Energy
extern double UPotential;    // Potential Energy
extern double UTotal;        // Total Energy

extern unsigned int output_restart_steps;		//how often restart files are written 
extern unsigned int RunStep;				//recording the time step for a next restart simulation
extern unsigned int NumberOfSteps;                      //number of MD (or MD-MPC coupling) steps
extern unsigned int NumberOfCells;                      ////for MPC, number of collision cells
extern unsigned int NumberOfParticles;                  ////number of SOLVENT particles
extern unsigned int NumberHeatColloids;			///////old code, not used anymore
extern unsigned int NumberOfColloids;			//
extern unsigned int Nx;					//number of colloids along x-direction
extern unsigned int Ny;					// *y*
extern unsigned int Nz;					// *z*
extern unsigned int Nobs;				//////old code, not used anymore
extern unsigned int NumberChainColloids;		//////number of chains, used for old code
//number of collision cells along x-y-z directions

extern int NumberOfBins;
extern double BinSize;	
//some subroutines
void WritePdb(FILE *FilePtr);			//write a .xyz movie
void Readdat(void);				//read parameters from outside data file
void Init(void);				////initialization
void Stream(void);				////streaming step in MPC step
void Collision(void);				////collision step in MPC step
void BDloop(void);				//Langevin Dynamics loop
void MSD(void);                                 //////old code
void Flow(void);				////flow field measurement in MPC code,has not been used yet
void GridCells(int);                            ////for MPC code
void StoreFlow(void);				////flow field measurement in MPC code,has not been used yet
void SlabTemp(void);				//////old code	
void TempK(void);				//////old code
void PrintTemp(void);				//////print temperature for old code
void Exchange(int EXPart);			////for mpc resorting particles
void Force(void);				//MD force calculation
//void Integrate(void);
void Integrate(int Switch);			//MD-LD integration
void Neighbours(void);				////for MPC code
//void CountNeighbours(void);
void LocalHeat(void);
void MeasureHeat(void);				//////old code
void Thermostats(void);				//////old code
//void ReSequence(void);
//void ReScale(void);

void SampleRDF(int Switch);			//function for measuring g(r) and S(q)
void SampleDiff(int Switch);			//function for measuring spatio-temperal correlations
void output_pos_vel_restart(int run, int t);	//write restart files
void initial_pos_vel_restart(int run, int t);	//read restart files
