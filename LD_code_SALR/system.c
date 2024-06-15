#include "system.h"

double Temperature;
unsigned long int JZ,JSR;
int HZ,IZ,KN[128];
float WN[128],FN[128];

unsigned long int seed[MAXTHREAD];
unsigned int N_threads;
unsigned int NUM_THREADS;


int SampleSteps;
double xif;//for the solvent friction
double c0;
double c1;
double c2;
double c3;
double c2c3;

double Cx;
double Gx;
double Ex;

unsigned int PrintStep;
unsigned int Step;
unsigned int RecordStep;

int SampleCounter[MAXT];
double Vacf[MAXT];
double Wacf[MAXT];
double R2[MAXT];
double U2[MAXT];
double S_qt[MAX_q][MAXT];
//int mark[MAXCOLL];//if Philic Colloid mark 0, if Phobic Colloid mark 1

int Counter;
int Needle;
double beta;


int CountBox_x[MAXCOLL];
int CountBox_y[MAXCOLL];
int CountBox_z[MAXCOLL];


//linked-list for MD calculation
double OutOffi;//the Radius of the outer shell for Verlet-list or Cell-list approaches
double OutOffo;
double OutOffcc;
double CellSize;
int Nix[MAXCOLL];
int Niy[MAXCOLL];
int Niz[MAXCOLL];

int iNeighbour;//number of neighbour without overlapping
int IntegrateNeighbour[MAXCOLL*100000];
int Neighbour[MAXCOLL*100000];
int StartOfNeighbour[MAXCOLL];
int NumberOfNeighbour[MAXCOLL];

int icell;
int Rcelli;
int Rcello;

//for md part
double Wall0;
double Wallz;

double UWalls;
double USpring;
double UBending;

double cCon;
double ColloidBending;//elastic constants
double ColloidSpring;
double RodLength;
double BondLength;

VECTOR CNP;
VECTOR COP;

VECTOR ColloidNewPositions[MAXCOLL];
VECTOR ColloidOldPositions[MAXCOLL];
VECTOR Colloidt0Positions[MAXCOLL];

VECTOR ColloidPositions[MAXCOLL];
VECTOR ColloidVelocities[MAXCOLL];
VECTOR ColloidForces[MAXCOLL];

VECTOR ColloidPositionsPriv[MAXCOLL][MAXTHREAD];
VECTOR ColloidVelocitiesPriv[MAXCOLL][MAXTHREAD];
VECTOR ColloidForcesPriv[MAXCOLL][MAXTHREAD];
double dq;
double Mcol;
double invMcol;
double CutOffi;        // Cut-Off Radius
double HeatOn;
double CutOffo;
double Eshif;         // Cut-Off energy Sigma
double Ecut;
double EcutYK;
double FcutYK;
double CutOffcc;
//double CutOffwy;
//double CutOffwz;
double CutOffw;
double CutOffYK;
double DebyeLength;
double qq;
double epsilon;

//extern VECTOR TempPosition;//
double Rcol;          //colloid radius
double Rcc;//colloid-colloid interaction radius
double Inertia;
double invInertia;
VECTOR COM;
VECTOR OmegaDt0;//angular acceleration

/*
VECTOR Omega[MAXCOLL];//angular velocity
VECTOR OmegaDt[MAXCOLL];//angular
VECTOR TorqueAll[MAXCOLL];//torque respect to the center of mass of rods
double Arms[2*MAXCOLL];
VECTOR Torque[MAXCOLL];
VECTOR Os[MAXCOLL];//unit vector along bond
VECTOR OsDot[MAXCOLL];//d(Os)/dt
VECTOR Os2Dot[MAXCOLL];//d2(Os)/dt2
VECTOR SignedForce[MAXCOLL];
VECTOR OmegaxOs[MAXCOLL];
VECTOR Os2Dot1[MAXCOLL];
VECTOR Os2Dot2[MAXCOLL];
*/
double Rho;                //density

double MDt;            // MD Timestep

double BoxSize_x;           // Boxlengths, cubic simulation boxes
double BoxSize_y;
double BoxSize_z;
//double CellSize;          //collision cubic collision cells

double UKinetic;      // Kinetic Energy
double UPotential;    // Potential Energy
double UTotal;        // Total Energy

unsigned int output_restart_steps;
unsigned int RunStep;
unsigned int NumberOfSteps;
unsigned int NumberOfCells;              
unsigned int NumberOfParticles;  
unsigned int NumberHeatColloids;
unsigned int NumberOfColloids;
unsigned int Nx;
unsigned int Ny;
unsigned int Nz;
unsigned int Nobs;
unsigned int NumberChainColloids;

//trasport phenomena properties
int NumberOfBins;
double BinSize;
