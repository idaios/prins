/*! \file 
*\author Pavlos Pavlidis
*/ 

#ifndef _STAT_H
#define _STAT_H

#include "qsort.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>

#define MAXLN 100
#define MAXFILE 1000


#define MATDIM1 20
#define MATDIM2 60
#define MAXENV 3
#define MAXAA 20


#ifdef _USE_PTHREADS
#include <pthread.h>
#include <unistd.h>

#define EXIT 127
#define BUSYWAIT 0
#define SORT 1
#define DIST 2

#endif



/* char iAA[20] = {'A', 'R', 'N', 'D', */
/* 		'C', 'E', 'Q', 'G',  */
/* 		'H', 'I', 'L', 'K',  */
/* 		'M', 'F', 'P', 'S',  */
/* 		'T', 'W', 'Y', 'V'}; */


/* int polarity[20] = { 0, 1, 1, 1,  */
/* 		     0, 1, 1, 0, */
/* 		     1, 0, 0, 1,  */
/* 		     0, 0, 0, 1, */
/* 		     1, 0, 1, 0 }; */

/* double hydropathy[20] = {1.8, -4.5, -3.5, -3.5,  */
/* 			 2.5, -3.5, -3.5, -0.4, */
/* 			 -3.2, 4.5, 3.8, -3.9,  */
/* 			 1.9, 2.8, -1.6, -0.8,  */
/* 			 -0.7, -0.9, -1.3, 4.2}; */

/* int charge[20] = {0, 1, 0, -1,  */
/* 		  0, -1, 0, 0,  */
/* 		  0, 0, 0, 1,  */
/* 		  0, 0, 0, 0,  */
/* 		  0, 0, 0, 0}; */


typedef struct 
{

  float x, y, z;
  char aa[4];
  int pos;
  char index[6];
  int counter;
  int neighbors;
  int size;
  int env;
  int file;
  char chain;
  double score;
}atom;

typedef struct
{

  double **distanceMatrix;

  int **pDistanceMatrix;

  int size;

  int *dim2;
  
} dist_struct;

typedef struct
{
  int polar;
  
  int charge;
  
  double hydropathy_index;

  int acidity;

  char name;
  double distance;
  char aa[4];

}features_struct;


typedef struct
{
  //double **neighborMatrix;
  //double ***pNeighborMatrix;
  //char **aminoacids;
  //double **distances;
  features_struct **features;
  
  int* numberOfNeighbors;
  
  int size;
  
  double distanceThreshold;
  
  int neighborsMin;
  
  int neighborsMax;
  
} neighbor_struct;

typedef struct
{
  int neighboraa;
  
  /* membrane */
  int counter;
  double z1;
  double z2;
} extraFeatures;

typedef struct
{
  /* stores the atoms for each file (we only need pointers)*/
  atom **atoms;

  /* the first dimension of the atoms*/
  int nofiles;
  
  /* the dimension of atoms */
  int y;
  /* the number of atoms */
  int natoms;
  int totalAtoms;
  int neighboraa;
  /* the M matrix of the paper */
  double **matrix;
  /* the first dimension of the matrix M*/
  int dim1;
  /* the second dimension of the matrix M*/
  int dim2;
  
  int *natomsEnv;

  int** envContacts;

  int envs;
  
} env;




typedef struct
{
  int activeEnvIndex;
  env *e;
} environments;


typedef struct
{
  char name[MAXFILE];
  
} file_struct;


typedef struct
{
  int aaCounts[MAXAA];
  int total;

} total_aa;

typedef struct
{
  atom **atoms;			/**< stores atoms info from all files. First dim: files, second dim: atoms per file */
  int *size;			/**< number of atoms per file */

  dist_struct *dis;

  int **indexUsedAtoms;

  int *nUsedAtoms;

  double *scores;

  double *membraneZ1;

  double *membraneZ2;
  
} data_struct;


file_struct * allFiles;

dist_struct * dist;
//atom  *atoms;
extraFeatures *ef;
neighbor_struct *allNeighbors;
total_aa *taa;
data_struct *data;
//env *tip;

int membraneLim;


#ifdef _USE_PTHREADS

typedef struct
{
  int threadID;
  int threadTOTAL;
  
  int threadBARRIER;
  int threadOPERATION;

}threadData_t;


pthread_t * workerThreadL;
threadData_t *threadDataL;

void sortIndexes_thread (int tid, int threads, dist_struct *dist);

void sortIndexes_parallel(dist_struct *dist);

void calcDistances_thread( int tid, int threads, dist_struct *dist);

void calcDistances_parallel( dist_struct  *dist);

void syncThreadsBARRIER();

#endif


void commandLineParser(int argc, char **argv,
		       char *infile,
		       char *ofile, 
		       char *listfile,
		       double *threshold,
		       char *matrixOutputFileName,
		       int *input,
		       int *output, 
		       int *nrand,
		       char *randomOutputFileName);



int parseInputFile(char *infile);

int parseInputFile2(char *infile, atom **atoms, int fileindex);

int parseListFile(char *listfile);

void calculateDistances(dist_struct *dist, env *e, double threshold);

void sortIndexes ( dist_struct *dist);

void printResults();

/**
 * \fn
 * Converts the string of AA name to integer (0..19)
 * 
 * @param name a char array that describes the AA name (char[4] )
 * 
 * @return integer between 0 and 19
 */
int convertAAToInt( char * name );

void copyMatrix(double **source, double **destination, int size);

int generateDistanceMatrix(dist_struct *dist, env *e, double threshold);

int totalEnvironmentMatrix(env *e);

void initializeEnvDegree(env *e);

void classifyAtoms(env *e, dist_struct *dist);

int incrementEnvironmentMatrix(env* e, dist_struct *dist);

int incrementEnvironmentMatrixUniform(env* e, dist_struct *dist, int envMode, double threshold);

void scoreData_i(env* e, data_struct *d, int index, FILE *outputfile, char *fn, double *scores);

void scoreRandomData_i(env *e, data_struct *d, int index, double **randomScores, double *mean, double *sd, int nrand, int *nsums);

void print2DArrayDouble(double **a, int x, int *y);

void print2DArrayInt(int **a, int x, int *y);

void printEnvironMatrix(env *e, FILE *output);

void readEnvironmentMatrixFromFile(FILE *matrixInputFile, env *e);

void freeGlobalMemory();

void freeData();

int getFileNames(char *listfile, char ***files, int currentX, int currentY);

int getFileNamesFromCMD(int argc, char **argv, char ***files, int currentX, int currentY);

#endif
