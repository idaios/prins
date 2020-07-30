/*! \file 
*\author Pavlos Pavlidis
*/ 
#include "prins.h"


double gettime(void)
{
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
}


#ifdef _USE_PTHREADS


void initializeThreadData(threadData_t * cur, int i, int threads){
  cur->threadID = i;
  
  cur->threadTOTAL = threads;
  
  cur->threadBARRIER = 0;
  
  cur->threadOPERATION = BUSYWAIT;
}


void *thread(void * x)
{
  threadData_t *currentThread = (threadData_t*)x;
  
  int tid = currentThread->threadID, threads = currentThread->threadTOTAL;

  while(1)
    {
      sleep(0);
      
      if(threadDataL[tid].threadOPERATION == EXIT )
	return NULL;

      if(threadDataL[tid].threadOPERATION == SORT )
	{
	  sortIndexes_thread (tid, threads, dist);
	  
	  threadDataL[tid].threadBARRIER = 1;
	  
	  while(threadDataL[tid].threadBARRIER == 1)
	    sleep(0);
	  
	}
      else if(threadDataL[tid].threadOPERATION == DIST)
	{

	  calcDistances_thread(tid, threads, dist);

	  threadDataL[tid].threadBARRIER = 1;
	  
	  while(threadDataL[tid].threadBARRIER == 1)
	    sleep(0);

	}

    }
  return NULL;
}


void terminateWorkerThreads()
{
	int i, threads=threadDataL[0].threadTOTAL;
	
	for(i=0;i<threads;i++)
		threadDataL[i].threadOPERATION = EXIT;

	sleep(0);			

	for(i=1;i<threads;i++)
		pthread_join( workerThreadL[i-1], NULL);

	free(threadDataL);
}


void syncThreadsBARRIER()
{
	int i, threads = threadDataL[0].threadTOTAL, barrierS=0;

	while(barrierS!=threads)
	{
		barrierS=0;
		for(i=0;i<threads;i++)
			barrierS += threadDataL[i].threadBARRIER;
	}

	for(i=0;i<threads;i++)
	{
		threadDataL[i].threadOPERATION=BUSYWAIT;
		threadDataL[i].threadBARRIER=0;
	}
}


#endif 




 /** 
  * Stores all the atom information from all .pdb files in the list
  * 
  * @param listfile a file that contains all the file names that will be used in the analysis. Each file is a pdb file
  * @param f a file_struct
  * @param d a data_struct
  * 
  * @return 
  */
 static int parseAllFiles(char* listfile, data_struct **d)
 {
   int i,nofiles = 0, n = 0;

   nofiles = parseListFile(listfile);

   (*d) = malloc(sizeof(data_struct));

   (*d)->atoms = malloc(nofiles * sizeof(atom*) );

   (*d)->size = calloc(nofiles, sizeof(int));

   (*d)->dis = malloc(nofiles *sizeof(dist_struct) );

   (*d)->membraneZ1 = calloc(nofiles, sizeof(double));

   (*d)->membraneZ2 = calloc(nofiles, sizeof(double));

   for(i = 0; i <nofiles; ++i)
     {

       (*d)->membraneZ1[i] = -99999999.;

       (*d)->membraneZ2[i] = 99999999.;


       n = parseInputFile2(allFiles[i].name, &((*d)->atoms[i]), i);
       (*d)->size[i] = n;

       (*d)->membraneZ1[i] = ef->z1;

       (*d)->membraneZ2[i] = ef->z2;
     }

   return nofiles;
 }


int getFilename(char *f, int sz, char t[MAXFILE] )
{
  int i=0, j = -1, k = 0, end = 0;
  while( t[i] != '\0' && i < MAXFILE)
     {
       end++;
       if(t[i] == '/')
	 j = i;
       i++;
     }

   

   k = 0;
   for(i = j+1; i < end; ++i,j++)
     {

       f[k++] = t[i];
       if(k > sz)
	 return 0;
     }

   f[k] = '\0';
   return 1;
 }

void freeEnv(env *e)
{
  int i;

  for(i = 0; i < e->dim1; ++i)
    free(e->matrix[i]);
  free(e->matrix);

  free(e->natomsEnv);

  for(i=0; i < e->envs; ++i)
    free(e->envContacts[i]);
  free(e->envContacts);

  free(e);
  
}




int runList(double threshold, char* listfile, int envMode, char* matrixFileName, int input, int output, double **scores, double ***randomScores, int nrand, char *randomOutputFileName)
{
  assert(ef->counter == 0 || ef->counter == 2);
  
  int i, n, j, jj = 0,   nofiles,  activeEnvIndex = 0,  envi = 0, nsums = 0;
  
  char *outputfilename = calloc(50, sizeof(char)), *fname = calloc(50, sizeof(char));
  
  FILE *outputfile, *randomOutputFile = NULL;

  double mean = 0, sd = 0;
  

  nofiles = parseAllFiles(listfile, &data);

  *scores = calloc( nofiles, sizeof(double) );

  if( nrand > 0 )
    {

      *randomScores = calloc( nofiles, sizeof(double*) );
      
      for( i = 0; i < nofiles; ++i)
	{
	  (*randomScores)[i] = calloc( nrand, sizeof(double) );
	}
    }
  
  /* allocate memory for the used atoms */
  data->indexUsedAtoms = malloc( nofiles * sizeof(int*) );
  
  data->nUsedAtoms = calloc( nofiles, sizeof(int) );

  
  
  /* not sure if local or global variable is better here */
  environments *envir = malloc(MAXENV * sizeof(environments) );
  
  /* environment 0 : uniform
     environemnt 1 : membrane
     environment 2 : degree of connectivity  */
  for(i = 0; i < MAXENV; ++i)
    {
      activeEnvIndex = i;
      
      //printf("active environment index: %d %d\n", activeEnvIndex, envMode);
      
      /* inactive environment... continue */
      if( (activeEnvIndex != envMode) )
	continue;
      
      envir[i].e = malloc( sizeof(env) );
      
      //fprintf(stderr, "sizeof env[i]->e: %lu\n", sizeof(env));

      (envir[i].e)->nofiles = nofiles;
      
      //(envir[i].e)->y = 10;

      //(envir[i].e)->atoms = malloc(envir[i].e->y * sizeof(atom*) );
      
      (envir[i].e)->neighboraa = 0;
      
      if(i == 0)
	{
	  (envir[i].e)->dim1 = envir[i].e->dim2 = MATDIM1;
	  
	  (envir[i].e)->matrix = malloc(MATDIM1 * sizeof(double*) );
	  
	  for( j = 0; j < MATDIM1; ++j )
	    (envir[i].e)->matrix[j] = calloc(MATDIM1, sizeof(double) );

	  envir[i].e->dim1 = envir[i].e->dim2 = MATDIM1;
	}
      else if( i == 1 || i == 2)
	{
	  envir[i].e->matrix = malloc(MATDIM2 * sizeof(double*) );
	  for( j = 0; j < MATDIM2; ++j )
	    envir[i].e->matrix[j] = calloc(MATDIM2, sizeof(double) );
	  
	  envir[i].e->dim1 = envir[i].e->dim2 = MATDIM2;
	}
      
    }
  
  /* see which environments are active */
  for(envi = 0; envi < MAXENV; ++envi)
    {
      activeEnvIndex = envi;

      //printf("MODE Size: %d %d\n", activeEnvIndex, envMode);
      
      /* inactive environment... continue */
      if( (activeEnvIndex != envMode ) )
	continue;
      
      if(envi == 2)
	initializeEnvDegree(envir[envi].e);
      
      /* read the data from the file */
      for(i = 0; i < nofiles; ++i)
	{
	  n = data->size[i];

	  /* realloc space for the atoms */
	  
	  envir[envi].e->atoms = malloc( n * sizeof(atom*)); //realloc(envir[envi].e->atoms, n * sizeof(atom*));
	  envir[envi].e->y = n;
	      
	      
	      //}
	  
	  /* counter == 2 means that the z coordinates have been set */
	  jj = 0;
	  for(j = 0; j < n; ++j)
	    {
	      /* if the atom is between the membranes */
	      if(ef->counter == 0 || (data->atoms[i][j].z > data->membraneZ1[i] &&
				      data->atoms[i][j].z < data->membraneZ2[i]  ) )
		{
		  // (envir[envi].e->atoms[jj]) is a pointer
		  envir[envi].e->atoms[jj] = &(data->atoms[i][j]); //atoms[j];
		  
		  ++jj;
		  
		  data->nUsedAtoms[i] = jj;
		  if(jj == 1)
		    data->indexUsedAtoms[i] = calloc(jj, sizeof(int) );
		  else
		    data->indexUsedAtoms[i] = realloc( data->indexUsedAtoms[i], jj * sizeof(int) );		      
		  data->indexUsedAtoms[i][jj-1] = j;
		}
	      
	      
	    }
	  
	  envir[envi].e->natoms = jj;
	  /* generate distance matrix for the specific file and environment */
	  if(generateDistanceMatrix( &(data->dis[i]), envir[envi].e, threshold) == 0 )
	    continue;
	  	      
	  classifyAtoms(envir[envi].e, &(data->dis[i]));
	      
	  if(input == 0)
	    {
	      incrementEnvironmentMatrix(envir[envi].e, &(data->dis[i]));
	    }
	  

	    free(envir[envi].e->atoms);
	  
	}

      if(input == 1)
	{
	  
	  FILE *matrixInputFile = fopen(matrixFileName, "r");
	  
	  readEnvironmentMatrixFromFile(matrixInputFile, envir[envi].e);
	  
	}
      else
	{
	  totalEnvironmentMatrix(envir[envi].e);
	}
    	  
      
      if(output == 1)
	{

	  FILE *matrixOutputFile = fopen(matrixFileName, "w");
      

	  
	  printEnvironMatrix(envir[envi].e, matrixOutputFile);
      
	  fclose(matrixOutputFile);

	}
      
      if(envi == 2)
	{
	  
	  
	  if(nrand > 0)
	    {
	      
	      if( randomOutputFileName[0] != '\0' )
		randomOutputFile = fopen(randomOutputFileName, "w");
	      else
		randomOutputFile = stdout;
	      
	      fprintf( randomOutputFile, "score\tz-score");
	      
	      for( j = 0; j <nrand; ++j)  
		fprintf( randomOutputFile, "\tr%d", j);
	      
	      fprintf( randomOutputFile, "\n");
	      
	    }
	  
	  
	  /* scoring each aa for each file */
	  for( i = 0; i < nofiles; ++i)
	    {
	      getFilename(fname, 50, allFiles[i].name);
	      
	      sprintf(outputfilename, "score_%s.txt", fname);

	      outputfile = fopen(outputfilename, "w");
	      
	      fprintf(outputfile, "aa\tindex\tpos\tfileIndex\tfile\tscore\taverageScore\tmaxScore\tchain\tNeighbors\n");
	      
	      /* score file_i */
	      scoreData_i(envir[envi].e, data, i, outputfile, fname, *scores);

	      if(nrand > 0)
		{
		  
		  scoreRandomData_i( envir[envi].e, data, i, *randomScores, &mean, &sd, nrand, &nsums);

		  fprintf(randomOutputFile, "%f\t%f\t", (*scores)[i] / (double)nsums, (mean - (*scores)[i] / (double)nsums )/sd );

		  for( j = 0; j <nrand; ++j)  
		    fprintf(randomOutputFile, "%f\t", (*randomScores)[i][j]); 
		  		  
		  fprintf(randomOutputFile, "\n");
		  
		}
	      
	      fclose(outputfile);
	    }

	  
	  if( randomOutputFileName && nrand > 0)
	    fclose(randomOutputFile);

	}
            
    }
  
  for(i=0; i < nofiles; ++i)
    {
      free(data->atoms[i]);
    }

  free(data->atoms);

  free(data->size);
  
  freeGlobalMemory();
  
  for(i=0; i < MAXENV; ++i)
    {
      if( i != envMode )
	continue;
      /* fprintf(stderr, "natoms: %d\n", envir[i].e->natoms); */
      freeEnv(envir[i].e);
    }
  
  free(envir);

  free(outputfilename);

  free(fname);
  
  return 1;
}


int getData(char **files, data_struct **d, int k)
{

  int i, n;

  (*d) = malloc( sizeof(data_struct) );

  (*d)->atoms = malloc(k * sizeof(atom*) );
  
  (*d)->size = calloc(k, sizeof(int) );

  (*d)->dis = malloc(k * sizeof(dist_struct) );

  (*d)->membraneZ1 = calloc( k, sizeof(double) );

  (*d)->membraneZ2 = calloc(k, sizeof(double) );

  

  for(i = 0; i < k; ++i)
    {

      (*d)->membraneZ1[i] = -99999999.;

      (*d)->membraneZ2[i] = 99999999.;
      
      n = parseInputFile2(files[i], &((*d)->atoms[i]), i);
       
       (*d)->size[i] = n;

       (*d)->membraneZ1[i] = ef->z1;

       (*d)->membraneZ2[i] = ef->z2;
     }
  
  return k;

}


int runData(double threshold, data_struct *data, char **files, int nofiles, int envMode, char* matrixFileName, int input, int output, double **scores, double ***randomScores, int nrand, char *randomOutputFileName)
{
  
  assert(ef->counter == 0 || ef->counter == 2);
  
  int i, n, j, jj = 0, activeEnvIndex = 0,  envi = 0, nsums = 0;
  
  char *outputfilename = calloc(50, sizeof(char)), *fname = calloc(50, sizeof(char));
  
  FILE *outputfile, *randomOutputFile = NULL;

  double mean = 0, sd = 0;
  

  *scores = calloc( nofiles, sizeof(double) );

  if( nrand > 0 )
    {

      *randomScores = calloc( nofiles, sizeof(double*) );
      
      for( i = 0; i < nofiles; ++i)
	{
	  (*randomScores)[i] = calloc( nrand, sizeof(double) );
	}
    }
  
  /* allocate memory for the used atoms */
  data->indexUsedAtoms = malloc( nofiles * sizeof(int*) );
  
  data->nUsedAtoms = calloc( nofiles, sizeof(int) );

  
  
  /* not sure if local or global variable is better here */
  environments *envir = malloc(MAXENV * sizeof(environments) );
  
  /* environment 0 : uniform
     environemnt 1 : membrane
     environment 2 : degree of connectivity  */
  for(i = 0; i < MAXENV; ++i)
    {
      activeEnvIndex = i;
      
      //printf("active environment index: %d %d\n", activeEnvIndex, envMode);
      
      /* inactive environment... continue */
      if( (activeEnvIndex != envMode) )
	continue;
      
      envir[i].e = malloc( sizeof(env) );
      
      //fprintf(stderr, "sizeof env[i]->e: %lu\n", sizeof(env));

      (envir[i].e)->nofiles = nofiles;
      
      //(envir[i].e)->y = 10;
      
      //(envir[i].e)->atoms = malloc(envir[i].e->y * sizeof(atom*) );
      
      (envir[i].e)->neighboraa = 0;
      
      if(i == 0)
	{
	  (envir[i].e)->dim1 = envir[i].e->dim2 = MATDIM1;
	  
	  (envir[i].e)->matrix = malloc(MATDIM1 * sizeof(double*) );
	  
	  for( j = 0; j < MATDIM1; ++j )
	    (envir[i].e)->matrix[j] = calloc(MATDIM1, sizeof(double) );
	  
	  envir[i].e->dim1 = envir[i].e->dim2 = MATDIM1;
	}
      else if( i == 1 || i == 2)
	{
	  envir[i].e->matrix = malloc(MATDIM2 * sizeof(double*) );
	  for( j = 0; j < MATDIM2; ++j )
	    envir[i].e->matrix[j] = calloc(MATDIM2, sizeof(double) );
	  
	  envir[i].e->dim1 = envir[i].e->dim2 = MATDIM2;
	}
      
    }
  
  /* see which environments are active */
  for(envi = 0; envi < MAXENV; ++envi)
    {
      activeEnvIndex = envi;

      //printf("MODE Size: %d %d\n", activeEnvIndex, envMode);
      
      /* inactive environment... continue */
      if( (activeEnvIndex != envMode ) )
	continue;
      
      if(envi == 2)
	initializeEnvDegree(envir[envi].e);
      
      /* read the data from the file */
      for(i = 0; i < nofiles; ++i)
	{
	  n = data->size[i];

	  /* realloc space for the atoms */
	  
	  envir[envi].e->atoms = malloc( n * sizeof(atom*)); //realloc(envir[envi].e->atoms, n * sizeof(atom*));
	  envir[envi].e->y = n;
	      
	      
	      //}
	  
	  /* counter == 2 means that the z coordinates have been set */
	  jj = 0;
	  for(j = 0; j < n; ++j)
	    {
	      /* if the atom is between the membranes */
	      if(ef->counter == 0 || (data->atoms[i][j].z > data->membraneZ1[i] &&
				      data->atoms[i][j].z < data->membraneZ2[i]  ) )
		{
		  // (envir[envi].e->atoms[jj]) is a pointer
		  envir[envi].e->atoms[jj] = &(data->atoms[i][j]); //atoms[j];
		  
		  ++jj;
		  
		  data->nUsedAtoms[i] = jj;
		  if(jj == 1)
		    data->indexUsedAtoms[i] = calloc(jj, sizeof(int) );
		  else
		    data->indexUsedAtoms[i] = realloc( data->indexUsedAtoms[i], jj * sizeof(int) );		      
		  data->indexUsedAtoms[i][jj-1] = j;
		}
	      
	      
	    }
	  
	  envir[envi].e->natoms = jj;
	  /* generate distance matrix for the specific file and environment */
	  if(generateDistanceMatrix( &(data->dis[i]), envir[envi].e, threshold) == 0 )
	    continue;
	  	      
	  classifyAtoms(envir[envi].e, &(data->dis[i]));
	      
	  if(input == 0)
	    {
	      incrementEnvironmentMatrix(envir[envi].e, &(data->dis[i]));
	    }
	  

	    free(envir[envi].e->atoms);
	  
	}

      if(input == 1)
	{
	  
	  FILE *matrixInputFile = fopen(matrixFileName, "r");
	  
	  readEnvironmentMatrixFromFile(matrixInputFile, envir[envi].e);
	  
	}
      else
	{
	  totalEnvironmentMatrix(envir[envi].e);
	}
    	  
      
      if(output == 1)
	{

	  FILE *matrixOutputFile = fopen(matrixFileName, "w");
      

	  
	  printEnvironMatrix(envir[envi].e, matrixOutputFile);
      
	  fclose(matrixOutputFile);

	}
      
      if(envi == 2)
	{
	  
	  
	  if(nrand > 0)
	    {
	      
	      if( randomOutputFileName[0] != '\0' )
		randomOutputFile = fopen(randomOutputFileName, "w");
	      else
		randomOutputFile = stdout;
	      
	      fprintf( randomOutputFile, "score\tz-score");
	      
	      for( j = 0; j <nrand; ++j)  
		fprintf( randomOutputFile, "\tr%d", j);
	      
	      fprintf( randomOutputFile, "\n");
	      
	    }
	  
	  
	  /* scoring each aa for each file */
	  for( i = 0; i < nofiles; ++i)
	    {
	      getFilename(fname, 50, files[i]);
	      
	      sprintf(outputfilename, "score_%s.txt", fname);

	      outputfile = fopen(outputfilename, "w");
	      
	      fprintf(outputfile, "aa\tindex\tpos\tfileIndex\tfile\tscore\taverageScore\tmaxScore\tchain\tNeighbors\n");
	      
	      /* score file_i */
	      scoreData_i(envir[envi].e, data, i, outputfile, fname, *scores);

	      if(nrand > 0)
		{
		  
		  scoreRandomData_i( envir[envi].e, data, i, *randomScores, &mean, &sd, nrand, &nsums);

		  fprintf(randomOutputFile, "%f\t%f\t", (*scores)[i] / (double)nsums, (mean - (*scores)[i] / (double)nsums )/sd );

		  for( j = 0; j <nrand; ++j)  
		    fprintf(randomOutputFile, "%f\t", (*randomScores)[i][j]); 
		  		  
		  fprintf(randomOutputFile, "\n");
		  
		}
	      
	      fclose(outputfile);
	    }

	  
	  if( randomOutputFileName && nrand > 0)
	    fclose(randomOutputFile);

	}
            
    }
  
  for(i=0; i < nofiles; ++i)
    {
      free(data->atoms[i]);
    }

  free(data->atoms);

  free(data->size);
  
  freeGlobalMemory();
  
  for(i=0; i < MAXENV; ++i)
    {
      if( i != envMode )
	continue;
      /* fprintf(stderr, "natoms: %d\n", envir[i].e->natoms); */
      freeEnv(envir[i].e);
    }
  
  free(envir);

  free(outputfilename);

  free(fname);
  
  return 1;
}





int main(int argc, char **argv)
{

   

#ifdef _USE_PTHREADS
   
   int threads = 2;
   
#endif

   double t0 = gettime();

   char infile[MAXFILE];

   char ofile[MAXFILE];

   char matrixOutputFileName[MAXFILE], randomOutputFileName[MAXFILE];

   randomOutputFileName[0] = '\0';

   strcpy(matrixOutputFileName, "AminoacidInteractionMatrix.txt");

   char listfile[MAXFILE]; 
   listfile[0] = '\0';

   //int nofiles = 0;

   //   int i;

   double threshold = 6.5;

   int input = 0, output = 0, nrand = 0, nfiles = 0, j = 0;

   double *scores, **randomScores;

   char **files = NULL;

   commandLineParser(argc, argv, infile, ofile, listfile, &threshold, matrixOutputFileName, &input, &output, &nrand, randomOutputFileName);

   ef = malloc( sizeof(extraFeatures) );

   /* initialize values for membrane*/
   ef->counter = 0;
   ef->z1 = -10000.;
   ef->z2 = 10000.;

   if( ( nfiles = getFileNamesFromCMD(argc, argv, &files, 0, 0) ) > 0 )
     {
       
       
       taa = malloc( sizeof(total_aa) );
       
       for(j = 0; j < MAXAA; ++j)
	 taa->aaCounts[j] = 0;
       taa->total = 0;
              
       getData ( files, &data, nfiles);
       
       assert( listfile[0] == '\0' );

       runData( threshold, data, files, nfiles, 2, matrixOutputFileName, input, output, &scores, &randomScores, nrand, randomOutputFileName);
       
     }
   else if(listfile[0] != '\0')
    {
      runList( threshold, listfile, 2, matrixOutputFileName, input, output, &scores, &randomScores, nrand, randomOutputFileName);
      
      /* double class1 = analyzeNeighborFeatures( allNeighbors, 1); */
      
      /* double class2 = analyzeNeighborFeatures( allNeighbors, 2); */

      /* double class3 = analyzeNeighborFeatures( allNeighbors, 3); */

      /* printf("class1: %f\tclass2: %f\tclass3: %f\n", class1, class2, class3); */
      
    }
  
  fprintf(stderr, "Total time: %f\n", gettime() - t0);

  
/*   exit(1); */

/*   int n = parseInputFile(infile); */



/*   dist = malloc( sizeof(dist_struct)); */
/*   dist->size = n; */
  
/*   dist-> distanceMatrix = malloc( n * sizeof(double* ) ); */

/*   dist->pDistanceMatrix = malloc ( n * sizeof(double**)); */

  
/*   for( i = 0; i < n; ++i ) */
/*     { */
/*       dist->distanceMatrix[i] = calloc (n, sizeof(double )); */
      
/*       dist->pDistanceMatrix[i] = malloc (n * sizeof(double * )); */
/*     } */





/* #ifdef _USE_PTHREADS */

/*   workerThreadL = malloc( sizeof(pthread_t) * (threads-1));	 */
/*   assert(workerThreadL != NULL); */

/*   threadDataL = malloc( threads * sizeof(threadData_t)); */
/*   assert(threadDataL != NULL); */

/*   for( i = 0; i<threads; ++i ) */
/*     { */
/*       initializeThreadData ( &threadDataL[i], i, threads); */
/*     } */

/*   for ( i = 1; i < threads; ++i) */
/*     { */
/*       pthread_create(&workerThreadL[i-1], NULL, thread, (void*) (&threadDataL[i])); */
/*     } */

/*   calcDistances_parallel(dist); */

/*   sortIndexes_parallel(dist); */
  
/*   terminateWorkerThreads(); */

/* #else */
/*   calculateDistances( dist ); */
/*   sortIndexes( dist ); */
/* #endif */


/*   /\* double **distanceMatrix = malloc( n * sizeof(double* ) ); *\/ */

/*   /\* double ***pDistanceMatrix = malloc( n * sizeof(double**) ); *\/ */
  
/*   /\* for( i = 0; i<n; ++i) *\/ */
/*   /\*   { *\/ */
/*   /\*     distanceMatrix[i] = malloc ( n * sizeof (double ) ); *\/ */

/*   /\*     pDistanceMatrix[i] = malloc (n * sizeof( double* ) ); *\/ */
/*   /\*   } *\/ */


/*   //printResults(); */


/*   fprintf(stderr, "Total time: %f\n", gettime() - t0); */

  free(ef);

  return 1;
}
 
