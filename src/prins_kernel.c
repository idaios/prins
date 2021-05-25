/*! \file 
*\author Pavlos Pavlidis
*/ 
#include "prins.h"


int convertAAToInt( char * name )
{
  if(!strcmp("ALA", name)) return 0;
  if(!strcmp("ARG", name)) return 1;
  if(!strcmp("ASN", name)) return 2;
  if(!strcmp("ASP", name)) return 3;
  if(!strcmp("CYS", name)) return 4;
  if(!strcmp("GLU", name)) return 5;
  if(!strcmp("GLN", name)) return 6;
  if(!strcmp("GLY", name)) return 7;
  if(!strcmp("HIS", name)) return 8;
  if(!strcmp("ILE", name)) return 9;
  if(!strcmp("LEU", name)) return 10;
  if(!strcmp("LYS", name)) return 11;
  if(!strcmp("MET", name)) return 12;
  if(!strcmp("PHE", name)) return 13;
  if(!strcmp("PRO", name)) return 14;
  if(!strcmp("SER", name)) return 15;
  if(!strcmp("THR", name)) return 16;
  if(!strcmp("TRP", name)) return 17;
  if(!strcmp("TYR", name)) return 18;
  if(!strcmp("VAL", name)) return 19;
  
  fprintf(stderr, "Aminoacid %s cannot be found in the list\n", name);
  exit(-1);
  //if(!strcmp("VAL", name)) return 21;
  return 100;
}


int convertIntToAA( int i, char * name )
{
  
  int j = i % 20;

  switch(j){
  case 0: strcpy(name, "ALA"); 
    break;
  case 1: strcpy(name, "ARG"); 
    break;  
  case 2: strcpy(name, "ASN"); 
    break;
  case 3: strcpy(name, "ASP"); 
    break;
  case 4: strcpy(name, "CYS"); 
    break;
  case 5: strcpy(name, "GLU"); 
    break;
  case 6: strcpy(name, "GLN"); 
    break;
  case 7: strcpy(name, "GLY"); 
    break;
  case 8: strcpy(name, "HIS"); 
    break;
  case 9: strcpy(name, "ILE"); 
    break;
  case 10: strcpy(name, "LEU"); 
    break;
  case 11: strcpy(name, "LYS"); 
    break;
  case 12: strcpy(name, "MET"); 
    break;
  case 13: strcpy(name, "PHE"); 
    break;
  case 14: strcpy(name, "PRO"); 
    break;
  case 15: strcpy(name, "SER"); 
    break;
  case 16: strcpy(name, "THR"); 
    break;
  case 17: strcpy(name, "TRP"); 
    break;
  case 18: strcpy(name, "TYR"); 
    break;
  case 19: strcpy(name, "VAL"); 
    break;
  default: 
    return 0;
    break;
  }
  return 1;
}



void printEnvironMatrix(env *e, FILE *outputfile)
{
  int i,j;
  char aa[4];


  
  for(j = 0; j < e->dim2; ++j)
    {      
      if( convertIntToAA(j, aa) == 1 && j != e->dim2 -1)
	fprintf(outputfile, "%s ", aa);
      else if( convertIntToAA(j, aa) == 1 && j == e->dim2 -1)
	fprintf(outputfile, "%s\n", aa);
      
      if(convertIntToAA(j, aa) == 0)
	{
	  fprintf(stderr, "Error int to AA assignemnt (%d)\n", j);
	  exit(-1);
	}
    }
      


  for( i = 0; i < e->dim1; ++i)
    {
      
      for(j = 0; j < e->dim2; ++j)
	{

	  if( j == 0)
	    {
	      convertIntToAA(i, aa);
	      fprintf(outputfile, "%s ", aa);
	    }
	  
	  if(j != e->dim2 - 1)
	    fprintf(outputfile, "%f ", e->matrix[i][j]);
	  else
	    fprintf(outputfile, "%f\n", e->matrix[i][j]);
	}
    }


}

void print2DArrayDouble(double **a, int x, int *y)
{
  int i,j;
  for(i = 0; i < x; ++i)
    {
      for(j = 0; j < y[i]; ++j)
	printf("%f ", a[i][j]);
      printf("\n");
    }
  
}


void print2DArrayInt(int **a, int x, int *y)
{
  int i,j;
  for(i = 0; i < x; ++i)
    {
      for(j = 0; j < y[i]; ++j)
	printf("%d ", a[i][j]);
      printf("\n");
    }
  
}



/**
 * Scores each aminoacid based on what neighbors it has.
 Using a set of pdb files a matrix M has been constructed that scores each interaction. Then, each aminoacid is scored based on its neighbors (interactions).
 *
 * @param e the environment structure. It stores the matrix M with the interactions
 * @param d the data structure. It stores the information for all aminoacids
 * @param index indexes the dataset (i.e. which pdb file we are working on)
 * @param outputfile handles the output
 * @param fn a string that gives the filename of the output
 *
 * @param scores keeps the total score for each sequence
 *
 */
void scoreData_i(env* e, data_struct *d, int index, FILE *outputfile, char *fn, double *scores)
{


  int i,j, n = d->nUsedAtoms[index], ind0, ind, ind_tmp, aaInd1, aaInd2;
  
  char neighbors[10000], conv[1000];

  

  char maxNeighbor[4]; //  = calloc(4, sizeof(char));
  

  double score = 0., averageScore = 0., maxScore = -999999999.9999;

  scores[ index ] = 0.;
  
  if(n == 0)
     return;

  
    
  for(i = 0; i < n; ++i)
    {
      /* map i to ind0. i refers to the atoms have been used after membrane filtering. However, the name of the aa has been stored in data_struct BEFORE filtering. Thus, we need to map the coordinates */
      ind0 = d->indexUsedAtoms[index][i]; 

      maxScore = -999999999.9999;
      
      aaInd1 = ( convertAAToInt(d->atoms[index][ind0].aa) + (d->atoms[index][ind0].env - 1) * 20 );
      
      score = 0.;
      
      neighbors[0] = '\0';



      for(j = 0; j < d->dis[index].dim2[i]; ++j)
	{
	  ind_tmp = d->dis[index].pDistanceMatrix[i][j];

	  ind = d->indexUsedAtoms[index][ind_tmp];
	  
	  aaInd2 = ( convertAAToInt(d->atoms[index][ind].aa) + (d->atoms[index][ind].env - 1) * 20 );
	  

	  if( e->matrix[aaInd1][aaInd2] > maxScore)
	    {
	      maxScore = e->matrix[aaInd1][aaInd2];
	      strcpy(maxNeighbor, d->atoms[index][ind].aa);
	    }
	  
	  score += e->matrix[aaInd1][aaInd2];
	  strcat(neighbors, d->atoms[index][ind].aa);
	  strcat(neighbors, ":");
	  sprintf(conv, "%d", d->atoms[index][ind].pos);
	  strcat(neighbors,conv);
	  strcat(neighbors, ":");
	  sprintf(conv, "%d", d->atoms[index][ind].env);
	  strcat(neighbors, conv);
	  strcat(neighbors, ":");
	  sprintf(conv, "%f", e->matrix[aaInd1][aaInd2]);
	  strcat(neighbors, conv);
	  if(j != d->dis[index].dim2[i] - 1)
	    strcat(neighbors, ",");
	}
      if(d->dis[index].dim2[i] > 0)
	averageScore = score/d->dis[index].dim2[i];            
      else
	averageScore = 0.;

      
      if(d->dis[index].dim2[i] == 0)
	{
	  neighbors[0]='-';
	  neighbors[1] = '\0';
	}
      
      d->atoms[index][ind0].score = score;

      scores[ index ] += score;
      

      fprintf(outputfile, "%s\t%s\t%d\t%d\t%s\t%f\t%f\t%f\t*%s*\t%c\t%s\n", d->atoms[index][ind0].aa, d->atoms[index][ind0].index, d->atoms[index][ind0].pos, d->atoms[index][ind0].file, fn, d->atoms[index][ind0].score, averageScore, maxScore, maxNeighbor,d->atoms[index][ind0].chain, neighbors);
    }

  
}

static int rand_int(int n) {
  int limit = RAND_MAX - RAND_MAX % n;
  int rnd;

  do {
    rnd = rand();
  } while (rnd >= limit);
  return rnd % n;
}


void shuffle(int *array, int *randomArray, int n) {
  int i, j, tmp;

  for( i = 0; i < n; ++i)
    {
      randomArray[i] = array[i];
    }

  for (i = n - 1; i > 0; i--) {
    j = rand_int(i + 1);
    tmp = array[j];
    randomArray[j] = randomArray[i];
    randomArray[i] = tmp;
  }
}


/** 
 * Creates random permutations of aminoacids for each sequence and scores them
 * 
 * @param e 
 * @param d 
 * @param index 
 * @param outputfile 
 * @param fn 
 */
void scoreRandomData_i(env *e, data_struct *d, int index, double **randomScores, double *mean, double *sd, int nrand, int *nsums)
{
  
  int i, j, r, n = d->nUsedAtoms[index], ind0, ind, ind_tmp, aaInd1, aaInd2, maxInd = 0;
  
  double score = 0.;
  
  *mean = *sd = 0;
  
  *nsums = 0;
  
  if(n == 0)
    return;
  
  int *intSequence, *randomSeq;
  
  intSequence = calloc( n, sizeof(int));
  
  randomSeq = calloc( n, sizeof(int));

  for(i = 0; i <n; ++i)
    {
      ind0 = d->indexUsedAtoms[index][i]; 
      
      if(ind0 > maxInd )
	maxInd = ind0;
    }
  
  int indArray[maxInd + 1];
  
  for( i = 0; i < n; ++i)
    {
      
      /* map i to ind0. i refers to the atoms have been used after
	 membrane filtering. However, the name of the aa has been
	 stored in data_struct BEFORE filtering. Thus, we need to map
	 the coordinates */

      ind0 = d->indexUsedAtoms[index][i]; 
      
      intSequence[i] = convertAAToInt( d -> atoms[index][ind0].aa );
      
      
      /* It would be nice if C would have a hash function.
      	 Since there is no one, I use the indArray to keep the indexes
      */
      indArray[ ind0 ] = i;
    }
  
  
  *nsums = 0;
  
  for(r = 0; r < nrand; ++r)
    {
      
      shuffle(intSequence, randomSeq, n);
      
      for(i = 0; i < n; ++i)
	{
	  
	  ind0 = d->indexUsedAtoms[index][i];
	  
	  aaInd1 = randomSeq[i] + ( d->atoms[index][ind0].env - 1) * 20 ;
	  
	  score = 0.;

	  for(j = 0; j < d->dis[index].dim2[i]; ++j)
	    {
	      
	      ind_tmp = d->dis[index].pDistanceMatrix[i][j];
	      
	      ind = d->indexUsedAtoms[index][ind_tmp];

	      aaInd2 = randomSeq[ indArray[ ind ] ] + (d->atoms[index][ind].env - 1) * 20 ;
	      	      
	      score += e->matrix[aaInd1][aaInd2];
	      
	      if( r == 0)
		(*nsums)++;
	    }

	  randomScores[index][ r ] += score;
	  
	}
      
      randomScores[index][ r ] /= (double)(*nsums);
      
      *mean += randomScores[index][r];
      
      *sd  += randomScores[index][r] * randomScores[index][r];
      
    }

  *mean /= nrand;

  *sd /= nrand;

  *sd -= (*mean) * (*mean);

  *sd = sqrt( (*sd) );
}




/** 
 * Calculates the distance between two coordinates in the 3d space. Coordinates are attributes of the atoms
 * 
 * @param m first atom
 * @param n second atom
 * 
 * @return euclidean distance between the two atoms
 */
static double getDistance(atom m, atom n) 
{
  /* euclidean distance */
  double distance = ( (m.x - n.x) * (m.x - n.x) ) + 
    ( (m.y - n.y) * (m.y - n.y) ) +
    ( (m.z - n.z) * (m.z - n.z) );
  return  ( sqrt ( distance ) );
  
}


void freeGlobalMemory()
{
  free(allFiles);
}


/* void freeData(int n, data_struct *d) */
/* { */
/*   int i; */

/*   free(d->atoms); */
  
/* } */

/** 
 * Calculate distances of all atoms and stores them in a distance matrix
 * 
 * @param dist structure that holds the distance matrix and a matrix of pointers to it
 */
void calculateDistances(dist_struct* dist, env *e, double threshold)
{
  /* n = dist->size = e->natoms */
  int i=0, j = 0, n = dist->size, sz = 0, cnt = 0;
  
  double d = 0.;
  
  for(i = 0; i < n ; ++i)
    {
      cnt = 0;
      sz = 0;
      for ( j = 0; j < n; ++j )
	{
	  d = getDistance( *(e->atoms[i]), *(e->atoms[j]) );
	  
	  //dist->distanceMatrix[j][i] = dist->distanceMatrix[i][j];

	  

	  if( d < threshold &&
	      ( ( /* same chain */ (*e->atoms[i]).chain == (*e->atoms[j]).chain &&
	  	  /* not neigbhors */ abs( (*e->atoms[i]).pos - (*e->atoms[j]).pos ) > 1 ) || (/* different chains */(*e->atoms[i]).chain != (*e->atoms[j]).chain) ) )
	    ///abs(i-j) > 1 )
	    {
	      cnt++;
	      
	      if(sz == 0)
		{
		  sz = 1;
		  dist->distanceMatrix[i] = calloc( sz, sizeof(double) );
		  dist->pDistanceMatrix[i] = calloc( sz, sizeof(int) );
		  
		}
	      else if(sz < cnt)
		{
		  sz = cnt;
		  dist->distanceMatrix[i] = realloc( dist->distanceMatrix[i], cnt*sizeof(double) );
		  dist->pDistanceMatrix[i] = realloc( dist->pDistanceMatrix[i], cnt * sizeof(int) );
		}
	      
	      dist->distanceMatrix[i][cnt-1] = d;
	      dist->pDistanceMatrix[i][cnt-1] = j;
	      
	    }
	  dist->dim2[i] = cnt;
	  	  
	}
    }
  
}



/* /\**  */
/*  * The compare function used in the qsort algorithm */
/*  *  */
/*  * @param a  */
/*  * @param b  */
/*  *  */
/*  * @return 1 if the first elements is greater than the second, -1 if the second is greater than the first, or 0 otherwise */
/*  *\/ */
/* static int compare (const void *a, const void *b)  */
/* { */

/*   const double xx = **(double**)a;  */
/*   const double yy = **(double**)b; */
  
/*   if( xx > yy ) */
/*     return 1; */
  
/*   if( xx < yy ) */
/*     return -1; */
  
/*   return 0; */
/* } */


static inline void double_qsort( double **arr, unsigned int n)
{
#define double_lt(a,b)((**a) < (**b) )
  QSORT(double*, arr, n, double_lt);
}


/* /\**  */
/*  * q-sorts the pointers to the distance matrix */
/*  *  */
/*  * @param dist structure that holds the distance matrix and a matrix of pointers to the distance matrix */
/*  *\/ */
/* void sortIndexes ( dist_struct *dist) */
/* { */
/*   int i = 0, n = dist->size; */
/*   for(i=0; i< n; ++i ) */
/*     { */
/*       double_qsort(dist->pDistanceMatrix[i], n); */
/*       //qsort(dist->pDistanceMatrix[i], n, sizeof(dist->pDistanceMatrix[i][0]), compare); */
/*     } */
  
/* } */



/* /\**  */
/*  * Prints the distance matrix (sorted) */
/*  *  */
/*  *\/ */
/* void printResults(){ */
  
/*   int i = 0, j = 0, ind = 0, n = dist->size; */

/*   for(i=0; i<n; ++i) */
/*     { */
/*       for( j = 0; j < n;  ++j ) */
/* 	{ */
/* 	  ind = dist->pDistanceMatrix[i][j] - dist->distanceMatrix[i]; */
/* 	  printf("%s:%f ", atoms[ind].aa, ( *dist->pDistanceMatrix[i][j]) ); */
/* 	} */
/*       printf("\n\n"); */
/*     } */
/* } */


/** 
 * element-wise copy of the elements of a matrix to another matrix
 * Matrices should be square-matrices of the same size

 * @param source the source matrix
 * @param destination destination matrix
 * @param size dimension of the matrix
 */
void copyMatrix(double **source, double **destination, int size)
{
  int i,j;
  for(i = 0; i <size; ++i)
    for(j = 0; j < size; ++j)
      destination[i][j] = source[i][j];
}


/** 
 * It generates the M matrix (aminoacid 'interactions') for a specific environemnt. The matrix dimensions depend on the environment. For example if Uniformity is assumed, i.e. all atoms belong the same environment, then matrix dimensions are 20x20. If atoms are classified according to the number of their neighbors then according to Jha et al. 2011, we have three environments and the matrix dimensions are 60x60
 * 
 * @param dist structure that holds the distance matrix and a matrix of pointers to the distance matrix
 * @param e structure that holds the environment. Contains the aa interaction matrix, M
 * @return 1 if the size of the matrix is greater than 0, and 0 if the size of the matrix is 0
 */
int generateDistanceMatrix(dist_struct *dist, env *e, double threshold)
{
  int n = e->natoms; 
  
  // if there is no valid input then skip everything
  if( n == 0 )
    {
      dist->size = 0;
      return 0;
    }
    
  dist->size = n;
  
  dist->distanceMatrix = malloc( n * sizeof(double* ));
  
  dist->pDistanceMatrix = malloc( n * sizeof(int * ) );

  dist->dim2 = calloc( n, sizeof(int) );

  
  /*********************/
  /* Analyze the set of atoms in an environment */
  /*********************/
    
  calculateDistances( dist, e, threshold);

  //sortIndexes( dist );
  
  return 1;
}



/* /\**  */
/*  * Returns the index of the second dimension of the distance matrix, whose element is largest less than the threshold. Elements at each row of the matrix have to be sorted */
/*  *  */
/*  * @param dist  structure that holds the distance matrix */
/*  * @param index the index of the row */
/*  * @param threshold threshold value. Elements should not be larger than this threshold */
/*  * @param length the second dimension of the matrix */
/*  *  */
/*  * @return value of the index whose element is the greatest but lower than the threshold */
/*  *\/ */
/* int getIndexLess(dist_struct *dist, int index, double threshold, int length){ */
/*   int i = 0; */
    
/*   while(i < length && (*dist->pDistanceMatrix[index][i]) <= threshold ) */
/*     { */
/*       i++; */
/*     } */

/*   return i; */
/* } */





int totalEnvironmentMatrix(env *e)
{
  int i,j, env1, env2, aai1, aai2, g;
  for(i = 0; i < e->dim1; ++i)
    {
      env1 = (int)(i/20) + 1;
      aai1 = i - 20 * (env1 - 1);
      for(j = 0; j < e->dim2; ++j)
	{
	  env2 = (int)(j/20) + 1;
	  aai2 = j - 20 * ( env2 - 1 );
	  g = (i == j)?1:2;
	  
	  if( e->matrix[i][j]  < 1e-6 )
	    e->matrix[i][j] = 1e-6;
	  if( taa->aaCounts[aai2] == 0)
	    {
	      taa->aaCounts[aai2] = 1;
	      taa->total++;
	    }
	  if( taa->aaCounts[aai1] == 0)
	    {
	      taa->aaCounts[aai1] = 1;
	      taa->total++;
	    }
	  if(e->envContacts[env1 - 1][env2 -1]  == 0)
	    {
	      e->envContacts[env1 - 1][env2 -1] = 1;
	    }

	  
	  /* if(i == 0 && j == 0) */
	  /*   { */
	  /*     printf("%f %d %d %d %d %d\n", e->matrix[i][j], g, taa->aaCounts[aai1], taa->total, taa->aaCounts[aai2], e->envContacts[env2 - 1][env1-1]); */
	  /*     exit(-1); */
	  /*   } */

	  
	  e->matrix[i][j] = -log( e->matrix[i][j] / ( g * (1.*taa->aaCounts[aai1]/taa->total) * (1.* taa->aaCounts[aai2]/taa->total) * e->envContacts[env2 - 1][env1-1]  ) );

	}
    }
  /* for(i = 0; i < e->dim1; ++i) */
  /*   { */
  /*     for(j = 0; j < e->dim2; ++j) */
  /* 	printf("%f ", e->matrix[i][j]); */
  /*     printf("\n"); */
  /*   } */
  return 1;
  
}


void initializeEnvDegree(env *e)
{
  int i;
  
  e->natomsEnv = calloc( 3, sizeof(int) );
  e->totalAtoms = 0;
  e->envs = 3;
  
  e->envContacts = calloc(e->envs, sizeof(int*));
  for(i = 0; i < e->envs; ++i)
    e->envContacts[i] = calloc(e->envs, sizeof(int));
  
}

				
static int classifyToEnvironmentDegree(atom *a, int n)
{
  if(a->neighbors < n)
    return 1;
  if(a->neighbors > n)
    return 3;
  return 2;
}


void classifyAtoms(env *e, dist_struct *dist)
{
  int i = 0, n = dist->size, ind0;

  
  if(e->natoms != dist->size)
    {
      fprintf(stderr, "ERROR: e->natoms (%d) != dist->size (%d)\n", e->natoms, dist->size);
      assert(e->natoms == dist->size);
    }
  
  
  for(i = 0; i < n; ++i)
    {
      ind0 = i;
      
      (*e->atoms[ind0]).neighbors = dist->dim2[ind0];
      
      (*e->atoms[ind0]).env = classifyToEnvironmentDegree( e->atoms[ind0], 6);
      
      ++e->natomsEnv[ (*e->atoms[ind0]).env -1 ];
      
      ++e->totalAtoms;
    }
  
}	


int incrementEnvironmentMatrix(env* e, dist_struct *dist)
{
  
  int i = 0, j = 0, n = dist->size, aaInd1 = 0, aaInd2 = 0, ind, ind0;
  
  for(i = 0; i < n; ++i)
    {
      ind0 = i;
      
      aaInd1 = ( convertAAToInt( (*e->atoms[ind0]).aa) + ( (*e->atoms[ind0]).env - 1) * 20 );
      for(j = 0; j < dist->dim2[i]; ++j)
	{
	  ind = dist->pDistanceMatrix[i][j];
	  	      
	  ++(e->envContacts[ (*e->atoms[ind0]).env - 1][ (*e->atoms[ind]).env -1]);
	  e->envContacts[ (*e->atoms[ind]).env - 1][ (*e->atoms[ind0]).env -1]
	    = e->envContacts[ (*e->atoms[ind0]).env - 1][ (*e->atoms[ind]).env -1];
	      
	  aaInd2 = ( convertAAToInt( (*e->atoms[ind]).aa) + ( (*e->atoms[ind]).env - 1) * 20 );
	  
	  e->matrix[aaInd1][aaInd2] = e->matrix[aaInd1][aaInd2] + 1.;
	  e->matrix[aaInd2][aaInd1] = e->matrix[aaInd1][aaInd2];
	  
	}
    }
  
   /* for(i = 0; i < e->envs; ++i) */
   /*   { */
   /*     printf("Environment %d has %d elements out of %d\n", i, e->natomsEnv[i], e->totalAtoms); */
   /*   } */

 return 1;
 }


void readEnvironmentMatrixFromFile(FILE *matrixInputFile, env *e)
{
  int header = 1, x = e->dim1, y = e->dim2, i, j;
  char aa[4];
  
  if(header == 1)
    {
      for(j = 0; j < y; ++j)
	{
	  if( fscanf(matrixInputFile, "%s", aa) != 1)
	    {
	      fprintf(stderr, "Error in reading the matrix file\n");
	      exit(-1);
	    }
	}
      header = 0;
    }
    
  
  for(i = 0; i < x; ++i)
    {
      
      if(fscanf(matrixInputFile, "%s", aa) != 1)
	{
	  fprintf(stderr, "Error in reading line %d from the matrix\n", i);
	  exit(-1);
	}
      
      for(j = 0; j < y; ++j)
	{
	  if( fscanf(matrixInputFile, "%lf", &e->matrix[i][j]) != 1)
	    {
	      fprintf(stderr, "Error in reading line %d %d from the matrix\n", i, j);
	      exit(-1);
	    }
	}
    }
  
}



/* int incrementEnvironmentMatrixUniform(env* e, dist_struct *dist, int envMode, double threshold) */
/*  { */
/*    int i = 0, j = 0, n = dist->size, aaInd1 = 0, aaInd2 = 0, ind, ind0; */

/*    if(e->natoms != dist->size) */
/*      { */
/*        fprintf(stderr, "ERROR: e->natoms (%d) != dist->size (%d)\n", e->natoms, dist->size); */
/*        exit(-1); */
/*      } */

/*    if(envMode == 1) //Uniform */
/*      { */
/*        for(i = 0; i <n; ++i) */
/* 	 { */
/* 	   ind0 = dist->pDistanceMatrix[i][0] - dist->distanceMatrix[i]; */
/* 	   aaInd1 = convertAAToInt(e->atoms[ind0].aa); */
/* 	   e->atoms[ind0].env = 0; */
/* 	   for(j = 0; j < n; ++j) */
/* 	     { */


/* 	       if( *(dist->pDistanceMatrix[i][j]) < threshold ) */
/* 		 { */


/* 		   if( abs(e->atoms[ind].index - e->atoms[ind0].index) < 2 ) */
/* 		     continue; */


/* 		   ind = dist->pDistanceMatrix[i][j] - dist->distanceMatrix[i]; */
/* 		   aaInd2 = convertAAToInt( e->atoms[ind].aa ); */
/* 		   e->matrix[aaInd1][aaInd2] = e->matrix[aaInd1][aaInd2] + 1.; */

/* 		   printf("Distance between atom %s (%d) and %s (%d) is %f\n",  */
/* 			  e->atoms[ind0].aa, e->atoms[ind0].index, e->atoms[ind].aa, e->atoms[ind].index, *(dist->pDistanceMatrix[i][j]) ); */
/* 		 } */
/* 	       else */
/* 		 break; */

/* 	     } */
/* 	 } */
/*      } */

/*    for(i = 0; i < e->dim1; ++i) */
/*      { */
/*        for(j = 0; j < e->dim2; ++j) */
/* 	 printf("%f ", e->matrix[i][j]); */
/*        printf("\n"); */
/*      } */

/*    return 1; */
/*  } */



#ifdef _USE_PTHREADS



void sortIndexes_thread (int tid, int threads, dist_struct *dist)
{
  int i = 0;

  for(i = 0; i < dist->size; ++i)
    {
      if(i % threads == tid)
	{
	  qsort(dist->pDistanceMatrix[i], dist->size, sizeof(dist->pDistanceMatrix[i][0]), compare);
	}
    }
}

void sortIndexes_parallel(dist_struct *dist)
{

  int i, threads = threadDataL[0].threadTOTAL;
  
  for(i = 1; i<threads; ++i)
    {
      threadDataL[i].threadOPERATION = SORT;
    }

  sortIndexes_thread(0, threads, dist);

  threadDataL[0].threadBARRIER = 1;

  syncThreadsBARRIER();
}



void calcDistances_thread( int tid, int threads, dist_struct *dist)
{
  int i=0, j = 0, n = dist->size;

  for(i = 0; i<n-1; ++i)
    {
      for(j = i+1; j < n; ++j )
	{
	  if( (i+j) % threads == tid )
	    {
	      dist->distanceMatrix[i][j] = getDistance(atoms[i], atoms[j]);
	      dist->distanceMatrix[j][i] = dist->distanceMatrix[i][j];

	      dist->pDistanceMatrix[i][j] = &dist->distanceMatrix[i][j];
	      dist->pDistanceMatrix[j][i] = &dist->distanceMatrix[j][i];
	    }
	}
      dist->distanceMatrix[i][i] = 0.;
      dist->pDistanceMatrix[i][i] = &dist->distanceMatrix[i][i];
    }
  dist->distanceMatrix[n-1][n-1] = 0.;
  dist->pDistanceMatrix[n-1][n-1] = &dist->distanceMatrix[n-1][n-1];
}

void calcDistances_parallel( dist_struct  *dist)
{
  int i, threads = threadDataL[0].threadTOTAL;

  for(i = 1; i < threads; ++i)
    {
      threadDataL[i].threadOPERATION = DIST;
    }

  calcDistances_thread(0, threads, dist);

  threadDataL[0].threadBARRIER = 1;

  syncThreadsBARRIER();
}


#endif
