/*! \file 
*\author Pavlos Pavlidis
*/ 

#include "prins.h"

/** 
 * Checks if a command line argument is correct. In detail, it checks whether the flag is the last argument (i.e. value is missing) or if the next argument does not start with a '-'
 * 
 * @param i index of the command line argument
 * @param argc 
 * @param argv 
 * 
 * @return 
 */
static int check_arg(int i, int argc, char **argv)
{
  if( i == argc - 1)
    {
      fprintf(stderr, " Error in argument %s ( %s <STRING> )\n", argv[i], argv[i]);
      return 0;
    }
  
  if( argv[i+1][0] == '-' )
    {
      fprintf(stderr, " Error in argument %s ( Next argument starts with a '-' ) \n", argv[i]);
      return 0;
    }
  
  return 1;
}

/** 
 * \brief Prints out the header of PRInS
 * 
 */
void header()
{
  fprintf(stderr, "\n\nProtein Residuals INteraction Statistics\n---- PRInS ----\n\n");
  fprintf(stderr, "Current Version: 1.0\n");
  fprintf(stderr, "Release Data: May 2014\n");
  fprintf(stderr, "Author: Pavlos Pavlidis\n");
  fprintf(stderr, "Address: FORTH-IMBB, Crete, Greece\n");
  fprintf(stderr, "Web: pop-gen.eu\n");
  fprintf(stderr, "-----------------------------------------\n\n\n");	  
  fprintf(stderr, "For help with the command line, type prins -h\n\n");
}

/** 
 * \brief Prints out the help for PRInS
 * 
 */
void usage()
{
  fprintf(stderr, "\n\nprins -h | --help: Prints out this help\n\n");

  fprintf(stderr, "\t-iMatrix <FileName>: imports a matrix that describes\nstatistical interactions\n\n");

  fprintf(stderr, "\t-oMatrix <FileName>: outputs the statistics for the\ninteractions between residuals in a file, in a matrix form\n\n");

  fprintf(stderr, "\t-pdb <FileName>: reads in the pdb file\n\n");

  fprintf(stderr, "\t-o <FileName>: outputs result to a file\n\n");

  fprintf(stderr, "\t-list <FileName>: processes a list of files. Their filenames\nare written in a list, one name in each line\n\n");

  fprintf(stderr, "\t-threshold <float>: distance in Angstrom to consider two\nneighboring Ca atoms as interacting (default: 6.5)\n\n");

  fprintf(stderr, "\t-membrane: analyze atoms only within the membrane\nboundaries\n\n");

  fprintf(stderr, "\t-rfile <file name>: print out statistics from permutations of the aminoacids in protein structures\n\n");

  fprintf(stderr, "\t-nrand <int>: perform <int> randomizations (permutations) of the aminoacid chain (default: 0, don't do it)\n\n");

}


/** 
 * 
 * 
 * @param argc number of arguments in the command line
 * @param argv command line arguments
 * @param infile input file. It is a .pdb file
 * @param ofile  output file. Where the results will be printed
 * @param listfile  input file. Actually it is a file that contains a list of file names (each is a .pdb file)
 * @param threshold maximum distance between two C_a atoms in order to be considered as neighbors
 */
void commandLineParser(int argc, char **argv,
		       char *infile,
		       char *ofile,
		       char *listfile,
		       double *threshold,
		       char *matrixFileName, 
		       int *input,
		       int *output, 
		       int *nrand, 
		       char *randomOutputFileName)
{

  int i;
  membraneLim = 0;
  
  *input = 0;
  
  *output = 1;

  if(argc < 2 )
    {
      header();
      exit(-1);
    }

  for(i=1; i<argc; ++i)
    {

      if(strcmp( argv[i], "--help" ) == 0 || strcmp( argv[i], "-h") == 0)
	{
	  usage();
	  exit(-1);
	}

      

      if(!strcmp( argv[i], "-iMatrix") || !strcmp( argv[i], "-imatrix"))
	{
	  if( check_arg(i, argc, argv ) == 0 )
	    exit(1);
	  
	  *input = 1;

	  *output = 0;

	  strcpy(matrixFileName, argv[++i]);
	  continue;
	}
      
      
      if(!strcmp( argv[i], "-oMatrix") || !strcmp( argv[i], "-omatrix") ) 
	{
	  if ( check_arg( i, argc, argv) == 0 )
	    exit(1);
	 
	  	  
	  if(*input == 1)
	    {
	      fprintf(stderr, "\nIt is not possible to input and output the matrix\n");
	      exit(-1);
	    }

	  strcpy(matrixFileName, argv[++i]);
	  continue;
	}

      
      if( !strcmp( argv[i], "-pdb") )
	{
	  
	  if ( check_arg( i, argc, argv) == 0 )
	    exit(1);

	  strcpy(infile, argv[++i]);
	  continue;
	}

      if( !strcmp( argv[i], "-o") )
	{
	  if( check_arg(i, argc, argv ) == 0 )
	    exit(1);

	  strcpy(ofile, argv[++i]);
	  continue;
	}

      if( !strcmp( argv[i], "-list") )
	{
	  if( check_arg(i, argc, argv ) == 0 )
	    exit(1);

	  strcpy(listfile, argv[++i]);
	  continue;
	}

      if( !strcmp( argv[i], "-threshold") )
	{
	  if( check_arg( i, argc, argv ) == 0 )
	    exit(1);
	  
	  *threshold = atof( argv[++i] );
	  continue;
	}
      if( !strcmp( argv[i], "--membrane") )
	{
	  membraneLim = 1;
	  continue;
	}
      
      if( !strcmp( argv[i], "-nrand") )
	{
	  *nrand = atoi( argv[++i]);
	  continue;
	}

      if( !strcmp( argv[i], "-rfile") )
	{
	  strcpy( randomOutputFileName, argv[++i]);
	  continue;
	}
      
      usage();

      fprintf(stderr, "\n\nArgument %s is not supported\n", argv[i]);
      exit(-1);
    }
}



/** 
 * \brief Matches a string in a text (exact match)
 * 
 * @param str the text
 * @param flag the string to be matched on the text
 * @param flagLength length of the string
 * @param strLength length of the text
 * 
 * @return 1: match, 0: non-match
 */
static int flagMatchAnywhere(char * str, char *flag, int flagLength, int strLength)
{
  int i, match = 0, j = 0;

  
  if(flagLength > strLength)
    return 0;

  while(j < strLength || str[j] == '\0' )
    {
      for(i = 0; i < flagLength; ++i)
	{
	  if(str[j] != flag[i])
	    {
	      match = 0;
	      ++j;
	      break;
	    }
	  else
	    {
	      ++j;
	    }
	  if( (j < strLength && str[j] == 10) || j == strLength )
	    break;
	}
      if( (j < strLength && str[j] == 10) || j == strLength )
	break;

      if(i == flagLength)
	{
	  match = 1;
	  break;
	}
    }
  
  return match;
}


/** 
 * \brief Matches a string to a text at the beginning of the text
 * 
 * 
 * @param str the text
 * @param flag the string to be matched on the text
 * @param flagLength length of the string
 * @param strLength length of the text
 * 
 * @return 1: match, 0: non-match
 */
static int flagMatch(char * str, char *flag, int flagLength, int strLength)
{
  int i;
  
  if(flagLength > strLength)
    return 0;

  for(i = 0; i < flagLength; ++i)
    {
      
      if(str[i] != flag[i])
	return 0;
    }
  
  return 1;
}

/** 
 * \brief matches a substring to given positions (start, end) in a text
 * 
 * @param str the text
 * @param flag the string to be matched
 * @param returnFlag the matched string
 * @param start coordinate of the text to start matching
 * @param end coordinate of the text to end matching
 * @param strLength total length of the text
 * 
 * @return NULL if no match. The matched String otherwise.
 */
static char *substringMatch( char *str, char *flag, char * returnFlag, int start, int end, int strLength)
{
  assert( start >= 0 && end < strLength);
  
  int j = 0, i = 0;
  
  for(i=start; i< end; ++i)
    {
      if( str[i] != flag[j])
	break;
      returnFlag[j] = str[i];
      
      j++;
    }

  if(j < end - start )
    {
      returnFlag = NULL;
      return NULL;
    }

  return returnFlag;
  
}


/** 
 * \brief assigns a substring of a text to a char array
 * 
 * Substring ends with a \0
 * @param str the text
 * @param returnFlag substring (returned by reference)
 * @param start matching starts at this position
 * @param end matching ends at this position
 * @param strLength the length of the text
 */
static void substringAssignment( char *str, char * returnFlag, int start, int end, int strLength)
{
  assert( start >= 0 && end < strLength);
  
  int j = 0, i = 0;
  
  for(i=start; i< end; ++i)
    {
      
      returnFlag[j] = str[i];
      
      j++;
    }
  returnFlag[j] = '\0';
  
}

/** 
 * Skips blank lines
 * \r is skipped
 * \n is considered the line feeder

 * @param ln line to be checked
 * @param length the length of the line
 * 
 * @return 1 if the line consists only of space and tabls and ENTER as a new line. 
 */
static int skipBlankLines(char *ln, int length)
{
  int i = 0;
  
  char *p = ln;
  
  while( *p == 32 || *p == 9 )
    {
      i++;
      p++;
      
      if( i >= length )
	{
	  fprintf(stderr, "Error reading file line %s \n", ln);
	  exit(0);
	}
    }

  if ( *p == '\r' ){
    p++;
    i++;
    if( i >= length )
      {
	fprintf(stderr, "Error reading file line %s \n", ln);
	exit(0);
      }
  }

  
  if( *p == '\n' )
    return 1;

  return 0;
}

/** 
 * Replaces \n with \0
 * 
 * @param ln the line that replacement will take place
 * @param length length of the line
 * 
 * @return 1 if replacements occurs. 0 otherwise
 */
static int chopLn( char *ln, int length )
{
  int i = 0;

  while( i < length)
    {
      
      if( *ln == '\n')
	{
	  *ln = '\0';
	  return 1;
	}
      
      ln++;
      i++;
    }

  return 0;
}


int parseListFile(char *listfile)
{
  int numberOfFiles = 0, y = 1, i = 0, j = 0;
  
  char ln[MAXFILE];
  
  FILE *fp = fopen(listfile, "r");

  if( fp == NULL)
    return 0;
  
  allFiles = calloc ( y, sizeof(file_struct) );
  
  taa = malloc( sizeof(total_aa) );

  for(j = 0; j < MAXAA; ++j)
    taa->aaCounts[j] = 0;
  taa->total = 0;
  
  while( fgets( ln, MAXFILE-1, fp) != NULL )
    {
      
      if( skipBlankLines( ln, MAXFILE ) == 1)
	continue;
      
      if ( chopLn (ln, MAXFILE) == 0 )
	assert( 0 );
      
      if( i == y )
	{
	  ++y;
	  allFiles = realloc( allFiles, y * sizeof( file_struct ) );
	}
      
      strcpy ( allFiles[i].name, ln );

      /* printf("File: %s\n", allFiles[i].name); */
      
      ++i;
    }
  numberOfFiles = i;

  return numberOfFiles;
}

int getFileNamesFromCMD(int argc, char **argv, char ***files, int currentX, int currentY)
{
  int i, j = 0, k = 0, x = 10, y = 1024;

  
  if(currentX != 0)
    x = currentX;

  if(currentY != 0)
    y = currentY;

  if(*files == NULL)
    {
  
      (*files) = calloc( x, sizeof(char*) );
      
      for( i = 0; i < x; ++i)
	{
	  (*files)[i] = calloc( y, sizeof(char) );
	}
    }
  
  for( i = 0; i < argc; ++i)
    {

      if(j == x)
	{
	  *files = realloc( *files, (x + 5)*sizeof(char*));
	  
	  for( k = x; k < x + 5; ++k)
	    (*files)[k] = calloc( y, sizeof(char));

	  x += 5;
	  
	}

      if( strcmp(argv[i], "-pdb") == 0)
	{
	  (*files)[j++] = argv[++i];
	}
      
    }
  
  return j;
}


int getFileNames(char *listfile, char ***files, int currentX, int currentY)
{
  
  int i = 0,j = 0,x = 10, y = 1024;

  char ln[MAXLN];

  if(currentX != 0)
    x = currentX;
  
  if(currentY != 0)
    y = currentY;
  
  if( *files == NULL)
    {
      *files = calloc(x, sizeof(char*));
      
      for( i = 0; i < x; ++i)
	(*files)[i] = calloc( y, sizeof(char));
    }

  FILE *fp = fopen(listfile, "r");
  
  if( fp == NULL)
    return 0;

  while( fgets( ln, MAXLN-1, fp) != NULL)
    {
      if( skipBlankLines ( ln, MAXFILE ) == 1)
	continue;
      
      if(chopLn( ln, MAXFILE) == 0 )
	assert(0);

      if( i == x)
	{
	  
	  *files = realloc( *files, x+5 * sizeof( char* ));

	  for(j = x; j < x + 5; ++j)
	    {
	      (*files)[j] = calloc( y, sizeof(char) );
	    }

	  x+=5;
	  
	}
      
      strcpy( (*files)[i], ln);
      
      ++i;
    }
  
  return i;
}


static void getMembraneZCoordinate(char *str, extraFeatures *ef)
{
  char record[10];

  substringAssignment(str, record, 46, 54, MAXLN);
  if(ef->counter == 0)
    {
      ef->z1 = atof(record);
      ef->counter = 1;
    }
  else if(ef->counter == 1)
    {
      ef->z2 = atof(record);
      ef->counter = 2;
    }
  else
    {
      fprintf(stderr, "ef->counter should be 0|1. Now it is %d\n", ef->counter);
      exit(-1);
    }
}




/* int parseInputFile(char *infile) */
/* { */
  
/*   int y = 10, i=0, flagLength = 4, countLine = 0; */
/*   //    x1start = 30, x1end = 37, x2start = 38, x2end = 45, x3start = 46, x3end = 53; */
/*   char ln[MAXLN], flag[] = "ATOM", record[10]; // *pch */
  

/*   FILE *fp = fopen(infile, "r"); */

/*   if(fp == NULL) */
/*     return 0; */

/*   atoms = malloc( y * sizeof(atom) ); */

  
    

/*   /\* first read all lines and... */
/*      1. store the information in the atom structure */
/*      2. check if there is membrane information */
/*   *\/ */
/*   while( fgets( ln, MAXLN-1, fp) != NULL ) */
/*     { */
      
/*       if( ef->counter < 2 && */
/* 	  (flagMatch(ln, "HETATM", 6, MAXLN) == 1) &&  */
/* 	  (flagMatchAnywhere(ln, "DUM", 3, MAXLN) == 1 ) */
/* 	  ) */
/* 	{ */
/* 	  /\* */
/* 	    Membrane has been found  */
/* 	  *\/ */
/* 	  //getMembraneZCoordinate(ln, ef); */
/* 	  printf("counter: %d, z1: %f, z2: %f\n", ef->counter, ef->z1, ef->z2); */
/* 	} */
      
/*       if( flagMatch(ln, flag, flagLength, MAXLN) == 0) */
/* 	continue; */
  
/*       /\* read the C atom information *\/ */
/*       if( substringMatch(ln, " C  ", record, 12, 16, MAXLN) != NULL ) */
/* 	{ */
/* 	  countLine++; */
	  
/* 	  if(countLine > y ) */
/* 	    { */
/* 	      y += 10; */
/* 	      atoms = realloc ( atoms, y * sizeof(atom) ); */
/* 	    } */
	  
/* 	  substringAssignment( ln, atoms[i].aa, 17, 20, MAXLN); */

/* 	  substringAssignment(ln, atoms[i].index, 6, 11, MAXLN); */
	  

/* 	  substringAssignment(ln, record, 22, 26, MAXLN); */
/* 	  atoms[i].pos = atoi(record); */
	  
/* 	  substringAssignment (ln, record, 30, 38, MAXLN); */
/* 	  atoms[i].x = atof(record); */
	  
/* 	  substringAssignment (ln, record, 38, 46, MAXLN); */
/* 	  atoms[i].y = atof(record); */
	  
/* 	  substringAssignment(ln, record, 46, 54, MAXLN); */
/* 	  atoms[i].z = atof(record);	   */

/* 	  atoms[i].counter = countLine-1; */

/* 	  taa->total++; */
/* 	  taa->aaCounts[ convertAAToInt(atoms[i].aa) ]++; */

/* 	  ++i; */
/* 	} */
/*     } */
/*   return countLine; */
/* } */


int parseInputFile2(char *infile, atom** a, int fileindex)
{
  
  int y = 1, i=0, flagLength = 4, countLine = 0;
  //    x1start = 30, x1end = 37, x2start = 38, x2end = 45, x3start = 46, x3end = 53;
  char ln[MAXLN], flag[] = "ATOM", record[10]; // *pch
  
  ef->counter = 0;

  FILE *fp = fopen(infile, "r");

  if(fp == NULL)
    return 0;

  *a = malloc( y * sizeof(atom) );

    /* first read all lines and...
     1. store the information in the atom structure
     2. check if there is membrane information
  */
  while( fgets( ln, MAXLN-1, fp) != NULL )
    {

      
      if(flagMatch(ln, "END", 3, MAXLN) == 1)
      	break;
            
      if( membraneLim == 1 && ef->counter < 2 &&
      	  (flagMatch(ln, "HETATM", 6, MAXLN) == 1) &&
      	  (flagMatchAnywhere(ln, "DUM", 3, MAXLN) == 1 )
      	  )
      	{
      	  /*
      	    Membrane has been found
      	  */
      	  getMembraneZCoordinate(ln, ef);

      	}
      
      if( flagMatch(ln, flag, flagLength, MAXLN) == 0)
	continue;

      if(flagMatchAnywhere(ln, "UNK", 3, MAXLN) == 1)
	continue;
  
      /* read the C atom information */
      if( substringMatch(ln, " CA ", record, 12, 16, MAXLN) != NULL )
	{
	  countLine++;
	  
	  if(countLine > y )
	    {
	      ++y;
	      (*a) = realloc ( (*a), y * sizeof(atom) );
	    }
	  
	  substringAssignment( ln, (*a)[i].aa, 17, 20, MAXLN);
	  (*a)[i].aa[3] = '\0';
	  
	  substringAssignment(ln, record, 21, 22, MAXLN);
	  (*a)[i].chain = record[0];

	  substringAssignment(ln, (*a)[i].index, 6, 11, MAXLN);
	  //(*a)[i].index = atoi(record);

	  substringAssignment(ln, record, 22, 26, MAXLN);
	  (*a)[i].pos = atoi(record);

	  
	  
	  substringAssignment (ln, record, 30, 38, MAXLN);
	  (*a)[i].x = atof(record);
	  
	  substringAssignment (ln, record, 38, 46, MAXLN);
	  (*a)[i].y = atof(record);
	  
	  substringAssignment(ln, record, 46, 54, MAXLN);
	  (*a)[i].z = atof(record);	  

	  (*a)[i].counter = countLine-1;
	  
	  (*a)[i].file = fileindex;

	  (*a)[i].score = 0.;

	  taa->total++;
	  
	  taa->aaCounts[ convertAAToInt((*a)[i].aa) ]++;

	  ++i;
	}
    }
  return countLine;
}


