/*
 * harness.c
 *
 *  Created on: 8 Jul 2010
 *      Author: Ben Smithers (bs8959 AT brist DOT ac DOT uk)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "sp_mwdc_harness.h"

int main(int argc, char *argv[]){

	//defaults
	char wildCardChar = '?';
	int repeats = 3;
	short printResults = 0;

	if(argc < 2){
		printf("Usage: filename [-flags for match to perform]+ [options] \n");
		printf("       At least one algorithm must be specified \n\n");
		printf("Flags: \n");
		printf("  -naive              Perform naive match in O(nm) time\n");
		printf("  -fft                Perform standard FFT in O(n log n) time using FFTW\n");
		printf("  -fft-fast           Use standard trick for O(n log m) time using FFWT\n");
		printf("Variations: both -fft and -fft may have any of the following suffixes for additional options (e.g. -fft-pad)\n\n");
		printf("  -measure            Use FFTW_MEASURE option. Not recommended - very slow \n");
		printf("  -pad                Pad transform sizes to a power of two \n");
		printf("  -wisdom             Pads to a power of 2, and uses system 'Wisdom'\n");
		printf("  -minsize            Force a minimum transform size of 2048(fft-fast only)'\n");
		printf("  -r2r                Use FFTW R2HC and HC2R transforms'\n");
		printf("  -NWIT               There are no wild cards in the text, be more efficient'\n");
		printf("Valid combinations: \n");
		printf("  -measure-pad\n");
		printf("  -NWIT-wisdom\n");
		printf("  -NWIT-minsize\n");
		printf("  -r2r-wisdom\n");
		printf("  -r2r-wisdom-minsize\n");
		printf("  -r2r-minsize\n");

		printf("Options: \n\n");
		printf("  -data file          Output a data file readable by gnuplot \n");
		printf("  -repeats            Number of repeats to perform. Default is 3\n");
		printf("  -wildcard           Set the wildcard character. Defaults to ?\n");
		printf("  -print              Print the matches found, rather than just the number of matches\n");
		exit(0);
	}

	FILE *file;

	file = fopen(argv[1],"r");
	if(file == NULL){
		printf("Could not open file %s for reading\n\n",argv[1]);
		exit(0);
	}

	//determine lengths for n and m
	int n, m;
	n = getLineLength(file);
	//read new line
	//fgetc(file);
	m = getLineLength(file);

	if(m > n || m <= 0 || n <= 0){
		printf("Invalid lengths for text or pattern. File must contain text\\npatterh\\n only\n");
		printf("n: %d; m: %d\n",n,m);
		exit(0);
	}

	char *text,*pattern;

	//Plus one to add space for \0
	text = (char *) malloc(sizeof(char) * (n+1));
	pattern = (char *) malloc(sizeof(char) * (m+1));

	//reset file pointer and read the text and pattern

	rewind(file);

	//first ensure we only read the correct amount

	char *returned;
	returned = fgets(text,n+1,file);//plus one because of \0
	if(returned==NULL){
		fprintf(stderr,"Unable to read text from the file");
		exit(1);
	}
	//read the new line
	fgetc(file);
	returned = fgets(pattern,m+1,file);
	if(returned==NULL){
		fprintf(stderr,"Unable to read pattern from the file");
		exit(1);
	}

	fclose(file);

	//determine tests to perform

	int i;
	FILE *dataFile = NULL;
	//Rough idea of the the number of different matches to perform. May be
	//A slight over-estimate with options, but it shouldn't matter.
	int *algorithms = (int *) calloc(argc -2,sizeof(int));
	int numAlgorithms = 0;//num algorithms being used

	if(argc > 2){
		for(i=2;i<argc;i++){
			/*Check for n log n FFT methods */
			if(strcmp("-fft",argv[i])==0){
				algorithms[numAlgorithms++] = FFT;
			}else if(strcmp("-fft-measure",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_MEASURE;
			}else if(strcmp("-fft-pad",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_PAD;
			}else if(strcmp("-fft-measure-pad",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_MEASURE_PAD;
			}else if(strcmp("-fft-wisdom",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_WISDOM;
			}else if(strcmp("-fft-r2r",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_R2R;
			}else if(strcmp("-fft-r2r-wisdom",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_R2R_WISDOM;
			}else if(strcmp("-fft-r2r-NWIT-wisdom",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_R2R_NWIT_WISDOM;
			}else if(strcmp("-fft-NWIT",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_NWIT;
			}else if(strcmp("-fft-NWIT-wisdom",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_NWIT_WISDOM;
			}else if(strcmp("-fft-fast",argv[i])==0){ /*Check for n log m FFT methods */
				algorithms[numAlgorithms++] = FFT_FAST;
			}else if(strcmp("-fft-fast-measure",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_MEASURE;
			}else if(strcmp("-fft-fast-pad",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_PAD;
			}else if(strcmp("-fft-fast-measure-pad",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_MEASURE_PAD;
			}else if(strcmp("-fft-fast-wisdom",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_WISDOM;
			}else if(strcmp("-fft-fast-minsize",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_MINSIZE;
			}else if(strcmp("-fft-fast-r2r",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_R2R;
			}else if(strcmp("-fft-fast-r2r-wisdom",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_R2R_WISDOM;
			}else if(strcmp("-fft-fast-r2r-minsize",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_R2R_MINSIZE;
			}else if(strcmp("-fft-fast-r2r-NWIT-minsize",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_R2R_NWIT_MINSIZE;
			}else if(strcmp("-fft-fast-NWIT-wisdom",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_NWIT_WISDOM;
			}else if(strcmp("-fft-fast-NWIT-minsize",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_FAST_NWIT_MINSIZE;
			}else if(strcmp("-naive",argv[i])==0){ /*Check for other methods */
				algorithms[numAlgorithms++] = NAIVE;
			}else if(strcmp("-naive-convolutions",argv[i])==0){ /*Check for other methods */
				algorithms[numAlgorithms++] = NAIVE_CONVOLUTIONS;
			}else if(strcmp("-flint",argv[i])==0){
				algorithms[numAlgorithms++] = FLINT;
			}else if(strcmp("-fft-random",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_RANDOM;
			}else if(strcmp("-fft-random-check",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_RANDOM_CHECK;
			}else if(strcmp("-fft-random-r2r",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_RANDOM_R2R;
			}else if(strcmp("-fft-random-r2r-check",argv[i])==0){
				algorithms[numAlgorithms++] = FFT_RANDOM_R2R_CHECK;
			}else if(strcmp("-kalai",argv[i])==0){
				algorithms[numAlgorithms++] = KALAI;
			}else if(strcmp("-kalai-check",argv[i])==0){
				algorithms[numAlgorithms++] = KALAI_CHECK;
			}else if(strcmp("-kalai-r2r",argv[i])==0){
				algorithms[numAlgorithms++] = KALAI_R2R;
			}else if(strcmp("-kalai-r2r-check",argv[i])==0){
				algorithms[numAlgorithms++] = KALAI_R2R_CHECK;
				/*Check for options */
			}else if(strcmp("-data",argv[i])==0){
				i++;
				if(i < argc){
					dataFile = fopen(argv[i],"a");
					if(file == NULL){
						printf("Could not open file %s for writing graph data to\n\n",argv[i]);
						exit(0);
					}
				}else{
					printf("No file specified\n\n");
					exit(0);
				}
			}else if(strcmp("-repeats",argv[i])==0){
				i++;
				if(i < argc){
					repeats = atoi(argv[i]);
				}else{
					printf("No repeat number specified\n\n");
					exit(0);
				}
			}else if(strcmp("-wildcard",argv[i])==0){
				i++;
				if(i < argc){
					wildCardChar = argv[i][0];
				}else{
					printf("No character given for wild card\n\n");
					exit(0);
				}
			}else if(strcmp("-print",argv[i])==0){
				printResults = 1;
			}else{
				printf("Warning: unknown flag: %s\n",argv[i]);
			}
		}
	}

	printf("+==============================================+\n");
	printf("|             Problem Information              |\n");
	printf("|             Text size   : %-10d         |\n",n);
	printf("|             Pattern size: %-10d         |\n",m);
	printf("+==============================================+\n");

	//print pattern size to datafile
	if(dataFile != NULL){
		fprintf(dataFile,"%d ",m);
	}

	for(i=0;i<numAlgorithms;i++){
		switch(algorithms[i]){
		case NAIVE:
			runTest(text,pattern,wildCardChar,n,m,NAIVE_START,NULL,repeats,dataFile,printResults,NAIVE_FP);
			break;
		case NAIVE_CONVOLUTIONS:
			runTest(text,pattern,wildCardChar,n,m,NAIVE_CONVOLUTIONS_START,NULL,repeats,dataFile,printResults,NAIVE_CONVOLUTIONS_FP);
			break;
		case FFT:
			runTest(text,pattern,wildCardChar,n,m,FFT_START,NULL,repeats,dataFile,printResults,FFT_FP);
			break;
		case FFT_MEASURE:
			runTest(text,pattern,wildCardChar,n,m,FFT_MEASURE_START,NULL,repeats,dataFile,printResults,FFT_MEASURE_FP);
			break;
		case FFT_PAD:
			runTest(text,pattern,wildCardChar,n,m,FFT_PAD_START,NULL,repeats,dataFile,printResults,FFT_PAD_FP);
			break;
		case FFT_MEASURE_PAD:
			runTest(text,pattern,wildCardChar,n,m,FFT_MEASURE_PAD_START,NULL,repeats,dataFile,printResults,FFT_MEASURE_PAD_FP);
			break;
		case FFT_WISDOM:
			runTest(text,pattern,wildCardChar,n,m,FFT_WISDOM_START,NULL,repeats,dataFile,printResults,FFT_WISDOM_FP);
			break;
		case FFT_R2R:
			runTest(text,pattern,wildCardChar,n,m,FFT_R2R_START,NULL,repeats,dataFile,printResults,FFT_R2R_FP);
			break;
		case FFT_R2R_WISDOM:
			runTest(text,pattern,wildCardChar,n,m,FFT_R2R_WISDOM_START,NULL,repeats,dataFile,printResults,FFT_R2R_WISDOM_FP);
			break;
		case FFT_R2R_NWIT_WISDOM:
			runTest(text,pattern,wildCardChar,n,m,FFT_R2R_NWIT_WISDOM_START,NULL,repeats,dataFile,printResults,FFT_R2R_NWIT_WISDOM_FP);
			break;
		case FFT_NWIT:
			runTest(text,pattern,wildCardChar,n,m,FFT_NWIT_START,NULL,repeats,dataFile,printResults,FFT_NWIT_FP);
			break;
		case FFT_NWIT_WISDOM:
			runTest(text,pattern,wildCardChar,n,m,FFT_NWIT_WISDOM_START,NULL,repeats,dataFile,printResults,FFT_NWIT_WISDOM_FP);
			break;
		case FFT_FAST:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_START,NULL,repeats,dataFile,printResults,FFT_FAST_FP);
			break;
		case FFT_FAST_MEASURE:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_MEASURE_START,NULL,repeats,dataFile,printResults,FFT_FAST_MEASURE_FP);
			break;
		case FFT_FAST_PAD:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_PAD_START,NULL,repeats,dataFile,printResults,FFT_FAST_PAD_FP);
			break;
		case FFT_FAST_MEASURE_PAD:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_MEASURE_PAD_START,NULL,repeats,dataFile,printResults,FFT_FAST_MEASURE_PAD_FP);
			break;
		case FFT_FAST_WISDOM:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_WISDOM_START,NULL,repeats,dataFile,printResults,FFT_FAST_WISDOM_FP);
			break;
		case FFT_FAST_MINSIZE:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_MINSIZE_START,NULL,repeats,dataFile,printResults,FFT_FAST_MINSIZE_FP);
			break;
		case FFT_FAST_NWIT_WISDOM:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_NWIT_WISDOM_START,NULL,repeats,dataFile,printResults,FFT_FAST_NWIT_WISDOM_FP);
			break;
		case FFT_FAST_NWIT_MINSIZE:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_NWIT_MINSIZE_START,NULL,repeats,dataFile,printResults,FFT_FAST_NWIT_MINSIZE_FP);
			break;
		case FFT_FAST_R2R:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_R2R_START,NULL,repeats,dataFile,printResults,FFT_FAST_R2R_FP);
			break;
		case FFT_FAST_R2R_WISDOM:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_R2R_WISDOM_START,NULL,repeats,dataFile,printResults,FFT_FAST_R2R_WISDOM_FP);
			break;
		case FFT_FAST_R2R_MINSIZE:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_R2R_MINSIZE_START,NULL,repeats,dataFile,printResults,FFT_FAST_R2R_MINSIZE_FP);
			break;
		case FFT_FAST_R2R_NWIT_MINSIZE:
			runTest(text,pattern,wildCardChar,n,m,FFT_FAST_R2R_NWIT_MINSIZE_START,NULL,repeats,dataFile,printResults,FFT_FAST_R2R_NWIT_MINSIZE_FP);
			break;
		case FLINT:
			runTest(text,pattern,wildCardChar,n,m,FLINT_START,NULL,repeats,dataFile,printResults,FLINT_FP);
			break;
		case FFT_RANDOM:
			runTest(text,pattern,wildCardChar,n,m,FFT_RANDOM_START,NULL,repeats,dataFile,printResults,FFT_RANDOM_FP);
			break;
		case FFT_RANDOM_CHECK:
			runTest(text,pattern,wildCardChar,n,m,FFT_RANDOM_CHECK_START,NULL,repeats,dataFile,printResults,FFT_RANDOM_CHECK_FP);
			break;
		case FFT_RANDOM_R2R:
			runTest(text,pattern,wildCardChar,n,m,FFT_RANDOM_R2R_START,NULL,repeats,dataFile,printResults,FFT_RANDOM_R2R_FP);
			break;
		case FFT_RANDOM_R2R_CHECK:
			runTest(text,pattern,wildCardChar,n,m,FFT_RANDOM_R2R_CHECK_START,NULL,repeats,dataFile,printResults,FFT_RANDOM_R2R_CHECK_FP);
			break;
		case KALAI:
			runTest(text,pattern,wildCardChar,n,m,KALAI_START,NULL,repeats,dataFile,printResults,KALAI_FP);
			break;
		case KALAI_CHECK:
			runTest(text,pattern,wildCardChar,n,m,KALAI_CHECK_START,NULL,repeats,dataFile,printResults,KALAI_CHECK_FP);
			break;
		case KALAI_R2R:
			runTest(text,pattern,wildCardChar,n,m,KALAI_R2R_START,NULL,repeats,dataFile,printResults,KALAI_R2R_FP);
			break;
		case KALAI_R2R_CHECK:
			runTest(text,pattern,wildCardChar,n,m,KALAI_R2R_CHECK_START,NULL,repeats,dataFile,printResults,KALAI_R2R_CHECK_FP);
			break;
		default:
			printf("Unknown algorithm!\n\n");
			exit(0);
		}
	}

	if(dataFile != NULL){
		fprintf(dataFile,"\n");
	}

	//free text and pattern
	free(text);
	free(pattern);
}

void runTest(char *text, char *pattern, char wildCardChar, int n, int m, char *startText, char *endText,int repeats,FILE *dataFile,short printResults,
		void (*functionPointer)(char *,char *,char,int,int,int *,struct SP_MWDC_MATCHING_POSITIONS*)){
	clock_t start;
	printf("%s\n\n",startText);
	//use the same text for the end as for the start if it's null
	if(endText == NULL){
		endText = startText;
	}
	struct SP_MWDC_MATCHING_POSITIONS *listOfMatches = NULL;
	int numMatches;
	if(!printResults){
		listOfMatches = sp_mwdc_create_new_list_of_matches();
	}
	int i;
	double totalTime = 0.0;
	for(i=0;i<repeats;i++){
		numMatches = 0;
		start = clock();
		functionPointer(text,pattern,wildCardChar,n,m,&numMatches,listOfMatches);
		double result = ((double)clock()-start)/(double)CLOCKS_PER_SEC;
		printf("Number of matches found: %d",numMatches);
		printf("(TIMER)Search #%d took %.2f\n",i,result);
		totalTime += result;
	}
	double avgTime = totalTime/(double)repeats;
	printf("\n(TIMER)Average Time For %s was %.2f\n\n",endText,avgTime);
	printf("------------------------------------------\n");
	if(dataFile != NULL){
		fprintf(dataFile,"%.2f ",avgTime);
	}
	if(!printResults){
		sp_mwdc_freeListOfMatches(listOfMatches);
	}

}


int getLineLength(FILE *file){
	int length;
	int c;
	length = 0;
	while((c=fgetc(file))!='\n'){
		if(c==EOF){
			break;
		}
		length++;
	}
	return length;
}
