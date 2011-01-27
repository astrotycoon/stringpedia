/*
 * harness.c
 *
 *  Created on: 8 Jul 2010
 *      Author: ben
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "sp_km_harness.h"

int main(int argc, char *argv[]){

	//Default
	int repeats = 3;
	int k = -1;
	short printResults = 0;

	if(argc < 3){
		printf("Usage: filename [-flags for match to perform]+ [options] \n");
		printf("       At least one algorithm must be specified \n\n");
		printf("Flags: \n");
		printf("  -naive              Perform naive match in O(nm) time\n");
		printf("  -abrahamson         Use the Abrahamson/Kosaraju method to find the hamming distance at all locations\n");
		printf("  -kosaraju           Use the Abrahamson/Kosaraju method to find the hamming distance at all locations\n");
		printf("Options: \n\n");
		printf("  -k number           The k in k-mismatches. Defaults to m\n");
		printf("  -data filename      Output a data file readable by gnuplot \n");
		printf("  -repeats number     Number of times to repeat each algorithm \n");
		printf("                      Defaults to 3\n");
		printf("  -print              Print the k-mismatches, rather than just the number of matches\n");
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
		printf("Invalid lengths for text or pattern. File must contain text\\npattern\\n only\n");
		printf("n: %d; m: %d\n",n,m);
		exit(0);
	}

	char *text,*pattern;

	//Plus one to add space for \0
	text = (char *) malloc(sizeof(char) * (n+1));
	pattern = (char *) malloc(sizeof(char) * (m+1));

	//reset file pointer and read the text and pattern

	rewind(file);

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

	for(i=2;i<argc;i++){
		if(strcmp("-naive",argv[i])==0){
			algorithms[numAlgorithms++] = NAIVE;
		}else if(strcmp("-unbounded",argv[i])==0){
			algorithms[numAlgorithms++] = UNBOUNDED;
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
		}else if(strcmp("-k",argv[i])==0){
			i++;
			if(i < argc){
				k = atoi(argv[i]);
			}else{
				printf("No number specified for k\n\n");
				exit(0);
			}
		}else if(strcmp("-print",argv[i])==0){
			printResults = 1;
		}else{
			printf("Warning: unknown flag: %s\n",argv[i]);
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
			runTest(text,pattern,n,m,k,NAIVE_START,NAIVE_END,repeats,dataFile,printResults,NAIVE_FP);
			break;
		case UNBOUNDED:
			runTest(text,pattern,n,m,k,UNBOUNDED_START,UNBOUNDED_END,repeats,dataFile,printResults,UNBOUNDED_FP);
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
	free(algorithms);
}

void runTest(char *text, char *pattern, int n, int m, int k, char *startText,
		char *endText,int repeats,FILE *dataFile,short print,
		void (*FP)(char *,char *,int,int,int,int *,struct SP_KM_MATCHING_POSITIONS *, unsigned int)){

	clock_t start;


	printf("%s\n\n",startText);
	int i;
	struct SP_KM_MATCHING_POSITIONS *listOfMatches;
	if(print){
		listOfMatches = NULL;
	}else{
		listOfMatches = sp_km_create_new_list_of_matches();
	}
	int numResults = 0;
	double totalTime = 0.0;
	//Compute the LCE table for the bounded algorithm
	start = clock();
	for(i=0;i<repeats;i++){
		start = clock();
		FP(text,pattern,n,m,k,&numResults,listOfMatches,0);
		printf("Number of k-mismatches found: %d (out of a possible %d alignments)\n",numResults,n-m+1);
		double result = ((double)clock()-start)/(double)CLOCKS_PER_SEC;
		printf("(TIMER)Search #%d took %.2f\n",i,result);
		totalTime += result;
	}
	double avgTime = totalTime/(double)repeats;
	printf("\n(TIMER)Average Time For %s was %.2f\n",endText,avgTime);
	printf("------------------------------------------\n");
	if(dataFile != NULL){
		fprintf(dataFile,"%.2f ",avgTime);
	}
	if(!print){
		sp_km_freeListOfMatches(listOfMatches);
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



