/*
 * inputFromCorpus.c
 *
 *  Created on: 13 Aug 2010
 *      Author: ben
 *      Creates pattern and text by parsing the supplied file, stripping new
 *      lines and control characters the randomly taking a chunk of text to
 *      use as the pattern and garbling it such that an exact match does not
 *      occur, but a k-mismatch does
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>

void printUsage(){
	printf("About: Takes a file and generates a text and pattern from that ");
    printf("(strips new lines and other control characters in the algorithms) before writing them to a file");
	printf("Usage: ./inputFromCorpus inputfile m k [options]\n\n");
	printf("Options: -n value -- specify a maximum length on the created text");
	printf("         -s seed  -- use the given seed, rather than the current time");
	printf("         -o file  -- file to write to. Default is generated.txt");
}

/**
 * Comparator for qsort.
 */
int comparator(const void *a, const void *b){
	return *((int *)a) - *((int  *)b);
}

int main (int argc, char *argv[]){
	if(argc < 4 || (argc==2 && strcmp(argv[1],"-help")==0)){
		printUsage();
		exit(1);
	}
	int n,m,k,i,c,j;
	n = INT_MAX;
	int seedSpecified = 0;
	unsigned seed = 0;
	FILE *in, *out;
	out = NULL;

	in = fopen(argv[1],"r");
	if(in==NULL){
		printf("Could not open input file: %s for reading\n",argv[1]);
		exit(1);
	}

	m = atoi(argv[2]);
	k = atoi(argv[3]);

	if(k >= m){
		printf("k must be less than m!\n");
		exit(1);
	}

	for(i=4;i<argc;i++){
		if(strcmp(argv[i],"-n")==0){
			i++;
			if(i>=argc){
				printf("No value specified for option -n\n");
				exit(1);
			}
			n = atoi(argv[i]);
		}else if(strcmp(argv[i],"-s")==0){
			i++;
			if(i>=argc){
				printf("No value specified for option -s\n");
				exit(1);
			}
			seed = atoi(argv[i]);
			seedSpecified = 1;
		}else if(strcmp(argv[i],"-o")==0){
			i++;
			if(i>=argc){
				printf("No value specified for option -o\n");
				exit(1);
			}
			out = fopen(argv[i],"w+");
			if(out==NULL){
				printf("Could not open input file: %s for writing\n",argv[i]);
				exit(1);
			}
		}else{
			printf("Warning: unknown option %s\n",argv[i]);
		}
	}

	//If out still not opened, use default
	if(out==NULL){
		out = fopen("generated.txt","w+");
		if(out==NULL){
			printf("Could not open input file: generated.txt for writing\n");
			exit(1);
		}
	}

	//Write text to file
	i = 0;
	while((c=fgetc(in))!=EOF){
		if(c != '\n' && c != 7 && c != '$'){
			fputc(c,out);
			if(++i >= n){
				break;
			}
		}
	}
	fclose(in);
	fputc('\n',out);
	n = i;
	if(i <= m){
		printf("Text is only %d characters long; cannot produce a pattern of length %d!\n",i,m);
	}

	if(seedSpecified==0){
		seed = time(NULL);
	}
	srand(seed);

	//Create an array of locations to create mis-matches at
	int *mismatches = (int *) malloc(sizeof(int) * m);

	for(i=0;i<m;i++){
		mismatches[i] = i;
	}

	//Now shuffle them
	for(i=m-1;i>=0;i--){
		int random = rand() % (i+1); //scale numbers between 0 and i
		int temp = mismatches[random];
		mismatches[random] = mismatches[i];
		mismatches[i] = temp;
	}

	//Reduce array size to k and sort
	mismatches = (int *) realloc(mismatches,sizeof(int) * k);
	qsort(mismatches,k,sizeof(int),comparator);

	//Find a position in the text to extract a pattern from
	j = rand() % (n-m);
	//Plus one for \0 which is appended
	char *pattern = (char *)malloc(sizeof(char) * (m+1));
	fseek(out,j,SEEK_SET);
	//printf("Reading pattern from location: %d\n",j);
	fgets(pattern,(m+1),out);
	if(pattern==NULL){
		printf("Failed to read pattern!\n");
		exit(1);
	}else{
		//printf("pattern: %s\n",pattern);
	}

	//Do some slight garballingmismatches[i]
	for(i=0;i<k;i++){
		//printf("garbling pattern[%d]=%c\n",mismatches[i],pattern[mismatches[i]]);
		c = pattern[mismatches[i]];
		c++;
		if(c != '\n' && c != 7 && c != '$'){
			pattern[mismatches[i]] = c;
		}
	}

	fseek(out,0,SEEK_END);
	fputs(pattern,out);
	fputc('\n',out);
	fclose(out);
	free(mismatches);
	free(pattern);

	//Write information to seed-log
	FILE *seedlog = fopen("seeds.log","a");
	fprintf(seedlog,"Seed: %u | m: %d | k: %d | n: %d\n",seed,m,k,n);
	fclose(seedlog);


}
