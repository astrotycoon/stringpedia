/*
 * fft-random.c
 *
 *  Created on: 26 Jul 2010
 *      Author: ben
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <fftw3.h>
#include <time.h>
#include "sp_mwdc_fft_random_matcher.h"


void sp_mwdc_match_with_fftw_randomized(char *text, char *pattern,char wildCardChar,int n, int m,
		int *numMatches, struct SP_MWDC_MATCHING_POSITIONS *listOfMatches, unsigned int flags){

	//Set variables for use in "internal" functions
	unsigned int FFTW_FLAGS = (flags & SP_MWDC_FFTW_MEASURE) ? FFTW_MEASURE : FFTW_ESTIMATE;
	short pad = (flags & SP_MWDC_NO_PADDING) ? 0 : 1;
	short import = (flags & SP_MWDC_NO_WISDOM) ? 0 : 1;
	short dontClean = (flags & SP_MWDC_DONT_CLEAN_WISDOM) ? 0 : 1;
	short minSize = (flags & SP_MWDC_NO_MINSIZE) ? 0 : SP_MWDC_MINSIZE;
	short firstMatchOnly = (flags & SP_MWDC_FIRST_MATCH_ONLY) ? 1 : 0;
	short check = (flags & SP_MWDC_VERIFY_NAIVELY) ? 1 : 0;

	//Reset the counts
	*numMatches = 0;

	if(flags & SP_MWDC_NO_WILDS_IN_TEXT){
		if(flags & SP_MWDC_REAL2REAL){
			matchWithRandomizedFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_FLAGS,pad,import,minSize,check,dontClean,firstMatchOnly,numMatches,listOfMatches);
		}else{
			matchWithRandomizedFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_FLAGS,pad,import,minSize,check,dontClean,firstMatchOnly,numMatches,listOfMatches);
		}
	}else if(flags & SP_MWDC_REAL2REAL){
		matchWithKalai_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_FLAGS,pad,import,minSize,check,dontClean,firstMatchOnly,numMatches,listOfMatches);
	}else{
		matchWithKalai_Internal(text,pattern,wildCardChar,n,m,FFTW_FLAGS,pad,import,minSize,check,dontClean,firstMatchOnly,numMatches,listOfMatches);
	}
}

void matchWithRandomizedFFT_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short checkNaively,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){


	int i;
	int transformSize = 2*m;
	int characterShift = 96;

	if(padToPowerOf2){
		transformSize = roundUpToPowerOf2(transformSize);
	}

	if(minTransformSize < n && transformSize < minTransformSize){
		transformSize = minTransformSize;
	}
	if(importSysWisdom){
		if(!fftw_import_system_wisdom()){
			printf("Failed to read wisdom!\n");
		}
	}


	//Create an array to map values in text and pattern to.

	int M = (int) sqrt( (double)UINT_MAX/(double)m  );
	if (M > RAND_MAX){
		M = RAND_MAX;
	}

	int *map = malloc(sizeof(int) * M);

	for(i=0;i<M;i++){
		map[i] = i;
	}
	//Now shuffle them

	//Try to read seed from text file -- allows tracking of the seed number so
	//we can re-test any false positives
	FILE *seedfile = fopen("randseed.txt","r");
	unsigned seed;
	if(seedfile==NULL || fscanf(seedfile,"%u",&seed) !=1){
		seed = time(NULL);
	}

	srand(seed);
	fclose(seedfile);

	for(i=M-1;i>=0;i--){
		int random = rand() % (i+1); //scale numbers between 0 and i
		int temp = map[random];
		map[random] = map[i];
		map[i] = temp;
	}


	//start by performing the work on the pattern -- this is unchanged in each chunk

	double *p,sumOfPSquared;
	p = fftw_malloc(sizeof(double) * transformSize);
	sumOfPSquared = 0.0;

	for(i=0; i<m; i++){
		if(pattern[m-i-1]==wildCardChar){
			p[i] = 0.0;
		}else{
			if((int)pattern[m-i-1]-characterShift >= M){
				printf("(SEED) Using seed: %u\n",seed);
				printf("(ERROR) alphabet too large. Cannot map character %c (value: %d) to value in range[0...%d]\n",pattern[m-i-1],(int)pattern[m-i-1]-characterShift,M);
				exit(0);
			}
			p[i] = (double) map[(int)pattern[m-i-1]-characterShift];
			sumOfPSquared += p[i] * p[i];
		}
		//printf("p[%d] set to %f\n",i,p[i]);
	}

	for(;i<transformSize;i++){
		p[i] = 0.0;
	}

	//printf("Sum of p squared: %f\n",sumOfPSquared);
	/*Calculate the DFT of the pattern */

	fftw_complex *DFTofPattern = fftw_malloc(sizeof(fftw_complex) * transformSize);

	fftw_plan forward,inverse;
	forward = fftw_plan_dft_r2c_1d(transformSize,p,DFTofPattern,FFTW_FLAGS);
	inverse = 0;

	fftw_execute(forward);
	//Free the original pattern values
 	fftw_free(p);

	//pre-compute values for t, tsquared and tcubed

	//also add transformsize-m high order 0 values for potential "overflow" of the chunks.
	int overflow = transformSize - m;
	double *t_PRECOMPUTE = (double *) fftw_malloc(sizeof(double) * (n+overflow));

	for(i=0; i<n; i++){
		if((int)text[i]-characterShift >= M){
			printf("(SEED) Using seed: %u\n",seed);
			printf("(ERROR) alphabet too large. Cannot map character %c (value: %d) to value in range[0...%d]\n",text[i],(int)text[i]-characterShift,M);
			exit(0);
		}
		t_PRECOMPUTE[i] = (double) map[(int)text[i]-characterShift];
		//printf("t_PRECOMPUTE[%d] set to %f\n",i,t_PRECOMPUTE[i]);
	}

	for(; i<n+overflow; i++){
		t_PRECOMPUTE[i] = 0.0;
	}

 	//Allocate space for text values and DFTs of text
	double *t = fftw_malloc(sizeof(double) * transformSize);

	fftw_complex *DFTofText = fftw_malloc(sizeof(fftw_complex) * transformSize);

	//Also allocate space for the products

	fftw_complex *pTimesT_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
	double *pTimesT = fftw_malloc(sizeof(double) * transformSize);

	//Iteratively divide the text into overlapping chunks of size 2m, and search these.
	//char *textPortion = (char *) malloc(sizeof(char) * chunkSize);
	int start = 0;
	//strncpy(textPortion,text,m*2);

	/* Unless n =m the index n-m will have been tested on the previous loop, so use less than, not lte */
	while(start < n-m || start==0){
		memcpy(t,t_PRECOMPUTE+start,sizeof(double)*transformSize);
		fftw_execute_dft_r2c(forward,t,DFTofText);

		/* Multiply the point representations*/
		for(i=0;i<transformSize;i++){
			pTimesT_PR[i][0] = DFTofPattern[i][0] * DFTofText[i][0] - DFTofPattern[i][1] * DFTofText[i][1];
			pTimesT_PR[i][1] = DFTofPattern[i][0] * DFTofText[i][1] + DFTofPattern[i][1] * DFTofText[i][0];
			//printf("pTimesT_PR[%d] = %f %+f i\n",i,pTimesT_PR[i][0],pTimesT_PR[i][1]);
		}

		/* Convert back to a coefficient representation */
		if(start==0){
			inverse = fftw_plan_dft_c2r_1d(transformSize,pTimesT_PR,pTimesT,FFTW_FLAGS);
			fftw_execute(inverse);
		}else{
			fftw_execute_dft_c2r(inverse,pTimesT_PR,pTimesT);
		}
		printRandomizedFFTMatches(pTimesT,sumOfPSquared,transformSize,n,m,start,checkNaively,text,pattern,seed,wildCardChar,firstMatchOnly,numMatches,listOfMatches);
		start +=transformSize-m+1;

	}
	//free allocated memory
	fftw_free(t);
	fftw_free(t_PRECOMPUTE);

	fftw_free(DFTofText);
	fftw_free(DFTofPattern);

	fftw_free(pTimesT);
	fftw_free(pTimesT_PR);
	free(map);

	fftw_destroy_plan(forward);
	fftw_destroy_plan(inverse);
	if(!dontCleanWisdom){
		fftw_forget_wisdom();
	}
	fftw_cleanup();

}

void matchWithRandomizedFFT_R2R_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short checkNaively,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;
	int transformSize = 2*m;
	int characterShift = 96;

	if(padToPowerOf2){
		transformSize = roundUpToPowerOf2(transformSize);
	}

	if(minTransformSize < n && transformSize < minTransformSize){
		transformSize = minTransformSize;
	}
	if(importSysWisdom){
		if(!fftw_import_system_wisdom()){
			printf("Failed to read wisdom!\n");
		}
	}


	//Create an array to map values in text and pattern to.

	int M = (int) sqrt( (double)UINT_MAX/(double)m  );
	if (M > RAND_MAX){
		M = RAND_MAX;
	}

	int *map = malloc(sizeof(int) * M);

	for(i=0;i<M;i++){
		map[i] = i;
	}
	//Now shuffle them

	//Try to read seed from text file -- allows tracking of the seed number so
	//we can re-test any false positives
	FILE *seedfile = fopen("randseed.txt","r");
	unsigned seed;
	if(seedfile==NULL || fscanf(seedfile,"%u",&seed) !=1){
		seed = time(NULL);
	}

	srand(seed);
	fclose(seedfile);

	for(i=M-1;i>=0;i--){
		int random = rand() % (i+1); //scale numbers between 0 and i
		int temp = map[random];
		map[random] = map[i];
		map[i] = temp;
	}


	//start by performing the work on the pattern -- this is unchanged in each chunk

	double *p,sumOfPSquared;
	p = fftw_malloc(sizeof(double) * transformSize);
	sumOfPSquared = 0.0;

	for(i=0; i<m; i++){
		if(pattern[m-i-1]==wildCardChar){
			p[i] = 0.0;
		}else{
			if((int)pattern[m-i-1]-characterShift >= M){
				printf("(SEED) Using seed: %u\n",seed);
				printf("(ERROR) alphabet too large. Cannot map character %c (value: %d) to value in range[0...%d]\n",pattern[m-i-1],(int)pattern[m-i-1]-characterShift,M);
				exit(0);
			}
			p[i] = (double) map[(int)pattern[m-i-1]-characterShift];
			sumOfPSquared += p[i] * p[i];
		}
		//printf("p[%d] set to %f\n",i,p[i]);
	}

	for(;i<transformSize;i++){
		p[i] = 0.0;
	}

	//printf("Sum of p squared: %f\n",sumOfPSquared);
	/*Calculate the DFT of the pattern */

	double *DFTofPattern = fftw_malloc(sizeof(double) * transformSize);

	fftw_plan forward,inverse;
	forward = fftw_plan_r2r_1d(transformSize,p,DFTofPattern,FFTW_R2HC,FFTW_FLAGS);
	inverse = 0;

	fftw_execute(forward);
	//Free the original pattern values
 	fftw_free(p);

	//pre-compute values for t, tsquared and tcubed

	//also add transformsize-m high order 0 values for potential "overflow" of the chunks.
	int overflow = transformSize - m;
	double *t_PRECOMPUTE = (double *) fftw_malloc(sizeof(double) * (n+overflow));

	for(i=0; i<n; i++){
		if((int)text[i]-characterShift >= M){
			printf("(SEED) Using seed: %u\n",seed);
			printf("(ERROR) alphabet too large. Cannot map character %c (value: %d) to value in range[0...%d]\n",text[i],(int)text[i]-characterShift,M);
			exit(0);
		}
		t_PRECOMPUTE[i] = (double) map[(int)text[i]-characterShift];
		//printf("t_PRECOMPUTE[%d] set to %f\n",i,t_PRECOMPUTE[i]);
	}

	for(; i<n+overflow; i++){
		t_PRECOMPUTE[i] = 0.0;
	}

 	//Allocate space for text values and DFTs of text
	double *t = fftw_malloc(sizeof(double) * transformSize);

	double *DFTofText = fftw_malloc(sizeof(double) * transformSize);

	//Also allocate space for the products

	double *pTimesT_PR = fftw_malloc(sizeof(double) * transformSize);
	double *pTimesT = fftw_malloc(sizeof(double) * transformSize);

	//Iteratively divide the text into overlapping chunks of size 2m, and search these.
	//char *textPortion = (char *) malloc(sizeof(char) * chunkSize);
	int start = 0;
	//strncpy(textPortion,text,m*2);

	/* Unless n =m the index n-m will have been tested on the previous loop, so use less than, not lte */
	while(start < n-m || start==0){
		memcpy(t,t_PRECOMPUTE+start,sizeof(double)*transformSize);
		fftw_execute_r2r(forward,t,DFTofText);
		/* Multiply the point representations*/
		pTimesT_PR[0] = DFTofPattern[0] * DFTofText[0];
 		if(transformSize % 2==0){
 			for(i=1;i<transformSize/2;i++){
 				pTimesT_PR[i] = DFTofPattern[i] * DFTofText[i] - DFTofPattern[transformSize-i] * DFTofText[transformSize-i];
 				pTimesT_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofText[i] + DFTofPattern[i] * DFTofText[transformSize-i];

 			}
 			pTimesT_PR[i] = DFTofPattern[i] * DFTofText[i];
 		}else{
 			for(i=1;i<=transformSize/2;i++){
 				pTimesT_PR[i] = DFTofPattern[i] * DFTofText[i] - DFTofPattern[transformSize-i] * DFTofText[transformSize-i];
 				pTimesT_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofText[i] + DFTofPattern[i] * DFTofText[transformSize-i];

 			}
 		}

		/* Convert back to a coefficient representation */
		if(start==0){
			inverse = fftw_plan_r2r_1d(transformSize,pTimesT_PR,pTimesT,FFTW_HC2R,FFTW_FLAGS);
			fftw_execute(inverse);
		}else{
			fftw_execute_r2r(inverse,pTimesT_PR,pTimesT);
		}
		printRandomizedFFTMatches(pTimesT,sumOfPSquared,transformSize,n,m,start,checkNaively,text,pattern,seed,wildCardChar,firstMatchOnly,numMatches,listOfMatches);
		start +=transformSize-m+1;

	}
	//free allocated memory
	fftw_free(t);
	fftw_free(t_PRECOMPUTE);

	fftw_free(DFTofText);
	fftw_free(DFTofPattern);

	fftw_free(pTimesT);
	fftw_free(pTimesT_PR);
	free(map);

	fftw_destroy_plan(forward);
	fftw_destroy_plan(inverse);
	if(!dontCleanWisdom){
		fftw_forget_wisdom();
	}
	fftw_cleanup();

}


void printRandomizedFFTMatches(double *pTimesT,double sumOfPSquared,int chunkSize,int n,int m, int indexOffset,
		short checkNaively,char *text, char *pattern,unsigned seed,char wildCardChar,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;

	for(i=0;i<=chunkSize-m && (i+indexOffset)<=n-m;i++){
		double result = sumOfPSquared - pTimesT[m+i-1]/(chunkSize);
		if(fabs(result) <= 0.1){
			if(!checkNaively){
				(*numMatches)++;
				if(listOfMatches != NULL){
					sp_mwdc_addToListOfMatches(listOfMatches,i+indexOffset);
				}else{
					printf("Match at index %d\n",i+indexOffset);
				}
				if(firstMatchOnly){
					return;
				}
			}else{
				//double check result
				short matches = 1;
				int j;
				int shift = i+indexOffset;
				for(j=0;j<m;j++){
					if(pattern[j]!=wildCardChar && pattern[j] != text[shift+j]){
						matches = 0;
						break;
					}
				}
				if(matches){
					(*numMatches)++;
					if(listOfMatches != NULL){
						sp_mwdc_addToListOfMatches(listOfMatches,i+indexOffset);
					}else{
						printf("Match at index %d verified naively\n",i+indexOffset);
					}
					if(firstMatchOnly){
						return;
					}
				}else{
					printf("INCORRECT(ORIGINAL) MATCH IDENTIFIED AT INDEX %d (using seed: %u)\n",i+indexOffset,seed);
				}

			}
		}
	}
}

void matchWithKalai_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short checkNaively,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;
	int transformSize = 2*m;

	if(padToPowerOf2){
		transformSize = roundUpToPowerOf2(transformSize);
	}

	if(minTransformSize < n && transformSize < minTransformSize){
		transformSize = minTransformSize;
	}
	if(importSysWisdom){
		if(!fftw_import_system_wisdom()){
			printf("Failed to read wisdom!\n");
		}
	}
	//Try to read seed from text file -- allows tracking of the seed number so
	//we can re-test any false positives
	FILE *seedfile = fopen("randseed.txt","r");
	unsigned seed;
	if(seedfile==NULL || fscanf(seedfile,"%u",&seed) !=1){
		seed = time(NULL);
	}

	srand(seed);
	fclose(seedfile);

	//M=max random value
	unsigned M = ( sqrt(DBL_MAX/(double)m)  > RAND_MAX/(m/2) ) ? RAND_MAX/(m/2) : ((unsigned)sqrt(DBL_MAX/(double)m));
	//int M = (int) sqrt( (double)UINT_MAX/(double)m  );
	//printf("M: %d\n",M);

	double *r, *rTimesP;
	r = fftw_malloc(sizeof(double) * transformSize);
	rTimesP = fftw_malloc(sizeof(double) * transformSize);

	for(i=0;i<m;i++){
		if(pattern[i]==wildCardChar){
			r[m-i-1] = 0.0;
			rTimesP[m-i-1] = 0.0;
		}else{
			r[m-i-1] = (double) (rand()%(M-1) + 1);
			rTimesP[m-i-1] = pattern[i] * r[m-i-1];
		}
	}

	for(;i<transformSize;i++){
		r[i] = 0.0;
		rTimesP[i] = 0.0;
	}
	/*Calculate the DFT of r and r times p */

	fftw_complex *DFTofR = fftw_malloc(sizeof(fftw_complex) * transformSize);
	fftw_complex *DFTofRTimesP = fftw_malloc(sizeof(fftw_complex) * transformSize);

	fftw_plan forward,inverse;
	forward = fftw_plan_dft_r2c_1d(transformSize,r,DFTofR,FFTW_FLAGS);
	inverse = 0;
	fftw_execute(forward);
	fftw_execute_dft_r2c(forward,rTimesP,DFTofRTimesP);

	//Free the original pattern values
 	fftw_free(r);
 	fftw_free(rTimesP);

	//pre-compute values for t, tsquared and tcubed

	//also add transformsize-m high order 0 values for potential "overflow" of the chunks.
	int overflow = transformSize - m;
	double *t_PRECOMPUTE = (double *) fftw_malloc(sizeof(double) * (n+overflow));
	double *tMask_PRECOMPUTE = (double *) fftw_malloc(sizeof(double) * (n+overflow));

	for(i=0; i<n; i++){
		if(text[i]==wildCardChar){
			t_PRECOMPUTE[i] = 0.0;
			tMask_PRECOMPUTE[i] = 0.0;
		}else{
			t_PRECOMPUTE[i] = (double)text[i];
			tMask_PRECOMPUTE[i] = 1.0;
		}
	}

	for(; i<n+overflow; i++){
		t_PRECOMPUTE[i] = 0.0;
		tMask_PRECOMPUTE[i] =0.0;
	}

 	//Allocate space for text values and DFTs of text
	double *t = fftw_malloc(sizeof(double) * transformSize);
	double *tmask = fftw_malloc(sizeof(double) * transformSize);
	fftw_complex *DFTofText = fftw_malloc(sizeof(fftw_complex) * transformSize);
	fftw_complex *DFTofTextMask = fftw_malloc(sizeof(fftw_complex) * transformSize);

	//Also allocate space for the products

	fftw_complex *tTimesR_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
	fftw_complex *tMaskTimesPTimesR_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
	double *tTimesR = fftw_malloc(sizeof(double) * transformSize);
	double *tMaskTimesPTimesR = fftw_malloc(sizeof(double) * transformSize);

	//Iteratively divide the text into overlapping chunks of size 2m, and search these.
	//char *textPortion = (char *) malloc(sizeof(char) * chunkSize);
	int start = 0;
	//strncpy(textPortion,text,m*2);

	/* Unless n =m the index n-m will have been tested on the previous loop, so use less than, not lte */
	while(start <= n-m){
		memcpy(t,t_PRECOMPUTE+start,sizeof(double)*transformSize);
		memcpy(tmask,tMask_PRECOMPUTE+start,sizeof(double)*transformSize);
		fftw_execute_dft_r2c(forward,t,DFTofText);
		fftw_execute_dft_r2c(forward,tmask,DFTofTextMask);

		/* Multiply the point representations*/
		for(i=0;i<transformSize;i++){
			tTimesR_PR[i][0] = DFTofText[i][0] * DFTofR[i][0] - DFTofText[i][1] * DFTofR[i][1];
			tTimesR_PR[i][1] = DFTofText[i][0] * DFTofR[i][1] + DFTofText[i][1] * DFTofR[i][0];

			tMaskTimesPTimesR_PR[i][0] = DFTofTextMask[i][0] * DFTofRTimesP[i][0] - DFTofTextMask[i][1] * DFTofRTimesP[i][1];
			tMaskTimesPTimesR_PR[i][1] = DFTofTextMask[i][0] * DFTofRTimesP[i][1] + DFTofTextMask[i][1] * DFTofRTimesP[i][0];

		}

		/* Convert back to a coefficient representation */
		if(start==0){
			inverse = fftw_plan_dft_c2r_1d(transformSize,tTimesR_PR,tTimesR,FFTW_FLAGS);
			fftw_execute(inverse);
		}else{
			fftw_execute_dft_c2r(inverse,tTimesR_PR,tTimesR);
		}
		fftw_execute_dft_c2r(inverse,tMaskTimesPTimesR_PR,tMaskTimesPTimesR);

		printKalaiFFTMatches(tTimesR,tMaskTimesPTimesR,transformSize,n,m,start,checkNaively,text,pattern,seed,wildCardChar,firstMatchOnly,numMatches,listOfMatches);
		start +=transformSize-m+1;

	}
	//free allocated memory
	fftw_free(t);
	fftw_free(tmask);
	fftw_free(t_PRECOMPUTE);
	fftw_free(tMask_PRECOMPUTE);


	fftw_free(DFTofR);
	fftw_free(DFTofRTimesP);
	fftw_free(DFTofText);
	fftw_free(DFTofTextMask);

	fftw_free(tTimesR);
	fftw_free(tTimesR_PR);
	fftw_free(tMaskTimesPTimesR);
	fftw_free(tMaskTimesPTimesR_PR);

	fftw_destroy_plan(forward);
	fftw_destroy_plan(inverse);
	if(!dontCleanWisdom){
		fftw_forget_wisdom();
	}
	fftw_cleanup();

}

void matchWithKalai_R2R_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short checkNaively,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;
	int transformSize = 2*m;

	if(padToPowerOf2){
		transformSize = roundUpToPowerOf2(transformSize);
	}

	if(minTransformSize < n && transformSize < minTransformSize){
		transformSize = minTransformSize;
	}
	if(importSysWisdom){
		if(!fftw_import_system_wisdom()){
			printf("Failed to read wisdom!\n");
		}
	}

	//Try to read seed from text file -- allows tracking of the seed number so
	//we can re-test any false positives
	FILE *seedfile = fopen("randseed.txt","r");
	unsigned seed;
	if(seedfile==NULL || fscanf(seedfile,"%u",&seed) !=1){
		seed = time(NULL);
	}

	srand(seed);
	fclose(seedfile);

	//M=max random value
	unsigned M = ( sqrt(DBL_MAX/(double)m)  > RAND_MAX/(m/2) ) ? RAND_MAX/(m/2) : ((unsigned)sqrt(DBL_MAX/(double)m));
	//int M = (int) sqrt( (double)UINT_MAX/(double)m  );
	//printf("M: %d\n",M);

	double *r, *rTimesP;
	r = fftw_malloc(sizeof(double) * transformSize);
	rTimesP = fftw_malloc(sizeof(double) * transformSize);

	for(i=0;i<m;i++){
		if(pattern[i]==wildCardChar){
			r[m-i-1] = 0.0;
			rTimesP[m-i-1] = 0.0;
		}else{
			r[m-i-1] = (double) (rand()%(M-1) + 1);
			rTimesP[m-i-1] = pattern[i] * r[m-i-1];
		}
	}

	for(;i<transformSize;i++){
		r[i] = 0.0;
		rTimesP[i] = 0.0;
	}
	/*Calculate the DFT of r and r times p */

	double *DFTofR = fftw_malloc(sizeof(double) * transformSize);
	double *DFTofRTimesP = fftw_malloc(sizeof(double) * transformSize);

	fftw_plan forward,inverse;
	inverse = 0;
	forward = fftw_plan_r2r_1d(transformSize,r,DFTofR,FFTW_R2HC,FFTW_FLAGS);
	fftw_execute(forward);
	fftw_execute_r2r(forward,rTimesP,DFTofRTimesP);

	//Free the original pattern values
 	fftw_free(r);
 	fftw_free(rTimesP);

	//pre-compute values for t, tsquared and tcubed

	//also add transformsize-m high order 0 values for potential "overflow" of the chunks.
	int overflow = transformSize - m;
	double *t_PRECOMPUTE = (double *) fftw_malloc(sizeof(double) * (n+overflow));
	double *tMask_PRECOMPUTE = (double *) fftw_malloc(sizeof(double) * (n+overflow));

	for(i=0; i<n; i++){
		if(text[i]==wildCardChar){
			t_PRECOMPUTE[i] = 0.0;
			tMask_PRECOMPUTE[i] = 0.0;
		}else{
			t_PRECOMPUTE[i] = (double)text[i];
			tMask_PRECOMPUTE[i] = 1.0;
		}
	}

	for(; i<n+overflow; i++){
		t_PRECOMPUTE[i] = 0.0;
		tMask_PRECOMPUTE[i] =0.0;
	}

 	//Allocate space for text values and DFTs of text
	double *t = fftw_malloc(sizeof(double) * transformSize);
	double *tmask = fftw_malloc(sizeof(double) * transformSize);
	double *DFTofText = fftw_malloc(sizeof(double) * transformSize);
	double *DFTofTextMask = fftw_malloc(sizeof(double) * transformSize);

	//Also allocate space for the products

	double *tTimesR_PR = fftw_malloc(sizeof(double) * transformSize);
	double *tMaskTimesPTimesR_PR = fftw_malloc(sizeof(double) * transformSize);
	double *tTimesR = fftw_malloc(sizeof(double) * transformSize);
	double *tMaskTimesPTimesR = fftw_malloc(sizeof(double) * transformSize);

	//Iteratively divide the text into overlapping chunks of size 2m, and search these.
	//char *textPortion = (char *) malloc(sizeof(char) * chunkSize);
	int start = 0;
	//strncpy(textPortion,text,m*2);

	/* Unless n =m the index n-m will have been tested on the previous loop, so use less than, not lte */
	while(start <= n-m){
		memcpy(t,t_PRECOMPUTE+start,sizeof(double)*transformSize);
		memcpy(tmask,tMask_PRECOMPUTE+start,sizeof(double)*transformSize);
		fftw_execute_r2r(forward,t,DFTofText);
		fftw_execute_r2r(forward,tmask,DFTofTextMask);

		/* Multiply the point representations*/
		tTimesR_PR[0] = DFTofText[0] * DFTofR[0];
		tMaskTimesPTimesR_PR[0] = DFTofTextMask[0] * DFTofRTimesP[0];

 		if(transformSize % 2==0){
 			for(i=1;i<transformSize/2;i++){
 				tTimesR_PR[i] = DFTofText[i] * DFTofR[i] - DFTofText[transformSize-i] * DFTofR[transformSize-i];
 				tTimesR_PR[transformSize-i] = DFTofText[transformSize-i] * DFTofR[i] + DFTofText[i] * DFTofR[transformSize-i];

 				tMaskTimesPTimesR_PR[i] = DFTofTextMask[i] * DFTofRTimesP[i] - DFTofTextMask[transformSize-i] * DFTofRTimesP[transformSize-i];
 				tMaskTimesPTimesR_PR[transformSize-i] = DFTofTextMask[transformSize-i] * DFTofRTimesP[i] + DFTofTextMask[i] * DFTofRTimesP[transformSize-i];

 			}
 			tTimesR_PR[i] = DFTofText[i] * DFTofR[i];
 			tMaskTimesPTimesR_PR[i] = DFTofTextMask[i] * DFTofRTimesP[i];
 		}else{
 			for(i=1;i<=transformSize/2;i++){
 				tTimesR_PR[i] = DFTofText[i] * DFTofR[i] - DFTofText[transformSize-i] * DFTofR[transformSize-i];
 				tTimesR_PR[transformSize-i] = DFTofText[transformSize-i] * DFTofR[i] + DFTofText[i] * DFTofR[transformSize-i];

 				tMaskTimesPTimesR_PR[i] = DFTofTextMask[i] * DFTofRTimesP[i] - DFTofTextMask[transformSize-i] * DFTofRTimesP[transformSize-i];
 				tMaskTimesPTimesR_PR[transformSize-i] = DFTofTextMask[transformSize-i] * DFTofRTimesP[i] + DFTofTextMask[i] * DFTofRTimesP[transformSize-i];

 			}
 		}

		/* Convert back to a coefficient representation */
		if(start==0){
			inverse = fftw_plan_r2r_1d(transformSize,tTimesR_PR,tTimesR,FFTW_HC2R,FFTW_FLAGS);
			fftw_execute(inverse);
		}else{
			fftw_execute_r2r(inverse,tTimesR_PR,tTimesR);
		}
		fftw_execute_r2r(inverse,tMaskTimesPTimesR_PR,tMaskTimesPTimesR);

		printKalaiFFTMatches(tTimesR,tMaskTimesPTimesR,transformSize,n,m,start,checkNaively,text,pattern,seed,wildCardChar,firstMatchOnly,numMatches,listOfMatches);
		start +=transformSize-m+1;

	}
	//free allocated memory
	fftw_free(t);
	fftw_free(tmask);
	fftw_free(t_PRECOMPUTE);
	fftw_free(tMask_PRECOMPUTE);


	fftw_free(DFTofR);
	fftw_free(DFTofRTimesP);
	fftw_free(DFTofText);
	fftw_free(DFTofTextMask);

	fftw_free(tTimesR);
	fftw_free(tTimesR_PR);
	fftw_free(tMaskTimesPTimesR);
	fftw_free(tMaskTimesPTimesR_PR);

	fftw_destroy_plan(forward);
	fftw_destroy_plan(inverse);
	if(!dontCleanWisdom){
		fftw_forget_wisdom();
	}
	fftw_cleanup();

}


void printKalaiFFTMatches(double *tTimesR,double *tMaskTimesPTimesR,int chunkSize,int n,int m, int indexOffset,
		short checkNaively,char *text, char *pattern,unsigned seed,char wildCardChar,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;

	for(i=0;i<=chunkSize-m && (i+indexOffset)<=n-m;i++){
		double result = tTimesR[m+i-1] - tMaskTimesPTimesR[m+i-1];
		//printf("Index %d; result=%f\n",i+indexOffset,result);
		if(fabs(result) <= 0.1){
			if(!checkNaively){
				(*numMatches)++;
				if(listOfMatches != NULL){
					sp_mwdc_addToListOfMatches(listOfMatches,i+indexOffset);
				}else{
					printf("Match at index %d\n",i+indexOffset);
				}
				if(firstMatchOnly){
					return;
				}
			}else{
				//double check result
				short matches = 1;
				int j;
				int shift = i+indexOffset;
				for(j=0;j<m;j++){
					//printf("Comparing pattern[%d] with text[%d]. Values are: %c and %c\n",j,shift+j,pattern[j],text[shift+j]);
					if(pattern[j]!=wildCardChar && pattern[j] != text[shift+j]){
						//printf("Mismatch found at pattern[%d] with text[%d]. Values are: %c and %c\n",j,shift+j,pattern[j],text[shift+j]);
						matches = 0;
						break;
					}
				}
				if(matches){
					(*numMatches)++;
					if(listOfMatches != NULL){
						sp_mwdc_addToListOfMatches(listOfMatches,i+indexOffset);
					}else{
						printf("Match at index %d verified naively\n",i+indexOffset);
					}
					if(firstMatchOnly){
						return;
					}
				}else{
					printf("INCORRECT(KALAI) MATCH IDENTIFIED AT INDEX %d (using seed: %u)\n",i+indexOffset,seed);
				}

			}
		}
	}
}

void matchWithKalai(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithKalai_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,0,0,0,numMatches,listOfMatches);
}

void matchWithKalai_CheckNaively(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithKalai_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}

void matchWithKalai_PrintOnlyFailures(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithKalai_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}

void matchWithKalai_R2R(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithKalai_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,0,0,0,numMatches,listOfMatches);
}

void matchWithKalai_R2R_CheckNaively(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithKalai_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}

void matchWithKalai_R2R_PrintOnlyFailures(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithKalai_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}




void matchWithRandomizedFFT(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithRandomizedFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,0,0,0,numMatches,listOfMatches);
}

void matchWithRandomizedFFT_CheckNaively(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithRandomizedFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}

void matchWithRandomizedFFT_PrintOnlyFailures(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithRandomizedFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}

void matchWithRandomizedFFT_R2R(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithRandomizedFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,0,0,0,numMatches,listOfMatches);
}

void matchWithRandomizedFFT_R2R_CheckNaively(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithRandomizedFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}

void matchWithRandomizedFFT_R2R_PrintOnlyFailures(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithRandomizedFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}




