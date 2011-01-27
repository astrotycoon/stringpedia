/*
 * fft.c
 *
 *  Created on: 8 Jul 2010
 *      Author: ben
 */

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "sp_mwdc_fft_matcher.h"

#define TIME_MEASURE 0
#define TIME_ESTIMATE 0
#define TIME_WISDOM 0


void sp_mwdc_match_with_fftw(char *text, char *pattern,char wildCardChar,int n, int m,
		int *numMatches, struct SP_MWDC_MATCHING_POSITIONS *listOfMatches, unsigned int flags){

	//Set variables for use in "internal" functions
	unsigned int FFTW_FLAGS = (flags & SP_MWDC_FFTW_MEASURE) ? FFTW_MEASURE : FFTW_ESTIMATE;
	short pad = (flags & SP_MWDC_NO_PADDING) ? 0 : 1;
	short import = (flags & SP_MWDC_NO_WISDOM) ? 0 : 1;
	short noWilds = (flags & SP_MWDC_NO_WILDS_IN_TEXT) ? 0 : 1;
	short dontClean = (flags & SP_MWDC_DONT_CLEAN_WISDOM) ? 0 : 1;
	short minSize = (flags & SP_MWDC_NO_MINSIZE) ? 0 : SP_MWDC_MINSIZE;
	short firstMatchOnly = (flags & SP_MWDC_FIRST_MATCH_ONLY) ? 1 : 0;

	//Reset the counts
	*numMatches = 0;

	if(flags & SP_MWDC_NLOGM){
		if(flags & SP_MWDC_REAL2REAL){
			matchWithFasterFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_FLAGS,
					pad,import,minSize,noWilds,dontClean,firstMatchOnly,numMatches,listOfMatches);
		}else{
			matchWithFasterFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_FLAGS,
					pad,import,minSize,noWilds,dontClean,firstMatchOnly,numMatches,listOfMatches);
		}
	}else if(flags & SP_MWDC_REAL2REAL){
		matchWithFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_FLAGS,
				pad,import,noWilds,dontClean,firstMatchOnly,numMatches,listOfMatches);
	}else{
		matchWithFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_FLAGS,
				pad,import,noWilds,dontClean,firstMatchOnly,numMatches,listOfMatches);
	}
}

/**
 * In order to allow function pointers in the harness, this "internal" function is used by other functions
 * which set additional parameters
 */
void matchWithFFT_Internal(char *text, char *pattern,char wildCardChar,int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2,short importSysWisdom,short noWildCardsInText,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;
	clock_t startTime;

	if(importSysWisdom){
		if(!fftw_import_system_wisdom()){
			printf("Failed to read wisdom!\n");
		}
	}

	/* Convert text to an array of doubles; and double the degree bound.
	 * Compute squares and cubes of each element as we go */
	//unsigned transformSize = n*2;
	unsigned transformSize = n;
	if(padToPowerOf2){
		transformSize = roundUpToPowerOf2(transformSize);
	}
	double *t, *tsquared, *tcubed;
	t = fftw_malloc(sizeof(double) * transformSize);
	tsquared = fftw_malloc(sizeof(double) * transformSize);
	tcubed = NULL;
	//If no wild cards in text, we don't need tcubed (or it's transform etc)
	if(!noWildCardsInText){
		tcubed = fftw_malloc(sizeof(double) * transformSize);
	}

	//if we reuse plans, we'll need these defined
	fftw_plan forward, inverse;

	for(i=0; i<n; i++){
		t[i] = (text[i]==wildCardChar) ? 0.0 : (double) text[i];
		tsquared[i] = t[i] * t[i];
		if (!noWildCardsInText) tcubed[i] = t[i] * tsquared[i];
	}

	for(;i<transformSize;i++){
		t[i] = 0.0;
		tsquared[i] = 0.0;
		if (!noWildCardsInText) tcubed[i] = 0.0;
	}

	/* Convert pattern to an array of doubles and pad to the same size as the text
	 * Reverse pattern as we go
	 *Compute squares and cubes of each element as we go */

	//If no wild cards in text, we an also compute the fixed sum of p^3

	double *p,*psquared, *pcubed,sumOfPCubed;
	p = fftw_malloc(sizeof(double) * transformSize);
	psquared = fftw_malloc(sizeof(double) * transformSize);
	pcubed = NULL;
	if (noWildCardsInText){
		sumOfPCubed = 0;
	}else{
		pcubed= fftw_malloc(sizeof(double) * transformSize);
	}


	for(i=0; i<m; i++){
		p[i] = (pattern[m-i-1]==wildCardChar) ? 0.0 : (double) pattern[m-i-1];
		psquared[i] = p[i] * p[i];
		if (noWildCardsInText){
			sumOfPCubed += p[i] * psquared[i];
		}else{
			pcubed[i] = p[i] * psquared[i];
		}

	}

	for(;i<transformSize;i++){
		p[i] = 0.0;
		psquared[i] = 0.0;
		if (!noWildCardsInText) pcubed[i] = 0.0;
	}

	/*Calculate the DFT of the text, square of text and cube of text*/
	fftw_complex *DFTofText, *DFTofTextSquared, *DFTofTextCubed;
	DFTofText = fftw_malloc(sizeof(fftw_complex) * transformSize);
	DFTofTextSquared = fftw_malloc(sizeof(fftw_complex) * transformSize);
	if (!noWildCardsInText) DFTofTextCubed = fftw_malloc(sizeof(fftw_complex) * transformSize);

	if( (TIME_MEASURE && FFTW_FLAGS==FFTW_MEASURE) || (TIME_ESTIMATE && FFTW_FLAGS==FFTW_ESTIMATE) || (importSysWisdom && TIME_WISDOM))
		startTime = clock();
	forward = fftw_plan_dft_r2c_1d(transformSize,t,DFTofText,FFTW_FLAGS);
	if( (TIME_MEASURE && FFTW_FLAGS==FFTW_MEASURE) || (TIME_ESTIMATE && FFTW_FLAGS==FFTW_ESTIMATE) || (importSysWisdom && TIME_WISDOM))
		printf("(INTERNAL TIMER)FFTW_MEASURE for transform of size %d took %.2f\n",transformSize,((double)clock()-startTime)/(double)CLOCKS_PER_SEC);
	fftw_execute(forward);
	fftw_execute_dft_r2c(forward,tsquared,DFTofTextSquared);
	if (!noWildCardsInText) fftw_execute_dft_r2c(forward,tcubed,DFTofTextCubed);

	//Free the original text values
	fftw_free(t);
	fftw_free(tsquared);
	if (!noWildCardsInText) fftw_free(tcubed);

	/*Calculate the DFT of the pattern, square of pattern and cube of pattern */

	fftw_complex *DFTofPattern, *DFTofPatternSquared, *DFTofPatternCubed;
	DFTofPattern = fftw_malloc(sizeof(fftw_complex) * transformSize);
	DFTofPatternSquared = fftw_malloc(sizeof(fftw_complex) * transformSize);
	if (!noWildCardsInText) DFTofPatternCubed = fftw_malloc(sizeof(fftw_complex) * transformSize);

	fftw_execute_dft_r2c(forward,p,DFTofPattern);
	fftw_execute_dft_r2c(forward,psquared,DFTofPatternSquared);
	if (!noWildCardsInText) fftw_execute_dft_r2c(forward,pcubed,DFTofPatternCubed);

	//Free the original pattern values
 	fftw_free(p);
 	fftw_free(psquared);
 	if (!noWildCardsInText) fftw_free(pcubed);

	/* Multiply the point representations*/

 	fftw_complex *pCubedTimesT_PR, *pSquaredTimesTSquared_PR, *pTimesTCubed_PR, *pTimesTSquared_PR, *pSquaredTimesT_PR;
 	if(noWildCardsInText){
 		pTimesTSquared_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
 		pSquaredTimesT_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
 	}else{
		pCubedTimesT_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
		pSquaredTimesTSquared_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
		pTimesTCubed_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
 	}

 	if(noWildCardsInText){
		for(i=0;i<transformSize;i++){

			pSquaredTimesT_PR[i][0] = DFTofPatternSquared[i][0] * DFTofText[i][0] - DFTofPatternSquared[i][1] * DFTofText[i][1];
			pSquaredTimesT_PR[i][1] = DFTofPatternSquared[i][0] * DFTofText[i][1] + DFTofPatternSquared[i][1] * DFTofText[i][0];

			pTimesTSquared_PR[i][0] = DFTofPattern[i][0] * DFTofTextSquared[i][0] - DFTofPattern[i][1] * DFTofTextSquared[i][1];
			pTimesTSquared_PR[i][1] = DFTofPattern[i][0] * DFTofTextSquared[i][1] + DFTofPattern[i][1] * DFTofTextSquared[i][0];
		}
 	}else{
		for(i=0;i<transformSize;i++){
			pCubedTimesT_PR[i][0] = DFTofPatternCubed[i][0] * DFTofText[i][0] - DFTofPatternCubed[i][1] * DFTofText[i][1];
			pCubedTimesT_PR[i][1] = DFTofPatternCubed[i][0] * DFTofText[i][1] + DFTofPatternCubed[i][1] * DFTofText[i][0];

			pSquaredTimesTSquared_PR[i][0] = DFTofPatternSquared[i][0] * DFTofTextSquared[i][0] - DFTofPatternSquared[i][1] * DFTofTextSquared[i][1];
			pSquaredTimesTSquared_PR[i][1] = DFTofPatternSquared[i][0] * DFTofTextSquared[i][1] + DFTofPatternSquared[i][1] * DFTofTextSquared[i][0];

			pTimesTCubed_PR[i][0] = DFTofPattern[i][0] * DFTofTextCubed[i][0] - DFTofPattern[i][1] * DFTofTextCubed[i][1];
			pTimesTCubed_PR[i][1] = DFTofPattern[i][0] * DFTofTextCubed[i][1] + DFTofPattern[i][1] * DFTofTextCubed[i][0];

		}
 	}

	//Free the DFTs
 	fftw_free(DFTofText);
 	fftw_free(DFTofTextSquared);
 	if (!noWildCardsInText) fftw_free(DFTofTextCubed);
 	fftw_free(DFTofPattern);
 	fftw_free(DFTofPatternSquared);
 	if (!noWildCardsInText) fftw_free(DFTofPatternCubed);

	/* Convert back to a coefficient representation */
	double *pCubedTimesT, *pSquaredTimesTSquared, *pTimesTCubed,*pTimesTSquared,*pSquaredTimesT;
	if(noWildCardsInText){
		pTimesTSquared = fftw_malloc(sizeof(double) * transformSize);
		pSquaredTimesT = fftw_malloc(sizeof(double) * transformSize);
	}else{
		pCubedTimesT = fftw_malloc(sizeof(double) * transformSize);
		pSquaredTimesTSquared = fftw_malloc(sizeof(double) * transformSize);
		pTimesTCubed= fftw_malloc(sizeof(double) * transformSize);
	}

	if(noWildCardsInText){
		inverse = fftw_plan_dft_c2r_1d(transformSize,pTimesTSquared_PR,pTimesTSquared,FFTW_FLAGS);
		fftw_execute(inverse);
		fftw_execute_dft_c2r(inverse,pSquaredTimesT_PR,pSquaredTimesT);

		printFFTMatches_NoWildInText(sumOfPCubed,pSquaredTimesT,pTimesTSquared,transformSize,n,m,0,firstMatchOnly,numMatches,listOfMatches);

		//Free the rest of the allocated values
		fftw_free(pTimesTSquared_PR);
		fftw_free(pSquaredTimesT_PR);
		fftw_free(pTimesTSquared);
		fftw_free(pSquaredTimesT);

	}else{
		if( (TIME_MEASURE && FFTW_FLAGS==FFTW_MEASURE) || (TIME_ESTIMATE && FFTW_FLAGS==FFTW_ESTIMATE))
			startTime = clock();
		inverse = fftw_plan_dft_c2r_1d(transformSize,pCubedTimesT_PR,pCubedTimesT,FFTW_FLAGS);
		if( (TIME_MEASURE && FFTW_FLAGS==FFTW_MEASURE) || (TIME_ESTIMATE && FFTW_FLAGS==FFTW_ESTIMATE))
			printf("(INTERNAL TIMER)FFTW_MEASURE for transform of size %d took %.2f\n",transformSize,((double)clock()-startTime)/(double)CLOCKS_PER_SEC);
		fftw_execute(inverse);
		fftw_execute_dft_c2r(inverse,pSquaredTimesTSquared_PR,pSquaredTimesTSquared);
		fftw_execute_dft_c2r(inverse,pTimesTCubed_PR,pTimesTCubed);

		printFFTMatches(pCubedTimesT,pSquaredTimesTSquared,pTimesTCubed,transformSize,n,m,0,firstMatchOnly,numMatches,listOfMatches);

		//Free the rest of the allocated values

		fftw_free(pCubedTimesT_PR);
		fftw_free(pSquaredTimesTSquared_PR);
		fftw_free(pTimesTCubed_PR);

		fftw_free(pCubedTimesT);
		fftw_free(pSquaredTimesTSquared);
		fftw_free(pTimesTCubed);
	}

	fftw_destroy_plan(forward);
	fftw_destroy_plan(inverse);
	if(!dontCleanWisdom){
		fftw_forget_wisdom();
	}
 	fftw_cleanup();
}


void matchWithFFT_R2R_Internal(char *text, char *pattern,char wildCardChar,int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2,short importSysWisdom,short noWildCardsInText,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;

	if(importSysWisdom){
		if(!fftw_import_system_wisdom()){
			printf("Failed to read wisdom!\n");
		}
	}

	/* Convert text to an array of doubles; and double the degree bound.
	 * Compute squares and cubes of each element as we go */
	unsigned transformSize = n;
	if(padToPowerOf2){
		transformSize = roundUpToPowerOf2(transformSize);
	}
	double *t, *tsquared, *tcubed;
	t = fftw_malloc(sizeof(double) * transformSize);
	tsquared = fftw_malloc(sizeof(double) * transformSize);
	tcubed = NULL;

	//If no wild cards in text, we don't need tcubed (or it's transform etc)
	if(!noWildCardsInText){
		tcubed = fftw_malloc(sizeof(double) * transformSize);
	}

	//if we reuse plans, we'll need these defined
	fftw_plan forward, inverse;

	for(i=0; i<n; i++){
		t[i] = (text[i]==wildCardChar) ? 0.0 : (double) text[i];
		tsquared[i] = t[i] * t[i];
		if (!noWildCardsInText) tcubed[i] = t[i] * tsquared[i];
	}

	for(;i<transformSize;i++){
		t[i] = 0.0;
		tsquared[i] = 0.0;
		if (!noWildCardsInText) tcubed[i] = 0.0;
	}

	/* Convert pattern to an array of doubles and pad to the same size as the text
	 * Reverse pattern as we go
	 *Compute squares and cubes of each element as we go */

	//If no wild cards in text, we an also compute the fixed sum of p^3

	double *p,*psquared, *pcubed,sumOfPCubed;
	p = fftw_malloc(sizeof(double) * transformSize);
	psquared = fftw_malloc(sizeof(double) * transformSize);
	pcubed = NULL;
	if (noWildCardsInText){
		sumOfPCubed = 0;
	}else{
		pcubed= fftw_malloc(sizeof(double) * transformSize);
	}

	for(i=0; i<m; i++){
		p[i] = (pattern[m-i-1]==wildCardChar) ? 0.0 : (double) pattern[m-i-1];
		psquared[i] = p[i] * p[i];
		if (noWildCardsInText){
			sumOfPCubed += p[i] * psquared[i];
		}else{
			pcubed[i] = p[i] * psquared[i];
		}

	}

	for(;i<transformSize;i++){
		p[i] = 0.0;
		psquared[i] = 0.0;
		if (!noWildCardsInText) pcubed[i] = 0.0;
	}

	/*Calculate the DFT of the text, square of text and cube of text*/
	double *DFTofText, *DFTofTextSquared, *DFTofTextCubed;
	DFTofText = fftw_malloc(sizeof(double) * transformSize);
	DFTofTextSquared = fftw_malloc(sizeof(double) * transformSize);
	if (!noWildCardsInText) DFTofTextCubed = fftw_malloc(sizeof(double) * transformSize);

	forward = fftw_plan_r2r_1d(transformSize,t,DFTofText,FFTW_R2HC,FFTW_FLAGS);
	fftw_execute(forward);
	fftw_execute_r2r(forward,tsquared,DFTofTextSquared);
	if (!noWildCardsInText) fftw_execute_r2r(forward,tcubed,DFTofTextCubed);

	//Free the original text values
	fftw_free(t);
	fftw_free(tsquared);
	if (!noWildCardsInText) fftw_free(tcubed);

	/*Calculate the DFT of the pattern, square of pattern and cube of pattern */

	double *DFTofPattern, *DFTofPatternSquared, *DFTofPatternCubed;
	DFTofPattern = fftw_malloc(sizeof(double) * transformSize);
	DFTofPatternSquared = fftw_malloc(sizeof(double) * transformSize);
	if (!noWildCardsInText) DFTofPatternCubed = fftw_malloc(sizeof(double) * transformSize);

	fftw_execute_r2r(forward,p,DFTofPattern);
	fftw_execute_r2r(forward,psquared,DFTofPatternSquared);
	if (!noWildCardsInText) fftw_execute_r2r(forward,pcubed,DFTofPatternCubed);

	//Free the original pattern values
 	fftw_free(p);
 	fftw_free(psquared);
 	if (!noWildCardsInText) fftw_free(pcubed);

	/* Multiply the point representations*/

 	double *pCubedTimesT_PR, *pSquaredTimesTSquared_PR, *pTimesTCubed_PR, *pTimesTSquared_PR, *pSquaredTimesT_PR;
 	if(noWildCardsInText){
 		pTimesTSquared_PR = fftw_malloc(sizeof(double) * transformSize);
 		pSquaredTimesT_PR = fftw_malloc(sizeof(double) * transformSize);
 	}else{
		pCubedTimesT_PR = fftw_malloc(sizeof(double) * transformSize);
		pSquaredTimesTSquared_PR = fftw_malloc(sizeof(double) * transformSize);
		pTimesTCubed_PR = fftw_malloc(sizeof(double) * transformSize);
 	}

 	if(noWildCardsInText){

		pSquaredTimesT_PR[0] = DFTofPatternSquared[0] * DFTofText[0];
		pTimesTSquared_PR[0] = DFTofPattern[0] * DFTofTextSquared[0];
		if(transformSize % 2==0){

			for(i=1;i<transformSize/2;i++){
				pSquaredTimesT_PR[i] = DFTofPatternSquared[i] * DFTofText[i] - DFTofPatternSquared[transformSize-i] * DFTofText[transformSize-i];
				pSquaredTimesT_PR[transformSize-i] = DFTofPatternSquared[transformSize-i] * DFTofText[i] + DFTofPatternSquared[i] * DFTofText[transformSize-i];

				pTimesTSquared_PR[i] = DFTofPattern[i] * DFTofTextSquared[i] - DFTofPattern[transformSize-i] * DFTofTextSquared[transformSize-i];
				pTimesTSquared_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofTextSquared[i] + DFTofPattern[i] * DFTofTextSquared[transformSize-i];

			}

			pSquaredTimesT_PR[i] = DFTofPatternSquared[i] * DFTofText[i];
			pTimesTSquared_PR[i] = DFTofPattern[i] * DFTofTextSquared[i];

		}else{

			for(i=1;i<=transformSize/2;i++){
				pSquaredTimesT_PR[i] = DFTofPatternSquared[i] * DFTofText[i] - DFTofPatternSquared[transformSize-i] * DFTofText[transformSize-i];
				pSquaredTimesT_PR[transformSize-i] = DFTofPatternSquared[transformSize-i] * DFTofText[i] + DFTofPatternSquared[i] * DFTofText[transformSize-i];

				pTimesTSquared_PR[i] = DFTofPattern[i] * DFTofTextSquared[i] - DFTofPattern[transformSize-i] * DFTofTextSquared[transformSize-i];
				pTimesTSquared_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofTextSquared[i] + DFTofPattern[i] * DFTofTextSquared[transformSize-i];
			}

		}

 	}else{
 		pCubedTimesT_PR[0] = DFTofPatternCubed[0] * DFTofText[0];
 		pSquaredTimesTSquared_PR[0] = DFTofPatternSquared[0] * DFTofTextSquared[0];
 		pTimesTCubed_PR[0] = DFTofPattern[0] * DFTofTextCubed[0];

 		if(transformSize % 2==0){

 			for(i=1;i<transformSize/2;i++){
 				pCubedTimesT_PR[i] = DFTofPatternCubed[i] * DFTofText[i] - DFTofPatternCubed[transformSize-i] * DFTofText[transformSize-i];
 				pCubedTimesT_PR[transformSize-i] = DFTofPatternCubed[transformSize-i] * DFTofText[i] + DFTofPatternCubed[i] * DFTofText[transformSize-i];

 				pSquaredTimesTSquared_PR[i] = DFTofPatternSquared[i] * DFTofTextSquared[i] - DFTofPatternSquared[transformSize-i] * DFTofTextSquared[transformSize-i];
 				pSquaredTimesTSquared_PR[transformSize-i] = DFTofPatternSquared[transformSize-i] * DFTofTextSquared[i] + DFTofPatternSquared[i] * DFTofTextSquared[transformSize-i];

 				pTimesTCubed_PR[i] = DFTofPattern[i] * DFTofTextCubed[i] - DFTofPattern[transformSize-i] * DFTofTextCubed[transformSize-i];
 				pTimesTCubed_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofTextCubed[i] + DFTofPattern[i] * DFTofTextCubed[transformSize-i];

 			}

 			pCubedTimesT_PR[i] = DFTofPatternCubed[i] * DFTofText[i];
 			pSquaredTimesTSquared_PR[i] = DFTofPatternSquared[i] * DFTofTextSquared[i];
 			pTimesTCubed_PR[i] = DFTofPattern[i] * DFTofTextCubed[i];
 		}else{

 			for(i=1;i<=transformSize/2;i++){
 				pCubedTimesT_PR[i] = DFTofPatternCubed[i] * DFTofText[i] - DFTofPatternCubed[transformSize-i] * DFTofText[transformSize-i];
 				pCubedTimesT_PR[transformSize-i] = DFTofPatternCubed[transformSize-i] * DFTofText[i] + DFTofPatternCubed[i] * DFTofText[transformSize-i];

 				pSquaredTimesTSquared_PR[i] = DFTofPatternSquared[i] * DFTofTextSquared[i] - DFTofPatternSquared[transformSize-i] * DFTofTextSquared[transformSize-i];
 				pSquaredTimesTSquared_PR[transformSize-i] = DFTofPatternSquared[transformSize-i] * DFTofTextSquared[i] + DFTofPatternSquared[i] * DFTofTextSquared[transformSize-i];

 				pTimesTCubed_PR[i] = DFTofPattern[i] * DFTofTextCubed[i] - DFTofPattern[transformSize-i] * DFTofTextCubed[transformSize-i];
 				pTimesTCubed_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofTextCubed[i] + DFTofPattern[i] * DFTofTextCubed[transformSize-i];

 			}

 		}
 	}

	//Free the DFTs
 	fftw_free(DFTofText);
 	fftw_free(DFTofTextSquared);
 	if (!noWildCardsInText) fftw_free(DFTofTextCubed);
 	fftw_free(DFTofPattern);
 	fftw_free(DFTofPatternSquared);
 	if (!noWildCardsInText) fftw_free(DFTofPatternCubed);

	/* Convert back to a coefficient representation */
	double *pCubedTimesT, *pSquaredTimesTSquared, *pTimesTCubed,*pTimesTSquared,*pSquaredTimesT;
	if(noWildCardsInText){
		pTimesTSquared = fftw_malloc(sizeof(double) * transformSize);
		pSquaredTimesT = fftw_malloc(sizeof(double) * transformSize);
	}else{
		pCubedTimesT = fftw_malloc(sizeof(double) * transformSize);
		pSquaredTimesTSquared = fftw_malloc(sizeof(double) * transformSize);
		pTimesTCubed= fftw_malloc(sizeof(double) * transformSize);
	}

	if(noWildCardsInText){
		inverse = fftw_plan_r2r_1d(transformSize,pTimesTSquared_PR,pTimesTSquared,FFTW_HC2R,FFTW_FLAGS);
		fftw_execute(inverse);
		fftw_execute_r2r(inverse,pSquaredTimesT_PR,pSquaredTimesT);

		printFFTMatches_NoWildInText(sumOfPCubed,pSquaredTimesT,pTimesTSquared,transformSize,n,m,0,firstMatchOnly,numMatches,listOfMatches);

		//Free the rest of the allocated values
		fftw_free(pTimesTSquared_PR);
		fftw_free(pSquaredTimesT_PR);
		fftw_free(pTimesTSquared);
		fftw_free(pSquaredTimesT);

	}else{
		inverse = fftw_plan_r2r_1d(transformSize,pCubedTimesT_PR,pCubedTimesT,FFTW_HC2R,FFTW_FLAGS);
		fftw_execute(inverse);
		fftw_execute_r2r(inverse,pSquaredTimesTSquared_PR,pSquaredTimesTSquared);
		fftw_execute_r2r(inverse,pTimesTCubed_PR,pTimesTCubed);

		printFFTMatches(pCubedTimesT,pSquaredTimesTSquared,pTimesTCubed,transformSize,n,m,0,firstMatchOnly,numMatches,listOfMatches);

		//Free the rest of the allocated values

		fftw_free(pCubedTimesT_PR);
		fftw_free(pSquaredTimesTSquared_PR);
		fftw_free(pTimesTCubed_PR);

		fftw_free(pCubedTimesT);
		fftw_free(pSquaredTimesTSquared);
		fftw_free(pTimesTCubed);
	}

	fftw_destroy_plan(forward);
	fftw_destroy_plan(inverse);

	if(!dontCleanWisdom){
		fftw_forget_wisdom();
	}
 	fftw_cleanup();
}

void matchWithFasterFFT_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,
		short noWildsInText,short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;
	clock_t startTime;
	int chunkSize = 2*m;
	//int transformSize = 2*chunkSize;
	int transformSize = chunkSize;
	if(padToPowerOf2){
		transformSize = roundUpToPowerOf2(transformSize);
	}
	//fftw_export_wisdom_to_file(stdout);
	//Round up to the specified minimum transform size, so long as it's less than the text size
	if(minTransformSize < n && transformSize < minTransformSize){
		transformSize = minTransformSize;
	}
	//printf("Transform size: %d\n",transformSize);
	//start by performing the work on the pattern -- this is unchanged in each chunk

	if(importSysWisdom){
		if(!fftw_import_system_wisdom()){
			printf("Failed to read wisdom!\n");
		}
	}

	/* Reverse pattern as we go
	 * Compute squares and cubes of each element as we go
	 * If no wild cards in text, we an also compute the fixed sum of p^3*/

	double *p,*psquared, *pcubed,sumOfPCubed;
	p = fftw_malloc(sizeof(double) * transformSize);
	psquared = fftw_malloc(sizeof(double) * transformSize);
	pcubed = NULL;
	if (noWildsInText){
		sumOfPCubed = 0;
	}else{
		pcubed= fftw_malloc(sizeof(double) * transformSize);
	}

	for(i=0; i<m; i++){
		p[i] = (pattern[m-i-1]==wildCardChar) ? 0.0 : (double) pattern[m-i-1];
		psquared[i] = p[i] * p[i];
		if (noWildsInText){
			sumOfPCubed += p[i] * psquared[i];
		}else{
			pcubed[i] = p[i] * psquared[i];
		}
	}

	for(;i<transformSize;i++){
		p[i] = 0.0;
		psquared[i] = 0.0;
		if (!noWildsInText) pcubed[i] = 0.0;
	}

	/*Calculate the DFT of the pattern, square of pattern and cube of pattern */

	fftw_complex *DFTofPattern, *DFTofPatternSquared, *DFTofPatternCubed;
	DFTofPattern = fftw_malloc(sizeof(fftw_complex) * transformSize);
	DFTofPatternSquared = fftw_malloc(sizeof(fftw_complex) * transformSize);
	DFTofPatternCubed = NULL;
	if(!noWildsInText) DFTofPatternCubed = fftw_malloc(sizeof(fftw_complex) * transformSize);

	//If reusing plans, we'll need these defined
	fftw_plan forward, inverse;
	inverse = 0;//stop warning of potentially uninitialized variables.

	if(FFTW_FLAGS==FFTW_MEASURE){
		startTime = clock();
	}
	forward = fftw_plan_dft_r2c_1d(transformSize,p,DFTofPattern,FFTW_FLAGS);
	if( (TIME_MEASURE && FFTW_FLAGS==FFTW_MEASURE) || (TIME_ESTIMATE && FFTW_FLAGS==FFTW_ESTIMATE) || (importSysWisdom && TIME_WISDOM))
		printf("(INTERNAL TIMER)FFTW_MEASURE for transform of size %d took %.2f\n",transformSize,((double)clock()-startTime)/(double)CLOCKS_PER_SEC);
	fftw_execute(forward);
	fftw_execute_dft_r2c(forward,psquared,DFTofPatternSquared);
	if(!noWildsInText) fftw_execute_dft_r2c(forward,pcubed,DFTofPatternCubed);

	//Free the original pattern values
 	fftw_free(p);
 	fftw_free(psquared);
 	if(!noWildsInText) fftw_free(pcubed);

 	//Allocate space for text values and DFTs of text
	double *t, *tsquared, *tcubed;
	t = fftw_malloc(sizeof(double) * transformSize);
	tsquared = fftw_malloc(sizeof(double) * transformSize);
	tcubed = NULL;
	if(!noWildsInText) tcubed = fftw_malloc(sizeof(double) * transformSize);

	fftw_complex *DFTofText, *DFTofTextSquared, *DFTofTextCubed;
	DFTofText = fftw_malloc(sizeof(fftw_complex) * transformSize);
	DFTofTextSquared = fftw_malloc(sizeof(fftw_complex) * transformSize);
	DFTofTextCubed = NULL;
	if(!noWildsInText) DFTofTextCubed = fftw_malloc(sizeof(fftw_complex) * transformSize);

	//Also allocate space for the products

	fftw_complex *pCubedTimesT_PR,*pSquaredTimesTSquared_PR,*pTimesTCubed_PR,*pSquaredTimesT_PR,*pTimesTSquared_PR;
	double *pCubedTimesT, *pSquaredTimesTSquared, *pTimesTCubed,*pSquaredTimesT,*pTimesTSquared;

	pCubedTimesT_PR = pSquaredTimesTSquared_PR = pTimesTCubed_PR = pSquaredTimesT_PR = pTimesTSquared_PR = NULL;
	pCubedTimesT = pSquaredTimesTSquared = pTimesTCubed = pSquaredTimesT = pTimesTSquared = NULL;

	if(noWildsInText){
		pSquaredTimesT_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
		pTimesTSquared_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);

		pSquaredTimesT = fftw_malloc(sizeof(double) * transformSize);
		pTimesTSquared = fftw_malloc(sizeof(double) * transformSize);
	}else{
		pCubedTimesT_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
		pSquaredTimesTSquared_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);
		pTimesTCubed_PR = fftw_malloc(sizeof(fftw_complex) * transformSize);

		pCubedTimesT = fftw_malloc(sizeof(double) * transformSize);
		pSquaredTimesTSquared = fftw_malloc(sizeof(double) * transformSize);
		pTimesTCubed= fftw_malloc(sizeof(double) * transformSize);
	}

	//Iteratively divide the text into overlapping chunks of size 2m, and search these.
	//char *textPortion = (char *) malloc(sizeof(char) * chunkSize);
	int start = 0;
	//strncpy(textPortion,text,m*2);

	//pre-compute values for t, tsquared and tcubed

	//also add transformsize-m high order 0 values for potential "overflow" of the chunks.
	int overflow = transformSize - m;
	double *t_PRECOMPUTE = (double *) malloc(sizeof(double) * (n+overflow));
	double *tsquared_PRECOMPUTE = (double *) malloc(sizeof(double) * (n+overflow));
	double *tcubed_PRECOMPUTE = NULL;
	if(!noWildsInText) tcubed_PRECOMPUTE = (double *) malloc(sizeof(double) * (n+overflow));

	for(i=0; i<n; i++){
		t_PRECOMPUTE[i] = (text[i]==wildCardChar) ? 0.0 : (double) text[i];
		tsquared_PRECOMPUTE[i] = t_PRECOMPUTE[i] * t_PRECOMPUTE[i];
		if(!noWildsInText) tcubed_PRECOMPUTE[i] = t_PRECOMPUTE[i] * tsquared_PRECOMPUTE[i];
	}

	for(; i<n+overflow; i++){
		t_PRECOMPUTE[i] = 0.0;
		tsquared_PRECOMPUTE[i] = 0.0;
		if(!noWildsInText) tcubed_PRECOMPUTE[i] = 0.0;
	}

	if(noWildsInText){
		/* Unless n =m the index n-m will have been tested on the previous loop, so use less than, not lte */
		while(start <= n-m || start==0){

			memcpy(t,t_PRECOMPUTE+start,sizeof(double)*transformSize);
			memcpy(tsquared,tsquared_PRECOMPUTE+start,sizeof(double)*transformSize);
			fftw_execute_dft_r2c(forward,t,DFTofText);
			fftw_execute_dft_r2c(forward,tsquared,DFTofTextSquared);

			for(i=0;i<transformSize;i++){
				pSquaredTimesT_PR[i][0] = DFTofPatternSquared[i][0] * DFTofText[i][0] - DFTofPatternSquared[i][1] * DFTofText[i][1];
				pSquaredTimesT_PR[i][1] = DFTofPatternSquared[i][0] * DFTofText[i][1] + DFTofPatternSquared[i][1] * DFTofText[i][0];

				pTimesTSquared_PR[i][0] = DFTofPattern[i][0] * DFTofTextSquared[i][0] - DFTofPattern[i][1] * DFTofTextSquared[i][1];
				pTimesTSquared_PR[i][1] = DFTofPattern[i][0] * DFTofTextSquared[i][1] + DFTofPattern[i][1] * DFTofTextSquared[i][0];
			}

			/* Convert back to a coefficient representation */

			//On first iteration, need to create the inverse plan
			if(start==0){
				inverse = fftw_plan_dft_c2r_1d(transformSize,pSquaredTimesT_PR,pSquaredTimesT,FFTW_FLAGS);
				fftw_execute(inverse);
			}else{
				fftw_execute_dft_c2r(inverse,pSquaredTimesT_PR,pSquaredTimesT);
			}
			fftw_execute_dft_c2r(inverse,pTimesTSquared_PR,pTimesTSquared);

			printFFTMatches_NoWildInText(sumOfPCubed,pSquaredTimesT,pTimesTSquared,transformSize,n,m,start,firstMatchOnly,numMatches,listOfMatches);
			start +=transformSize-m+1;

		}

		fftw_free(pSquaredTimesT_PR);
		fftw_free(pTimesTSquared_PR);
		fftw_free(pSquaredTimesT);
		fftw_free(pTimesTSquared);
	}else{

		/* Unless n =m the index n-m will have been tested on the previous loop, so use less than, not lte */
		while(start <= n-m || start==0){

			memcpy(t,t_PRECOMPUTE+start,sizeof(double)*transformSize);
			memcpy(tsquared,tsquared_PRECOMPUTE+start,sizeof(double)*transformSize);
			memcpy(tcubed,tcubed_PRECOMPUTE+start,sizeof(double)*transformSize);

			fftw_execute_dft_r2c(forward,t,DFTofText);
			fftw_execute_dft_r2c(forward,tsquared,DFTofTextSquared);
			fftw_execute_dft_r2c(forward,tcubed,DFTofTextCubed);


			/* Multiply the point representations*/


			for(i=0;i<transformSize;i++){
				pCubedTimesT_PR[i][0] = DFTofPatternCubed[i][0] * DFTofText[i][0] - DFTofPatternCubed[i][1] * DFTofText[i][1];
				pCubedTimesT_PR[i][1] = DFTofPatternCubed[i][0] * DFTofText[i][1] + DFTofPatternCubed[i][1] * DFTofText[i][0];

				pSquaredTimesTSquared_PR[i][0] = DFTofPatternSquared[i][0] * DFTofTextSquared[i][0] - DFTofPatternSquared[i][1] * DFTofTextSquared[i][1];
				pSquaredTimesTSquared_PR[i][1] = DFTofPatternSquared[i][0] * DFTofTextSquared[i][1] + DFTofPatternSquared[i][1] * DFTofTextSquared[i][0];

				pTimesTCubed_PR[i][0] = DFTofPattern[i][0] * DFTofTextCubed[i][0] - DFTofPattern[i][1] * DFTofTextCubed[i][1];
				pTimesTCubed_PR[i][1] = DFTofPattern[i][0] * DFTofTextCubed[i][1] + DFTofPattern[i][1] * DFTofTextCubed[i][0];

			}

			/* Convert back to a coefficient representation */

			//On first iteration, need to create the inverse plan
			if(start==0){
				if(FFTW_FLAGS==FFTW_MEASURE){
					startTime = clock();
				}
				inverse = fftw_plan_dft_c2r_1d(transformSize,pCubedTimesT_PR,pCubedTimesT,FFTW_FLAGS);
				if( (TIME_MEASURE && FFTW_FLAGS==FFTW_MEASURE) || (TIME_ESTIMATE && FFTW_FLAGS==FFTW_ESTIMATE) || (importSysWisdom && TIME_WISDOM))
					printf("(INTERNAL TIMER)FFTW_MEASURE for transform of size %d took %.2f\n",transformSize,((double)clock()-startTime)/(double)CLOCKS_PER_SEC);
				fftw_execute(inverse);
			}else{
				fftw_execute_dft_c2r(inverse,pCubedTimesT_PR,pCubedTimesT);
			}
			fftw_execute_dft_c2r(inverse,pSquaredTimesTSquared_PR,pSquaredTimesTSquared);
			fftw_execute_dft_c2r(inverse,pTimesTCubed_PR,pTimesTCubed);

			printFFTMatches(pCubedTimesT,pSquaredTimesTSquared,pTimesTCubed,transformSize,n,m,start,firstMatchOnly,numMatches,listOfMatches);

			start +=transformSize-m+1;

		}
		//free allocated memory
		fftw_free(tcubed);
		fftw_free(DFTofTextCubed);
		fftw_free(DFTofPatternCubed);
		fftw_free(tcubed_PRECOMPUTE);
		fftw_free(pCubedTimesT_PR);
		fftw_free(pSquaredTimesTSquared_PR);
		fftw_free(pTimesTCubed_PR);
		fftw_free(pCubedTimesT);
		fftw_free(pSquaredTimesTSquared);
		fftw_free(pTimesTCubed);
	}

	//free the rest of the allocated memory
	fftw_free(t);
	fftw_free(tsquared);
	fftw_free(DFTofText);
	fftw_free(DFTofTextSquared);
	fftw_free(DFTofPattern);
	fftw_free(DFTofPatternSquared);
	fftw_free(t_PRECOMPUTE);
	fftw_free(tsquared_PRECOMPUTE);


	fftw_destroy_plan(forward);
	fftw_destroy_plan(inverse);
	if(!dontCleanWisdom){
		fftw_forget_wisdom();
	}
	fftw_cleanup();

}


void matchWithFasterFFT_R2R_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short noWildsInText,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;
	clock_t startTime;
	int chunkSize = 2*m;
	//int transformSize = 2*chunkSize;
	int transformSize = chunkSize;
	if(padToPowerOf2){
		transformSize = roundUpToPowerOf2(transformSize);
	}
	//fftw_export_wisdom_to_file(stdout);
	//Round up to the specified minimum transform size, so long as it's less than the text size
	if(minTransformSize < n && transformSize < minTransformSize){
		transformSize = minTransformSize;
	}
	//printf("Transform size: %d\n",transformSize);
	//start by performing the work on the pattern -- this is unchanged in each chunk

	if(importSysWisdom){
		if(!fftw_import_system_wisdom()){
			printf("Failed to read wisdom!\n");
		}
	}

	/* Reverse pattern as we go
	 * Compute squares and cubes of each element as we go
	 * If no wild cards in text, we an also compute the fixed sum of p^3*/

	double *p,*psquared, *pcubed,sumOfPCubed;
	p = fftw_malloc(sizeof(double) * transformSize);
	psquared = fftw_malloc(sizeof(double) * transformSize);
	pcubed = NULL;
	if (noWildsInText){
		sumOfPCubed = 0;
	}else{
		pcubed= fftw_malloc(sizeof(double) * transformSize);
	}

	for(i=0; i<m; i++){
		p[i] = (pattern[m-i-1]==wildCardChar) ? 0.0 : (double) pattern[m-i-1];
		psquared[i] = p[i] * p[i];
		if (noWildsInText){
			sumOfPCubed += p[i] * psquared[i];
		}else{
			pcubed[i] = p[i] * psquared[i];
		}
	}

	for(;i<transformSize;i++){
		p[i] = 0.0;
		psquared[i] = 0.0;
		if (!noWildsInText) pcubed[i] = 0.0;
	}

	/*Calculate the DFT of the pattern, square of pattern and cube of pattern */

	double *DFTofPattern, *DFTofPatternSquared, *DFTofPatternCubed;
	DFTofPattern = fftw_malloc(sizeof(double) * transformSize);
	DFTofPatternSquared = fftw_malloc(sizeof(double) * transformSize);
	DFTofPatternCubed = NULL;
	if(!noWildsInText) DFTofPatternCubed = fftw_malloc(sizeof(double) * transformSize);

	//If reusing plans, we'll need these defined
	fftw_plan forward, inverse;
	inverse = 0;//stop warning of potentially uninitialized variables.

	if(FFTW_FLAGS==FFTW_MEASURE){
		startTime = clock();
	}
	forward = fftw_plan_r2r_1d(transformSize,p,DFTofPattern,FFTW_R2HC,FFTW_FLAGS);
	fftw_execute(forward);
	fftw_execute_r2r(forward,psquared,DFTofPatternSquared);
	if(!noWildsInText) fftw_execute_r2r(forward,pcubed,DFTofPatternCubed);

	//Free the original pattern values
 	fftw_free(p);
 	fftw_free(psquared);
 	if(!noWildsInText) fftw_free(pcubed);

 	//Allocate space for text values and DFTs of text
	double *t, *tsquared, *tcubed;
	t = fftw_malloc(sizeof(double) * transformSize);
	tsquared = fftw_malloc(sizeof(double) * transformSize);
	tcubed = NULL;
	if(!noWildsInText) tcubed = fftw_malloc(sizeof(double) * transformSize);
	double *DFTofText, *DFTofTextSquared, *DFTofTextCubed;
	DFTofText = fftw_malloc(sizeof(double) * transformSize);
	DFTofTextSquared = fftw_malloc(sizeof(double) * transformSize);
	DFTofTextCubed = NULL;
	if(!noWildsInText) DFTofTextCubed = fftw_malloc(sizeof(double) * transformSize);

	//Also allocate space for the products

	double *pCubedTimesT_PR,*pSquaredTimesTSquared_PR,*pTimesTCubed_PR,*pSquaredTimesT_PR,*pTimesTSquared_PR;
	double *pCubedTimesT, *pSquaredTimesTSquared, *pTimesTCubed,*pSquaredTimesT,*pTimesTSquared;

	pCubedTimesT_PR = pSquaredTimesTSquared_PR = pTimesTCubed_PR = pSquaredTimesT_PR = pTimesTSquared_PR = NULL;
	pCubedTimesT = pSquaredTimesTSquared = pTimesTCubed = pSquaredTimesT = pTimesTSquared = NULL;

	if(noWildsInText){
		pSquaredTimesT_PR = fftw_malloc(sizeof(double) * transformSize);
		pTimesTSquared_PR = fftw_malloc(sizeof(double) * transformSize);

		pSquaredTimesT = fftw_malloc(sizeof(double) * transformSize);
		pTimesTSquared = fftw_malloc(sizeof(double) * transformSize);
	}else{
		pCubedTimesT_PR = fftw_malloc(sizeof(double) * transformSize);
		pSquaredTimesTSquared_PR = fftw_malloc(sizeof(double) * transformSize);
		pTimesTCubed_PR = fftw_malloc(sizeof(double) * transformSize);

		pCubedTimesT = fftw_malloc(sizeof(double) * transformSize);
		pSquaredTimesTSquared = fftw_malloc(sizeof(double) * transformSize);
		pTimesTCubed= fftw_malloc(sizeof(double) * transformSize);
	}

	//Iteratively divide the text into overlapping chunks of size 2m, and search these.
	//char *textPortion = (char *) malloc(sizeof(char) * chunkSize);
	int start = 0;
	//strncpy(textPortion,text,m*2);

	//pre-compute values for t, tsquared and tcubed

	//also add transformsize-m high order 0 values for potential "overflow" of the chunks.
	int overflow = transformSize - m;
	double *t_PRECOMPUTE = (double *) malloc(sizeof(double) * (n+overflow));
	double *tsquared_PRECOMPUTE = (double *) malloc(sizeof(double) * (n+overflow));
	double *tcubed_PRECOMPUTE = NULL;
	if(!noWildsInText) tcubed_PRECOMPUTE = (double *) malloc(sizeof(double) * (n+overflow));

	for(i=0; i<n; i++){
		t_PRECOMPUTE[i] = (text[i]==wildCardChar) ? 0.0 : (double) text[i];
		tsquared_PRECOMPUTE[i] = t_PRECOMPUTE[i] * t_PRECOMPUTE[i];
		if(!noWildsInText) tcubed_PRECOMPUTE[i] = t_PRECOMPUTE[i] * tsquared_PRECOMPUTE[i];
	}

	for(; i<n+overflow; i++){
		t_PRECOMPUTE[i] = 0.0;
		tsquared_PRECOMPUTE[i] = 0.0;
		if(!noWildsInText) tcubed_PRECOMPUTE[i] = 0.0;
	}

	if(noWildsInText){
		/* Unless n =m the index n-m will have been tested on the previous loop, so use less than, not lte */
		while(start <= n-m || start==0){

			memcpy(t,t_PRECOMPUTE+start,sizeof(double)*transformSize);
			memcpy(tsquared,tsquared_PRECOMPUTE+start,sizeof(double)*transformSize);
			fftw_execute_r2r(forward,t,DFTofText);
			fftw_execute_r2r(forward,tsquared,DFTofTextSquared);


			pSquaredTimesT_PR[0] = DFTofPatternSquared[0] * DFTofText[0];
			pTimesTSquared_PR[0] = DFTofPattern[0] * DFTofTextSquared[0];

			if(transformSize % 2==0){
				for(i=1;i<transformSize/2;i++){
					pSquaredTimesT_PR[i] = DFTofPatternSquared[i] * DFTofText[i] - DFTofPatternSquared[transformSize-i] * DFTofText[transformSize-i];
					pSquaredTimesT_PR[transformSize-i] = DFTofPatternSquared[transformSize-i] * DFTofText[i] + DFTofPatternSquared[i] * DFTofText[transformSize-i];

					pTimesTSquared_PR[i] = DFTofPattern[i] * DFTofTextSquared[i] - DFTofPattern[transformSize-i] * DFTofTextSquared[transformSize-i];
					pTimesTSquared_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofTextSquared[i] + DFTofPattern[i] * DFTofTextSquared[transformSize-i];
				}
				pSquaredTimesT_PR[i] = DFTofPatternSquared[i] * DFTofText[i];
				pTimesTSquared_PR[i] = DFTofPattern[i] * DFTofTextSquared[i];
			}else{
				for(i=1;i<=transformSize/2;i++){
					pSquaredTimesT_PR[i] = DFTofPatternSquared[i] * DFTofText[i] - DFTofPatternSquared[transformSize-i] * DFTofText[transformSize-i];
					pSquaredTimesT_PR[transformSize-i] = DFTofPatternSquared[transformSize-i] * DFTofText[i] + DFTofPatternSquared[i] * DFTofText[transformSize-i];

					pTimesTSquared_PR[i] = DFTofPattern[i] * DFTofTextSquared[i] - DFTofPattern[transformSize-i] * DFTofTextSquared[transformSize-i];
					pTimesTSquared_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofTextSquared[i] + DFTofPattern[i] * DFTofTextSquared[transformSize-i];

				}
			}
			/* Convert back to a coefficient representation */

			//On first iteration, need to create the inverse plan
			if(start==0){
				inverse = fftw_plan_r2r_1d(transformSize,pSquaredTimesT_PR,pSquaredTimesT,FFTW_HC2R,FFTW_FLAGS);
				fftw_execute(inverse);
			}else{
				fftw_execute_r2r(inverse,pSquaredTimesT_PR,pSquaredTimesT);
			}
			fftw_execute_r2r(inverse,pTimesTSquared_PR,pTimesTSquared);

			printFFTMatches_NoWildInText(sumOfPCubed,pSquaredTimesT,pTimesTSquared,transformSize,n,m,start,firstMatchOnly,numMatches,listOfMatches);
			start +=transformSize-m+1;

		}

		fftw_free(pSquaredTimesT_PR);
		fftw_free(pTimesTSquared_PR);
		fftw_free(pSquaredTimesT);
		fftw_free(pTimesTSquared);
	}else{

		/* Unless n =m the index n-m will have been tested on the previous loop, so use less than, not lte */
		while(start <= n-m || start==0){

			memcpy(t,t_PRECOMPUTE+start,sizeof(double)*transformSize);
			memcpy(tsquared,tsquared_PRECOMPUTE+start,sizeof(double)*transformSize);
			memcpy(tcubed,tcubed_PRECOMPUTE+start,sizeof(double)*transformSize);

			fftw_execute_r2r(forward,t,DFTofText);
			fftw_execute_r2r(forward,tsquared,DFTofTextSquared);
			fftw_execute_r2r(forward,tcubed,DFTofTextCubed);


			/* Multiply the point representations*/

			pCubedTimesT_PR[0] = DFTofPatternCubed[0] * DFTofText[0];
			pSquaredTimesTSquared_PR[0] = DFTofPatternSquared[0] * DFTofTextSquared[0];
			pTimesTCubed_PR[0] = DFTofPattern[0] * DFTofTextCubed[0];

			if(transformSize % 2==0){

				for(i=1;i<transformSize/2;i++){
					pCubedTimesT_PR[i] = DFTofPatternCubed[i] * DFTofText[i] - DFTofPatternCubed[transformSize-i] * DFTofText[transformSize-i];
					pCubedTimesT_PR[transformSize-i] = DFTofPatternCubed[transformSize-i] * DFTofText[i] + DFTofPatternCubed[i] * DFTofText[transformSize-i];

					pSquaredTimesTSquared_PR[i] = DFTofPatternSquared[i] * DFTofTextSquared[i] - DFTofPatternSquared[transformSize-i] * DFTofTextSquared[transformSize-i];
					pSquaredTimesTSquared_PR[transformSize-i] = DFTofPatternSquared[transformSize-i] * DFTofTextSquared[i] + DFTofPatternSquared[i] * DFTofTextSquared[transformSize-i];

					pTimesTCubed_PR[i] = DFTofPattern[i] * DFTofTextCubed[i] - DFTofPattern[transformSize-i] * DFTofTextCubed[transformSize-i];
					pTimesTCubed_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofTextCubed[i] + DFTofPattern[i] * DFTofTextCubed[transformSize-i];

				}

				pCubedTimesT_PR[i] = DFTofPatternCubed[i] * DFTofText[i];
				pSquaredTimesTSquared_PR[i] = DFTofPatternSquared[i] * DFTofTextSquared[i];
				pTimesTCubed_PR[i] = DFTofPattern[i] * DFTofTextCubed[i];
			}else{

				for(i=1;i<=transformSize/2;i++){
					pCubedTimesT_PR[i] = DFTofPatternCubed[i] * DFTofText[i] - DFTofPatternCubed[transformSize-i] * DFTofText[transformSize-i];
					pCubedTimesT_PR[transformSize-i] = DFTofPatternCubed[transformSize-i] * DFTofText[i] + DFTofPatternCubed[i] * DFTofText[transformSize-i];

					pSquaredTimesTSquared_PR[i] = DFTofPatternSquared[i] * DFTofTextSquared[i] - DFTofPatternSquared[transformSize-i] * DFTofTextSquared[transformSize-i];
					pSquaredTimesTSquared_PR[transformSize-i] = DFTofPatternSquared[transformSize-i] * DFTofTextSquared[i] + DFTofPatternSquared[i] * DFTofTextSquared[transformSize-i];

					pTimesTCubed_PR[i] = DFTofPattern[i] * DFTofTextCubed[i] - DFTofPattern[transformSize-i] * DFTofTextCubed[transformSize-i];
					pTimesTCubed_PR[transformSize-i] = DFTofPattern[transformSize-i] * DFTofTextCubed[i] + DFTofPattern[i] * DFTofTextCubed[transformSize-i];

				}

			}

			/* Convert back to a coefficient representation */

			//On first iteration, need to create the inverse plan
			if(start==0){
				inverse = fftw_plan_r2r_1d(transformSize,pCubedTimesT_PR,pCubedTimesT,FFTW_HC2R,FFTW_FLAGS);
				fftw_execute(inverse);
			}else{
				fftw_execute_r2r(inverse,pCubedTimesT_PR,pCubedTimesT);
			}
			fftw_execute_r2r(inverse,pSquaredTimesTSquared_PR,pSquaredTimesTSquared);
			fftw_execute_r2r(inverse,pTimesTCubed_PR,pTimesTCubed);

			printFFTMatches(pCubedTimesT,pSquaredTimesTSquared,pTimesTCubed,transformSize,n,m,start,firstMatchOnly,numMatches,listOfMatches);

			start +=transformSize-m +1;

		}
		//free allocated memory
		fftw_free(tcubed);
		fftw_free(DFTofTextCubed);
		fftw_free(DFTofPatternCubed);
		fftw_free(tcubed_PRECOMPUTE);
		fftw_free(pCubedTimesT_PR);
		fftw_free(pSquaredTimesTSquared_PR);
		fftw_free(pTimesTCubed_PR);
		fftw_free(pCubedTimesT);
		fftw_free(pSquaredTimesTSquared);
		fftw_free(pTimesTCubed);
	}

	//free the rest of the allocated memory
	fftw_free(t);
	fftw_free(tsquared);
	fftw_free(DFTofText);
	fftw_free(DFTofTextSquared);
	fftw_free(DFTofPattern);
	fftw_free(DFTofPatternSquared);
	fftw_free(t_PRECOMPUTE);
	fftw_free(tsquared_PRECOMPUTE);
	fftw_destroy_plan(forward);
	fftw_destroy_plan(inverse);
	if(!dontCleanWisdom){
		fftw_forget_wisdom();
	}
	fftw_cleanup();

}


void printFFTMatches_NoWildInText(double sumOfPCubed,double *pSquaredTimesT,double *pTimesTSquared,int chunkSize,
		int n,int m, int indexOffset,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;

	//divide each by n as fftw does an un-normalised transform
	//Must ensure we don't test past the last possible match in this chunk OR the last possible match overall (same for nlog n flavour)
	//for(i=0;i<n;i++){
	for(i=0;i<=chunkSize-m && (i+indexOffset)<=n-m;i++){
		//printf("testing index %d, offset value %d\n",i,i+indexOffset);
		double result = sumOfPCubed - 2*(pSquaredTimesT[m+i-1]/chunkSize) + pTimesTSquared[m+i-1]/(chunkSize);
//		printf("index: %d\n",i);
//		printf("First part : %f\n",sumOfPCubed);
//		printf("Second part: %f\n",2*(pSquaredTimesT[m+i-1]/chunkSize));
//		printf("Third part : %f\n",pTimesTSquared[m+i-1]/chunkSize);
//		printf("Result: %f\n",result);
		//TODO: Work required to ascertain a suitable threshold to accept input
		if(fabs(result) <= 0.1){
			(*numMatches)++;
			if(listOfMatches != NULL){
				sp_mwdc_addToListOfMatches(listOfMatches,i+indexOffset);
			}else{
				printf("Match at index %d\n",i+indexOffset);
			}
			if(firstMatchOnly){
				return;
			}
		}
//			printf("Match at index %d\n",i+indexOffset);
//
//		}else if(fabs(result) <= 0.1){
//			printf("Near match at index %d with value %f\n",i+indexOffset,result);
//		}
	}
}

void printFFTMatches(double *pCubedTimesT,double *pSquaredTimesTSquared,double *pTimesTCubed,int chunkSize,
		int n,int m, int indexOffset,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	int i;

	//divide each by n as fftw does an un-normalised transform
	//Must ensure we don't test past the last possible match in this chunk OR the last possible match overall (same for nlog n flavour)
	//for(i=0;i<n;i++){
	for(i=0;i<=chunkSize-m && (i+indexOffset)<=n-m;i++){
		//printf("testing index %d, offset value %d\n",i,i+indexOffset);
		double result = pCubedTimesT[m+i-1]/(chunkSize) - 2*(pSquaredTimesTSquared[m+i-1]/chunkSize) + pTimesTCubed[m+i-1]/(chunkSize);
//		printf("index: %d\n",i);
//		printf("First part : %f\n",pCubedTimesT[m+i-1]/chunkSize);
//		printf("Second part: %f\n",2*(pSquaredTimesTSquared[m+i-1])/chunkSize);
//		printf("Third part : %f\n",pTimesTCubed[m+i-1]/chunkSize);
		if(fabs(result) <= 0.1){
			(*numMatches)++;
			if(listOfMatches != NULL){
				sp_mwdc_addToListOfMatches(listOfMatches,i+indexOffset);
			}else{
				printf("Match at index %d\n",i+indexOffset);
			}
			if(firstMatchOnly){
				return;
			}
		}
	}
}


/********************************************************/
/*Convenience methods/methods to allow function pointers*/
/********************************************************/


/*Match using the standard n log n FFT algorithm, but re-use the plans. Don't measure, Don't pad transform sizes to a power of 2 */
void matchWithFFT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,0,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the standard n log n FFT algorithm, but re-use the plans and measure FFT parameters. Don't pad transform sizes to a power of 2 */
void matchWithFFT_Measure(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_MEASURE,0,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the standard n log n FFT algorithm, but re-use the plans and pad the transform size to a power of 2. Don't measure. */
void matchWithFFT_Pad(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the standard n log n FFT algorithm, but re-use the plans, measure FFT parameters and pad the transform size to a power of 2.*/
void matchWithFFT_Measure_Pad(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_MEASURE,1,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the standard n log n FFT algorithm, but re-use the plans and pad the transform size to a power of 2. Don't measure, but use
 * prebuilt wisdom */
void matchWithFFT_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,0,0,0,numMatches,listOfMatches);
}

/*Match using the standard n log n FFT algorithm, but re-use the plans. Don't measure, Don't pad.
 * Use an R2R transform */
void matchWithFFT_R2R(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,0,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the standard n log n FFT algorithm, but re-use the plans and pad. Don't measure.
 * Use an R2R transform,  and pre-built wisdom */
void matchWithFFT_R2R_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,0,0,0,numMatches,listOfMatches);
}


/*Match using the standard n log n FFT algorithm, but re-use the plans. Don't measure, Don't pad.
 * Use an R2R transform */
void matchWithFFT_R2R_NWIT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,0,0,1,0,0,numMatches,listOfMatches);
}

/*Match using the standard n log n FFT algorithm, but re-use the plans and pad. Don't measure.
 * Use an R2R transform,  and pre-built wisdom */
void matchWithFFT_R2R_Wisdom_NWIT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,1,0,0,numMatches,listOfMatches);
}

/*Match using the standard n log n FFT algorithm, but re-use the plans and pad the transform size to a power of 2. Don't measure, but use
 * no wilds in text*/
void matchWithFFT_NWIT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,0,0,1,0,0,numMatches,listOfMatches);
}

/*Match using the standard n log n FFT algorithm, but re-use the plans and pad the transform size to a power of 2. Don't measure, but use
 * prebuilt wisdom. Assume no wild cards in text */
void matchWithFFT_NWIT_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,1,0,0,numMatches,listOfMatches);
}

/*Match using the 'faster' n log m FFT algorithm. Don't reuse plans, Don't measure, Don't pad transform sizes to a power of 2 */
/*Don't set a minimum transform size */
void matchWithFasterFFT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,0,0,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the 'faster' n log m FFT algorithm, but re-use the plans. Don't measure, Don't pad transform sizes to a power of 2 */
/*Don't set a minimum transform size */
void matchWithFasterFFT_Measure(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_MEASURE,0,0,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the 'faster' n log m FFT algorithm, but re-use the plans and pad the transform size to a power of 2. Don't measure. */
/*Don't set a minimum transform size */
void matchWithFasterFFT_Pad(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,0,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the 'faster' n log m FFT algorithm, but re-use the plans, measure FFT parameters and pad the transform size to a power of 2.*/
/*Don't set a minimum transform size */
void matchWithFasterFFT_Measure_Pad(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_MEASURE,1,0,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the 'faster' n log m FFT algorithm, but re-use the plans and pad the transform size to a power of 2. Use system wisdom */
/*Don't set a minimum transform size */
void matchWithFasterFFT_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,0,0,0,0,numMatches,listOfMatches);
}

/*Match using the 'faster' n log m FFT algorithm, but re-use the plans and pad the transform size to a power of 2. Use system wisdom */
/*Ensure a minimum transform size of 2048 to ensure sane run-time graphs */
void matchWithFasterFFT_MinSize(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,0,0,0,numMatches,listOfMatches);
}

/*Match using the 'faster' n log m FFT algorithm, but re-use the plans and pad the transform size to a power of 2. Use system wisdom */
/*Assume no wild cards in text */
void matchWithFasterFFT_NWIT_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,0,1,0,0,numMatches,listOfMatches);
}

/*Match using the 'faster' n log m FFT algorithm, but re-use the plans and pad the transform size to a power of 2. Use system wisdom */
/*Ensure a minimum transform size of 2048 to ensure sane run-time graphs ; Assume no wild cards in text */
void matchWithFasterFFT_NWIT_MinSize(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}

void matchWithFasterFFT_R2R(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,0,0,0,0,0,0,numMatches,listOfMatches);
}

void matchWithFasterFFT_R2R_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,0,0,0,0,numMatches,listOfMatches);
}

void matchWithFasterFFT_R2R_minsize(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,0,0,0,numMatches,listOfMatches);
}

void matchWithFasterFFT_R2R_minsize_NWIT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchWithFasterFFT_R2R_Internal(text,pattern,wildCardChar,n,m,FFTW_ESTIMATE,1,1,2048,1,0,0,numMatches,listOfMatches);
}
