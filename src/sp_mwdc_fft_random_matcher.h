/*
 * fft-random.h
 *
 *  Created on: 26 Jul 2010
 *      Author: ben
 */

#ifndef FFTRANDOM_H_
#define FFTRANDOM_H_

#include "sp_mwdc.h"

/* Flags */
/* Note that most of the flags in fft_matcher.h also apply (exception: n log m -- no n log n randomized method was implemented) */

#define SP_MWDC_VERIFY_NAIVELY (1U << 10) //Verify locations suggested naively

void sp_mwdc_match_with_fftw_randomized(char *text, char *pattern,char wildCardChar,int n, int m,
		int *numMatches, struct SP_MWDC_MATCHING_POSITIONS *listOfMatches, unsigned int flags);

void matchWithRandomizedFFT_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short checkNaively,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithRandomizedFFT_R2R_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short checkNaively,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);


void printRandomizedFFTMatches(double *pTimesT,double sumOfPSquared,int chunkSize,int n,int m, int indexOffset,
		short checkNaively,char *text, char *pattern,unsigned seed,char wildCardChar,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithKalai_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short checkNaively,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithKalai_R2R_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short checkNaively,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);


void printKalaiFFTMatches(double *tTimesR,double *tMaskTimesPTimesR,int chunkSize,int n,int m, int indexOffset,
		short checkNaively,char *text, char *pattern,unsigned seed,char wildCardChar,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithKalai(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithKalai_CheckNaively(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithKalai_PrintOnlyFailures(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithKalai_R2R(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithKalai_R2R_CheckNaively(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithKalai_R2R_PrintOnlyFailures(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithRandomizedFFT(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithRandomizedFFT_CheckNaively(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithRandomizedFFT_PrintOnlyFailures(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithRandomizedFFT_R2R(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithRandomizedFFT_R2R_CheckNaively(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithRandomizedFFT_R2R_PrintOnlyFailures(char *text, char *pattern, char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

#endif /* FFTRANDOM_H_ */
