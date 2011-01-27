/*
 * fft.h
 *
 *  Created on: 8 Jul 2010
 *      Author: ben
 */

#ifndef SP_MWDC_FFT_MATCHER_H_
#define SP_MWDC_FFT_MATCHER_H_

/* Function operation bit field constants */


/* Time Complexity */
#define SP_MWDC_NLOGM (1U << 1)
#define SP_MWDC_NLOGN (1U << 2)

/* Types*/
#define SP_MWDC_REAL2REAL (1U << 3)

/* Flags */
/* most of these are turned on by default, so these mostly turn off behaviour */
#define SP_MWDC_NO_MINSIZE (1U << 4) //Don't set a minimum size for the transform in the n log m case
#define SP_MWDC_NO_PADDING (1U << 5) //Don't pad to a power to two
#define SP_MWDC_NO_WISDOM (1U << 6) //Don't try to load FFTW wisdom from system
#define SP_MWDC_DONT_CLEAN_WISDOM (1U << 7) //Don't clean FFTW wisdom after running the match
#define SP_MWDC_FFTW_MEASURE (1U << 8) //Get FFTW to measure plans.
#define SP_MWDC_NO_WILDS_IN_TEXT (1U << 9) //No wilds in the text -- be more efficient


/******************************************/

/* Specified as a constant as investigation into varying its size
   has not yet been made */
#ifndef SP_MWDC_MINSIZE
#define SP_MWDC_MINSIZE (2048)
#endif


#include <fftw3.h>
#include "sp_mwdc.h"

void sp_mwdc_match_with_fftw(char *text, char *pattern,char wildCardChar,int n, int m,
		int *numMatches, struct SP_MWDC_MATCHING_POSITIONS *listOfMatches, unsigned int flags);

void matchWithFFT_Internal(char *text, char *pattern,char wildCardChar,int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2,short importSysWisdom,short noWildCardsInText,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short noWildsInText,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_R2R_Internal(char *text, char *pattern,char wildCardChar,int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2,short importSysWisdom,short noWildCardsInText,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_R2R_Internal(char *text, char *pattern, char wildCardChar, int n, int m,
		unsigned int FFTW_FLAGS,short padToPowerOf2, short importSysWisdom, int minTransformSize,short noWildsInText,
		short dontCleanWisdom,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);



void doDFT(int n,double *in, fftw_complex * out);

void doInverseDFT(int n,fftw_complex *in, double * out);

void printFFTMatches(double *pCubedTimesT,double *pSquaredTimesTSquared,double *pTimesTCubed,int chunkSize,
		int n,int m, int indexOffset,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void printFFTMatches_NoWildInText(double sumOfPCubed,double *pSquaredTimesT,double *pTimesTSquared,int chunkSize,
		int n,int m, int indexOffset,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_Measure(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_Pad(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_Measure_Pad(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_NWIT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_NWIT_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_Measure(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_Pad(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_Measure_Pad(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_MinSize(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_NWIT_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_NWIT_MinSize(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_R2R(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_R2R_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_R2R_NWIT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFFT_R2R_Wisdom_NWIT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_R2R(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_R2R_Wisdom(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_R2R_minsize(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFasterFFT_R2R_minsize_NWIT(char *text, char *pattern,char wildCardChar,int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

#endif /* SP_MWDC_FFT_MATCHER_H_ */

