/*
 * naive.h
 *
 *  Created on: 8 Jul 2010
 *      Author: ben
 */

#ifndef NAIVE_H_
#define NAIVE_H_

/* Flags */
#define SP_MWDC_DO_NAIVE_CONVOLUTIONS (1U << 1) //Academic interest only. Naively perform 3 convolutions between text and pattern, rather than the standard O(nm) method

#include "sp_mwdc.h"

void sp_mwdc_match_naively(char *text, char *pattern,char wildCardChar,int n, int m,
		int *numMatches, struct SP_MWDC_MATCHING_POSITIONS *listOfMatches, unsigned int flags);

void matchNaively_Internal(char *text, char *pattern,char wildCardChar, int n, int m,short firstMatchOnly,
		int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchNaivelyWithConvolutions_Internal(char *text, char *pattern,char wildCardChar,
		int n, int m,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchNaively(char *text, char *pattern,char wildCardChar, int n, int m,int *numMatches,
		struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchNaivelyWithConvolutions(char *text, char *pattern,char wildCardChar, int n, int m,
		int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void multiplyPolynomialTruncated(double *t, double *p, double *result,int n, int m);


#endif /* NAIVE_H_ */
