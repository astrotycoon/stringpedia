/*
 * flint.h
 *
 *  Created on: 14 Jul 2010
 *      Author: ben
 *      Don't include any FLINT header files 
 */

#ifndef FLINT_H_
#define FLINT_H_

#include "sp_mwdc.h"
#include <stdlib.h>

void sp_mwdc_match_with_flint(char *text, char *pattern,char wildCardChar,int n, int m,
		int *numMatches, struct SP_MWDC_MATCHING_POSITIONS *listOfMatches, unsigned int flags);

void matchWithFlint_Internal(char *text, char *pattern,char wildCardChar,int n, int m,
		short firstMatchOnly, int *numMatches, struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

void matchWithFlint(char *text, char *pattern,char wildCardChar,int n, int m,
		int *numMatches, struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

#endif /* FLINT_H_ */
