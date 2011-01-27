/*
 * naive.h
 *
 *  Created on: 29 Jul 2010
 *      Author: ben
 */

#ifndef NAIVE_H_
#define NAIVE_H_

#include "sp_km.h"

void sp_km_naive_kmismatch(char *text, char *pattern, int n, int m,int k,
		int *numMatches, struct SP_KM_MATCHING_POSITIONS *listOfMatches, unsigned int flag);

#endif /* NAIVE_H_ */

