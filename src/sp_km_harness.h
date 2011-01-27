/*
 * harness.h
 *
 *  Created on: 29 Jul 2010
 *      Author: ben
 */

#ifndef HARNESS_H_
#define HARNESS_H_

#include "sp_km_naive_matcher.h"
#include "sp_km_unbounded_matcher.h"
#include "sp_km.h"

#define NAIVE 1
#define NAIVE_START "Naive"
#define NAIVE_END "Naive"
#define NAIVE_FP &sp_km_naive_kmismatch

#define UNBOUNDED 2
#define UNBOUNDED_START "Unbounded by k, O(n SQRT(m log m))"
#define UNBOUNDED_END "Unbounded by k, O(n SQRT(m log m))"
#define UNBOUNDED_FP &sp_km_unbounded_kmismatch

int getLineLength(FILE *file);

void runTest(char *text, char *pattern, int n, int m, int k, char *startText,
		char *endText,int repeats,FILE *dataFile,short print,
		void (*FP)(char *,char *,int,int,int,int *,struct SP_KM_MATCHING_POSITIONS *, unsigned int));

#endif /* HARNESS_H_ */
