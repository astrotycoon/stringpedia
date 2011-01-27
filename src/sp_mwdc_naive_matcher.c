/*
 * naive.c
 *
 *  Created on: 8 Jul 2010
 *      Author: ben
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sp_mwdc_naive_matcher.h"

void sp_mwdc_match_naively(char *text, char *pattern,char wildCardChar,int n, int m,
		int *numMatches, struct SP_MWDC_MATCHING_POSITIONS *listOfMatches, unsigned int flags){

	//set variable for flags
	short firstMatchOnly = (flags & SP_MWDC_FIRST_MATCH_ONLY) ? 1 : 0;

	//Reset counts
	*numMatches = 0;

	if(flags & SP_MWDC_DO_NAIVE_CONVOLUTIONS){
		matchNaivelyWithConvolutions_Internal(text,pattern,wildCardChar,n,m,firstMatchOnly,numMatches,listOfMatches);
	}else{
		matchNaively_Internal(text,pattern,wildCardChar,n,m,firstMatchOnly,numMatches,listOfMatches);
	}
}

void matchNaively_Internal(char *text, char *pattern,char wildCardChar, int n, int m,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	//loop though each possible index and check each value
	int i,j,matchFound;

	for(i=0;i<=n-m;i++){
		matchFound = 1;
		j = 0;
		for(j=0;j<m;j++){
			if(text[i+j] != pattern[j] && text[i+j] != wildCardChar && pattern[j] != wildCardChar){
				matchFound = 0;
				break;
			}
		}
		if(matchFound){
			(*numMatches)++;
			if(listOfMatches != NULL){
				sp_mwdc_addToListOfMatches(listOfMatches,i);
			}else{
				printf("Match at index %d\n",i);
			}
			if(firstMatchOnly){
				return;
			}
		}
	}
}

void matchNaivelyWithConvolutions_Internal(char *text, char *pattern,char wildCardChar, int n, int m,short firstMatchOnly,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	//loop though each possible index and check each value
	int i;

	//Allocate space and assign to double arrays
	//Reverse pattern.

	double *t, *tsquared, *tcubed,*p,*psquared, *pcubed;
	t = (double *)  malloc(sizeof(double) * n);
	tsquared = (double *)  malloc(sizeof(double) * n);
	tcubed = (double *)  malloc(sizeof(double) * n);
	p = (double *)  malloc(sizeof(double) * m);
	psquared = (double *)  malloc(sizeof(double) * m);
	pcubed = (double *)  malloc(sizeof(double) * m);

	for(i=0; i<n; i++){
		t[i] = (text[i]==wildCardChar) ? 0.0 : (double) text[i];
		tsquared[i] = t[i] * t[i];
		tcubed[i] = t[i] * tsquared[i];
	}

	for(i=0; i<m; i++){
		p[i] = (pattern[m-i-1]==wildCardChar) ? 0.0 : (double) pattern[m-i-1];
		psquared[i] = p[i] * p[i];
		pcubed[i] = p[i] * psquared[i];
	}

	double *pCubedTimesT,*pSquaredTimesTSquared,*pTimesTCubed;
	pCubedTimesT = (double *)  malloc(sizeof(double) * n);
	pSquaredTimesTSquared = (double *)  malloc(sizeof(double) * n);
	pTimesTCubed = (double *)  malloc(sizeof(double) * n);

	multiplyPolynomialTruncated(t,pcubed,pCubedTimesT,n,m);
	multiplyPolynomialTruncated(tsquared,psquared,pSquaredTimesTSquared,n,m);
	multiplyPolynomialTruncated(tcubed,p,pTimesTCubed,n,m);


//	for(i=0;i<n;i++){
//		//printf("testing index %d, offset value %d\n",i,i+indexOffset);
//		double result = pCubedTimesT[i] - 2*pSquaredTimesTSquared[i] + pTimesTCubed[i];
//		printf("index: %d\n",i);
//		printf("First part : %f\n",pCubedTimesT[i]);
//		printf("Second part: %f\n",2*pSquaredTimesTSquared[i]);
//		printf("Third part : %f\n",pTimesTCubed[i]);
//		if(fabs(result) == 0.0){
//			printf("Match at index %d\n",i);
//			//Work required to ascertain a suitable threshold to accept input
//		}else if(fabs(result) <= 0.1){
//			printf("Near match at index %d with value %f\n",i,result);
//		}
//	}

	for(i=0;i<=n-m;i++){
		//printf("testing index %d, offset value %d\n",i,i+indexOffset);
		double result = pCubedTimesT[m+i-1] - 2*pSquaredTimesTSquared[m+i-1] + pTimesTCubed[m+i-1];
//		printf("index: %d\n",i);
//		printf("First part : %f\n",pCubedTimesT[m+i-1]);
//		printf("Second part: %f\n",2*pSquaredTimesTSquared[m+i-1]);
//		printf("Third part : %f\n",pTimesTCubed[m+i-1]);
		if(fabs(result) <= 0.1){
			(*numMatches)++;
			if(listOfMatches != NULL){
				sp_mwdc_addToListOfMatches(listOfMatches,i);
			}else{
				printf("Match at index %d\n",i);
			}
			if(firstMatchOnly){
				return;
			}
		}
	}

	free(t);
	free(tsquared);
	free(tcubed);

	free(p);
	free(psquared);
	free(pcubed);

	free(pTimesTCubed);
	free(pSquaredTimesTSquared);
	free(pCubedTimesT);

}

void matchNaively(char *text, char *pattern,char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchNaively_Internal(text,pattern,wildCardChar,n,m,0,numMatches,listOfMatches);
}

void matchNaivelyWithConvolutions(char *text, char *pattern,char wildCardChar, int n, int m,int *numMatches,struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	matchNaivelyWithConvolutions_Internal(text,pattern,wildCardChar,n,m,0,numMatches,listOfMatches);
}

inline void multiplyPolynomialTruncated(double *t, double *p, double *result,int n, int m){
	int i,j,start,end;
	start = 0;

	for(i=0;i<n;i++){
		//start = (i-n+1 > 0)? i-n+1 : 0;
		end = (i < m-1)? i : m-1;
		//result[i] = p[start] * t[i-start];
		result[i] = p[start] * t[i];
		//printf("result[%d]= %f * %f = %f\n",i,p[start],t[i-start],p[start] * t[i-start]);
		for(j=1;j<=end;j++){
			result[i] += p[j] * t[i-j];
			//printf("result[%d]+= %f * %f = %f\n",i,p[j],t[i-j],result[i]);
		}
		//printf("start: %d, end: %d \n",start,end);
	}
}
