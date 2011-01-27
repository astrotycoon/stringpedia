/*
 * naive.c
 *
 *  Created on: 29 Jul 2010
 *      Author: ben
 *      Perform k-mismatches completely naively.
 *      Time-Complexity: O(nm)
 */
#include <stdio.h>
#include "sp_km_naive_matcher.h"
/*
 * Find k-mismatches naively.
 * @param text A character array containing the text
 * @param pattern A character array containing the pattern
 * @param n the length of the text
 * @param m the length of the pattern
 * @param k The threshold, k. If k < 0 or >= m hamming distance at all
 * alignments will be found
 * @param numMatches an unsigned pointer to a location to store the number
 * of k-mismatches found. Note that if k < 0 or >= m this is just n-m.
 */
void sp_km_naive_kmismatch(char *text, char *pattern, int n, int m,int k,
		int *numMatches, struct SP_KM_MATCHING_POSITIONS *listOfMatches, unsigned int flags){
	int i,j;
	if( k < 0){
		k = m+1;
	}
	*numMatches = 0;
	for(i=0;i<n-m+1;i++){
		int hamDistance = 0;
		for(j=0;j<m;j++){
			if(text[i+j] != pattern[j]){
				hamDistance++;
				//If hamming distance is greater than k, we dont need to bother
				//computing this anymore
				if(hamDistance > k){
					break;
				}
			}
		}

		if(hamDistance <= k){
			(*numMatches)++;
			if(listOfMatches!=NULL){
				sp_km_addToListOfMatches(listOfMatches,i,hamDistance);
			}else{
				printf("Hamming distance at text[%d]: %d\n",i,hamDistance);
			}
			if(flags & SP_KM_FIRST_MATCH_ONLY){
				return;
			}
		}

	}

}
