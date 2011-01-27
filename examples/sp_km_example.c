/*
 * sp_km_example.c
 *
 *  Created on: 5 Sep 2010
 *      Author: ben
 * Example use of k-mismatches. 
 */

#include "../include/stringpedia.h"

int main(void){
	char *text = "she sells sea shells on the sea shore";
	char *pattern = "see ";
	int k = 2;

	struct SP_KM_MATCHING_POSITIONS *listOfMatches = sp_km_create_new_list_of_matches();
	int numMatches = 0;

	//SP_KM_FIRST_MATCH_ONLY - stop after finding the first match
	sp_km_naive_kmismatch(text,pattern,37,4,k,&numMatches,listOfMatches,SP_KM_FIRST_MATCH_ONLY);

	/* An 'iterator' as well as pointers to the start and end of the linked list is provided, so that you may use  *
	 * traverse the linked list, then reset the iterator to the start if you so wish. Start should not be modified *
	 * otherwise the list cannot be freed properly                                                                 */

	printf("Printing start index of matches: \n");
	while(listOfMatches->iterator != NULL){
		printf("    Match at index: %d with hamming distance: %d\n",listOfMatches->iterator->i,listOfMatches->iterator->hammingDistance);
		listOfMatches->iterator = listOfMatches->iterator->next;
	}


	sp_km_freeListOfMatches(listOfMatches);
}
