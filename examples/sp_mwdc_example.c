/*
 * sp_mwdc_example.c
 *
 *  Created on: 5 Sep 2010
 *      Author: ben
 *  Example use of matching with don't cares
 */

#include "../include/stringpedia.h"

int main(void){
	char *text = "she sells sea shells on the sea shore";
	char *pattern = "s??l";

	struct SP_MWDC_MATCHING_POSITIONS *listOfMatches = sp_mwdc_create_new_list_of_matches();
	int numMatches = 0;


	sp_mwdc_match_with_fftw(text,pattern,'?',37,4,&numMatches,listOfMatches,SP_MWDC_NLOGM);

	printf("Printing start index of matches: \n");

	/* An 'iterator' as well as pointers to the start and end of the linked list is provided, so that you may use  *
	 * traverse the linked list, then reset the iterator to the start if you so wish. Start should not be modified *
	 * otherwise the list cannot be freed properly                                                                 */

	while(listOfMatches->iterator != NULL){
		printf("    Match at index: %d\n",listOfMatches->iterator->i);
		listOfMatches->iterator = listOfMatches->iterator->next;
	}


	sp_mwdc_freeListOfMatches(listOfMatches);
}
