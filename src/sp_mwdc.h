/*
 * SP_MWDC.h
 *
 *  Created on: 1 Sep 2010
 *      Author: Ben Smithers (bs8959 AT bris DOT ac DOT uk)
 */

#ifndef SP_MWDC_H_
#define SP_MWDC_H_

/*Flags that apply to the all match-with-dont-cares problems */
#define SP_MWDC_FIRST_MATCH_ONLY (1U << 0) //Stop after finding the first match

struct SP_MWDC_INT_LIST{
	struct SP_MWDC_INT_LIST *next;
	int i;
};

struct SP_MWDC_MATCHING_POSITIONS{
	struct SP_MWDC_INT_LIST *start;
	struct SP_MWDC_INT_LIST *iterator;
	struct SP_MWDC_INT_LIST *end;
};

struct SP_MWDC_MATCHING_POSITIONS *sp_mwdc_create_new_list_of_matches();

struct SP_MWDC_INT_LIST *sp_mwdc_newListItem(int i);

void sp_mwdc_addToListOfMatches(struct SP_MWDC_MATCHING_POSITIONS *listOfMatches,int i);

void sp_mwdc_freeListOfMatches(struct SP_MWDC_MATCHING_POSITIONS *listOfMatches);

unsigned roundUpToPowerOf2(unsigned n);

#include "sp_mwdc_fft_matcher.h"
#include "sp_mwdc_fft_random_matcher.h"
#include "sp_mwdc_flint_matcher.h"
#include "sp_mwdc_naive_matcher.h"


#endif /* SP_MWDC_H_ */
