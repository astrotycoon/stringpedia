/*
 * sp_mwdc.c
 *
 *  Created on: 1 Sep 2010
 *      Author: Ben Smithers (bs8959 AT brist DOT ac DOT uk)
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "sp_mwdc.h"

struct SP_MWDC_MATCHING_POSITIONS *sp_mwdc_create_new_list_of_matches(){
	struct SP_MWDC_MATCHING_POSITIONS *list = malloc(sizeof(struct SP_MWDC_MATCHING_POSITIONS));
	list->start = list->iterator = list->end = NULL;
	return list;
}

struct SP_MWDC_INT_LIST *sp_mwdc_newListItem(int i){
	struct SP_MWDC_INT_LIST *listItem = malloc(sizeof(struct SP_MWDC_INT_LIST));
	listItem->next = NULL;
	listItem->i = i;
	return listItem;
}

void sp_mwdc_addToListOfMatches(struct SP_MWDC_MATCHING_POSITIONS *listOfMatches,int i){
	struct SP_MWDC_INT_LIST *listItem = sp_mwdc_newListItem(i);
	if(listOfMatches->start == NULL){
		listOfMatches->start = listOfMatches->end = listOfMatches->iterator = listItem;
	}else{
		listOfMatches->end->next = listItem;
		listOfMatches->end = listItem;
	}
}

void sp_mwdc_freeListOfMatches(struct SP_MWDC_MATCHING_POSITIONS *listOfMatches){
	struct SP_MWDC_INT_LIST *next;
	while(listOfMatches->start !=NULL){
		next = listOfMatches->start->next;
		free(listOfMatches->start);
		listOfMatches->start = next;
	}
	free(listOfMatches);
}

unsigned roundUpToPowerOf2(unsigned n){
	if(n==0){
		return 1;
	}
	n--;
	int bits = sizeof(unsigned) * CHAR_BIT;
	int i;
	for(i=1;i<bits;i<<=1){
		n |= (n >> i);
	}
	n++;
	return n;
}
