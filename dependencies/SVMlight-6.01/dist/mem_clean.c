/************************************************************************/
/*                                                                      */
/*   mexcommon.c                                                        */
/*                                                                      */
/*   Common definitions and functions in support of the MATLAB/MEX      */
/*                                                                      */
/*   Author: Tom Briggs                                                 */
/*   Date: January 1, 2005                                              */
/*                                                                      */
/*                                                                      */
/* Based on SVMLite:                                                    */
/*   Copyright (c) 2004  Thorsten Joachims - All rights reserved        */
/*                                                                      */
/*   This software is available for non-commercial use only. It must    */
/*   not be modified and distributed without prior permission of the    */
/*   author. The author is not responsible for implications from the    */
/*   use of this software.                                              */
/*                                                                      */
/*   January 1, 2005 - Modifications by Tom Briggs for MATLAB           */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mem_clean.h"

#ifdef MATLAB_MEX
#include "mex.h"

#define MC_MALLOC mxMalloc
#define MC_FREE mxFree

#else

#define MC_MALLOC malloc
#define MC_FREE free

#endif

/**
 * hash_init_array() - initialize the array 
 */
void hash_init_array(HASH_ARRAY *array)
{
	int i = 0;

	if (array == NULL) return;


	for (i = 0; i < HASH_NBLOCKS; i++) {
		HASH_BLOCK *block = &((*array).blocks[i]);

		block->depth = 0;
		block->root = NULL;

		array->n = 0;
	}

}

/** 
 * compute the hash bucket / block number that the given ptr 
 * will be located in.  Generally, this will be p % b, where p
 * is the integer numerical address, and b is the number of blocks.
 */
unsigned int hash_get_bucket(void *ptr)
{

	unsigned int number = (long ) ptr;
	unsigned int bucket = number % HASH_NBLOCKS;

	return bucket;
}

/**
 * hash_add() - add the given ptr to the given hash_array  */
void hash_add(HASH_ARRAY *array, void *ptr)
{

	int bucketNum = hash_get_bucket(ptr);

	HASH_BLOCK *bucket = &(array->blocks[bucketNum]);

	if (bucket->depth == 0)
	{
		bucket->root = MC_MALLOC(sizeof(LIST_NODE));
		bucket->root->ptr = ptr;
		bucket->root->next = NULL;
		bucket->depth++;
		array->n++;
	}
	else
	{
		LIST_NODE *bktptr = bucket->root;
		while (bktptr->next != NULL)
		{
			bktptr = bktptr->next;
		}
		bktptr->next = MC_MALLOC(sizeof(LIST_NODE));
		bktptr = bktptr->next;
		bktptr->ptr = ptr;
		bktptr->next = NULL;

		bucket->depth++;
		array->n++;
	}

}


/**
 * hash_delete - delete the ptr from the given array
 */
int hash_delete(HASH_ARRAY *array, void *ptr)
{

	int bucketNum = hash_get_bucket(ptr);

	HASH_BLOCK *bucket = &(array->blocks[bucketNum]);

	if (bucket->depth == 0)
	{
		return 0;
	}
	else 
	{
		LIST_NODE *bktptr = bucket->root;
		LIST_NODE *prev = NULL;
		if (bktptr->ptr == ptr) {	/* delete at root */
			LIST_NODE *oldroot = bucket->root;
			bucket->root = bucket->root->next;
			bucket->depth--;
			array->n--;
			MC_FREE(oldroot);
			return 1;
		} else {

			while((bktptr != ptr) && (bktptr != NULL)) {
				prev = bktptr;
				bktptr = bktptr->next;
			}

			if (bktptr == ptr) {
				LIST_NODE *oldnode = bktptr;
				prev->next = bktptr->next;
				MC_FREE(oldnode);
				bucket->depth--;
				array->n--;
				return 1;
			}
			else {
				/* fprintf(stderr,"Warning: could not delete ptr %X\n", (int) ptr); */
				return 0;
			}
		}
	}       
	return 0;
}


/**
 * retrieve a pointer to the indicated hash bucket */
void *hash_getptr(HASH_ARRAY *array, int bucketNum)
{
	HASH_BLOCK *bucket;

	if (array == NULL) 
		return NULL;

	if ((bucketNum < 0) || (bucketNum >= HASH_NBLOCKS)) 
		return NULL;


	bucket = &(array->blocks[bucketNum]);

	if (bucket == NULL)
		return NULL;

	if (bucket->depth == 0)
		return NULL;

	return bucket->root;
}

/**
 * retrieve the number of elements in the entire array */
int hash_nitems(HASH_ARRAY *array)
{
	if (array == NULL) 
		return -1;

	return array->n;
}

/**
 * retrieve the number of elements in the given hash bucket. */
int hash_get_depth(HASH_ARRAY *array, int bucketNum)
{
	HASH_BLOCK *bucket;

	if (array == NULL) 
		return -1;

	if ((bucketNum < 0) || (bucketNum >= HASH_NBLOCKS)) 
		return -1;

	bucket = &(array->blocks[bucketNum]);
	if (bucket == NULL) 
		return -1;

	return bucket->depth;
}

/**
 * destroy the hash array, freeing each allocated block as it continues.
 * this does not need to be a deep operation.  When this method is complete
 * all memory allocated by a malloc with a corresponding hash_add() will
 * be cleared by this method (if it hasn't been cleared before with a 
 * free(). */
void hash_destroy_array(HASH_ARRAY *array)
{
	int i = 0;

	if (array == NULL)
	{
		return ;
	}

	for(i = 0; i < HASH_NBLOCKS; i++)
	{
		HASH_BLOCK *bucket = &(array->blocks[i]);
		if ((bucket != NULL) && (bucket->depth != 0)) {
			while(bucket->root != NULL) {
				void *oldptr = bucket->root->ptr;
				hash_delete(array, oldptr);
				MC_FREE(oldptr);
			}
		}
	}
}
