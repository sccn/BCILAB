/************************************************************************/
/*                                                                      */
/*   mem_clean.h                                                        */
/*                                                                      */
/*   Header for malloc block tracker.  Uses a hash array of linked      */
/*   to store the allocated blocks.  Worst-case runtime is O(n) for     */
/*   any operation.  However, this would be a very degenerate case.     */
/*   Average runtime is log(n) / log(b), where b is the prime base      */
/*                                                                      */
/*   Author: Tom Briggs                                                 */
/*   Date: January 1, 2005                                              */
/*                                                                      */
/************************************************************************/

#ifndef _MEM_CLEAN_H
#define _MEM_CLEAN_H

/* defined to be a sufficiently large 
 * prime integer - this integer provides an efficient
 * store for about 30,000 malloc blocks, and probably
 * does not need to be increased. If it is changed
 * make sure the value is prime to ensure that the 
 * additive and multiplicative properties of a finite
 * field/group are maintained.
 */
#define HASH_NBLOCKS 2053


/* linked list node */
typedef struct _list_node {
  void *ptr;
  void *next;
} LIST_NODE;

/* hash block */
typedef struct _hash_block {
  LIST_NODE *root;  /* root of singly linked list */
  int depth;        /* number of elements stored at this loc. */
} HASH_BLOCK;

/* hash_array */
typedef struct _hash_array {
  HASH_BLOCK blocks[HASH_NBLOCKS];  /* array of blocks */
  int n;    /* total number of elements stored in array */
} HASH_ARRAY;


/* hash_arrray functions */
void hash_init_array(HASH_ARRAY *array);
void hash_add(HASH_ARRAY *array, void *ptr);
int hash_delete(HASH_ARRAY *array, void *ptr);
void hash_destroy_array(HASH_ARRAY *array);

#endif
