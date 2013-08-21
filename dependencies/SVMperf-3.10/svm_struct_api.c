/***********************************************************************/
/*                                                                     */
/*   svm_struct_api.c (instantiated for SVM-perform)                   */
/*                                                                     */
/*   Definition of API for attaching implementing SVM learning of      */
/*   structures (e.g. parsing, multi-label classification, HMM)        */ 
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 31.10.05                                                    */
/*                                                                     */
/*   Copyright (c) 2005  Thorsten Joachims - All rights reserved       */
/*                                                                     */
/*   This software is available for non-commercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.                                             */
/*                                                                     */
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "svm_struct_api.h"
#include "svm_common.h"
#include "svm_struct_common.h"
#include "svm_struct_learn.h"

#define MAX(x,y)      ((x) < (y) ? (y) : (x))
#define MIN(x,y)      ((x) > (y) ? (y) : (x))
#define SIGN(x)       ((x) > (0) ? (1) : (((x) < (0) ? (-1) : (0))))

int compareup(const void *a, const void *b) 
{
  double va,vb;
  va=((STRUCT_ID_SCORE *)a)->score;
  vb=((STRUCT_ID_SCORE *)b)->score;
  if(va == vb) {
    va=((STRUCT_ID_SCORE *)a)->tiebreak;
    vb=((STRUCT_ID_SCORE *)b)->tiebreak;
  }
  return((va > vb) - (va < vb));
}

int comparedown(const void *a, const void *b) 
{
  return(-compareup(a,b));
}

LABEL       find_most_violated_constraint_errorrate(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm,
						     int loss_type);
LABEL       find_most_violated_constraint_thresholdmetric(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm,
						     int loss_type);
LABEL       find_most_violated_constraint_rankmetric(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm,
						     int loss_type);
LABEL       find_most_violated_constraint_avgprec(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm,
						     int loss_type);

double zeroone(int a, int b, int c, int d);
double fone(int a, int b, int c, int d);
double errorrate(int a, int b, int c, int d);
double prec(int a, int b, int c, int d);
double rec(int a, int b, int c, int d);
double swappedpairs(LABEL y, LABEL ybar);
double rocarea(LABEL y, LABEL ybar);
double prbep(LABEL y, LABEL ybar);
double avgprec(LABEL y, LABEL ybar);

double zeroone_loss(int a, int b, int c, int d);
double fone_loss(int a, int b, int c, int d);
double errorrate_loss(int a, int b, int c, int d);
double prbep_loss(int a, int b, int c, int d);
double prec_k_loss(int a, int b, int c, int d);
double rec_k_loss(int a, int b, int c, int d);
double swappedpairs_loss(LABEL y, LABEL ybar);
double avgprec_loss(LABEL y, LABEL ybar);

void        svm_struct_learn_api_init(int argc, char* argv[])
{
  /* Called in learning part before anything else is done to allow
     any initializations that might be necessary. */
}

void        svm_struct_learn_api_exit()
{
  /* Called in learning part at the very end to allow any clean-up
     that might be necessary. */
}

void        svm_struct_classify_api_init(int argc, char* argv[])
{
  /* Called in prediction part before anything else is done to allow
     any initializations that might be necessary. */
}

void        svm_struct_classify_api_exit()
{
  /* Called in prediction part at the very end to allow any clean-up
     that might be necessary. */
}

#ifdef COMPILE_MEX_INTERFACE
SAMPLE      read_struct_examples(const mxArray *file[], STRUCT_LEARN_PARM *sparm)
#else
SAMPLE      read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm)
#endif
{
  /* Reads struct examples and returns them in sample. The number of
     examples must be written into sample.n */
  SAMPLE   sample;  /* sample */
  EXAMPLE  *examples;
  long     n;       /* number of examples */
  long     totwords, maxlength=0, length, i, j, nump=0, numn=0;
  WORD     *words,*w;

  /* we have only one big example */
  examples=(EXAMPLE *)my_malloc(sizeof(EXAMPLE));
  /* Using the read_documents function from SVM-light */
  read_documents(file,&examples[0].x.doc,&examples[0].y.class,&totwords,&n);
  examples[0].x.totdoc=n;
  examples[0].y.totdoc=n;
  sample.n=1;
  sample.examples=examples;

  if(sparm->preimage_method==9) {
    for(i=0;i<n;i++) {
      examples[0].x.doc[i]->fvec->next=copy_svector(examples[0].x.doc[i]->fvec);
      examples[0].x.doc[i]->fvec->kernel_id=0;
      examples[0].x.doc[i]->fvec->next->kernel_id=2;
    }
  }

  for(i=0;i<sample.examples[0].x.totdoc;i++) {
    length=1;
    if(sample.examples[0].y.class[i]>0) 
      nump++;
    else 
      numn++;
    for(w=sample.examples[0].x.doc[i]->fvec->words;w->wnum;w++) {
      length++;
      if(length > maxlength) /* find vector with max elements */
	maxlength=length;
    }
  }

  /* add feature for bias if necessary */
  /* WARNING: Currently this works correctly only for linear kernel! */
  sparm->bias_featurenum=0;
  if(sparm->bias != 0) {
    words = (WORD *)my_malloc(sizeof(WORD)*(maxlength+1));
    sparm->bias_featurenum=totwords+1;
    totwords++;
    for(i=0;i<sample.examples[0].x.totdoc;i++) {
      for(j=0;sample.examples[0].x.doc[i]->fvec->words[j].wnum;j++) 
	words[j]=sample.examples[0].x.doc[i]->fvec->words[j];
      words[j].wnum=sparm->bias_featurenum; /* bias */
      words[j].weight=sparm->bias;
      j++;
      words[j].wnum=0; /* end */
      words[j].weight=0;
      free_svector(sample.examples[0].x.doc[i]->fvec);
      sample.examples[0].x.doc[i]->fvec=create_svector(words,"",1.0);      
    }
    free(words);
  }

  /* Remove all features with numbers larger than num_features, if
     num_features is set to a positive value. This is important for
     svm_struct_classify. */
  if((sparm->num_features > 0) && sparm->truncate_fvec)
    for(i=0;i<sample.examples[0].x.totdoc;i++)
      for(j=0;sample.examples[0].x.doc[i]->fvec->words[j].wnum;j++) 
	if(sample.examples[0].x.doc[i]->fvec->words[j].wnum>sparm->num_features) {
	  sample.examples[0].x.doc[i]->fvec->words[j].wnum=0;
	  sample.examples[0].x.doc[i]->fvec->words[j].weight=0;
	}

  /* change label value for better scaling using thresholdmetrics */
  if((sparm->loss_function == ZEROONE) 
     || (sparm->loss_function == FONE) 
     || (sparm->loss_function == ERRORRATE)
     || (sparm->loss_function == PRBEP) 
     || (sparm->loss_function == PREC_K) 
     || (sparm->loss_function == REC_K)) {
    for(i=0;i<sample.examples[0].x.totdoc;i++) {
      if(sample.examples[0].y.class[i]>0)
	sample.examples[0].y.class[i]=0.5*100.0/(numn+nump);
      else
	sample.examples[0].y.class[i]=-0.5*100.0/(numn+nump);
    }
  }
  /* change label value for easy computation of rankmetrics (i.e. ROC-area) */
  if(sparm->loss_function == SWAPPEDPAIRS) {
    for(i=0;i<sample.examples[0].x.totdoc;i++) {
      /*      if(sample.examples[0].y.class[i]>0)
	sample.examples[0].y.class[i]=numn*0.5;
      else
      sample.examples[0].y.class[i]=-nump*0.5; */
      if(sample.examples[0].y.class[i]>0)
	sample.examples[0].y.class[i]=0.5*100.0/nump;
      else
	sample.examples[0].y.class[i]=-0.5*100.0/numn;
    }
  }
  if(sparm->loss_function == AVGPREC) {
    for(i=0;i<sample.examples[0].x.totdoc;i++) {
      if(sample.examples[0].y.class[i]>0)
	sample.examples[0].y.class[i]=numn;
      else
	sample.examples[0].y.class[i]=-nump;
    }
  }

  return(sample);
}

void        init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
			      STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
			      KERNEL_PARM *kparm)
{
  /* Initialize structmodel sm. The weight vector w does not need to be
     initialized, but you need to provide the maximum size of the
     feature space in sizePsi. This is the maximum number of different
     weights that can be learned. Later, the weight vector w will
     contain the learned weights for the model. */
  long   i,j,k,totwords=0,totdoc=0,totexp=0,nump=0,numn=0,new_size,*select;
  WORD   *words,*w;
  DOC    **orgdoc,**basis;
  double *dummy;
  MATRIX *G,*L;
  double *indep,ii,weight;

  totdoc=sample.examples[0].x.totdoc;
  if(sparm->sparse_kernel_method > 0) { /* use nystrom or inc cholesky */
    sparm->sparse_kernel_type=kparm->kernel_type;
  }
  else
    sparm->sparse_kernel_type=0;
  sm->sparse_kernel_type=sparm->sparse_kernel_type;
  sm->invL=NULL; 
  sm->expansion=NULL;

  /* When using sparse kernel approximation, this replaces the
     original feature vector with the kernel values of the original
     feature vector and the expansion. */
  if(sm->sparse_kernel_type > 0) {
    kparm->kernel_type=sm->sparse_kernel_type;

    if(strcmp(sparm->sparse_kernel_file,"")) {
      if(struct_verbosity > 0) 
	printf("Reading basis functions for sparse kernel expansion from '%s'...",sparm->sparse_kernel_file); fflush(stdout);
#ifdef COMPILE_MEX_INTERFACE
			mexErrMsgTxt("Sorry, basis functions from a file is not implemented in MEX interface...");
#else
      read_documents(sparm->sparse_kernel_file,&basis,&dummy,&totwords,&totexp);
      free(dummy);
#endif
      if(struct_verbosity > 0) 
	printf("done.\n");
    }
    else {
      basis=(DOC **)malloc(sizeof(DOC *)*totdoc);
      for(i=0;i<totdoc;i++) {
	basis[i]=(DOC *)malloc(sizeof(DOC));
	(*(basis[i]))=(*(sample.examples[0].x.doc[i]));
        basis[i]->fvec=copy_svector(sample.examples[0].x.doc[i]->fvec);
      }
      totexp=totdoc;
    }
    if(sparm->sparse_kernel_size > totexp) sparm->sparse_kernel_size=totexp;

    /* determine basis functions to use in expansion: B */
    if(sparm->sparse_kernel_method==1) { 
      /* select expansion via random sampling */
      if(struct_verbosity > 0) 
	printf("Selecting random sample of basis functions..."); fflush(stdout);
      sm->expansion=(DOC **)malloc(sizeof(DOC *)*totexp);
      sm->expansion_size=0;
      for(ii=0.5;ii<totexp;ii+=((double)totexp/sparm->sparse_kernel_size)) {
	sm->expansion[sm->expansion_size]=(DOC *)malloc(sizeof(DOC));
	(*(sm->expansion[sm->expansion_size]))=(*(basis[(long)ii]));
	sm->expansion[sm->expansion_size]->fvec=copy_svector(basis[(long)ii]->fvec);
	sm->expansion_size++;
      }
      if(struct_verbosity > 0) 
	printf("done.\n");

      /* Make sure they are all independent. If not, select independent
	 subset. */
      if(struct_verbosity > 0) 
	printf("Finding independent subset of vectors..."); fflush(stdout);
      G=create_matrix(sm->expansion_size,sm->expansion_size);
      for(i=0;i<sm->expansion_size;i++) { /* need only upper triangle */
	for(j=i;j<sm->expansion_size;j++) {
	  G->element[i][j]=kernel(kparm,sm->expansion[i],sm->expansion[j]);
	}
      }
      indep=find_indep_subset_of_matrix(G,0.000001);
      free_matrix(G);
      new_size=0;
      for(i=0;i<sm->expansion_size;i++) {
	if(indep[i] != 0) {
	  sm->expansion[new_size]=sm->expansion[i];
	  new_size++;
	}
	else {
	  free_example(sm->expansion[i],1);
	}
      }
      free_nvector(indep);
      if(struct_verbosity>0) 
	printf("found %ld of %ld...",new_size,sm->expansion_size);
      sm->expansion_size=new_size;
      if(struct_verbosity > 0) 
	printf("done.\n");
      
      /* compute matrix B^T*B */
      if(struct_verbosity > 0) 
	printf("Computing Gram matrix for kernel expansion...");fflush(stdout);
      G=create_matrix(sm->expansion_size,sm->expansion_size);
      for(i=0;i<sm->expansion_size;i++) {/* need upper triangle for cholesky */
	for(j=i;j<sm->expansion_size;j++) {
	  G->element[i][j]=kernel(kparm,sm->expansion[i],sm->expansion[j]);
	}
      }
      if(struct_verbosity > 0) 
	printf("done.\n");
      if(struct_verbosity > 0) 
	printf("Computing Cholesky decomposition and inverting..."); fflush(stdout);
      L=cholesky_matrix(G);
      sm->invL=invert_ltriangle_matrix(L);
      free_matrix(L);
      free_matrix(G);
      if(struct_verbosity > 0) 
	printf("done.\n");
    }
    else if(sparm->sparse_kernel_method==2) { 
      /* select expansion via incomplete cholesky */
      if(struct_verbosity > 0) 
	printf("Computing incomplete Cholesky decomposition..."); fflush(stdout);
      L=incomplete_cholesky(basis,totexp,sparm->sparse_kernel_size,
				 0.000001,kparm,&select);
      sm->invL=invert_ltriangle_matrix(L);
      free_matrix(L);
      sm->expansion=(DOC **)malloc(sizeof(DOC *)*totexp);
      sm->expansion_size=0;
      for(i=0;select[i]>=0;i++) {
	sm->expansion[sm->expansion_size]=(DOC *)malloc(sizeof(DOC));
	(*(sm->expansion[sm->expansion_size]))=(*(basis[select[i]]));
	sm->expansion[sm->expansion_size]->fvec=copy_svector(basis[select[i]]->fvec);
	sm->expansion_size++;
      }
      if(struct_verbosity > 0) 
	printf("done.\n");
    }
    for(i=0;i<totexp;i++) 
      free_example(basis[i],1);
    free(basis);

    /* compute new features for each example x: B^T*x */
    if(struct_verbosity > 0) 
      printf("Replacing feature vectors..."); fflush(stdout);
    orgdoc=sample.examples[0].x.doc;
    sample.examples[0].x.doc=(DOC **)malloc(sizeof(DOC *)*totdoc);
    words=(WORD *)malloc(sizeof(WORD)*(sm->expansion_size+1));
    for(i=0;i<totdoc;i++) {
      k=0;
      for(j=0;j<sm->expansion_size;j++) {
	weight=kernel(kparm,orgdoc[i],sm->expansion[j]);
	if(weight != 0) {
	  words[k].wnum=j+1;
	  words[k].weight=weight;
	  k++;
	}
      }
      words[k].wnum=0;
      sample.examples[0].x.doc[i]=create_example(orgdoc[i]->docnum,
					  orgdoc[i]->queryid,
					  orgdoc[i]->slackid,
					  orgdoc[i]->costfactor,
					  create_svector(words,
					     orgdoc[i]->fvec->userdefined,1.0));
    }
    free(words);
    for(i=0;i<totdoc;i++) 
      free_example(orgdoc[i],1);
    free(orgdoc);
    kparm->kernel_type=LINEAR;
    if(struct_verbosity > 0) 
      printf("done.\n");
  }

  /* count number of positive and negative examples */
  for(i=0;i<sample.examples[0].x.totdoc;i++) {
    if(sample.examples[0].y.class[i]>0) 
      nump++;
    else 
      numn++;
  }

  totwords=0;
  for(i=0;i<totdoc;i++)  /* find highest feature number */
    for(w=sample.examples[0].x.doc[i]->fvec->words;w->wnum;w++) 
      if(totwords < w->wnum) 
	totwords=w->wnum;
  sparm->num_features=totwords;
  if(struct_verbosity>0)
    printf("Training set properties: %d features, %ld examples (%ld pos / %ld neg)\n",
	   sparm->num_features,totdoc,nump,numn);
  sm->sizePsi=sparm->num_features;
  if(struct_verbosity>=2)
    printf("Size of Phi: %ld\n",sm->sizePsi);

  /* fill in default value for k when using Prec@k or Rec@k */
  if(sparm->loss_function == PREC_K) {
    if(sparm->prec_rec_k_frac == 0) {
      sparm->prec_rec_k_frac=0.5;
    }
    else if(sparm->prec_rec_k_frac > 1) {
      printf("\nERROR: The value of option --k for Prec@k must not be larger than 1.0!\n\n");
      exit(0);
    }
  }
  if(sparm->loss_function == REC_K) {
    if(sparm->prec_rec_k_frac == 0) {
      sparm->prec_rec_k_frac=2;
    }
    else if(sparm->prec_rec_k_frac < 1) {
      printf("\nERROR: The value of option --k for Rec@k must not be smaller than 1.0!\n\n");
      exit(0);
    }
  }

  /* make sure that the bias feature is not used with kernels */
  if(((kparm->kernel_type != LINEAR) || (sparm->sparse_kernel_type != LINEAR))
     && (sparm->bias != 0)) {
    printf("\nThe value of option --b must be zero when kernels are used.\n");
    printf("The option to use a bias for non-linear kernels is not implemented yet!\n\n");
    exit(0);
  }
}

CONSTSET    init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm)
{
  /* Initializes the optimization problem. Typically, you do not need
     to change this function, since you want to start with an empty
     set of constraints. However, if for example you have constraints
     that certain weights need to be positive, you might put that in
     here. The constraints are represented as lhs[i]*w >= rhs[i]. lhs
     is an array of feature vectors, rhs is an array of doubles. m is
     the number of constraints. The function returns the initial
     set of constraints. */
  CONSTSET c;
  long     sizePsi=sm->sizePsi;
  long     i;
  WORD     words[2];

  if(1) { /* normal case: start with empty set of constraints */
    c.lhs=NULL;
    c.rhs=NULL;
    c.m=0;
  }
  else { /* add constraints so that all learned weights are
            positive. WARNING: Currently, they are positive only up to
            precision epsilon set by -e. */
    c.lhs=my_malloc(sizeof(DOC *)*sizePsi);
    c.rhs=my_malloc(sizeof(double)*sizePsi);
    for(i=0; i<sizePsi; i++) {
      words[0].wnum=i+1;
      words[0].weight=1.0;
      words[1].wnum=0;
      /* the following slackid is a hack. we will run into problems,
         if we have move than 1000000 slack sets (i.e. examples) */
      c.lhs[i]=create_example(i,0,1000000+i,1,create_svector(words,"",1.0));
      c.rhs[i]=0.0;
    }
  }
  return(c);
}

LABEL       classify_struct_example(PATTERN x, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label yhat for pattern x that scores the highest
     according to the linear evaluation function in sm, especially the
     weights sm.w. The returned label is taken as the prediction of sm
     for the pattern x. The weights correspond to the features defined
     by psi() and range from index 1 to index sm->sizePsi. If the
     function cannot find a label, it shall return an empty label as
     recognized by the function empty_label(y). */
  LABEL y;
  int i;

  y.totdoc=x.totdoc;
  y.class=(double *)my_malloc(sizeof(double)*y.totdoc);
  /* simply classify by sign of inner product between example vector
     and weight vector */
  for(i=0;i<x.totdoc;i++) {
    y.class[i]=classify_example(sm->svm_model,x.doc[i]);
  }
  return(y);
}

LABEL       find_most_violated_constraint_slackrescaling(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label ybar for pattern x that that is responsible for
     the most violated constraint for the slack rescaling
     formulation. It has to take into account the scoring function in
     sm, especially the weights sm.w, as well as the loss
     function. The weights in sm.w correspond to the features defined
     by psi() and range from index 1 to index sm->sizePsi. Most simple
     is the case of the zero/one loss function. For the zero/one loss,
     this function should return the highest scoring label ybar, if
     ybar is unequal y; if it is equal to the correct label y, then
     the function shall return the second highest scoring label. If
     the function cannot find a label, it shall return an empty label
     as recognized by the function empty_label(y). */
  LABEL ybar;
  if((sparm->loss_function == ZEROONE) 
     || (sparm->loss_function == FONE) 
     || (sparm->loss_function == ERRORRATE)
     || (sparm->loss_function == PRBEP) 
     || (sparm->loss_function == PREC_K) 
     || (sparm->loss_function == REC_K)) {
    ybar=find_most_violated_constraint_thresholdmetric(x,y,sm,sparm,
						       sparm->loss_type);
  }
  else if((sparm->loss_function == SWAPPEDPAIRS)) {
    printf("ERROR: Slack-rescaling is not implemented for this loss function!\n");
    exit(1);
  }
  else {
    printf("ERROR: Unknown loss function '%d'.\n",sparm->loss_function);
    exit(1);
  }
  return(ybar);
}

LABEL       find_most_violated_constraint_marginrescaling(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label ybar for pattern x that that is responsible for
     the most violated constraint for the margin rescaling
     formulation. It has to take into account the scoring function in
     sm, especially the weights sm.w, as well as the loss
     function. The weights in sm.w correspond to the features defined
     by psi() and range from index 1 to index sm->sizePsi. Most simple
     is the case of the zero/one loss function. For the zero/one loss,
     this function should return the highest scoring label ybar, if
     ybar is unequal y; if it is equal to the correct label y, then
     the function shall return the second highest scoring label. If
     the function cannot find a label, it shall return an empty label
     as recognized by the function empty_label(y). */
  LABEL ybar;
  if(sparm->loss_function == ERRORRATE) 
    ybar=find_most_violated_constraint_errorrate(x,y,sm,sparm,
	  sparm->loss_type); 
  /* ybar=find_most_violated_constraint_thresholdmetric(x,y,sm,sparm,
     sparm->loss_type); */
  else if((sparm->loss_function == ZEROONE) 
     || (sparm->loss_function == FONE) 
     || (sparm->loss_function == PRBEP) 
     || (sparm->loss_function == PREC_K) 
     || (sparm->loss_function == REC_K)) 
    ybar=find_most_violated_constraint_thresholdmetric(x,y,sm,sparm,
						       sparm->loss_type);
  else if((sparm->loss_function == SWAPPEDPAIRS))
    ybar=find_most_violated_constraint_rankmetric(x,y,sm,sparm,
						  sparm->loss_type);
  else if((sparm->loss_function == AVGPREC))
    ybar=find_most_violated_constraint_avgprec(x,y,sm,sparm,
					       sparm->loss_type);
  else {
    printf("ERROR: Unknown loss function '%d'.\n",sparm->loss_function);
    exit(1);
  }
  return(ybar);
}

LABEL       find_most_violated_constraint_errorrate(PATTERN x, LABEL y, 
						    STRUCTMODEL *sm, 
						    STRUCT_LEARN_PARM *sparm,
						    int loss_type)
{
  /* Finds the most violated constraint for errorrate under
     margin-rescaling!!! */
  LABEL ybar;
  int i,totwords;
  double *score,valmax=0,diff,loss_step,*ortho_weights;
  MODEL svm_model;

  ybar.totdoc=x.totdoc;
  ybar.class=(double *)my_malloc(sizeof(double)*x.totdoc);
  score=(double *)my_malloc(sizeof(double)*(ybar.totdoc+1));

  totwords=sm->svm_model->totwords;
  svm_model=(*sm->svm_model);  

  /* For sparse kernel, replace weight vector with beta=gamma^T*L^-1 */
  if(sm->sparse_kernel_type>0) {
    svm_model.lin_weights=(double *)my_malloc(sizeof(double)*(totwords+1));
    ortho_weights=prod_nvector_ltmatrix(sm->svm_model->lin_weights+1,sm->invL);
    for(i=0;i<sm->invL->m;i++)
      svm_model.lin_weights[i+1]=ortho_weights[i];
    svm_model.lin_weights[0]=0;
    free_nvector(ortho_weights);
  }

  loss_step=100.0/x.totdoc;
  for(i=0;i<x.totdoc;i++) {
    score[i]=classify_example(&svm_model,x.doc[i]);
    diff=loss_step-2.0*y.class[i]*score[i];
    if(diff > 0) {
      ybar.class[i]=-y.class[i];
      valmax+=loss_step;
    }
    else {
      ybar.class[i]=y.class[i];
    }
    valmax+=ybar.class[i]*score[i];
  }

  /* Restore weight vector that was modified above */
  if(sm->sparse_kernel_type>0) {
    free(svm_model.lin_weights);
  }

  if(struct_verbosity >= 2) {
    printf("\n max_ybar {loss(y_i,ybar)+w*Psi(x,ybar)}=%f\n",valmax);
    SVECTOR *fy=psi(x,y,sm,sparm);
    SVECTOR *fybar=psi(x,ybar,sm,sparm);
    DOC *exy=create_example(0,0,1,1,fy);
    DOC *exybar=create_example(0,0,1,1,fybar);
    printf(" -> w*Psi(x,y_i)=%f, w*Psi(x,ybar)=%f\n",
	   classify_example(sm->svm_model,exy),
	   classify_example(sm->svm_model,exybar));
    free_example(exy,1);
    free_example(exybar,1);
  }

  free(score);

  return(ybar);
}
   


LABEL       find_most_violated_constraint_thresholdmetric(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm,
						     int loss_type)
{
  /* Finds the most violated constraint for metrics that are based on
     a threshold. */
  LABEL ybar;
  int i,nump,numn,start,prec_rec_k,totwords;
  double *score,*sump,*sumn;
  STRUCT_ID_SCORE *scorep,*scoren;
  int threshp=0,threshn=0;
  int a,d;
  double val,valmax,loss,score_y,*ortho_weights;
  MODEL svm_model;

  ybar.totdoc=x.totdoc;
  ybar.class=(double *)my_malloc(sizeof(double)*x.totdoc);
  score=(double *)my_malloc(sizeof(double)*(ybar.totdoc+1));
  scorep=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));
  scoren=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));
  sump=(double *)my_malloc(sizeof(double)*(ybar.totdoc+1));
  sumn=(double *)my_malloc(sizeof(double)*(ybar.totdoc+1));

  totwords=sm->svm_model->totwords;
  svm_model=(*sm->svm_model);  

  /* For sparse kernel, replace weight vector with beta=gamma^T*L^-1 */
  if(sm->sparse_kernel_type>0) {
    svm_model.lin_weights=(double *)my_malloc(sizeof(double)*(totwords+1));
    ortho_weights=prod_nvector_ltmatrix(sm->svm_model->lin_weights+1,sm->invL);
    for(i=0;i<sm->invL->m;i++)
      svm_model.lin_weights[i+1]=ortho_weights[i];
    svm_model.lin_weights[0]=0;
    free_nvector(ortho_weights);
  }

  nump=0;
  numn=0;
  for(i=0;i<x.totdoc;i++) {
    score[i]=fabs(y.class[i])*classify_example(&svm_model,x.doc[i]);
    if(y.class[i] > 0) {
      scorep[nump].score=score[i];
      scorep[nump].tiebreak=0;
      scorep[nump].id=i;
      nump++;
    }
    else {
      scoren[numn].score=score[i];
      scoren[numn].tiebreak=0;
      scoren[numn].id=i;
      numn++;
    }
  }

  /* Restore weight vector that was modified above */
  if(sm->sparse_kernel_type>0) {
    free(svm_model.lin_weights);
  }

  /* compute score of target label */
  score_y=0;
  if(loss_type==SLACK_RESCALING) {
    for(i=0;i<x.totdoc;i++) 
      /*      score_y+=y.class[i]*score[i]; */
      score_y+=score[i]; 
 }
 
  if(nump)
    qsort(scorep,nump,sizeof(STRUCT_ID_SCORE),comparedown);
  sump[0]=0;
  for(i=0;i<nump;i++) {
    sump[i+1]=sump[i]+scorep[i].score;
  }
  if(numn)
    qsort(scoren,numn,sizeof(STRUCT_ID_SCORE),compareup);
  sumn[0]=0;
  for(i=0;i<numn;i++) {
    sumn[i+1]=sumn[i]+scoren[i].score;
  }

  /* find max of loss(ybar,y)+score(ybar) for margin rescaling or max
     of loss(ybar,y)+loss*(score(ybar)-score(y)) for slack
     rescaling */
  valmax=0;
  start=1;
  prec_rec_k=(int)(nump*sparm->prec_rec_k_frac);
  if(prec_rec_k<1) prec_rec_k=1;
  for(a=0;a<=nump;a++) {
    for(d=0;d<=numn;d++) {
      if(sparm->loss_function == ZEROONE)
	loss=zeroone_loss(a,numn-d,nump-a,d);
      else if(sparm->loss_function == FONE)
	loss=fone_loss(a,numn-d,nump-a,d);
      else if(sparm->loss_function == ERRORRATE)
	loss=errorrate_loss(a,numn-d,nump-a,d);
      else if((sparm->loss_function == PRBEP) && (a+numn-d == nump))
	loss=prbep_loss(a,numn-d,nump-a,d);
      else if((sparm->loss_function == PREC_K) && (a+numn-d >= prec_rec_k))
	loss=prec_k_loss(a,numn-d,nump-a,d);
      else if((sparm->loss_function == REC_K) && (a+numn-d <= prec_rec_k)) 
	loss=rec_k_loss(a,numn-d,nump-a,d);
      else {
	loss=0;
      }
      if(loss > 0) {
	if(loss_type==SLACK_RESCALING) {
	  val=loss+loss*(sump[a]-(sump[nump]-sump[a])-sumn[d]+(sumn[numn]-sumn[d] - score_y));
	}
	else if(loss_type==MARGIN_RESCALING) {
	  val=loss+sump[a]-(sump[nump]-sump[a])-sumn[d]+(sumn[numn]-sumn[d]);
	}
	else {
	  printf("ERROR: Unknown loss type '%d'.\n",loss_type);
	  exit(1);
	}
	if((val > valmax) || (start)) {
	  start=0;
	  valmax=val;
	  threshp=a;
	  threshn=d;
	}
      }
    }
  }

  /* assign labels that maximize score */
  /*  for(i=0;i<nump;i++) {
    if(i<threshp) 
      ybar.class[scorep[i].id]=1;
    else 
      ybar.class[scorep[i].id]=-1;
  }
  for(i=0;i<numn;i++) {
    if(i<threshn) 
      ybar.class[scoren[i].id]=-1;
    else 
      ybar.class[scoren[i].id]=1;
      } */
  for(i=0;i<nump;i++) {
    if(i<threshp) 
      ybar.class[scorep[i].id]=y.class[scorep[i].id];
    else 
      ybar.class[scorep[i].id]=-y.class[scorep[i].id];
  }
  for(i=0;i<numn;i++) {
    if(i<threshn) 
      ybar.class[scoren[i].id]=y.class[scoren[i].id];
    else 
      ybar.class[scoren[i].id]=-y.class[scoren[i].id];
  }

  if(struct_verbosity >= 2) {
    if(loss_type==SLACK_RESCALING) 
      printf("\n max_ybar {loss(y_i,ybar)+loss(y_i,ybar)[w*Psi(x,ybar)-w*Psi(x,y)]}=%f\n",valmax);
    else
      printf("\n max_ybar {loss(y_i,ybar)+w*Psi(x,ybar)}=%f\n",valmax);
    SVECTOR *fy=psi(x,y,sm,sparm);
    SVECTOR *fybar=psi(x,ybar,sm,sparm);
    DOC *exy=create_example(0,0,1,1,fy);
    DOC *exybar=create_example(0,0,1,1,fybar);
    printf(" -> w*Psi(x,y_i)=%f, w*Psi(x,ybar)=%f\n",
	   classify_example(sm->svm_model,exy),
	   classify_example(sm->svm_model,exybar));
    free_example(exy,1);
    free_example(exybar,1);
  }
  free(score);
  free(scorep);
  free(scoren);
  free(sump);
  free(sumn);

  return(ybar);
}

LABEL       find_most_violated_constraint_rankmetric(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm,
						     int loss_type)
{
  /* Finds the most violated constraint for metrics that are based on
     a threshold. 
     WARNING: Currently only implemented for margin-rescaling!!! */
  LABEL ybar;
  long i,nump,numn,sump,sumn;
  double *score,*ortho_weights;
  STRUCT_ID_SCORE *scorep,*scoren,*predset;
  int totwords;
  MODEL svm_model;

  ybar.totdoc=x.totdoc;
  ybar.class=(double *)my_malloc(sizeof(double)*x.totdoc);
  score=(double *)my_malloc(sizeof(double)*(ybar.totdoc+1));
  scorep=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));
  scoren=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));
  predset=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));
  totwords=sm->svm_model->totwords;
  svm_model=(*sm->svm_model);  

  /* For sparse kernel, replace weight vector with beta=gamma^T*L^-1 */
  if(sm->sparse_kernel_type>0) {
    svm_model.lin_weights=(double *)my_malloc(sizeof(double)*(totwords+1));
    ortho_weights=prod_nvector_ltmatrix(sm->svm_model->lin_weights+1,sm->invL);
    for(i=0;i<sm->invL->m;i++)
      svm_model.lin_weights[i+1]=ortho_weights[i];
    svm_model.lin_weights[0]=0;
    free_nvector(ortho_weights);
  }

  nump=0;
  numn=0;
  for(i=0;i<x.totdoc;i++) {
    score[i]=classify_example(&svm_model,x.doc[i]);
    if(y.class[i] > 0) {
      scorep[nump].score=score[i];
      scorep[nump].tiebreak=0;
      scorep[nump].id=i;
      nump++;
    }
    else {
      scoren[numn].score=score[i];
      scorep[numn].tiebreak=0;
      scoren[numn].id=i;
      numn++;
    }
  }
  if(nump)
    qsort(scorep,nump,sizeof(STRUCT_ID_SCORE),comparedown);
  if(numn)
    qsort(scoren,numn,sizeof(STRUCT_ID_SCORE),comparedown);

  /* Restore weight vector that was modified above */
  if(sm->sparse_kernel_type>0) {
    free(svm_model.lin_weights);
  }

  /* find max of loss(ybar,y)+score(ybar) */
  if(sparm->loss_function == SWAPPEDPAIRS) { /* number of swapped pairs (ie. ROC Area) */
    for(i=0;i<nump;i++) {
      predset[i]=scorep[i];
      predset[i].score-=(0.5); 
    }
    for(i=0;i<numn;i++) {
      predset[nump+i]=scoren[i];
      predset[nump+i].score+=(0.5); 
    }
    qsort(predset,nump+numn,sizeof(STRUCT_ID_SCORE),comparedown);
    sump=0;
    sumn=0;
    for(i=0;i<numn+nump;i++) {
      if(y.class[predset[i].id] > 0) {
	ybar.class[predset[i].id]=((numn-2.0*sumn)*0.5*100.0/nump)/numn;
	sump++;
      }
      else {
	ybar.class[predset[i].id]=-((nump-2.0*(nump-sump))*0.5*100.0/nump)/numn;
	sumn++;
      }
    }
  }
  else {
    printf("ERROR: Invalid loss function '%d'.\n",sparm->loss_function);
    exit(1);
  }

  if(struct_verbosity >= 1) {
    SVECTOR *fy=psi(x,y,sm,sparm);
    SVECTOR *fybar=psi(x,ybar,sm,sparm);
    DOC *exy=create_example(0,0,1,1,fy);
    DOC *exybar=create_example(0,0,1,1,fybar);
    printf("\n max_y (loss(y_i,y)+ w*Psi(x,y)=%f\n",
	   loss(y,ybar,sparm)+classify_example(sm->svm_model,exybar));
    printf("w*Psi(x,y)=%f, w*Psi(x,ybar)=%f\n",
	   classify_example(sm->svm_model,exy),
	   classify_example(sm->svm_model,exybar));
    free_example(exy,1);
    free_example(exybar,1);
  }
  free(score);
  free(scorep);
  free(scoren);
  free(predset);

  return(ybar);
}

double objective_avgprec(int *negabove,STRUCT_ID_SCORE *scorep,int nump,
			 STRUCT_ID_SCORE *scoren,int numn)
{
  double val;
  int i,j;
  val=0;
  for(i=0;i<nump;i++) { /* contribution from loss */
    val+=(double)(i+1)/(double)(i+1+negabove[i]);
    /*    printf("na[%d]=%d,",i,negabove[i]); */
  }
  val=100-100.0*val/(double)(nump);

  printf("argmaxloss=%f\n",val);

  for(i=0;i<nump;i++) { /* contribution from pos scores */
    val+=(numn-negabove[i])*scorep[i].score-negabove[i]*scorep[i].score;
    for(j=0;j<negabove[i];j++) { /* contribution from neg scores */
      val+=scoren[j].score;
    }
    for(j=negabove[i];j<numn;j++) { /* contribution from neg scores */
      val-=scoren[j].score;
    }
  }
  return(val);
}

LABEL       find_most_violated_constraint_avgprec(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm,
						     int loss_type)
{
  /* Finds the most violated constraint for metrics that are based on
     a threshold. 
     WARNING: Currently only implemented for margin-rescaling!!! */
  LABEL ybar;
  int i,ii,iii,j,nump,numn,*negabove,changed,totwords;
  double *score,v,vv,vmax,*ortho_weights;
  STRUCT_ID_SCORE *scorep,*scoren;
  MODEL svm_model;

  ybar.totdoc=x.totdoc;
  ybar.class=(double *)my_malloc(sizeof(double)*x.totdoc);
  score=(double *)my_malloc(sizeof(double)*(ybar.totdoc+1));
  scorep=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));
  scoren=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));

  totwords=sm->svm_model->totwords;
  svm_model=(*sm->svm_model);  

  /* For sparse kernel, replace weight vector with beta=gamma^T*L^-1 */
  if(sm->sparse_kernel_type>0) {
    svm_model.lin_weights=(double *)my_malloc(sizeof(double)*(totwords+1));
    ortho_weights=prod_nvector_ltmatrix(sm->svm_model->lin_weights+1,sm->invL);
    for(i=0;i<sm->invL->m;i++)
      svm_model.lin_weights[i+1]=ortho_weights[i];
    svm_model.lin_weights[0]=0;
    free_nvector(ortho_weights);
  }

  nump=0;
  numn=0;
  for(i=0;i<x.totdoc;i++) {
    score[i]=classify_example(&svm_model,x.doc[i]);
    if(y.class[i] > 0) {
      scorep[nump].score=score[i];
      scorep[nump].tiebreak=0;
      scorep[nump].id=i;
      nump++;
    }
    else {
      scoren[numn].score=score[i];
      scorep[numn].tiebreak=0;
      scoren[numn].id=i;
      numn++;
    }
  }
  if(nump)
    qsort(scorep,nump,sizeof(STRUCT_ID_SCORE),comparedown);
  if(numn)
    qsort(scoren,numn,sizeof(STRUCT_ID_SCORE),comparedown);

  /* Restore weight vector that was modified above */
  if(sm->sparse_kernel_type>0) {
    free(svm_model.lin_weights);
  }

  /* find max of loss(ybar,y)+score(ybar) */
  if(sparm->loss_function == AVGPREC) { /* Average Precision */
    negabove=(int *)my_malloc(sizeof(int)*(nump+1));
    v=0;
    for(i=0;i<nump;i++) {
      negabove[i]=0;
      v+=numn*scorep[i].score;
    }
    for(j=0;j<numn;j++) {
      v-=nump*scoren[j].score;
    }
    negabove[nump]=numn;
    vv=objective_avgprec(negabove,scorep,nump,scoren,numn);
    printf("v=%f == vv=%f\n",v,vv);
    vmax=v;
    changed=1;
    while(changed) {
      changed=0;
      for(i=nump-1;i>=0;i--) {
	if(negabove[i+1] > negabove[i]) {/* there is a neg between the 2 pos */
	  for(ii=i;ii>=0;ii--) { /* try to move neg up */
	    for(iii=i;iii>=ii;iii--) {
	      negabove[iii]++;
	      v=v+(-100.0/(double)(nump)*(double)(iii+1)/(double)(iii+1+negabove[iii])+100.0/(double)(nump)*(double)(iii+1)/(double)(iii+1+negabove[iii]-1)-2*(scorep[iii].score-scoren[negabove[iii]-1].score));
	    }
	    /* vv=objective_avgprec(negabove,scorep,nump,scoren,numn); 
  	       printf("v=%f == vv=%f\n",v,vv); */
	    if(v<=vmax) {
	      for(iii=i;iii>=ii;iii--) { /* undo step */
		v=v-(-100.0/(double)(nump)*(double)(iii+1)/(double)(iii+1+negabove[iii])+100.0/(double)(nump)*(double)(iii+1)/(double)(iii+1+negabove[iii]-1)-2*(scorep[iii].score-scoren[negabove[iii]-1].score));
		negabove[iii]--;
	      }
	    }
	    else {
	      ii=-1; /* stop */
	      changed=1;
	      vmax=v;
	      /*
	      for(i=0;i<nump;i++) {
		printf("na[%d]=%d,",i,negabove[i]); 
	      }
	      printf("vmax=%f\n",vmax);
	      */
	    }
	    if((ii>0) && (negabove[ii-1]!=negabove[ii])) 
	      ii=-1; 
	  }
	}
      }
    }
    vv=objective_avgprec(negabove,scorep,nump,scoren,numn);
    printf("v=%f == vv=%f\n",v,vv);
    printf("vmax=%f\n",vmax);
    /*
    for(i=0;i<nump;i++) {
      printf("na[%d]=%d,",i,negabove[i]); 
    }
    */

    for(i=0;i<ybar.totdoc;i++) { /* create Psi with maximum objective */
      ybar.class[i]=0;
    }
    for(i=0;i<nump;i++) { 
      ybar.class[scorep[i].id]=numn-2*negabove[i];
      for(j=0;j<negabove[i];j++) {
	ybar.class[scoren[j].id]++;
      }
      for(j=negabove[i];j<numn;j++) {
	ybar.class[scoren[j].id]--;
      }
    }
  }
  else {
    printf("ERROR: Invalid loss function '%d'.\n",sparm->loss_function);
    exit(1);
  }

  if(struct_verbosity >= 1) {
    SVECTOR *fy=psi(x,y,sm,sparm);
    SVECTOR *fybar=psi(x,ybar,sm,sparm);
    DOC *exy=create_example(0,0,1,1,fy);
    DOC *exybar=create_example(0,0,1,1,fybar);
    printf("\n max_y (loss(y_i,y)+ w*Psi(x,y)=%f\n",
	   loss(y,ybar,sparm)+classify_example(sm->svm_model,exybar));
    printf("w*Psi(x,y)=%f, w*Psi(x,ybar)=%f\n",
	   classify_example(sm->svm_model,exy),
	   classify_example(sm->svm_model,exybar));
    free_example(exy,1);
    free_example(exybar,1);
  }
  free(score);
  free(scorep);
  free(scoren);

  return(ybar);
}


int         empty_label(LABEL y)
{
  /* Returns true, if y is an empty label. An empty label might be
     returned by find_most_violated_constraint_???(x, y, sm) if there
     is no incorrect label that can be found for x, or if it is unable
     to label x at all */
  return(0);
}

SVECTOR     *psi(PATTERN x, LABEL y, STRUCTMODEL *sm,
		 STRUCT_LEARN_PARM *sparm)
{
  /* Returns a feature vector describing the match between pattern x
     and label y. The feature vector is returned as a list of
     SVECTOR's. Each SVECTOR is in a sparse representation of pairs
     <featurenumber:featurevalue>, where the last pair has
     featurenumber 0 as a terminator. Featurenumbers start with 1 and
     end with sizePsi. Featuresnumbers that are not specified default
     to value 0. As mentioned before, psi() actually returns a list of
     SVECTOR's. Each SVECTOR has a field 'factor' and 'next'. 'next'
     specifies the next element in the list, terminated by a NULL
     pointer. The list can be though of as a linear combination of
     vectors, where each vector is weighted by its 'factor'. This
     linear combination of feature vectors is multiplied with the
     learned (kernelized) weight vector to score label y for pattern
     x. Without kernels, there will be one weight in sm.w for each
     feature. Note that psi has to match
     find_most_violated_constraint_???(x, y, sm) and vice versa. In
     particular, find_most_violated_constraint_???(x, y, sm) finds
     that ybar!=y that maximizes psi(x,ybar,sm)*sm.w (where * is the
     inner vector product) and the appropriate function of the
     loss + margin/slack rescaling method. See that paper for details. */
  SVECTOR *fvec=NULL,*fvec2;
  double *sum,*sum2;
  long i,totwords;

  /* The following special case speeds up computation for the linear
     kernel. The lines add the feature vectors for all examples into
     a single vector. This is more efficient for the linear kernel,
     but is invalid for all other kernels. */
  if(sm->svm_model->kernel_parm.kernel_type == LINEAR) {
    totwords=sparm->num_features;
    sum=(double *)my_malloc(sizeof(double)*(totwords+1));
    clear_nvector(sum,totwords);
    for(i=0;i<y.totdoc;i++) 
      add_vector_ns(sum,x.doc[i]->fvec,y.class[i]);

    /* For sparse kernel, replace feature vector Psi with k=L^-1*Psi */
    if(sm->sparse_kernel_type>0) { 
      sum2=prod_ltmatrix_nvector(sm->invL,sum+1);
      for(i=0;i<totwords;i++)
	sum[i+1]=sum2[i];
      free_nvector(sum2);
    }

    fvec=create_svector_n(sum,totwords,"",1);
    free(sum);
  }
  else { /* general kernel */
    for(i=0;i<y.totdoc;i++) {
      fvec2=copy_svector(x.doc[i]->fvec);
      fvec2->next=fvec;
      fvec2->factor=y.class[i];
      /* printf("class %d: %f\n",i,y.class[i]); */
      fvec=fvec2;
    }
  }

  return(fvec);
}

double      loss(LABEL y, LABEL ybar, STRUCT_LEARN_PARM *sparm)
{
  /* loss for correct label y and predicted label ybar. The loss for
     y==ybar has to be zero. sparm->loss_function is set with the -l option. */
  int a=0,b=0,c=0,d=0,i;
  double loss=1;

  /* compute contingency table */
  for(i=0;i<y.totdoc;i++) {
    if((y.class[i] > 0) && (ybar.class[i] > 0)) {
      a++;
    }
    if((y.class[i] > 0) && (ybar.class[i] <= 0)) {
      c++;
    }
    if((y.class[i] < 0) && (ybar.class[i] > 0)) {
      b++;
    }
    if((y.class[i] < 0) && (ybar.class[i] <= 0)) {
      d++;
    }
    /* printf("%f %f\n",y.class[i],ybar.class[i]); */
  }
  /* Return the loss according to the selected loss function. */
  if(sparm->loss_function == ZEROONE) { /* type 0 loss: 0/1 loss */
                                  /* return 0, if y==ybar. return 1 else */
    loss=zeroone_loss(a,b,c,d);
  }
  else if(sparm->loss_function == FONE) {
    loss=fone_loss(a,b,c,d);
  }
  else if(sparm->loss_function == ERRORRATE) {
    loss=errorrate_loss(a,b,c,d);
  }
  else if(sparm->loss_function == PRBEP) {
    /* WARNING: only valid if called for a labeling that is at PRBEP */
    loss=prbep_loss(a,b,c,d);
  }
  else if(sparm->loss_function == PREC_K) {
    /* WARNING: only valid if for a labeling that predicts k positives */
    loss=prec_k_loss(a,b,c,d);
  }
  else if(sparm->loss_function == REC_K) {
    /* WARNING: only valid if for a labeling that predicts k positives */
    loss=rec_k_loss(a,b,c,d);
  }
  else if(sparm->loss_function == SWAPPEDPAIRS) {
    loss=swappedpairs_loss(y,ybar);
  }
  else if(sparm->loss_function == AVGPREC) {
    loss=avgprec_loss(y,ybar);
  }
  else {
    /* Put your code for different loss functions here. But then
       find_most_violated_constraint_???(x, y, sm) has to return the
       highest scoring label with the largest loss. */
    printf("Unkown loss function type: %d\n",sparm->loss_function);
    exit(1);
  }
  if(struct_verbosity >= 2)
    printf("    loss=%f; contingency_table=(%d, %d, %d, %d)\n",loss,a,b,c,d);
  return(loss);
}

int         finalize_iteration(double ceps, int cached_constraint,
			       SAMPLE sample, STRUCTMODEL *sm,
			       CONSTSET cset, double *alpha, 
			       STRUCT_LEARN_PARM *sparm)
{
  /* This function is called just before the end of each cutting plane iteration. ceps is the amount by which the most violated constraint found in the current iteration was violated. cached_constraint is true if the added constraint was constructed from the cache. If the return value is FALSE, then the algorithm is allowed to terminate. If it is TRUE, the algorithm will keep iterating even if the desired precision sparm->epsilon is already reached. */
  return(0);
}

void        print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
					CONSTSET cset, double *alpha, 
					STRUCT_LEARN_PARM *sparm)
{
  /* This function is called after training and allows final touches to
     the model sm. But primarly it allows computing and printing any
     kind of statistic (e.g. training error) you might want. */
  int i,j;
  MODEL *model=sm->svm_model;

  /* Move weight from bias feature to threshold parameter in svm model */
  /* WARNING: This following is correct only for the linear kernel! */
  if(sparm->bias_featurenum) { 
    sm->svm_model->b=-sparm->bias*sm->svm_model->lin_weights[sparm->bias_featurenum];
    /* Set the bias feature to zero in all examples */
    for(i=1;i<sm->svm_model->sv_num;i++) {
      for(j=0;sm->svm_model->supvec[i]->fvec->words[j].wnum;j++)
	if(sm->svm_model->supvec[i]->fvec->words[j].wnum==sparm->bias_featurenum) {
	  sm->svm_model->supvec[i]->fvec->words[j].wnum=0;
	  sm->svm_model->supvec[i]->fvec->words[j].weight=0;
	}
    }
  }
  /* Replace SV with single weight vector */
  if(model->kernel_parm.kernel_type == LINEAR) {
    if(struct_verbosity>=1) {
      printf("Compacting linear model..."); fflush(stdout);
    }
    sm->svm_model=compact_linear_model(model);
    sm->w=sm->svm_model->lin_weights; /* short cut to weight vector */
    free_model(model,1);
    if(struct_verbosity>=1) {
      printf("done\n"); fflush(stdout);
    }
  }  
}

void        print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
				       STRUCT_LEARN_PARM *sparm, 
				       STRUCT_TEST_STATS *teststats)
{
  /* This function is called after making all test predictions in
     svm_struct_classify and allows computing and printing any kind of
     evaluation (e.g. precision/recall) you might want. You can use
     the function eval_prediction to accumulate the necessary
     statistics for each prediction. */

  if(!teststats->test_data_unlabeled) {
    printf("NOTE: The loss reported above is the percentage of errors. The zero/one-error\n");
    printf("      is the multivariate zero/one-error regarding the whole prediction vector!\n");
    printf("Accuracy : %5.2f\n",100-teststats->errorrate);
    printf("Precision: %5.2f\n",teststats->precision);
    printf("Recall   : %5.2f\n",teststats->recall);
    printf("F1       : %5.2f\n",teststats->fone);
    printf("PRBEP    : %5.2f\n",teststats->prbep);
    printf("ROCArea  : %5.2f\n",teststats->rocarea);
    printf("AvgPrec  : %5.2f\n",teststats->avgprec);
  }
  else {
    printf("NOTE: %ld test examples are unlabeled, so performance cannot be computed. The\n",teststats->test_data_unlabeled);
    printf("      loss and the error reported above may be inaccurate!\n");
  }
}

void        eval_prediction(long exnum, EXAMPLE ex, LABEL ypred, 
			    STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, 
			    STRUCT_TEST_STATS *teststats)
{
  /* This function allows you to accumlate statistic for how well the
     predicition matches the labeled example. It is called from
     svm_struct_classify. See also the function
     print_struct_testing_stats. */
  long i;
  if(exnum == 0) { /* this is the first time the function is
		      called. So initialize the teststats */
  }
  teststats->test_data_unlabeled=0;
  for(i=0;i<ex.y.totdoc;i++) 
    if(ex.y.class[i]==0) 
      teststats->test_data_unlabeled++;
  if(!teststats->test_data_unlabeled) {
    sparm->loss_function=ERRORRATE;
    teststats->errorrate=loss(ex.y,ypred,sparm);
    sparm->loss_function=PREC_K;
    teststats->precision=100.0-loss(ex.y,ypred,sparm);
    sparm->loss_function=REC_K;
    teststats->recall=100.0-loss(ex.y,ypred,sparm);
    sparm->loss_function=FONE;
    teststats->fone=100.0-loss(ex.y,ypred,sparm);
    teststats->prbep=prbep(ex.y,ypred);
    teststats->rocarea=rocarea(ex.y,ypred);
    teststats->avgprec=avgprec(ex.y,ypred);
  }
}
#ifdef COMPILE_MEX_INTERFACE
void        write_struct_model(mxArray *file[], STRUCTMODEL *sm, 
							   STRUCT_LEARN_PARM *sparm)
#else
void        write_struct_model(char *file, STRUCTMODEL *sm, 
			       STRUCT_LEARN_PARM *sparm)
#endif
{
  /* Writes structural model sm to file file. */
  MODEL *svm_model;
  long i,totwords;
  WORD *w;
  double *ortho_weights;

  /* Store model in normal svm-light format */
  if(sm->sparse_kernel_type > 0) {
    add_weight_vector_to_linear_model(sm->svm_model);
    svm_model=(MODEL *)malloc(sizeof(MODEL));
    (*svm_model)=(*sm->svm_model);
    svm_model->sv_num=sm->expansion_size+1;
    svm_model->kernel_parm.kernel_type=sm->sparse_kernel_type;
    svm_model->supvec=sm->expansion-1;
    svm_model->alpha=(double *)my_malloc(sizeof(double)*(svm_model->sv_num));
    ortho_weights=prod_nvector_ltmatrix(sm->svm_model->lin_weights+1,sm->invL);
    for(i=0;i<sm->expansion_size;i++)
      svm_model->alpha[i+1]=ortho_weights[i];
    free_nvector(ortho_weights);
    totwords=0;
    for(i=1;i<svm_model->sv_num;i++)  /* find highest feature number */
      for(w=svm_model->supvec[i]->fvec->words;w->wnum;w++) 
	if(totwords < w->wnum) 
	  totwords=w->wnum;
    svm_model->totwords=totwords;
    write_model(file,svm_model);
    free(svm_model->alpha);
    free(svm_model);
  }
  else {
    write_model(file,sm->svm_model);
  }

}

#ifndef COMPILE_MEX_INTERFACE
STRUCTMODEL read_struct_model(char *file, STRUCT_LEARN_PARM *sparm)
#else
STRUCTMODEL read_struct_model(const mxArray *file, STRUCT_LEARN_PARM *sparm)
#endif
{
  /* Reads structural model sm from file file. This function is used
     only in the prediction module, not in the learning module. */
  STRUCTMODEL sm;
  
  sm.svm_model=read_model(file);
  sparm->loss_function=ERRORRATE;
  sparm->bias=0;
  sparm->bias_featurenum=0;
  sparm->num_features=sm.svm_model->totwords;
  sparm->truncate_fvec=(sm.svm_model->kernel_parm.kernel_type==LINEAR);
  if(sm.svm_model->kernel_parm.kernel_type == CUSTOM) /* double kernel */
    sparm->preimage_method=9;
  sm.invL=NULL;
  sm.expansion=NULL;
  sm.expansion_size=0;
  sm.sparse_kernel_type=0;
  sm.w=sm.svm_model->lin_weights;
  sm.sizePsi=sm.svm_model->totwords;
  if((sm.svm_model->kernel_parm.kernel_type!=LINEAR) && sparm->classify_dense)
    add_dense_vectors_to_model(sm.svm_model);
  return(sm);
}

#ifndef COMPILE_MEX_INTERFACE
void        write_label(FILE *fp, LABEL y)
{
  /* Writes label y to file handle fp. */
  int i;
  for(i=0;i<y.totdoc;i++) {
    fprintf(fp,"%.8f\n",y.class[i]);
  }
}
#else
void        write_label(double **ptr, LABEL y)
{
	/* Writes label y to file handle fp. */
	double *pos=(*ptr);
	int i;
	for(i=0;i<y.totdoc;i++,pos++) {
		(*pos)=y.class[i];
	}
	(*ptr)=pos;
}
#endif

void        free_pattern(PATTERN x) {
  /* Frees the memory of x. */
  int i;
  for(i=0;i<x.totdoc;i++) 
    free_example(x.doc[i],1);
  free(x.doc);
}

void        free_label(LABEL y) {
  /* Frees the memory of y. */
  free(y.class);
}

void        free_struct_model(STRUCTMODEL sm) 
{
  /* Frees the memory of model. */
  int i;
  /* if(sm.w) free(sm.w); */ /* this is free'd in free_model */
  if(sm.svm_model) free_model(sm.svm_model,1);
  /* add free calls for user defined data here */
  if(sm.invL) free_matrix(sm.invL);
  if(sm.expansion) {
    for(i=0;i<sm.expansion_size;i++) 
      free_example(sm.expansion[i],1);
    free(sm.expansion);
  }
}

void        free_struct_sample(SAMPLE s)
{
  /* Frees the memory of sample s. */
  int i;
  for(i=0;i<s.n;i++) { 
    free_pattern(s.examples[i].x);
    free_label(s.examples[i].y);
  }
  free(s.examples);
}

void        print_struct_help()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_learn. */
  printf("         --b float   -> value of L2-bias feature. A value of 0 implies not\n");
  printf("                        having a bias feature. (default 1)\n");
  printf("                        WARNING: This is implemented only for linear kernel!\n");
  printf("         --p [0..]   -> fraction of positive examples to use as value of k for\n");
  printf("                        Prec@k and Rec@k. 0 indicates to use (0.5 * #pos) for\n");
  printf("                        Prec@k and (2 * #pos) for Rec@k. #pos is the number of\n");
  printf("                        positive examples in the training set. (default 0)\n");
  printf("         --t [0..]   -> Use a priori selection of basis for sparse kernel\n");
  printf("                        training (must not use '-w 9'):\n");
  printf("                        0: do not use a priori selection. (default)\n");
  printf("                        1: sample basis vectors randomly from --f file.\n");
  printf("                        2: incomplete Cholesky using vectors from --f option.\n");
  printf("         --i [0..]   -> Use CPSP algorithm for sparse kernel training\n");
  printf("                        (must use '-w 9'):\n");
  printf("                        0: do not use CPSP algorithm from [3a].\n");
  printf("                        1: CPSP using subset selection for preimages via\n");
  printf("                           59/95 heuristic (see [3a])\n");
  printf("                        2: CPSP using fixed point search (see [3a], RBF kernel\n");
  printf("                           only) (default)\n");
  /*  printf("                        3: stochastic fixed point search (RBF kernel only)\n"); */
  printf("                        4: CPSP using fixed point search with starting point\n");
  printf("                           via 59/95 heuristic (see [3a], RBF kernel only)\n");
  printf("         --k [0..]   -> Specifies the number of basis functions to use\n");
  printf("                        for sparse kernel approximation (both --t and --i).\n");
  printf("                        (default 500)\n");
  printf("         --r [0..]   -> Number of times to recompute the sparse kernel\n");
  printf("                        approximation and restart the optimizer. (default 0)\n");
  printf("         --s [0,1]   -> Selects whether shrinking heuristic is used in the\n");
  printf("                        custom algorithms for linear kernel and errorrate loss.\n");
  printf("                        (default 1)\n");
  printf("         --f string  -> Specifies file that contains basis functions to use\n");
  printf("                        for sparse kernel approximation. (default training\n");
  printf("                        file)\n");
  printf("\nThe following loss functions can be selected with the -l option:\n");
  printf("    %2d  Zero/one loss: 1 if vector of predictions contains error, 0 otherwise.\n",ZEROONE);
  printf("    %2d  F1: 100 minus the F1-score in percent.\n",FONE);
  printf("    %2d  Errorrate: Percentage of errors in prediction vector.\n",ERRORRATE);
  printf("    %2d  Prec/Rec Breakeven: 100 minus PRBEP in percent.\n",PRBEP);
  printf("    %2d  Prec@k: 100 minus precision at k in percent.\n",PREC_K);
  printf("    %2d  Rec@k: 100 minus recall at k in percent.\n",REC_K);
  printf("    %2d  ROCArea: Percentage of swapped pos/neg pairs (i.e. 100 - ROCArea).\n\n",SWAPPEDPAIRS);
  printf("NOTE: The '-c' parameters in SVM-light and SVM-perf are related as\n");
  printf("      c_light = c_perf*100/n for the 'Errorrate' loss function, where n is the\n");
  printf("      number of training examples.\n\n");
  printf("The algorithms implemented in SVM-perf are described in:\n");
  printf("[1a] T. Joachims, A Support Vector Method for Multivariate Performance\n");
  printf("     Measures, Proceedings of the International Conference on Machine Learning\n");
  printf("     (ICML), 2005.\n");
  printf("[2a] T. Joachims, Training Linear SVMs in Linear Time, Proceedings of the\n");
  printf("     ACM Conference on Knowledge Discovery and Data Mining (KDD), 2006.\n");
  printf("[3a] T. Joachims, Chun-Nam Yu, Sparse Kernel SVMs via Cutting-Plane Training,\n");
  printf("     Proceedings of the European Conference on Machine Learning (ECML), 2009.\n");
  printf("-> Papers are available at http://www.joachims.org/\n\n");
}

void         parse_struct_parameters(STRUCT_LEARN_PARM *sparm)
{
  /* Parses the command line parameters that start with -- */
  int i;

  sparm->bias=1;
  sparm->prec_rec_k_frac=0.0;
  sparm->sparse_kernel_type=LINEAR;
  sparm->sparse_kernel_method=0;
  sparm->sparse_kernel_size=500;
  sparm->preimage_method=2;
  sparm->recompute_rset=0;
  sparm->shrinking=1;
  strcpy(sparm->sparse_kernel_file,"");
  /* set number of features to -1, indicating that it will be computed
     in init_struct_model() */
  sparm->num_features=-1;
  sparm->truncate_fvec=0;

  for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
    switch ((sparm->custom_argv[i])[2]) 
      { 
      case 'b': i++; sparm->bias=atof(sparm->custom_argv[i]); break;
      case 'p': i++; sparm->prec_rec_k_frac=atof(sparm->custom_argv[i]); break;
      case 'k': i++; sparm->sparse_kernel_size=atol(sparm->custom_argv[i]); break;
      case 't': i++; sparm->sparse_kernel_method=atol(sparm->custom_argv[i]); break;
      case 'i': i++; sparm->preimage_method=atol(sparm->custom_argv[i]); break;
      case 'r': i++; sparm->recompute_rset=atol(sparm->custom_argv[i]); break;
      case 's': i++; sparm->shrinking=atol(sparm->custom_argv[i]); break;
      case 'f': i++; strcpy(sparm->sparse_kernel_file,sparm->custom_argv[i]); break;
      default: printf("\nUnrecognized option %s!\n\n",sparm->custom_argv[i]);
	       exit(0);
      }
  }
  /* Note that the validity of the value for sparm->prec_rec_k_frac in
     relation to #pos is checked in read_struct_examples() */
  if(sparm->prec_rec_k_frac < 0) {
    printf("\nThe value of option --k must be greater then zero!\n\n");
    exit(0);
  }
  if((sparm->sparse_kernel_method > 0) && (sparm->preimage_method > 0)) {
    printf("\nYou cannot set both --i and --t to a non-zero value!\n\n");
    exit(0);
  }
}

void        print_struct_help_classify()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_classify. */
  printf("         --d [0,1]  -> use a dense vector representation when classifying\n");
  printf("                       new examples with svm_perf_classify. This is faster,\n");
  printf("                       but uses more memory. (default 1)\n");
}

void         parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm)
{
  /* Parses the command line parameters that start with -- for the
     classification module */
  int i;

  sparm->classify_dense=1;

  for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
    switch ((sparm->custom_argv[i])[2]) 
      { 
      case 'd': i++; sparm->classify_dense=atof(sparm->custom_argv[i]); break;
      default: printf("\nUnrecognized option %s!\n\n",sparm->custom_argv[i]);
	       exit(0);
      }
  }
}


/*------- Performance measures --------*/

double zeroone(int a, int b, int c, int d) 
{
  if((a+d) == (a+b+c+d)) 
    return(0.0);
  else
    return(1.0);
}

double fone(int a, int b, int c, int d) 
{
  if((a == 0) || (a+b == 0) || (a+c == 0)) return(0.0);
  double precision=prec(a,b,c,d);
  double recall=rec(a,b,c,d);
  return(2.0*precision*recall/(precision+recall));
}

double prec(int a, int b, int c, int d) 
{
  /* Returns precision as fractional value. */
  if((a+b) == 0) return(0.0);
  return((double)a/(double)(a+b));
}

double rec(int a, int b, int c, int d) 
{
  /* Returns recall as fractional value. */
  if((a+c) == 0) return(0.0);
  return((double)a/(double)(a+c));
}

double errorrate(int a, int b, int c, int d) 
{
  /* Returns number of errors. */
  if((a+b+c+d) == 0) return(0.0);
  return(((double)(b+c))/(double)(a+b+c+d));
}

double swappedpairs(LABEL y, LABEL ybar)
{
  /* Returns percentage of swapped pos/neg pairs (i.e. 100 - ROC Area) for
     prediction vectors that encode the number of misranked examples
     for each particular example. */
  /* WARNING: Works only for labels in the compressed representation */
  int i;
  double sum=0;
  for(i=0;i<y.totdoc;i++) 
    sum+=fabs(y.class[i]-ybar.class[i]);
  return(sum/2.0);
}

double rocarea(LABEL y, LABEL ybar)
{
  /* Returns ROC Area for ybar containing scores that define a ranking
     of examples. Breaks ties in ranking pessimistically. */
  long i,nump,numn;
  double swappedpairs;
  STRUCT_ID_SCORE *predset;

  predset=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));
  for(i=0;i<ybar.totdoc;i++) {
    predset[i].score=ybar.class[i];
    predset[i].tiebreak=-y.class[i];
    predset[i].id=i;
  }
  qsort(predset,ybar.totdoc,sizeof(STRUCT_ID_SCORE),comparedown);
  numn=0;
  nump=0;
  swappedpairs=0;
  for(i=0;i<ybar.totdoc;i++) {
    if(y.class[predset[i].id] > 0) {
      swappedpairs+=numn;
      nump++;
    }
    else {
      numn++;
    }
  }
  free(predset);
  return(100.0-100.0*swappedpairs/((double)numn)/((double)nump));
}

double prbep(LABEL y, LABEL ybar)
{
  /* Returns PRBEP for ybar containing scores that define a ranking
     of examples. Breaks ties in ranking pessimistically. */
  long i,nump,a;
  STRUCT_ID_SCORE *predset;

  predset=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));
  nump=0;
  for(i=0;i<ybar.totdoc;i++) {
    predset[i].score=ybar.class[i];
    predset[i].tiebreak=-y.class[i];
    predset[i].id=i;
    if(y.class[i] > 0) 
      nump++;
 }
  qsort(predset,ybar.totdoc,sizeof(STRUCT_ID_SCORE),comparedown);
  a=0;
  for(i=0;i<nump;i++) {
    if(y.class[predset[i].id] > 0) {
      a++;
    }
  }
  free(predset);
  return(100.0*prec(a,nump-a,0,0));
}

double avgprec_compressed(LABEL y, LABEL ybar)
{
  /* Returns Average Precision for y and ybar in compressed
     representation (also see avgprec()). Breaks ties in ranking
     pessimistically. */
  int i,ii,nump,numn,a,b;
  double apr;
  STRUCT_ID_SCORE *predset;

  nump=0;
  numn=0;
  for(i=0;i<ybar.totdoc;i++) {
    if(y.class[i] > 0) 
      nump++;
    else 
      numn++;
  }
  /*  printf("nump=%d, numn=%d\n", nump, numn); */

  ii=0;
  predset=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(nump+1));
  for(i=0;i<ybar.totdoc;i++) {
    if(y.class[i] > 0) {
      predset[ii].score=ybar.class[i];
      predset[ii].tiebreak=-y.class[i];
      predset[ii].id=i;
      ii++;
    }
  }
  qsort(predset,nump,sizeof(STRUCT_ID_SCORE),comparedown);

  apr=0;
  for(a=1;a<=nump;a++) {
    b=(numn-predset[a-1].score)/2;
    /* printf("negabove[%d]=%d,",a,b); */
    apr+=prec(a,b,0,0);
  }

  free(predset);

  return(100.0*(apr/(double)(nump)));
}

double avgprec(LABEL y, LABEL ybar)
{
  /* Returns Average Precision for ybar containing scores that define a ranking
     of examples. Breaks ties in ranking pessimistically. */
  long i,nump,numn;
  double apr;
  STRUCT_ID_SCORE *predset;

  predset=(STRUCT_ID_SCORE *)my_malloc(sizeof(STRUCT_ID_SCORE)*(ybar.totdoc+1));
  for(i=0;i<ybar.totdoc;i++) {
    predset[i].score=ybar.class[i];
    predset[i].tiebreak=-y.class[i];
    predset[i].id=i;
  }
  qsort(predset,ybar.totdoc,sizeof(STRUCT_ID_SCORE),comparedown);
  numn=0;
  nump=0;
  apr=0;
  for(i=0;i<ybar.totdoc;i++) {
    if(y.class[predset[i].id] > 0) {
      nump++;
      apr+=prec(nump,numn,0,0);
    }
    else {
      numn++;
    }
  }
  free(predset);
  return(100.0*(apr/(double)(nump)));
}

/*------- Loss functions based on performance measures --------*/

double zeroone_loss(int a, int b, int c, int d) 
{
  return(zeroone(a,b,c,d));
}

double fone_loss(int a, int b, int c, int d) 
{
  return(100.0*(1.0-fone(a,b,c,d)));
}

double errorrate_loss(int a, int b, int c, int d) 
{
  return(100.0*errorrate(a,b,c,d));
}

double prbep_loss(int a, int b, int c, int d) 
{
  /* WARNING: Returns lower bound on PRBEP, if b!=c. */
  double precision=prec(a,b,c,d);
  double recall=rec(a,b,c,d);
  if(precision < recall) 
    return(100.0*(1.0-precision));
  else
    return(100.0*(1.0-recall));
}

double prec_k_loss(int a, int b, int c, int d) 
{
  /* WARNING: Only valid if called with a+c==k. */
  return(100.0*(1.0-prec(a,b,c,d)));
}

double rec_k_loss(int a, int b, int c, int d) 
{
  /* WARNING: Only valid if called with a+c==k. */
  return(100.0*(1.0-rec(a,b,c,d)));
}

double swappedpairs_loss(LABEL y, LABEL ybar)
{  
  double nump=0,numn=0;
  long i;
  for(i=0;i<y.totdoc;i++) {
    if(y.class[i] > 0) 
      nump++;
    else 
      numn++;
  }
  /*  return(100.0*swappedpairs(y,ybar)/(nump*numn)); */
  return(swappedpairs(y,ybar));
}

double avgprec_loss(LABEL y, LABEL ybar)
{
  return(100.0-avgprec_compressed(y,ybar));
}



