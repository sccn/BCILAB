/***********************************************************************/
/*                                                                     */
/*   svm_struct_learn_custom.c (instantiated for SVM-perform)          */
/*                                                                     */
/*   Allows implementing a custom/alternate algorithm for solving      */
/*   the structual SVM optimization problem. The algorithm can use     */ 
/*   full access to the SVM-struct API and to SVM-light.               */
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 09.01.08                                                    */
/*                                                                     */
/*   Copyright (c) 2008  Thorsten Joachims - All rights reserved       */
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

long    find_preimage_count;
double  kernel_eval_count;
double  kernel_eval_count_preimage;

void svm_learn_struct_joint_custom_linear(SAMPLE sample, 
				   STRUCT_LEARN_PARM *sparm,
				   LEARN_PARM *lparm, KERNEL_PARM *kparm, 
	 			   STRUCTMODEL *sm);
void svm_learn_struct_joint_custom_kernel(SAMPLE sample, 
				   STRUCT_LEARN_PARM *sparm,
				   LEARN_PARM *lparm, KERNEL_PARM *kparm, 
	 			   STRUCTMODEL *sm);
MATRIX *init_kernel_matrix_rset(CONSTSET *rcset,KERNEL_PARM *kparm,
				STRUCTMODEL *sm);
MATRIX *update_kernel_matrix_rset(MATRIX *matrix, int newpos, CONSTSET *rcset, 
				  KERNEL_PARM *kparm, STRUCTMODEL *sm); 
void remove_inactive_constraints_rset(CONSTSET *cset, CONSTSET *rcset, 
				      double *alpha, long currentiter, 
				      long *alphahist, long mininactive);
void free_cset(CONSTSET cset, long deep);
void clear_reduced_set(STRUCTMODEL *sm);
void shrink_reduced_set(long rset_maxsize,KERNEL_PARM *kparm,
			STRUCT_LEARN_PARM *sparm,STRUCTMODEL *sm);
CONSTSET find_reduced_set_for_cset(CONSTSET *cset,DOC **x,long n,
				   long rset_maxsize,
				   double precision,double *precision_avg,
				   KERNEL_PARM *kparm,STRUCT_LEARN_PARM *sparm,
				   STRUCTMODEL *sm);
SVECTOR *find_reduced_set_for_psi(SVECTOR *psi_factors,DOC **x,long n,
				  long maxadd,double precision,
				  double *precision_avg,
				  KERNEL_PARM *kparm,STRUCT_LEARN_PARM *sparm,
				  STRUCTMODEL *sm);
void add_vector_to_reduced_set(SVECTOR *preimage,DOC **x,long n,
			       KERNEL_PARM *kparm,STRUCTMODEL *sm);
void delete_vector_from_reduced_set(long position,STRUCTMODEL *sm);
CONSTSET compute_projection_for_cset(CONSTSET *cset,int factors,DOC **x,long n,
				     KERNEL_PARM *kparm,STRUCTMODEL *sm);
SVECTOR *compute_projection_for_psi(SVECTOR *psi_factors,int factors,DOC **x,
				    long n,KERNEL_PARM *kparm,STRUCTMODEL *sm,
				    double *fit);
SVECTOR *solve_least_squares_projection(MATRIX *cholA,MATRIX *A,double *b, 
					SVECTOR **expansion, long rsize, 
					double *fit);
SVECTOR *solve_least_squares_projection_factors(MATRIX *cholA,MATRIX *A,
						double *b, 
						SVECTOR **expansion,long rsize,
						double *fit);
SVECTOR *find_preimage(SVECTOR *psi,SVECTOR *psihat,KERNEL_PARM *kparm,
		       STRUCT_LEARN_PARM *sparm,double *objective);
SVECTOR *find_preimage_fixedpoint(SVECTOR *psi,SVECTOR *psihat,
				  SVECTOR *start,KERNEL_PARM *kparm,
				  long maxfeatnum,double *objective);
SVECTOR *find_preimage_fixedpoint_sample(SVECTOR *psi,SVECTOR *psihat,
					 SVECTOR *start,KERNEL_PARM *kparm,
					 long maxfeatnum,
					 long samplesize,double *objective);
SVECTOR *find_preimage_subset(SVECTOR *psi,SVECTOR *psihat,KERNEL_PARM *kparm,
			      long samplesize,double *objective);
SVECTOR *find_preimage_average(SVECTOR *psi,SVECTOR *psihat,
			       KERNEL_PARM *kparm, long maxfeatnum,
			       double *objective);
CONSTSET optimize_reduced_set_for_cset(CONSTSET *cset_factors,DOC **x,long n,
				       KERNEL_PARM *kparm,STRUCTMODEL *sm);
SVECTOR *find_cset_preimage_fixedpoint(CONSTSET *cset_factors,
                          CONSTSET *rcset_factors,DOC **x,long n,
			  STRUCTMODEL *sm,
  		          SVECTOR *start,KERNEL_PARM *kparm,double *objective);

void svm_learn_struct_joint_custom(SAMPLE sample, STRUCT_LEARN_PARM *sparm,
								   LEARN_PARM *lparm, KERNEL_PARM *kparm, 
								   STRUCTMODEL *sm)
{
  if(sparm->sparse_kernel_method) {
    printf("ERROR: The algorithm '-w 9' is incompatible with the '--t' option.\n");
    printf("       Please us '-w 3' instead.\n");
    exit(1);
  }
  if((sparm->preimage_method == 0) && (kparm->kernel_type > 0)) {
    printf("ERROR: The algorithm '-w 9' is incompatible with the '--i 0' option when\n");
    printf("       using kernels. Please us '-w 3' instead.\n");
    exit(1);
  }

  if(kparm->kernel_type == LINEAR) 
    svm_learn_struct_joint_custom_linear(sample,sparm,lparm,kparm,sm);
  else
    svm_learn_struct_joint_custom_kernel(sample,sparm,lparm,kparm,sm);
}

void svm_learn_struct_joint_custom_linear(SAMPLE sample, 
				   STRUCT_LEARN_PARM *sparm,
				   LEARN_PARM *lparm, KERNEL_PARM *kparm, 
				   STRUCTMODEL *sm)
{
	int         i,j;
	int         numIt=0;
	long        argmax_count=0;
	long        totconstraints=0;
	double      epsilon,epsilon_shrink;
	double      lossval,factor,dist;
	double      margin=0;
	double      slack, slacksum, slacksum_full, ceps=99999;
	double      dualitygap,modellength,alphasum;
	long        sizePsi;
	double      *alpha=NULL;
	long        *alphahist=NULL,optcount=0;
	CONSTSET    cset;
	SVECTOR     *diff=NULL;
	double      *diff_n=NULL;
	SVECTOR     *lhs;
	MODEL       *svmModel=NULL;
	DOC         *doc;
	
	EXAMPLE     *ex=sample.examples;
	double      rt_total=0,rt_opt=0,rt_init=0,rt_psi=0,rt_viol=0,rt_kernel=0;
	double      rt1=0,rt2=0;
	
	double      score,diffscore,class,*class_old=NULL,*margin_old=NULL;
	double      *diffval=NULL,*diffvar=NULL,*diffmean=NULL,*diffconf=NULL;
	long        fullround=0,checkcount,flipcount;
	
	DOC         **x=(ex[0].x.doc);
	double      *y=(ex[0].y.class);
	long        n=ex[0].y.totdoc;
	
	rt1=get_runtime();
	
	init_struct_model(sample,sm,sparm,lparm,kparm); 
	sizePsi=sm->sizePsi+1;          /* sm must contain size of psi on return */
	
	if(sparm->slack_norm == 1) {
		lparm->svm_c=sparm->C;          /* set upper bound C */
		lparm->sharedslack=1;
	}
	else if(sparm->slack_norm == 2) {
		printf("ERROR: The joint algorithm does not apply to L2 slack norm!"); 
		fflush(stdout);
		exit(0); 
	}
	else {
		printf("ERROR: Slack norm must be L1 or L2!"); fflush(stdout);
		exit(0);
	}
	if(sparm->loss_function != ERRORRATE) {
		printf("ERROR: The custom algorithm can only optimize errorrate as the loss\n"); 
		printf("       function. Please use algorithm '-w 3' instead.\n"); 
		fflush(stdout);
		exit(1);
	}
	if(kparm->kernel_type != LINEAR) {
    printf("ERROR: Non-linear kernels are not supported in this algorithm!\n"); 
		printf("       Please use algorithm '-w 3' instead.\n"); 
		fflush(stdout);
		exit(1);
	}
	
	lparm->biased_hyperplane=0;     /* set threshold to zero */
	epsilon=100.0;                  /* start with low precision and
	 increase later */
	epsilon_shrink=epsilon;
	
	cset=init_struct_constraints(sample, sm, sparm);
	if(cset.m > 0) {
		alpha=(double *)realloc(alpha,sizeof(double)*cset.m);
		alphahist=(long *)realloc(alphahist,sizeof(long)*cset.m);
		for(i=0; i<cset.m; i++) {
			alpha[i]=0;
			alphahist[i]=-1; /* -1 makes sure these constraints are never removed */
		}
	}
	kparm->gram_matrix=init_kernel_matrix(&cset,kparm);
	
	/* set initial model and slack variables */
	svmModel=(MODEL *)my_malloc(sizeof(MODEL));
	lparm->epsilon_crit=epsilon;
	svm_learn_optimization(cset.lhs,cset.rhs,cset.m,sizePsi+n,
						   lparm,kparm,NULL,svmModel,alpha);
	add_weight_vector_to_linear_model(svmModel);
	sm->svm_model=svmModel;
	sm->w=svmModel->lin_weights; /* short cut to weight vector */
	
	diff_n=create_nvector(sm->sizePsi);
	clear_nvector(diff_n,sm->sizePsi);
	class_old=(double *)my_malloc(sizeof(double)*(n));
	margin_old=(double *)my_malloc(sizeof(double)*(n));
	diffval=(double *)my_malloc(sizeof(double)*(n));
	diffvar=(double *)my_malloc(sizeof(double)*(n));
	diffmean=(double *)my_malloc(sizeof(double)*(n));
	diffconf=(double *)my_malloc(sizeof(double)*(n));
	for(i=0; i<n; i++) {
		class_old[i]=y[i];
		margin_old[i]=0;
		diffval[i]=0;
		diffvar[i]=0;
		diffmean[i]=100.0/n;
		diffconf[i]=-1;
	}
	margin=0;
	flipcount=n;
	
	rt_init+=MAX(get_runtime()-rt1,0);
	rt_total+=rt_init;
	
    /*****************/
	/*** main loop ***/
	/*****************/
	do { /* iteratively find and add constraints to working set */
		
		if(struct_verbosity>=1) { 
			printf("Iter %i: ",++numIt); 
			fflush(stdout);
		}
		
		rt1=get_runtime();
		
		/* Should the algorithm test all examples in this iteration? */
		fullround=(epsilon_shrink < (epsilon/5.0))
		|| (numIt < 5) 
		|| (ceps < sparm->epsilon)
		|| (!sparm->shrinking);
		
		/**** compute current slack ****/
		alphasum=0;
		for(j=0;(j<cset.m);j++) 
			alphasum+=alpha[j];
		for(j=0,slack=0;(j<cset.m) && (slack==0);j++)  
			if(alpha[j] > alphasum/cset.m)
				slack=MAX(0,cset.rhs[j]-classify_example(svmModel,cset.lhs[j]));
		
		/**** find a violated joint constraint ****/
		lhs=NULL;
		dist=0;
		flipcount=0;
		checkcount=0;
		
		if(struct_verbosity >= 2) rt2=get_runtime();
		
		/**** find most violated joint constraint ***/
		for(i=0; i<n; i++) {
			
			if(fullround || (diffconf[i] < 0)) {
				argmax_count++;
				checkcount++;
				score=classify_example(svmModel,x[i]);
				diffscore=100.0/n-2.0*y[i]*score;
				if(diffscore > 0) 
					class=-y[i];
				else 
					class=y[i];
				
				/**** scale feature vector for margin rescaling ****/
				factor=-(y[i]-class_old[i])+(y[i]-class);
				
				/**** add current fy-fybar to constraint ****/
				if(factor != 0) {
					flipcount++;
					add_vector_ns(diff_n,x[i]->fvec,factor); 
				}
				
				/**** add loss to rhs ****/
				if(y[i] == class)                         
					lossval=0;
				else
					lossval=100;
				margin+=lossval/n-margin_old[i];           
				class_old[i]=class;
				margin_old[i]=lossval/n;
				
	  /**** keep statistics for mini-batches ****/
				diffvar[i]=0.5*diffvar[i]+0.5*pow(diffmean[i]-diffscore,2.0);
				diffmean[i]=0.5*diffmean[i]+0.5*diffscore;
				diffconf[i]=fabs(diffmean[i])-1.0*sqrt(diffvar[i]);
				diffval[i]=diffscore;
			}
			
		} /* end of example loop */
		
		if(struct_verbosity>=2) {
			printf("%5.1f%% flipped / %5.1f%% checked ",100.0*flipcount/n,
				   100.0*checkcount/n);
			rt_viol+=MAX(get_runtime()-rt2,0);
			rt2=get_runtime();
		}
		
		/* create sparse vector from dense sum and ignore rounding errors*/
		for(i=0;i<=sm->sizePsi;i++) 
			if(fabs(diff_n[i]) < 0.00000000000000001) 
				diff_n[i]=0;
		diff=create_svector_n(diff_n,sm->sizePsi,"",1.0);
		
		/**** if `error', then add constraint and recompute QP ****/
		doc=create_example(cset.m,0,1,1,diff);
		dist=classify_example(svmModel,doc);
		free_example(doc,0);
		ceps=MAX(0,margin-dist-slack);
		if(struct_verbosity>=2) rt_psi+=MAX(get_runtime()-rt2,0);	
		if(ceps > sparm->epsilon) { 
			/**** resize constraint matrix and add new constraint ****/
			if(struct_verbosity>=2) rt2=get_runtime();
			cset.lhs=(DOC **)realloc(cset.lhs,sizeof(DOC *)*(cset.m+1));
			cset.lhs[cset.m]=create_example(cset.m,0,1,1,diff);
			cset.rhs=(double *)realloc(cset.rhs,sizeof(double)*(cset.m+1));
			cset.rhs[cset.m]=margin;
			alpha=(double *)realloc(alpha,sizeof(double)*(cset.m+1));
			alpha[cset.m]=0;
			alphahist=(long *)realloc(alphahist,sizeof(long)*(cset.m+1));
			alphahist[cset.m]=optcount;
			cset.m++;
			totconstraints++;
			if(struct_verbosity>=2) rt_psi+=MAX(get_runtime()-rt2,0);	
			
			if(struct_verbosity>=1) {
				printf(":");fflush(stdout);
			}
			if(struct_verbosity>=2) rt2=get_runtime();
			kparm->gram_matrix=update_kernel_matrix(kparm->gram_matrix,cset.m-1,
													&cset,kparm);
			if(struct_verbosity>=2) rt_kernel+=MAX(get_runtime()-rt2,0);
			
			/**** get new QP solution ****/
			if(struct_verbosity>=1) {
				printf("*");fflush(stdout);
			}
			if(struct_verbosity>=2) rt2=get_runtime();
			/* set svm precision so that higher than eps of most violated constr */
			if(!fullround) {
				epsilon_shrink=MIN(epsilon_shrink,MAX(ceps,sparm->epsilon)); 
				lparm->epsilon_crit=epsilon_shrink/2.0; 
			}
			else {
				epsilon=MIN(epsilon,MAX(ceps,sparm->epsilon)); /* best eps so far */
				lparm->epsilon_crit=epsilon/2.0; 
				epsilon_shrink=epsilon;
			}
			/* Run the QP solver on cset. */
			free_model(svmModel,0);
			svmModel=(MODEL *)my_malloc(sizeof(MODEL));
			kparm->kernel_type=GRAM; /* use kernel stored in kparm */
			svm_learn_optimization(cset.lhs,cset.rhs,cset.m,sizePsi,
								   lparm,kparm,NULL,svmModel,alpha);
			kparm->kernel_type=LINEAR; 
			svmModel->kernel_parm.kernel_type=LINEAR;
			add_weight_vector_to_linear_model(svmModel);
			sm->svm_model=svmModel;
			sm->w=svmModel->lin_weights; /* short cut to weight vector */
			optcount++;
			/* keep track of when each constraint was last
			 active. constraints marked with -1 are not updated */
			for(j=0;j<cset.m;j++) 
				if((alphahist[j]>-1) && (alpha[j] != 0))  
					alphahist[j]=optcount;
			if(struct_verbosity>=2) rt_opt+=MAX(get_runtime()-rt2,0);
			
			/* Check if some of the linear constraints have not been
			 active in a while. Those constraints are then removed to
			 avoid bloating the working set beyond necessity. */
			if(struct_verbosity>=2)
				printf("-");fflush(stdout);
			remove_inactive_constraints(&cset,alpha,optcount,alphahist,50);
		}
		else {
			free_svector(diff);
		}
		
		if(struct_verbosity>=1)
			printf("(NumConst=%d, SV=%ld, CEps=%.4f, QPEps=%.4f)\n",cset.m,
				   svmModel->sv_num-1,ceps,svmModel->maxdiff);
		
		rt_total+=MAX(get_runtime()-rt1,0);
		
	} while((ceps > sparm->epsilon) || (!fullround) ||
			finalize_iteration(ceps,0,sample,sm,cset,alpha,sparm)
			);
	
	
	if(struct_verbosity>=1) {
		printf("Final epsilon on KKT-Conditions: %.5f\n",
			   MAX(svmModel->maxdiff,ceps));
		
		/**** compute sum of slacks ****/
		/**** WARNING: If positivity constraints are used, then the
		 maximum slack id is larger than what is allocated
		 below ****/
		slacksum=0;
		for(j=0;j<cset.m;j++) 
			slacksum=MAX(slacksum,
						 cset.rhs[j]-classify_example(svmModel,cset.lhs[j]));
		slacksum_full=0;
		for(j=0;j<n;j++) 
			slacksum_full+=MAX(0,diffval[j]);
		
		alphasum=0;
		for(i=0; i<cset.m; i++)  
			alphasum+=alpha[i]*cset.rhs[i];
    modellength=model_length_s(svmModel);
		dualitygap=(0.5*modellength*modellength+sparm->C*(slacksum+ceps))
		-(alphasum-0.5*modellength*modellength);
		
		printf("Upper bound on duality gap: %.5f\n", dualitygap);
		printf("Dual objective value: dval=%.5f\n",
			   alphasum-0.5*modellength*modellength);
		printf("Primal objective value: pval=%.5f\n",
			   0.5*modellength*modellength+sparm->C*slacksum_full);
		printf("Total number of constraints in final working set: %i (of %i)\n",(int)cset.m,(int)totconstraints);
		printf("Number of iterations: %d\n",numIt);
		printf("Number of calls to 'find_most_violated_constraint': %ld\n",argmax_count);
		if(sparm->slack_norm == 1) {
			printf("Number of SV: %ld \n",svmModel->sv_num-1);
			printf("Norm of weight vector: |w|=%.5f\n",
	     model_length_s(svmModel));
		}
		else if(sparm->slack_norm == 2){ 
			printf("Number of SV: %ld (including %ld at upper bound)\n",
				   svmModel->sv_num-1,svmModel->at_upper_bound);
			printf("Norm of weight vector (including L2-loss): |w|=%.5f\n",
	     model_length_s(svmModel));
		}
		printf("Value of slack variable (on working set): xi=%.5f\n",slacksum);
		printf("Norm of longest difference vector: ||Psi(x,y)-Psi(x,ybar)||=%.5f\n",
			   length_of_longest_document_vector(cset.lhs,cset.m,kparm));
		if(struct_verbosity>=2) 
			printf("Runtime in cpu-seconds: %.2f (%.2f%% for QP, %.2f%% for kernel, %.2f%% for Argmax, %.2f%% for Psi, %.2f%% for init)\n",
				   rt_total/100.0, (100.0*rt_opt)/rt_total, (100.0*rt_kernel)/rt_total,
				   (100.0*rt_viol)/rt_total, (100.0*rt_psi)/rt_total, 
				   (100.0*rt_init)/rt_total);
		else if(struct_verbosity==1) 
			printf("Runtime in cpu-seconds: %.2f\n",rt_total/100.0);
	}
	if(struct_verbosity>=4)
		printW(sm->w,sizePsi,n,lparm->svm_c);
	
	if(svmModel) {
		sm->svm_model=copy_model(svmModel);
		sm->w=sm->svm_model->lin_weights; /* short cut to weight vector */
	}
	
	print_struct_learning_stats(sample,sm,cset,alpha,sparm);
	
	free(diff_n);
	free(class_old);
	free(margin_old);
	free(diffval);
	free(diffvar);
	free(diffmean);
	free(diffconf);
	if(svmModel)
		free_model(svmModel,0);
	free(alpha); 
	free(alphahist); 
  free_cset(cset,1);
  if(kparm->gram_matrix)
    free_matrix(kparm->gram_matrix);
}


void svm_learn_struct_joint_custom_kernel(SAMPLE sample, 
				   STRUCT_LEARN_PARM *sparm,
				   LEARN_PARM *lparm, KERNEL_PARM *kparm, 
				   STRUCTMODEL *sm)
{
  int         i,j,k;
  int         numIt=0;
  long        totconstraints=0,kernel_type_org,const_prev=-1;
  double      epsilon;
  double      lossval,factor,dist;
  double      margin=0;
  double      slack, slacksum, slacksum_full, ceps=99999;
  double      maxslack=0, minslack=0;
  double      dualitygap,modellength,alphasum,sum;
  long        sizePsi;
  double      *alpha=NULL;
  long        *alphahist=NULL,optcount=0;
  CONSTSET    cset,rcset;
  SVECTOR     *diff=NULL,*rdiff=NULL;
  SVECTOR     *lhs;
  MODEL       *svmModel=NULL;

  EXAMPLE     *ex=sample.examples;
  double      rt_total=0,rt_opt=0,rt_init=0,rt_psi=0,rt_viol=0,rt_kernel=0;
  double      rt_rset=0,rt_proj=0,rt1=0,rt2=0,rt3=0;

  double      *scores,diffscore,class;
  long        projection_outdated=1;
  long        recompute_kernel_matrix=0;
  long        force_solve_qp=0;
  long        rset_size;

  DOC         **x=(ex[0].x.doc);
  double      *y=(ex[0].y.class);
  long        n=ex[0].y.totdoc;
  WORD        *factors,*ai;

  long rset_maxadd=1;
  long rset_not_yet_recomputed=sparm->recompute_rset;
  long rset_maxsize=sparm->sparse_kernel_size;
  double epsilon_rset=0;
  double rset_precision_avg=0.000001;

  rt1=get_runtime();

  kernel_eval_count=0;
  kernel_eval_count_preimage=0;
  kernel_cache_statistic=0;
  find_preimage_count=0;

  init_struct_model(sample,sm,sparm,lparm,kparm); 
  sizePsi=sm->sizePsi+1;          /* sm must contain size of psi on return */

  if(sparm->slack_norm == 1) {
    lparm->svm_c=sparm->C;          /* set upper bound C */
    lparm->sharedslack=1;
  }
  else if(sparm->slack_norm == 2) {
    printf("ERROR: The joint algorithm does not apply to L2 slack norm!"); 
    fflush(stdout);
    exit(0); 
  }
  else {
    printf("ERROR: Slack norm must be L1 or L2!"); fflush(stdout);
    exit(0);
  }
  if(sparm->loss_function != ERRORRATE) {
    printf("ERROR: The custom algorithm can only optimize errorrate as the loss\n"); 
    printf("       function. Please use algorithm '-w 3' instead.\n"); 
    fflush(stdout);
    exit(1);
  }

  lparm->biased_hyperplane=0;     /* set threshold to zero */
  epsilon=100.0;                  /* start with low precision and
				     increase later */

  cset=init_struct_constraints(sample, sm, sparm);
  if(cset.m > 0) {
    printf("ERROR: This algorithms cannot handle additional constraints that are added\n"); 
    printf("       during initialization!"); 
    fflush(stdout);
    exit(0); 
  }
  rcset=cset;
  sm->reducedset=NULL;
  sm->reducedset_kernel=NULL;
  sm->reducedset_gram=NULL;
  sm->reducedset_cholgram=NULL;
  sm->reducedset_size=0;
  kparm->gram_matrix=init_kernel_matrix(&cset,kparm);

  /* set initial model and slack variables */
  svmModel=(MODEL *)my_malloc(sizeof(MODEL));
  lparm->epsilon_crit=epsilon;
  svm_learn_optimization(cset.lhs,cset.rhs,cset.m,sizePsi+n,
			 lparm,kparm,NULL,svmModel,alpha);
  add_weight_vector_to_linear_model(svmModel);
  sm->svm_model=svmModel;
  sm->w=svmModel->lin_weights; /* short cut to weight vector */

  factors=(WORD *)my_malloc(sizeof(WORD)*(n+1));
  scores=(double *)my_malloc(sizeof(double)*(n));

  projection_outdated=1;

  rt_init+=MAX(get_runtime()-rt1,0);
  rt_total+=rt_init;

    /*****************/
   /*** main loop ***/
  /*****************/
  do { /* iteratively find and add constraints to working set */

      if(struct_verbosity>=1) { 
	printf("Iter %i: ",++numIt); 
	fflush(stdout);
      }
      
      rt1=get_runtime();

      /**** compute distance from hyperplane for each training example ****/
      if(struct_verbosity>=2) rt2=get_runtime();
      for(i=0; i<n; i++)	  	
	/* scores[i]=classify_example(svmModel,x[i]); */
	for(k=0,scores[i]=0;k<sm->reducedset_size;k++) 
	  scores[i]+=svmModel->lin_weights[k+1]*sm->reducedset_kernel[k][i];
      if(struct_verbosity>=2) rt_viol+=MAX(get_runtime()-rt2,0);

      /**** compute precision of rset approximation ****/
      const_prev=cset.m-1;
      if(const_prev >= 0) {
	alphasum=0;
	for(j=0;(j<cset.m);j++) 
	  alphasum+=alpha[j];
	maxslack=0;
	minslack=0;
	long first=1;
	for(i=0;i<cset.m;i++) {
	  if(alpha[i] >= 0.000001*alphasum/cset.m) {
	    for(ai=cset.lhs[i]->fvec->words,sum=0;ai->wnum;ai++) 
	      sum+=ai->weight*scores[ai->wnum-1];
	    slack=MAX(0,cset.rhs[i]-sum);
	    if(first) {
	      minslack=slack; 
	      maxslack=slack;
	      first=0;
	    }
	    maxslack=MAX(maxslack,slack);
	    minslack=MIN(minslack,slack);
	  }
	}
	epsilon_rset=fabs(minslack-maxslack);
      }

      /**** compute current slack ****/
      alphasum=0;
      for(j=0;(j<cset.m);j++) 
	alphasum+=alpha[j];
      for(j=0,slack=0;(j<cset.m) && (slack==0);j++)  
	if(alpha[j] > 0.01*alphasum/cset.m) {
	  for(ai=cset.lhs[j]->fvec->words,sum=0;ai->wnum;ai++) 
	    sum+=ai->weight*scores[ai->wnum-1];
	  slack=MAX(0,cset.rhs[j]-sum);
	}

      /**** find a violated joint constraint ****/
      lhs=NULL;
      dist=0;
      ai=factors;
      margin=0;

      if(struct_verbosity >= 2) rt2=get_runtime();
	
      /**** find most violated joint constraint ***/
      for(i=0; i<n; i++) {
	  	
	diffscore=100.0/n-2.0*y[i]*scores[i];
	if(diffscore > 0) 
	  class=-y[i];
	else 
	  class=y[i];
	
	/**** scale feature vector for margin rescaling ****/
	factor=(y[i]-class);
	
	/**** add current fy-fybar to constraint ****/
	if(factor != 0) {
	  ai->wnum=i+1;
	  ai->weight=factor;
	  ai++;
	}
	
	/**** add loss to rhs ****/
	if(y[i] == class)                         
	  lossval=0;
	else
	  lossval=100;
	margin+=lossval/n;    
	dist+=scores[i]*factor;       

      } /* end of example loop */
      ai->wnum=0;
      ceps=MAX(0,margin-dist-slack);
      force_solve_qp=0;
      
      if(struct_verbosity>=2) rt_viol+=MAX(get_runtime()-rt2,0);

      /**** if `error', then add constraint and recompute QP ****/
      if((ceps > sparm->epsilon) || (rset_maxsize>sm->reducedset_size) ||
	projection_outdated || rset_not_yet_recomputed) { 

	recompute_kernel_matrix=0;
	force_solve_qp=1;
	/* create sparse vector with kernel expansion factors */
	rt2=get_runtime();
	diff=create_svector(factors,NULL,1.0);
	if(struct_verbosity>=2) rt_psi+=MAX(get_runtime()-rt2,0);	
	/* create reduced set approximation of vector */
	rt2=get_runtime();
	rset_size=sm->reducedset_size;
	rdiff=find_reduced_set_for_psi(diff,x,n,
			    MIN(rset_maxadd,rset_maxsize-sm->reducedset_size),
			    rset_precision_avg/200,&rset_precision_avg,
			    kparm,sparm,sm);
	if(rset_size!=sm->reducedset_size)
	  projection_outdated++;
	if(struct_verbosity>=2) rt_rset+=MAX(get_runtime()-rt2,0);
	/**** resize constraint matrix and add new constraint ****/
	if(struct_verbosity>=2) rt2=get_runtime();
	cset.lhs=(DOC **)realloc(cset.lhs,sizeof(DOC *)*(cset.m+1));
	cset.lhs[cset.m]=create_example(cset.m,0,1,1,diff);
	cset.rhs=(double *)realloc(cset.rhs,sizeof(double)*(cset.m+1));
	cset.rhs[cset.m]=margin;
	rcset.lhs=(DOC **)realloc(rcset.lhs,sizeof(DOC *)*(rcset.m+1));
	rcset.lhs[cset.m]=create_example(rcset.m,0,1,1,rdiff);
	rcset.rhs=(double *)realloc(rcset.rhs,sizeof(double)*(rcset.m+1));
	rcset.rhs[cset.m]=margin;
	alpha=(double *)realloc(alpha,sizeof(double)*(cset.m+1));
	alpha[cset.m]=0;
	alphahist=(long *)realloc(alphahist,sizeof(long)*(cset.m+1));
	alphahist[cset.m]=optcount;
	cset.m++;
	rcset.m++;
	totconstraints++;
	if(struct_verbosity>=2) rt_psi+=MAX(get_runtime()-rt2,0);	

	if(struct_verbosity>=2) rt2=get_runtime();
	/* add more basis functions for all psi in working set */
	if((ceps <= sparm->epsilon) && 
	   ((projection_outdated<=1) || (epsilon_rset<=sparm->epsilon)) &&
	   (!rset_not_yet_recomputed) && (rset_maxsize>sm->reducedset_size)) {
	  if(struct_verbosity>=1) {
	    printf("(extending reduced set [%f] ",rset_precision_avg); 
	    fflush(stdout);	
	  }
	  projection_outdated=0;
	  recompute_kernel_matrix=1;
	  free_cset(rcset,1);
	  rcset=find_reduced_set_for_cset(&cset,x,n,
				  MIN(rset_maxsize,sm->reducedset_size+cset.m),
				  rset_precision_avg/400,
				  &rset_precision_avg,kparm,sparm,sm);
	  if(struct_verbosity>=1) { printf(")");fflush(stdout); }
	}
	/* recompute reduced set from scratch */
	else if(1 && rset_not_yet_recomputed && 
		(ceps <= sparm->epsilon) && 
		((projection_outdated<=1) || (epsilon_rset<=sparm->epsilon))) {
	  if(struct_verbosity>=1) {
	    printf("(recomputing reduced set [%f] ",rset_precision_avg);
	    fflush(stdout);
	  }
	  rset_not_yet_recomputed--;
	  projection_outdated=0;
	  recompute_kernel_matrix=1;
	  clear_reduced_set(sm);
	  free_cset(rcset,1);
	  rcset=find_reduced_set_for_cset(&cset,x,n,
					  MIN(0.75*rset_maxsize,5*cset.m),
					  rset_precision_avg/1000,
					  &rset_precision_avg,
					  kparm,sparm,sm);
	  if(0) {
	    if(struct_verbosity>=1) {
	      printf("optimizing globally"); fflush(stdout);
	    }
	    free_cset(rcset,1);
	    rcset=optimize_reduced_set_for_cset(&cset,x,n,kparm,sm);
	  }
	  if(struct_verbosity>=1) { printf(")");fflush(stdout); }
	}
	/* DEPRICATED: remove oldest part of reduced set and replace
	   in following iterations */
	else if(0 && rset_not_yet_recomputed && 
		(ceps <= sparm->epsilon) && 
		((projection_outdated<=1) || (epsilon_rset<=sparm->epsilon)) &&
		(rset_maxsize<=sm->reducedset_size)){
	  if(struct_verbosity>=1) {
	    printf("(shrinking reduced set "); 
	    fflush(stdout);	
	  }
	  rset_not_yet_recomputed--;
	  projection_outdated=9999999;
	  shrink_reduced_set(0.75*rset_maxsize,kparm,sparm,sm);
	  if(struct_verbosity>=1) {
	    printf(")");fflush(stdout);
	  }
	}
        /* DEPRICATED: try to optimize reduced set of iteratively
	   removing basis vectors and replacing them with newly
	   optimized ones */
	else if(0 && rset_not_yet_recomputed && 
		(ceps <= sparm->epsilon) && 
		((projection_outdated<=1) || (epsilon_rset<=sparm->epsilon)) &&
		(sm->reducedset_size>=rset_maxsize)) {
	  if(struct_verbosity>=1) {
	    printf("(optimizing reduced set [%f] ",rset_precision_avg);
	    fflush(stdout);
	  }
	  rset_not_yet_recomputed--;
	  projection_outdated=0;
	  recompute_kernel_matrix=1;
	  shrink_reduced_set(0.25*rset_maxsize,kparm,sparm,sm);
	  free_cset(rcset,1);
	  rcset=optimize_reduced_set_for_cset(&cset,x,n,kparm,sm);
	  if(struct_verbosity>=1) { printf(")");fflush(stdout); }
	}
	if(struct_verbosity>=2) rt_rset+=MAX(get_runtime()-rt2,0);

	/* periodically recompute projection of all vectors */
	if((projection_outdated > 10) ||
	   (projection_outdated && (epsilon_rset > 0.9*ceps)) ||
	   (projection_outdated && (rset_maxsize==sm->reducedset_size))) {
	  if(struct_verbosity>=2) {
	    rt2=get_runtime();
	    printf("(subspace projection"); fflush(stdout);
	  }	
	  projection_outdated=0;
	  recompute_kernel_matrix=1;
	  free_cset(rcset,1);
	  rcset=compute_projection_for_cset(&cset,1,x,n,kparm,sm);
	  if(struct_verbosity>=2) {
	    rt_proj+=MAX(get_runtime()-rt2,0);
	    printf(")"); fflush(stdout);
	  }
	}

	if(struct_verbosity>=1) { printf(":");fflush(stdout); }
	if(struct_verbosity>=2) rt3=get_runtime();
	if(recompute_kernel_matrix) {
	  /* compute gram matrix of one-slack dual qp from scratch */
	  free_matrix(kparm->gram_matrix);
	  kparm->gram_matrix=init_kernel_matrix_rset(&rcset,kparm,sm);
	}
        else {
	  /* only update the gram matrix for the newly added vector */
	  kparm->gram_matrix=update_kernel_matrix_rset(kparm->gram_matrix,
						       rcset.m-1,
						       &rcset,kparm,sm);
	}
	if(struct_verbosity>=2) rt_kernel+=MAX(get_runtime()-rt3,0);

	/**** get new QP solution ****/
	if(struct_verbosity>=1) {
	  printf("*");fflush(stdout);
	}
	if(struct_verbosity>=2) rt2=get_runtime();
	/* set svm precision so that higher than eps of most violated constr */
	epsilon=MIN(epsilon,MAX(ceps,sparm->epsilon)); /* best eps so far */
	lparm->epsilon_crit=epsilon/4.0; 
	/* Run the QP solver on cset. */
	free_model(svmModel,0);
	svmModel=(MODEL *)my_malloc(sizeof(MODEL));
	kernel_type_org=kparm->kernel_type; 
	kparm->kernel_type=GRAM; /* use kernel stored in kparm */
	kernel_eval_count+=kernel_cache_statistic;
	kernel_cache_statistic=0;
	svm_learn_optimization(rcset.lhs,rcset.rhs,rcset.m,sizePsi,
			       lparm,kparm,NULL,svmModel,alpha);
	/* linear weight vector contains factors in reduced set subspace */
	svmModel->totwords=sm->reducedset_size+2;	
	add_weight_vector_to_linear_model(svmModel);
	/* tidy up svm model */
	svmModel->totwords=sm->sizePsi;	
	kparm->kernel_type=kernel_type_org; 
	svmModel->kernel_parm.kernel_type=kernel_type_org;
	sm->svm_model=svmModel;
	sm->w=svmModel->lin_weights;
	optcount++;
	/* keep track of when each constraint was last
	   active. constraints marked with -1 are not updated */
	for(j=0;j<cset.m;j++) 
	  if((alphahist[j]>-1) && (alpha[j] != 0))  
	    alphahist[j]=optcount;
	if(struct_verbosity>=2) rt_opt+=MAX(get_runtime()-rt2,0);
	
	/* Check if some of the linear constraints have not been
	   active in a while. Those constraints are then removed to
	   avoid bloating the working set beyond necessity. */
	if(struct_verbosity>=2)
	  printf("-");fflush(stdout);
	remove_inactive_constraints_rset(&cset,&rcset,alpha,optcount,
					 alphahist,5);
      }

      if(struct_verbosity>=1)
	printf("(NumConst=%d, SV=%ld, RSet=%ld, CEps=%.4f, QPEps=%.4f, RSEps=%.4f)\n",
	       cset.m,svmModel->sv_num-1,sm->reducedset_size,ceps,
	       svmModel->maxdiff,epsilon_rset);
	
      rt_total+=MAX(get_runtime()-rt1,0);

  } while((ceps > sparm->epsilon) || 
	  (rset_maxsize > sm->reducedset_size) ||
	  projection_outdated ||
	  force_solve_qp ||
	  finalize_iteration(ceps,0,sample,sm,cset,alpha,sparm)
	 );
  

  if(struct_verbosity>=1) {
    printf("Final epsilon on KKT-Conditions: %.5f\n",
	   MAX(svmModel->maxdiff,ceps));

    /**** build normal SVM model from exansion ****/
    MODEL *m=(MODEL *)my_malloc(sizeof(MODEL));
    (*m)=(*svmModel);
    m->sv_num=sm->reducedset_size+1;
    m->supvec=my_malloc(sizeof(DOC *)*m->sv_num);
    m->alpha=my_malloc(sizeof(double)*m->sv_num);
    for(j=1;j<m->sv_num;j++) {
      m->alpha[j]=m->lin_weights[j];
      m->supvec[j]=create_example(j,0,1,1,copy_svector(sm->reducedset[j-1]));
    }
    m->lin_weights=NULL;
    m->index=NULL;
    free_model(svmModel,0);
    svmModel=m;

    /**** compute sum of slacks ****/
    /**** WARNING: If positivity constraints are used, then the
	  maximum slack id is larger than what is allocated
	  below ****/
    for(j=0,slacksum=0;j<cset.m;j++) {
      for(ai=cset.lhs[j]->fvec->words,sum=0;ai->wnum;ai++) {
	sum+=ai->weight*scores[ai->wnum-1];
      }
      slacksum=MAX(slacksum,cset.rhs[j]-sum);
      /* slacksum=MAX(slacksum,
	             rcset.rhs[j]-classify_example(svmModel,rcset.lhs[j]));*/
    }
    for(j=0,slacksum_full=0;j<n;j++) {
      slacksum_full+=MAX(0,100.0/n-2.0*y[j]*scores[j]);
    }

    alphasum=0;
	for(i=0;i<cset.m;i++) 
      alphasum+=alpha[i]*cset.rhs[i];
    modellength=model_length_s(svmModel);
    dualitygap=(0.5*modellength*modellength+sparm->C*(slacksum+ceps))
               -(alphasum-0.5*modellength*modellength);
    
    printf("Upper bound on duality gap: %.5f\n", dualitygap);
    printf("Dual objective value: dval=%.5f\n",
	    alphasum-0.5*modellength*modellength);
    printf("Primal objective value: pval=%.5f\n",
	    0.5*modellength*modellength+sparm->C*slacksum_full);
    printf("Total number of constraints in final working set: %i (of %i)\n",(int)cset.m,(int)totconstraints);
    printf("Number of iterations: %d\n",numIt);
    printf("Number of calls to find_preimage: %ld\n",find_preimage_count);
    kernel_eval_count-=kernel_eval_count_preimage;
    printf("Number of kernel evaluations in main algorithm : %.0f (%.3f%% of n^2)\n",kernel_eval_count,kernel_eval_count/(0.01*(double)n*(double)n));
    printf("Number of kernel evaluations in preimage method: %.0f (%.3f%% of n^2)\n",kernel_eval_count_preimage,kernel_eval_count_preimage/(0.01*(double)n*(double)n));
    printf("Number of SV: %ld \n",svmModel->sv_num-1);
    if(sparm->slack_norm == 1) 
      printf("Norm of weight vector: |w|=%.5f\n",model_length_s(svmModel));
    else if(sparm->slack_norm == 2)
      printf("Norm of weight vector (including L2-loss): |w|=%.5f\n",
	     model_length_s(svmModel));
    printf("Value of slack variable (on working set): xi=%.5f\n",slacksum);
    printf("Value of slack variable (global): xi=%.5f\n",slacksum_full);
    if(struct_verbosity>=2) 
      printf("Runtime in cpu-seconds: %.2f (%.2f%% for QP, %.2f%% for kernel, %.2f%% for Argmax, %.2f%% for Psi, %.2f%% for RSet, %.2f%% for proj, %.2f%% for init)\n",
	   rt_total/100.0, (100.0*rt_opt)/rt_total, (100.0*rt_kernel)/rt_total,
	   (100.0*rt_viol)/rt_total, (100.0*rt_psi)/rt_total, 
	   (100.0*rt_rset)/rt_total, (100.0*rt_proj)/rt_total, 
	   (100.0*rt_init)/rt_total);
    else if(struct_verbosity==1) 
      printf("Runtime in cpu-seconds: %.2f\n",rt_total/100.0);
  }
  if(struct_verbosity>=4)
    printW(sm->w,sizePsi,n,lparm->svm_c);

  double error=0;
  for(i=0;i<n;i++) {
    if(scores[i]*y[i]<=0) 
      error++;
  }
  printf("Training Accuracy: %.2f%%\n",100.0*(1.0-error/n));

  sm->svm_model=svmModel;
  sm->w=NULL;

  print_struct_learning_stats(sample,sm,cset,alpha,sparm);

  free(factors);
  free(scores);
  free(alpha); 
  free(alphahist); 
  free_cset(cset,1);
  free_cset(rcset,1);
	if(kparm->gram_matrix)
		free_matrix(kparm->gram_matrix);

  for(i=0;i<sm->reducedset_size;i++) 
    free_svector(sm->reducedset[i]);
  free(sm->reducedset);
  for(i=0;i<sm->reducedset_size;i++) 
    free(sm->reducedset_kernel[i]);
  free(sm->reducedset_kernel);
  free_matrix(sm->reducedset_gram);
  free_matrix(sm->reducedset_cholgram);
}

MATRIX *init_kernel_matrix_rset(CONSTSET *rcset,KERNEL_PARM *kparm,
				STRUCTMODEL *sm) 
     /* assigns a kernelid to each constraint in cset and creates the
	corresponding kernel matrix. */
{
  int i,j;
  double kval;
  MATRIX *matrix;
  WORD *ai,*bi;

  /* assign kernel id to each new constraint */
  for(i=0;i<rcset->m;i++)
    rcset->lhs[i]->kernelid=i;

  /* allocate kernel matrix as necessary */
  matrix=create_matrix(i+50,i+50);

  for(j=0;j<rcset->m;j++) {
    for(i=j;i<rcset->m;i++) {
      /* kval=kernel(kparm,rcset->lhs[j],rcset->lhs[i]); */
      for(ai=rcset->lhs[j]->fvec->words,kval=0;ai->wnum;ai++) 
	for(bi=rcset->lhs[i]->fvec->words;bi->wnum;bi++) 
	  kval+=ai->weight*bi->weight*sm->reducedset_gram->element[ai->wnum-1][bi->wnum-1];      
      matrix->element[j][i]=kval;
      matrix->element[i][j]=kval;
    }
  }
  return(matrix);
}

MATRIX *update_kernel_matrix_rset(MATRIX *matrix, int newpos, CONSTSET *rcset, 
				  KERNEL_PARM *kparm, STRUCTMODEL *sm) 
     /* assigns new kernelid to constraint in position newpos and
	fills the corresponding part of the kernel matrix */
{
  int i,maxkernelid=0,newid;
  double kval;
  double *used;
  WORD *ai,*bi;

  /* find free kernelid to assign to new constraint */
  for(i=0;i<rcset->m;i++) 
    if(i != newpos)
      maxkernelid=MAX(maxkernelid,rcset->lhs[i]->kernelid);
  used=create_nvector(maxkernelid+2);
  clear_nvector(used,maxkernelid+2);
  for(i=0;i<rcset->m;i++) 
    if(i != newpos)
      used[rcset->lhs[i]->kernelid]=1;
  for(newid=0;used[newid];newid++);
  free_nvector(used);
  rcset->lhs[newpos]->kernelid=newid;

  /* extend kernel matrix if necessary */
  maxkernelid=MAX(maxkernelid,newid);
  if((!matrix) || (maxkernelid>=matrix->m))
    matrix=realloc_matrix(matrix,maxkernelid+50,maxkernelid+50);

  for(i=0;i<rcset->m;i++) {
    /* kval=kernel(kparm,rcset->lhs[newpos],rcset->lhs[i]); */
    for(ai=rcset->lhs[newpos]->fvec->words,kval=0;ai->wnum;ai++) 
      for(bi=rcset->lhs[i]->fvec->words;bi->wnum;bi++) 
	kval+=ai->weight*bi->weight*sm->reducedset_gram->element[ai->wnum-1][bi->wnum-1];
    matrix->element[newid][rcset->lhs[i]->kernelid]=kval;
    matrix->element[rcset->lhs[i]->kernelid][newid]=kval;
  }
  return(matrix);
}

void free_cset(CONSTSET cset, long deep) {
  /* free memory used by cset */
  long i;
  free(cset.rhs); 
  for(i=0;i<cset.m;i++) 
    free_example(cset.lhs[i],deep);
  free(cset.lhs);
}

void remove_inactive_constraints_rset(CONSTSET *cset, CONSTSET *rcset, 
			         double *alpha, long currentiter, 
				 long *alphahist, long mininactive)
     /* removes the constraints from cset (and alpha) for which
	alphahist indicates that they have not been active for at
	least mininactive iterations */

{  
  long i,m;
  
  m=0;
  for(i=0;i<cset->m;i++) {
    if((alphahist[i]<0) || ((currentiter-alphahist[i]) < mininactive)) {
      /* keep constraints that are marked as -1 or which have recently
         been active */
      cset->lhs[m]=cset->lhs[i];      
      cset->lhs[m]->docnum=m;
      cset->rhs[m]=cset->rhs[i];
      rcset->lhs[m]=rcset->lhs[i];      
      rcset->lhs[m]->docnum=m;
      rcset->rhs[m]=rcset->rhs[i];
      alpha[m]=alpha[i];
      alphahist[m]=alphahist[i];
      m++;
    }
    else {
      free_example(cset->lhs[i],1);
      free_example(rcset->lhs[i],1);
    }
  }
  if(cset->m != m) {
    cset->m=m;
    cset->lhs=(DOC **)realloc(cset->lhs,sizeof(DOC *)*cset->m);
    cset->rhs=(double *)realloc(cset->rhs,sizeof(double)*cset->m);
    rcset->m=m;
    rcset->lhs=(DOC **)realloc(rcset->lhs,sizeof(DOC *)*rcset->m);
    rcset->rhs=(double *)realloc(rcset->rhs,sizeof(double)*rcset->m);
    /* alpha=realloc(alpha,sizeof(double)*cset->m); */
    /* alphahist=realloc(alphahist,sizeof(long)*cset->m); */
  }
}

void clear_reduced_set(STRUCTMODEL *sm)
     /* remove all basis vectors from reduced set */
{
  long     i;

  /* clean up existing reduced set */
  if(sm->reducedset_size>0) {
    for(i=0;i<sm->reducedset_size;i++) 
      free_svector(sm->reducedset[i]);
    free(sm->reducedset);
    for(i=0;i<sm->reducedset_size;i++) 
      free(sm->reducedset_kernel[i]);
    free(sm->reducedset_kernel);
    free_matrix(sm->reducedset_gram);
    free_matrix(sm->reducedset_cholgram);
  }
  /* initialize new reduced set data structures */
  sm->reducedset=NULL;
  sm->reducedset_kernel=NULL;
  sm->reducedset_gram=NULL;
  sm->reducedset_cholgram=NULL;
  sm->reducedset_size=0;
}



void shrink_reduced_set(long rset_maxsize,KERNEL_PARM *kparm,
			STRUCT_LEARN_PARM *sparm,STRUCTMODEL *sm)
     /* remove the oldest basis vectors from rset so that the number
	of remaining basis vectors is rset_maxsize */
{
  long osize,delnum,i,j;
  MATRIX  *A,*cholA;
  SVECTOR **expansion;
  float   **expansion_kernel;
  double   kval;

  osize=sm->reducedset_size;
  expansion=sm->reducedset;
  expansion_kernel=sm->reducedset_kernel;
  A=sm->reducedset_gram;
  cholA=sm->reducedset_cholgram;

  rset_maxsize=MAX(1,rset_maxsize);
  delnum=osize-rset_maxsize;

  if(rset_maxsize<osize) {
    /* remove expansion vectors */
    for(i=0;i<delnum;i++)
      free_svector(expansion[i]);
    for(i=delnum;i<osize;i++)
      expansion[i-delnum]=expansion[i];
    expansion=realloc(expansion,sizeof(SVECTOR *)*(rset_maxsize)); 
    /* shrink kernel matrix between expansion and training examples */
    for(i=0;i<delnum;i++)
      free(expansion_kernel[i]);
    for(i=delnum;i<osize;i++)
      expansion_kernel[i-delnum]=expansion_kernel[i];
    expansion_kernel=realloc(expansion_kernel,sizeof(float*)*(rset_maxsize));
    /* update expansion Gram matrix (just compute from scratch */
    free_matrix(A);
    A=create_matrix(rset_maxsize,rset_maxsize);
    for(i=0;i<rset_maxsize;i++) {
      for(j=i;j<rset_maxsize;j++) {
	if((!expansion[i]) || (!expansion[j])) exit(1);
	kval=kernel_s(kparm,expansion[j],expansion[i]);
	A->element[j][i]=kval;    
	A->element[i][j]=kval;
      }
    }
    /* update Cholesky decomposition */
    if(cholA)
      free_matrix(cholA);
    cholA=cholesky_matrix(A);
    /* extend linear part */

    sm->reducedset_size=rset_maxsize;
    sm->reducedset=expansion;
    sm->reducedset_kernel=expansion_kernel;
    sm->reducedset_gram=A;
    sm->reducedset_cholgram=cholA;
  }
}


CONSTSET find_reduced_set_for_cset(CONSTSET *cset,DOC **x,long n,
				   long rset_maxsize,double precision,
				   double *precision_avg,KERNEL_PARM *kparm,
				   STRUCT_LEARN_PARM *sparm,STRUCTMODEL *sm)
{
  CONSTSET rcset;
  long     i,reducedset_size_prev,iteradd;
  SVECTOR  *rpsi;

  /* search for new reduced set */
  iteradd=MIN(10,MAX(1,0.5*(rset_maxsize-sm->reducedset_size)/cset->m));
  do {
    reducedset_size_prev=sm->reducedset_size;
    for(i=0;(i<cset->m) && (sm->reducedset_size<rset_maxsize);i++) {
      rpsi=find_reduced_set_for_psi(cset->lhs[i]->fvec,x,n,
				 MIN(iteradd,rset_maxsize-sm->reducedset_size),
				 precision,precision_avg,kparm,sparm,sm);
      free_svector(rpsi);
    }
    if(struct_verbosity >= 2) {
      printf("|");fflush(stdout);
    }
  } while(reducedset_size_prev < sm->reducedset_size);

  rcset=compute_projection_for_cset(cset,1,x,n,kparm,sm);
  return(rcset);
}

SVECTOR *find_reduced_set_for_psi(SVECTOR *psi_factors,DOC **x,long n,
				  long maxadd,double precision,
				  double *precision_avg,KERNEL_PARM *kparm,
				  STRUCT_LEARN_PARM *sparm,STRUCTMODEL *sm)
{
  SVECTOR *rpsi_factors,*psi,*psihat,*next,*preimage,*from;
  WORD    *ai;
  long    i,rsize,msize,esize,iter,noprogress;
  double  *b,obj;
  long    trials;

  rsize=sm->reducedset_size;
  esize=MAX(0,maxadd);
  msize=rsize+esize;

  /* create linear part of least squares problem for psi */
  b=create_nvector(msize);
  clear_nvector(b,msize);
  for(i=0;i<rsize;i++) 
    for(ai=psi_factors->words;ai->wnum;ai++) 
      b[i]+=ai->weight*sm->reducedset_kernel[i][ai->wnum-1];

  /* create psi vector from factors and x */
  psi=NULL;
  for(ai=psi_factors->words;ai->wnum;ai++) {
    from=x[ai->wnum-1]->fvec;
    while(from) {
      next=create_svector_shallow(from->words,NULL,ai->weight);
      next->kernel_id=from->kernel_id;
      next->next=psi;
      psi=next;
      from=from->next;
    }
  }

  noprogress=0;
  for(iter=0;(iter<esize) && (noprogress<1);iter++) {
    /* compute best fit given current expansion */
    psihat=solve_least_squares_projection(sm->reducedset_cholgram,
					  sm->reducedset_gram,b,
					  sm->reducedset,rsize,NULL);

    /* compute next rset vector based on difference between psi and psihat */
    long kprev=kernel_cache_statistic;
    trials=3;
    preimage=NULL;
    do { 
      trials--;
      preimage=find_preimage(psi,psihat,kparm,sparm,&obj);
    } while((!preimage) && (trials>0)); /* restart search if not converged */
    kernel_eval_count_preimage+=kernel_cache_statistic-kprev;
    
    (*precision_avg)=0.9*(*precision_avg)+0.1*obj;
    if((preimage) && (obj>0) && 
       ((obj>precision) || (sm->reducedset_size==0))) {
      add_vector_to_reduced_set(preimage,x,n,kparm,sm);
      /* extend linear part */
      for(ai=psi_factors->words;ai->wnum;ai++)
	b[rsize]+=ai->weight*sm->reducedset_kernel[rsize][ai->wnum-1];
      rsize=sm->reducedset_size;
      if(struct_verbosity >= 2) {
	printf("~");fflush(stdout);
      }
      noprogress=0;
    }
    else {
      noprogress++;
      if(preimage)
	free_svector(preimage);
    }
    free_svector_shallow(psihat);
  }
  free_svector_shallow(psi);
  
  rpsi_factors=solve_least_squares_projection_factors(sm->reducedset_cholgram,
						      sm->reducedset_gram,
						      b,sm->reducedset,rsize,
						      NULL);
  free_nvector(b);

  return(rpsi_factors);
}

void add_vector_to_reduced_set(SVECTOR *preimage,DOC **x,long n,
			       KERNEL_PARM *kparm,STRUCTMODEL *sm)
     /* add preimage vector as basis vector to reduced set */
{
  long    i,j;
  double  kval;
  long    rsize=sm->reducedset_size;
  SVECTOR **expansion=sm->reducedset;
  float   **expansion_kernel=sm->reducedset_kernel;
  MATRIX  *A=sm->reducedset_gram;
  MATRIX  *cholA=sm->reducedset_cholgram;

  /* add new preimage to expansion */
  expansion=realloc(expansion,sizeof(SVECTOR *)*(rsize+1));
  expansion[rsize]=preimage;
  /* extend kernel matrix between expansion and training examples */
  expansion_kernel=realloc(expansion_kernel,sizeof(float*)*(rsize+1));
  expansion_kernel[rsize]=my_malloc(sizeof(float)*n);
  for(i=0;i<n;i++)
    expansion_kernel[rsize][i]=kernel_s(kparm,preimage,x[i]->fvec);
  /* extend expansion Gram matrix */
  A=realloc_matrix(A,rsize+1,rsize+1);
  for(j=0;j<rsize;j++) {
    if(expansion[j]) {
      kval=kernel_s(kparm,expansion[j],preimage);
      A->element[j][rsize]=kval;    
      A->element[rsize][j]=kval;
    }
  }
  A->element[rsize][rsize]=kernel_s(kparm,preimage,preimage);
  /* update Cholesky decomposition */
  if(cholA)
    cholA=cholesky_addcol_matrix(cholA,A->element[rsize]);
  else
    cholA=cholesky_matrix(A);
  rsize++;

  sm->reducedset_size=rsize;
  sm->reducedset=expansion;
  sm->reducedset_kernel=expansion_kernel;
  sm->reducedset_gram=A;
  sm->reducedset_cholgram=cholA;
}

void delete_vector_from_reduced_set(long position,STRUCTMODEL *sm)
     /* delete reduced set vector at position from reduced set */
{
  long    i,j;
  long    rsize=sm->reducedset_size;
  SVECTOR **expansion=sm->reducedset;
  float   **expansion_kernel=sm->reducedset_kernel;
  MATRIX  *A=sm->reducedset_gram;
  MATRIX  *cholA=sm->reducedset_cholgram;

  /* delete vector from expansion */
  free_svector(expansion[position]);
  for(i=position+1;i<rsize;i++) 
    expansion[i-1]=expansion[i];
  expansion=realloc(expansion,sizeof(SVECTOR *)*(rsize-1));
  /* shrink kernel matrix between expansion and training examples */
  free(expansion_kernel[position]);
  for(i=position+1;i<rsize;i++) 
    expansion_kernel[i-1]=expansion_kernel[i];
  expansion_kernel=realloc(expansion_kernel,sizeof(float*)*(rsize-1));
  /* shrink expansion Gram matrix */
  for(i=position+1;i<rsize;i++) 
    for(j=0;j<rsize;j++)
      A->element[i-1][j]=A->element[i][j];    
  for(i=position+1;i<rsize;i++) 
    for(j=0;j<rsize;j++)
      A->element[j][i-1]=A->element[j][i];    
  A=realloc_matrix(A,rsize-1,rsize-1);
  /* update Cholesky decomposition */
  if(cholA)
    free_matrix(cholA);
  cholA=cholesky_matrix(A);
  rsize--;

  sm->reducedset_size=rsize;
  sm->reducedset=expansion;
  sm->reducedset_kernel=expansion_kernel;
  sm->reducedset_gram=A;
  sm->reducedset_cholgram=cholA;
}


CONSTSET compute_projection_for_cset(CONSTSET *cset,int factors,DOC **x,long n,
				     KERNEL_PARM *kparm,STRUCTMODEL *sm)
     /* find least squares projection of all constraints in cset for
	the given reduced set expansion */
{
  CONSTSET rcset;
  long     i;
  SVECTOR  *rpsi;
  double   fit,fitsum=0;
  /* create new rcset data structure */
  rcset.m=cset->m;
  rcset.rhs=my_malloc(sizeof(double)*(cset->m));
  rcset.lhs=my_malloc(sizeof(SVECTOR *)*(cset->m));
  for(i=0;i<cset->m;i++) {
    if(struct_verbosity >= 3) {
      rpsi=compute_projection_for_psi(cset->lhs[i]->fvec,factors,x,n,kparm,sm,
				      &fit);
      fitsum+=fit;
    }
    else 
      rpsi=compute_projection_for_psi(cset->lhs[i]->fvec,factors,x,n,kparm,sm,
				      NULL);
    rcset.lhs[i]=create_example(cset->lhs[i]->docnum,0,1,1,rpsi);
    rcset.rhs[i]=cset->rhs[i];
  }
  if(struct_verbosity >= 3) {
    printf(" projfit=%lf",fitsum); fflush(stdout);
  }
  return(rcset);
}

SVECTOR *compute_projection_for_psi(SVECTOR *psi_factors,int factors,DOC **x,
				    long n,KERNEL_PARM *kparm,STRUCTMODEL *sm,
				    double *fit)
     /* find least squares projection of psi for the given reduced set
	expansion */
{
  MATRIX  *cholA=sm->reducedset_cholgram;
  MATRIX  *A=sm->reducedset_gram;
  double  *b;
  float   **expansion_kernel=sm->reducedset_kernel;
  SVECTOR *rpsi,**expansion=sm->reducedset;
  long    i,size=sm->reducedset_size;
  WORD    *ai;
  /* create linear part of least squares problem */
  b=create_nvector(size);
  clear_nvector(b,size);
  for(i=0;i<size;i++) 
    for(ai=psi_factors->words;ai->wnum;ai++) 
      b[i]+=ai->weight*expansion_kernel[i][ai->wnum-1];
      /* b[i]+=ai->weight*kernel_s(kparm,expansion[i],x[ai->wnum-1]->fvec);*/
  /* solve least squares problem */
  if(factors)
    rpsi=solve_least_squares_projection_factors(cholA,A,b,expansion,size,fit);
  else
    rpsi=solve_least_squares_projection(cholA,A,b,expansion,size,fit);
  free_nvector(b);
  return(rpsi);
}

SVECTOR *solve_least_squares_projection(MATRIX *cholA,MATRIX *A,double *b, 
					SVECTOR **expansion, long rsize, 
					double *fit)
     /* NOTE: the SVECTOR returned is only a shallow copy */
{
    SVECTOR *rpsi=NULL,*next,*from;
    double  *beta;
    long    j;

    if(rsize <=0 )
      return(NULL);

    beta=solve_psd_linear_system_cholesky(cholA,b);      
    if(fit)
      (*fit)=quad_nvector_matrix(beta,A)
	     -2.0*sprod_nvector_nvector(b,beta,rsize-1);
    for(j=0;j<rsize;j++) {
      from=expansion[j];
      while(from) {
	next=create_svector_shallow(from->words,NULL,beta[j]);
	next->kernel_id=from->kernel_id;
	next->next=rpsi;
	rpsi=next;
	from=from->next;
      }
    }
    free_nvector(beta);
    return(rpsi);
}

SVECTOR *solve_least_squares_projection_factors(MATRIX *cholA,MATRIX *A,
						double *b, 
						SVECTOR **expansion,long rsize,
						double *fit)
{
    SVECTOR *rpsi_factors=NULL;
    double  *beta;
    long    j;
    WORD    *factors;

    if(rsize <=0 )
      return(NULL);

    factors=(WORD *)my_malloc(sizeof(WORD)*(rsize+1));

    beta=solve_psd_linear_system_cholesky(cholA,b);     
    if(fit)
      (*fit)=quad_nvector_matrix(beta,A)
             -2.0*sprod_nvector_nvector(b,beta,rsize-1);
    for(j=0;j<rsize;j++) {
      factors[j].wnum=j+1;
      factors[j].weight=beta[j];
    }
    factors[j].wnum=0;
    free_nvector(beta);
    rpsi_factors=create_svector_shallow(factors,NULL,1.0);
    return(rpsi_factors);
}


SVECTOR *find_preimage(SVECTOR *psi,SVECTOR *psihat,KERNEL_PARM *kparm,
		       STRUCT_LEARN_PARM *sparm,double *obj) 
     /* selects between different methods for computing the best preimage approximation for psi-psihat */
{
  SVECTOR *preimage,*best_preimage,*start,*v;
  long sizePsi,pos,i,j;
  double best_obj;
  long fnum=sparm->num_features;

  find_preimage_count++;

  if(sparm->preimage_method == 1) {
    /* 59/95 sampling selection heuristic */
    preimage=find_preimage_subset(psi,psihat,kparm,59,obj);
  }
  else if(sparm->preimage_method == 2) {
    /* fixed point method for gaussian kernel */
    /* use random vector from psi expansion as starting point */
    sizePsi=listlength_svector(psi);
    pos=rand() % sizePsi;
    for(v=psi,j=0;j<pos;v=v->next,j++);
    start=create_svector(v->words,NULL,1.0); 
    preimage=find_preimage_fixedpoint(psi,psihat,start,kparm,fnum,obj);
    free_svector(start);
  }
  else if(sparm->preimage_method == 3) {
    /* stochastic version of fixed point method for gaussian kernel (BUGGY) */
    /* use random vector from psi expansion as starting point */
    sizePsi=listlength_svector(psi);
    pos=rand() % sizePsi;
    for(v=psi,j=0;j<pos;v=v->next,j++);
    start=create_svector(v->words,NULL,1.0); 
    preimage=find_preimage_fixedpoint_sample(psi,psihat,start,kparm,
					     fnum,20000,obj);
    free_svector(start);
  }
  else if(sparm->preimage_method == 4) {
    /* fixed point with subset selection heuristic as starting point */
    start=find_preimage_subset(psi,psihat,kparm,3,obj);
    preimage=find_preimage_fixedpoint(psi,psihat,start,kparm,fnum,obj);
    free_svector(start);
  }
  else if(sparm->preimage_method == 5) { 
    /* fixed point with 5 random restarts */
    /* use random vector from psi expansion as starting point */
    sizePsi=listlength_svector(psi);
    best_obj=0;
    best_preimage=NULL;
    for(i=0;i<5;i++) {
      pos=rand() % sizePsi;
      for(v=psi,j=0;j<pos;v=v->next,j++);
      start=create_svector(v->words,NULL,1.0); 
      preimage=find_preimage_fixedpoint(psi,psihat,start,kparm,fnum,obj);
      free_svector(start);
      if((*obj>best_obj) || (!best_preimage)) {
	if(best_preimage)
	  free_svector(best_preimage);
	best_preimage=preimage;
	best_obj=(*obj);
      }
    }
    preimage=best_preimage;
  }
  else if(sparm->preimage_method == 7) {
    /* pick a random vector from psi */
    preimage=find_preimage_subset(psi,psihat,kparm,1,obj);
  }
  else if(sparm->preimage_method == 8) {
    /* use mean vector in input space */
    preimage=find_preimage_average(psi,psihat,kparm,fnum,obj);
  }
  else if(sparm->preimage_method == 9) {
    /* NOT READY YET: version that handles kernels that are sums of individual kernels, especially linear + gaussian */
    SVECTOR *linpsi=NULL,*rbfpsi=NULL,*linpsihat=NULL,*rbfpsihat=NULL;
    SVECTOR *next,*psi_i;
    for(psi_i=psi;psi_i;psi_i=psi_i->next) {
      if(psi_i->kernel_id==0) {
	next=create_svector_shallow(psi_i->words,NULL,psi_i->factor);
	next->next=linpsi;
	linpsi=next;
      }
      else {
	next=create_svector_shallow(psi_i->words,NULL,psi_i->factor);
	next->next=rbfpsi;
	rbfpsi=next;
      }
    }
    for(psi_i=psihat;psi_i;psi_i=psi_i->next) {
      if(psi_i->kernel_id==0) {
	next=create_svector_shallow(psi_i->words,NULL,psi_i->factor);
	next->next=linpsihat;
	linpsihat=next;
      }
      else {
	next=create_svector_shallow(psi_i->words,NULL,psi_i->factor);
	next->next=rbfpsihat;
	rbfpsihat=next;
      }
    }
    preimage=find_preimage_average(linpsi,linpsihat,kparm,fnum,obj);
    preimage->kernel_id=0;

    sizePsi=listlength_svector(rbfpsi);
    pos=rand() % sizePsi;
    for(v=rbfpsi,j=0;j<pos;v=v->next,j++);
    start=create_svector(v->words,NULL,1.0); 
    long ktype=kparm->kernel_type;
    kparm->kernel_type=RBF;
    preimage->next=find_preimage_fixedpoint(rbfpsi,rbfpsihat,start,kparm,fnum,obj);
    kparm->kernel_type=ktype;
    if(preimage->next)
      preimage->next->kernel_id=2;
    free_svector(start);

    free_svector_shallow(linpsi);
    free_svector_shallow(linpsihat);
    free_svector_shallow(rbfpsi);
    free_svector_shallow(rbfpsihat);
  }
  else {
    printf("ERROR: Unknown preimage finding methods!\n");
    exit(1);
  }
  return(preimage);
}

SVECTOR *find_preimage_fixedpoint(SVECTOR *psi,SVECTOR *psihat,
				  SVECTOR *start,KERNEL_PARM *kparm, 
				  long maxfeatnum,
				  double *objective)
     /* works only for Gaussian kernel */
     /* uses fixed point method from S&S, page 548 */
     /* start is used as the starting point */
{
  SVECTOR *next,*prev,*psi_i;
  double  kval,sum,*next_n,delta;
  long    i;
  double  obj,obj_prev=0;

  double  stepsize=0.1;
  double  increase_stepsize=0.1;
  long    maxiter=30;

  if(kparm->kernel_type != RBF) {
    printf("ERROR: Fixed point method for preimage finding only applies to RBF kernel!\n");
    exit(1);
  }

  prev=create_svector(start->words,NULL,1.0);
  next_n=create_nvector(maxfeatnum);
  sum=1;
  for(i=0,delta=0;(delta<0.99999) && (i<maxiter) && (fabs(sum)>10E-90);i++) {
    sum=0;
    clear_nvector(next_n,maxfeatnum);
    for(psi_i=psi;psi_i;psi_i=psi_i->next) {
      kval=single_kernel_s(kparm,psi_i,prev);
      add_vector_ns(next_n,psi_i,kval);
      sum+=kval;
    }
    for(psi_i=psihat;psi_i;psi_i=psi_i->next) {
      kval=single_kernel_s(kparm,psi_i,prev);
      add_vector_ns(next_n,psi_i,-kval);
      sum-=kval;
    }
    obj=pow(sum,2.0);
    if(obj<obj_prev) {
      stepsize=0.1;
      increase_stepsize=0.0;
    }
    stepsize=MIN(1,stepsize+increase_stepsize);
    smult_nvector(next_n,maxfeatnum,stepsize/sum);
    add_vector_ns(next_n,prev,(1.0-stepsize));
    /* next=create_svector_n(next_n,maxfeatnum,NULL,1.0); */
    next=create_svector_nvector_n(next_n,maxfeatnum,NULL,1.0);
    delta=sprod_ss(prev,next)
          /(sqrt(sprod_ss(prev,prev))*sqrt(sprod_ss(next,next)));
    if(struct_verbosity>=4) {
      printf("(%ld) obj=%lf, delta=%lf, ||next||=%lf, ||prev||=%lf, sum=%lf\n",i,obj,delta,sqrt(sprod_ss(next,next)),sqrt(sprod_ss(prev,prev)),sum);
    }
    free_svector(prev);
    prev=next;
    obj_prev=obj;
  }
  free_nvector(next_n);
  if(struct_verbosity>=3) {
    obj=pow(kernel_s(kparm,next,psi)-kernel_s(kparm,next,psihat),2.0);
    printf("-> (%ld) obj=%lf %lf, delta=%lf, ||next||=%lf, sum=%lf\n",i,obj,pow(sum,2.0),delta,sqrt(sprod_ss(next,next)),sum);
  }
  (*objective)=obj;
  if(isnan(delta) || isnan(sum) || (fabs(sum)<10E-90)) { /* not converged */
    (*objective)=0;
    free_svector(next);
    next=NULL;
  }
  return(next);
}

SVECTOR *find_preimage_fixedpoint_sample(SVECTOR *psi,SVECTOR *psihat,
		       SVECTOR *start,KERNEL_PARM *kparm,long maxfeatnum,
		       long samplesize,double *objective)
     /* works only for Gaussian kernel */
     /* uses fixed point method from S&S, page 548 */
     /* start is used as the starting point */
     /* WARNING: THIS FUNCTION IS PROBABLY BUGGY */
{
  SVECTOR *next,*prev,*psi_i;
  double  kval,sum,*next_n,delta;
  long    i;
  double  obj,obj_prev=0;

  long    sizePsi,addnum;
  double  stepsize=0.2;
  double  increase_stepsize=0.1;
  long    maxiter=50;

  if(kparm->kernel_type != RBF) {
    printf("ERROR: Fixed point method for preimage finding only applies to RBF kernel!\n");
    exit(1);
  }

  prev=create_svector(start->words,NULL,1.0);
  sizePsi=listlength_svector(psi);

  next_n=create_nvector(maxfeatnum);
  sum=1;
  for(i=0,delta=0;(delta<0.99999) && (i<maxiter) && (fabs(sum)>10E-90);i++) {
    sum=0;
    clear_nvector(next_n,maxfeatnum);
    addnum=0;
    for(psi_i=psi;psi_i;psi_i=psi_i->next) {
      if((rand() % sizePsi) < samplesize) {
	kval=single_kernel_s(kparm,psi_i,prev);
	add_vector_ns(next_n,psi_i,kval);
	sum+=kval;
	addnum++;
      }
    }
    smult_nvector(next_n,maxfeatnum,(double)sizePsi/(double)addnum);
    sum*=(double)sizePsi/(double)addnum;
    for(psi_i=psihat;psi_i;psi_i=psi_i->next) {
      kval=single_kernel_s(kparm,psi_i,prev);
      add_vector_ns(next_n,psi_i,-kval);
      sum-=kval;
    }
    obj=pow(sum,2.0);
    if(obj<obj_prev) {
      stepsize=0.1;
      increase_stepsize=0.0;
    }
    stepsize=MIN(1,stepsize+increase_stepsize);
    smult_nvector(next_n,maxfeatnum,stepsize/sum);
    add_vector_ns(next_n,prev,(1.0-stepsize));
    next=create_svector_n(next_n,maxfeatnum,NULL,1.0);
    delta=sprod_ss(prev,next)
          /(sqrt(sprod_ss(prev,prev))*sqrt(sprod_ss(next,next)));
    if(struct_verbosity>=4) {
      printf("(%ld) obj=%lf, delta=%lf, ||next||=%lf, ||prev||=%lf, sum=%lf\n",
	     i,obj,delta,sqrt(sprod_ss(next,next)),
	     sqrt(sprod_ss(prev,prev)),sum);
    }
    free_svector(prev);
    prev=next;
    obj_prev=obj;
  }
  free_nvector(next_n);
  if(struct_verbosity>=3) {
    obj=pow(kernel_s(kparm,next,psi)-kernel_s(kparm,next,psihat),2.0);
    printf("-> (%ld) obj=%lf %lf, delta=%lf, ||next||=%lf, sum=%lf\n",i,obj,pow(sum,2.0),delta,sqrt(sprod_ss(next,next)),sum);
  }
  (*objective)=obj;
  if(isnan(delta) || isnan(sum) || (fabs(sum)<10E-90)) { /* not converged */
    (*objective)=0;
    free_svector(next);
    next=NULL;
  }
  return(next);
}

SVECTOR *find_preimage_subset(SVECTOR *psi,SVECTOR *psihat,
			      KERNEL_PARM *kparm,long samplesize,
			      double *objective)
     /* tries several (samplesize) vectors from psi and returns the
	one that maximizes cosine with psi-psihat */
     /* Works for any kernel. */
{
  SVECTOR *v,*curr,*best=NULL;
  double  kval;
  long    i,j,pos,sizePsi;
  double  obj,maxobj=-1;

  sizePsi=listlength_svector(psi);
  for(i=0;i<samplesize;i++) {
    pos=rand() % sizePsi;
    for(v=psi,j=0;j<pos;v=v->next,j++);
    curr=create_svector(v->words,NULL,1.0);
    kval=kernel_s(kparm,psi,curr);
    kval-=kernel_s(kparm,psihat,curr);
    kval/=kernel_s(kparm,curr,curr);
    obj=kval*kval;
    if(obj>maxobj) {
      if(best) 
	free_svector(best);
      best=copy_svector(curr);
      maxobj=obj;
    }
    if(struct_verbosity>=4) {
      printf("(%ld) obj=%lf\n",i,obj);
    }
    free_svector(curr);
  }
  if(struct_verbosity>=3) {
    obj=pow(kernel_s(kparm,best,psi)-kernel_s(kparm,best,psihat),2.0);
    printf("-> (%ld) obj=%lf, ||next||=%lf\n",i,obj,sqrt(sprod_ss(best,best)));
  }
  (*objective)=maxobj;
  return(best);
}

SVECTOR *find_preimage_average(SVECTOR *psi,SVECTOR *psihat,
			       KERNEL_PARM *kparm,long maxfeatnum,
			       double *objective)
     /* returns the mean vector in input space */
{
  SVECTOR *avg,*psi_i;
  double  kval,sum,*avg_n;

  avg_n=create_nvector(maxfeatnum);
  clear_nvector(avg_n,maxfeatnum);
  sum=0;
  for(psi_i=psi;psi_i;psi_i=psi_i->next) {
    add_vector_ns(avg_n,psi_i,psi_i->factor);
    sum++;
  }
  smult_nvector(avg_n,maxfeatnum,1/sum);
  avg=create_svector_nvector_n(avg_n,maxfeatnum,NULL,1.0);
  free_nvector(avg_n);
  kval=kernel_s(kparm,psi,avg);
  kval-=kernel_s(kparm,psihat,avg);
  kval/=kernel_s(kparm,avg,avg);
  (*objective)=kval*kval;
  if(struct_verbosity>=3) {
    printf("-> obj=%lf\n",(*objective));
  }
  return(avg);
}


CONSTSET optimize_reduced_set_for_cset(CONSTSET *cset_factors,DOC **x,long n,
				       KERNEL_PARM *kparm,STRUCTMODEL *sm)
     /* globally optimize the current cset by iteratively removing one
	basis vector, and replacing it with a freshly optimized
	preimage vector */
{
  CONSTSET rcset,rcset_factors;
  long     iteration;
  SVECTOR  *preimage,*start;
  double   obj;

  /* just to print objective */
  if(struct_verbosity >=3 ){
    rcset_factors=compute_projection_for_cset(cset_factors,1,x,n,kparm,sm);
    free_cset(rcset_factors,1);
  }

  for(iteration=0;iteration<2*sm->reducedset_size;iteration++) {
    /* remove reduced set vector to be optimized */
    start=copy_svector(sm->reducedset[0]);
    delete_vector_from_reduced_set(0,sm);
    /* find best projection of cset onto current reduced set */
    rcset_factors=compute_projection_for_cset(cset_factors,1,x,n,kparm,sm);
    preimage=find_cset_preimage_fixedpoint(cset_factors,&rcset_factors,x,n,sm,
					   start,kparm,&obj);
    add_vector_to_reduced_set(preimage,x,n,kparm,sm);
    free_cset(rcset_factors,1);
    free_svector(start);
  }

  rcset=compute_projection_for_cset(cset_factors,1,x,n,kparm,sm);
  return(rcset);
}

SVECTOR *find_cset_preimage_fixedpoint(CONSTSET *cset_factors,
                           CONSTSET *rcset_factors,DOC **x,long n,
			   STRUCTMODEL *sm,
			   SVECTOR *start,KERNEL_PARM *kparm,double *objective)
     /* works only for Gaussian kernel */
     /* uses fixed point method from S&S, page 548 */
     /* start is used as the starting point */
{
  SVECTOR *next,*prev;
  double  sum,*next_n,*zi,delta,beta=0;
  long    i,j,maxfeatnum;
  double  obj,obj_prev=0,obj_start=0;
  double  *kernel_cset,*kernel_rset,*factor_x,*factor_rset;
  double  stepsize=1.0;
  double  increase_stepsize=0.1;
  long    maxiter=5;
  long    rsize=sm->reducedset_size;
  WORD    *ai;

  if(kparm->kernel_type != RBF) {
    printf("ERROR: Fixed point method for preimage finding only applies to RBF kernel!\n");
    exit(1);
  }

  maxfeatnum=0; /* slow */
  for(i=0;i<n;i++)
    maxfeatnum=MAX(maxfeatnum,maxfeatnum_svector(x[i]->fvec));
  maxfeatnum=MAX(maxfeatnum,maxfeatnum_svector(start))+1;

  kernel_cset=my_malloc(sizeof(double)*(n+1));
  kernel_rset=my_malloc(sizeof(double)*(sm->reducedset_size+1));
  prev=create_svector(start->words,NULL,1.0);
  zi=create_nvector(maxfeatnum);
  next_n=create_nvector(maxfeatnum);
  factor_x=create_nvector(n);
  factor_rset=create_nvector(rsize);
  sum=1;
  for(i=0,delta=0;(delta<0.99999) && (i<maxiter) && (fabs(sum)>10E-90);i++) {
    clear_nvector(factor_x,n);
    clear_nvector(factor_rset,rsize);
    sum=0;
    for(j=0;j<n;j++)
      kernel_cset[j]=single_kernel(kparm,x[j]->fvec,prev);
    for(j=0;j<rsize;j++)
      kernel_rset[j]=single_kernel(kparm,sm->reducedset[j],prev);
    for(j=0;j<cset_factors->m;j++) {
      beta=0;
      for(ai=cset_factors->lhs[j]->fvec->words;ai->wnum;ai++) {
	beta+=ai->weight*kernel_cset[ai->wnum-1];
      }
      for(ai=rcset_factors->lhs[j]->fvec->words;ai->wnum;ai++) {
	beta-=ai->weight*kernel_rset[ai->wnum-1];
      }
      for(ai=cset_factors->lhs[j]->fvec->words;ai->wnum;ai++) {
	factor_x[ai->wnum-1]+=beta*ai->weight*kernel_cset[ai->wnum-1];
      }
      for(ai=rcset_factors->lhs[j]->fvec->words;ai->wnum;ai++) {
	factor_rset[ai->wnum-1]-=beta*ai->weight*kernel_rset[ai->wnum-1];
      }
      sum+=beta*beta;
    }
    clear_nvector(next_n,maxfeatnum);
    for(j=0;j<n;j++)
      if(factor_x[j] != 0)
	add_vector_ns(next_n,x[j]->fvec,factor_x[j]);
    for(j=0;j<rsize;j++)
      add_vector_ns(next_n,sm->reducedset[j],factor_rset[j]);
    obj=sum;
    if(i==0) obj_start=obj;
    if(obj<obj_prev) {
      stepsize=0.1;
      increase_stepsize=0.0;
    }
    stepsize=MIN(1,stepsize+increase_stepsize);
    smult_nvector(next_n,maxfeatnum,stepsize/sum);
    add_vector_ns(next_n,prev,(1.0-stepsize));
    next=create_svector_nvector_n(next_n,maxfeatnum,NULL,1.0);
    delta=sprod_ss(prev,next)
          /(sqrt(sprod_ss(prev,prev))*sqrt(sprod_ss(next,next)));
    if(struct_verbosity>=3) {
      printf("(%ld) obj=%lf, delta=%lf, ||next||=%lf, ||prev||=%lf, sum=%lf\n",
	     i,obj,delta,sqrt(sprod_ss(next,next)),sqrt(sprod_ss(prev,prev)),
	     sum);
    }
    free_svector(prev);
    prev=next;
    obj_prev=obj;
  }
  free_nvector(next_n);
  free_nvector(zi);
  free_nvector(factor_x);
  free_nvector(factor_rset);
  free(kernel_cset);
  free(kernel_rset);
  if(struct_verbosity>=3) {
    printf("-> (%ld) obj=%lf, obj_start=%lf, delta=%lf, ||next||=%lf, sum=%lf\n",i,obj,obj_start,delta,sqrt(sprod_ss(next,next)),sum);
  }
  (*objective)=obj;
  if(isnan(delta) || isnan(sum) || (fabs(sum)<10E-90)) { /* not converged */
    (*objective)=0;
    free_svector(next);
    next=NULL;
  }
  return(next);
}
