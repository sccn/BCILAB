/***********************************************************************/
/*                                                                     */
/*   svm_struct_main.c                                                 */
/*                                                                     */
/*   Command line interface to the alignment learning module of the    */
/*   Support Vector Machine.                                           */
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 03.07.04                                                    */
/*                                                                     */
/*   Copyright (c) 2004  Thorsten Joachims - All rights reserved       */
/*                                                                     */
/*   This software is available for non-commercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.                                             */
/*                                                                     */
/***********************************************************************/


/* the following enables you to use svm-learn out of C++ */
#ifdef __cplusplus
extern "C" {
#endif
#include "svm_common.h"
#include "svm_learn.h"
#ifdef __cplusplus
}
#endif
# include "svm_struct_learn.h"
# include "svm_struct_common.h"
# include "svm_struct_api.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
/* } */

#ifndef COMPILE_MEX_INTERFACE
char trainfile[200];           /* file with training examples */
char modelfile[200];           /* file for resulting classifier */

void   read_input_parameters(int, char **, char *, char *,long *, long *,
			     STRUCT_LEARN_PARM *, LEARN_PARM *, KERNEL_PARM *,
			     int *);
#else
#include "mex_interface.h"
void   read_input_parameters(int, char **, long *, long *,
							 STRUCT_LEARN_PARM *, LEARN_PARM *, KERNEL_PARM *,
							 int *);
#endif

void   wait_any_key();
void   print_help();

#ifdef COMPILE_MEX_INTERFACE
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
#else
int main (int argc, char* argv[])
{  
#endif
  SAMPLE sample;  /* training sample */
  LEARN_PARM learn_parm;
  KERNEL_PARM kernel_parm;
  STRUCT_LEARN_PARM struct_parm;
  STRUCTMODEL structmodel;
  int alg_type;
#ifdef COMPILE_MEX_INTERFACE
	int argc;
	char **argv;
	if (nrhs < 2) {
		print_help();
		return;
	}
	else if (nrhs==2) {
		argc=1;
		argv=(char **)my_malloc(MAX_ARGVS*sizeof(char *));
		argv[0]="OLR";
	}
	else
		create_argc_argv(prhs[2],&argc,&argv);
#endif
  svm_struct_learn_api_init(argc,argv);
#ifndef COMPILE_MEX_INTERFACE	
  read_input_parameters(argc,argv,trainfile,modelfile,&verbosity,
			&struct_verbosity,&struct_parm,&learn_parm,
			&kernel_parm,&alg_type);
#else
	read_input_parameters(argc,argv,&verbosity,
						  &struct_verbosity,&struct_parm,&learn_parm,
						  &kernel_parm,&alg_type);
#endif
   if(struct_verbosity>=1) {
    printf("Reading training examples..."); fflush(stdout);
  }
  /* read the training examples */
#ifndef COMPILE_MEX_INTERFACE
  sample=read_struct_examples(trainfile,&struct_parm);
#else
	sample=read_struct_examples(prhs,&struct_parm);
#endif
  if(struct_verbosity>=1) {
    printf("done\n"); fflush(stdout);
  }
  
  /* Do the learning and return structmodel. */
  if(alg_type == 0)
    svm_learn_struct(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,NSLACK_ALG);
  else if(alg_type == 1)
    svm_learn_struct(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,NSLACK_SHRINK_ALG);
  else if(alg_type == 2)
    svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_PRIMAL_ALG);
  else if(alg_type == 3)
    svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_DUAL_ALG);
  else if(alg_type == 4)
    svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_DUAL_CACHE_ALG);
  else if(alg_type == 9)
    svm_learn_struct_joint_custom(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel);
  else
    exit(1);

  /* Warning: The model contains references to the original data 'docs'.
     If you want to free the original data, and only keep the model, you 
     have to make a deep copy of 'model'. */
  if(struct_verbosity>=1) {
    printf("Writing learned model...");fflush(stdout);
  }
#ifndef COMPILE_MEX_INTERFACE
  write_struct_model(modelfile,&structmodel,&struct_parm);
#else
	write_struct_model(plhs,&structmodel,&struct_parm);
#endif
  if(struct_verbosity>=1) {
    printf("done\n");fflush(stdout);
  }

  free_struct_sample(sample);
  free_struct_model(structmodel);

  svm_struct_learn_api_exit();
#ifndef COMPILE_MEX_INTERFACE
  return 0;
#endif
}

/*---------------------------------------------------------------------------*/

#ifndef COMPILE_MEX_INTERFACE
void read_input_parameters(int argc,char *argv[],char *trainfile,
			   char *modelfile,
			   long *verbosity,long *struct_verbosity, 
			   STRUCT_LEARN_PARM *struct_parm,
			   LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm,
			   int *alg_type)
#else
void read_input_parameters(int argc,char *argv[],
							   long *verbosity,long *struct_verbosity, 
							   STRUCT_LEARN_PARM *struct_parm,
							   LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm,
							   int *alg_type)
#endif
{
  long i;
  char type[100];
  
  /* set default */
  (*alg_type)=DEFAULT_ALG_TYPE;
  struct_parm->C=-0.01;
  struct_parm->slack_norm=1;
  struct_parm->epsilon=DEFAULT_EPS;
  struct_parm->custom_argc=0;
  struct_parm->loss_function=DEFAULT_LOSS_FCT;
  struct_parm->loss_type=DEFAULT_RESCALING;
  struct_parm->newconstretrain=100;
  struct_parm->ccache_size=5;
  struct_parm->batch_size=100;

#ifndef COMPILE_MEX_INTERFACE
  strcpy (modelfile, "svm_struct_model");
#endif
  strcpy (learn_parm->predfile, "trans_predictions");
  strcpy (learn_parm->alphafile, "");
  (*verbosity)=0;/*verbosity for svm_light*/
  (*struct_verbosity)=1; /*verbosity for struct learning portion*/
  learn_parm->biased_hyperplane=1;
  learn_parm->remove_inconsistent=0;
  learn_parm->skip_final_opt_check=0;
  learn_parm->svm_maxqpsize=10;
  learn_parm->svm_newvarsinqp=0;
  learn_parm->svm_iter_to_shrink=-9999;
  learn_parm->maxiter=100000;
  learn_parm->kernel_cache_size=40;
  learn_parm->svm_c=99999999;  /* overridden by struct_parm->C */
  learn_parm->eps=0.001;       /* overridden by struct_parm->epsilon */
  learn_parm->transduction_posratio=-1.0;
  learn_parm->svm_costratio=1.0;
  learn_parm->svm_costratio_unlab=1.0;
  learn_parm->svm_unlabbound=1E-5;
  learn_parm->epsilon_crit=0.001;
  learn_parm->epsilon_a=1E-10;  /* changed from 1e-15 */
  learn_parm->compute_loo=0;
  learn_parm->rho=1.0;
  learn_parm->xa_depth=0;
  kernel_parm->kernel_type=0;
  kernel_parm->poly_degree=3;
  kernel_parm->rbf_gamma=1.0;
  kernel_parm->coef_lin=1;
  kernel_parm->coef_const=1;
  strcpy(kernel_parm->custom,"empty");
  strcpy(type,"c");

  for(i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) 
      { 
      case '?': print_help(); exit(0);
      case 'a': i++; strcpy(learn_parm->alphafile,argv[i]); break;
      case 'c': i++; struct_parm->C=atof(argv[i]); break;
      case 'p': i++; struct_parm->slack_norm=atol(argv[i]); break;
      case 'e': i++; struct_parm->epsilon=atof(argv[i]); break;
      case 'k': i++; struct_parm->newconstretrain=atol(argv[i]); break;
      case 'h': i++; learn_parm->svm_iter_to_shrink=atol(argv[i]); break;
      case '#': i++; learn_parm->maxiter=atol(argv[i]); break;
      case 'm': i++; learn_parm->kernel_cache_size=atol(argv[i]); break;
      case 'w': i++; (*alg_type)=atol(argv[i]); break;
      case 'o': i++; struct_parm->loss_type=atol(argv[i]); break;
      case 'n': i++; learn_parm->svm_newvarsinqp=atol(argv[i]); break;
      case 'q': i++; learn_parm->svm_maxqpsize=atol(argv[i]); break;
      case 'l': i++; struct_parm->loss_function=atol(argv[i]); break;
      case 'f': i++; struct_parm->ccache_size=atol(argv[i]); break;
      case 'b': i++; struct_parm->batch_size=atof(argv[i]); break;
      case 't': i++; kernel_parm->kernel_type=atol(argv[i]); break;
      case 'd': i++; kernel_parm->poly_degree=atol(argv[i]); break;
      case 'g': i++; kernel_parm->rbf_gamma=atof(argv[i]); break;
      case 's': i++; kernel_parm->coef_lin=atof(argv[i]); break;
      case 'r': i++; kernel_parm->coef_const=atof(argv[i]); break;
      case 'u': i++; strcpy(kernel_parm->custom,argv[i]); break;
      case '-': strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);i++; strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);break; 
      case 'v': i++; (*struct_verbosity)=atol(argv[i]); break;
      case 'y': i++; (*verbosity)=atol(argv[i]); break;
      default: printf("\nUnrecognized option (svm_struct_main) %s!\n\n",argv[i]);
	       print_help();
	       exit(0);
      }
  }
#ifndef COMPILE_MEX_INTERFACE
  if(i>=argc) {
    printf("\nNot enough input parameters!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  strcpy (trainfile, argv[i]);
  if((i+1)<argc) {
    strcpy (modelfile, argv[i+1]);
  }
#endif
  if(learn_parm->svm_iter_to_shrink == -9999) {
    learn_parm->svm_iter_to_shrink=100;
  }

  if((learn_parm->skip_final_opt_check) 
     && (kernel_parm->kernel_type == LINEAR)) {
    printf("\nIt does not make sense to skip the final optimality check for linear kernels.\n\n");
    learn_parm->skip_final_opt_check=0;
  }    
  if((learn_parm->skip_final_opt_check) 
     && (learn_parm->remove_inconsistent)) {
    printf("\nIt is necessary to do the final optimality check when removing inconsistent \nexamples.\n");
    wait_any_key();
    print_help();
    exit(0);
  }    
  if((learn_parm->svm_maxqpsize<2)) {
    printf("\nMaximum size of QP-subproblems not in valid range: %ld [2..]\n",learn_parm->svm_maxqpsize); 
    wait_any_key();
    print_help();
    exit(0);
  }
  if((learn_parm->svm_maxqpsize<learn_parm->svm_newvarsinqp)) {
    printf("\nMaximum size of QP-subproblems [%ld] must be larger than the number of\n",learn_parm->svm_maxqpsize); 
    printf("new variables [%ld] entering the working set in each iteration.\n",learn_parm->svm_newvarsinqp); 
    wait_any_key();
    print_help();
    exit(0);
  }
  if(learn_parm->svm_iter_to_shrink<1) {
    printf("\nMaximum number of iterations for shrinking not in valid range: %ld [1,..]\n",learn_parm->svm_iter_to_shrink);
    wait_any_key();
    print_help();
    exit(0);
  }
  if(struct_parm->C<0) {
    printf("\nYou have to specify a value for the parameter '-c' (C>0)!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(((*alg_type) < 0) || (((*alg_type) > 5) && ((*alg_type) != 9))) {
    printf("\nAlgorithm type must be either '0', '1', '2', '3', '4', or '9'!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(learn_parm->transduction_posratio>1) {
    printf("\nThe fraction of unlabeled examples to classify as positives must\n");
    printf("be less than 1.0 !!!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(learn_parm->svm_costratio<=0) {
    printf("\nThe COSTRATIO parameter must be greater than zero!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(struct_parm->epsilon<=0) {
    printf("\nThe epsilon parameter must be greater than zero!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if((struct_parm->ccache_size<=0) && ((*alg_type) == 4)) {
    printf("\nThe cache size must be at least 1!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(((struct_parm->batch_size<=0) || (struct_parm->batch_size>100))  
     && ((*alg_type) == 4)) {
    printf("\nThe batch size must be in the interval ]0,100]!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if((struct_parm->slack_norm<1) || (struct_parm->slack_norm>2)) {
    printf("\nThe norm of the slacks must be either 1 (L1-norm) or 2 (L2-norm)!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if((struct_parm->loss_type != SLACK_RESCALING) 
     && (struct_parm->loss_type != MARGIN_RESCALING)) {
    printf("\nThe loss type must be either 1 (slack rescaling) or 2 (margin rescaling)!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(learn_parm->rho<0) {
    printf("\nThe parameter rho for xi/alpha-estimates and leave-one-out pruning must\n");
    printf("be greater than zero (typically 1.0 or 2.0, see T. Joachims, Estimating the\n");
    printf("Generalization Performance of an SVM Efficiently, ICML, 2000.)!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if((learn_parm->xa_depth<0) || (learn_parm->xa_depth>100)) {
    printf("\nThe parameter depth for ext. xi/alpha-estimates must be in [0..100] (zero\n");
    printf("for switching to the conventional xa/estimates described in T. Joachims,\n");
    printf("Estimating the Generalization Performance of an SVM Efficiently, ICML, 2000.)\n");
    wait_any_key();
    print_help();
    exit(0);
  }

  parse_struct_parameters(struct_parm);
}

void wait_any_key()
{
#ifndef COMPILE_MEX_INTERFACE
  printf("\n(more)\n");
  (void)getc(stdin);
#endif
}

void print_help()
{
  printf("\nSVM-struct learning module: %s, %s, %s\n",INST_NAME,INST_VERSION,INST_VERSION_DATE);
  printf("   includes SVM-struct %s for learning complex outputs, %s\n",STRUCT_VERSION,STRUCT_VERSION_DATE);
  printf("   includes SVM-light %s quadratic optimizer, %s\n",VERSION,VERSION_DATE);
  copyright_notice();
  printf("   usage: svm_struct_learn [options] example_file model_file\n\n");
  printf("Arguments:\n");
  printf("         example_file-> file with training data\n");
  printf("         model_file  -> file to store learned decision rule in\n");

  printf("General Options:\n");
  printf("         -?          -> this help\n");
  printf("         -v [0..3]   -> verbosity level (default 1)\n");
  printf("         -y [0..3]   -> verbosity level for svm_light (default 0)\n");
  printf("Learning Options:\n");
  printf("         -c float    -> C: trade-off between training error\n");
  printf("                        and margin (default 0.01)\n");
  printf("         -p [1,2]    -> L-norm to use for slack variables. Use 1 for L1-norm,\n");
  printf("                        use 2 for squared slacks. (default 1)\n");
  printf("         -o [1,2]    -> Rescaling method to use for loss.\n");
  printf("                        1: slack rescaling\n");
  printf("                        2: margin rescaling\n");
  printf("                        (default %d)\n",DEFAULT_RESCALING);
  printf("         -l [0..]    -> Loss function to use.\n");
  printf("                        0: zero/one loss\n");
  printf("                        ?: see below in application specific options\n");
  printf("                        (default %d)\n",DEFAULT_LOSS_FCT);
  printf("Optimization Options (see [2][5]):\n");
  printf("         -w [0,..,9] -> choice of structural learning algorithm (default %d):\n",(int)DEFAULT_ALG_TYPE);
  printf("                        0: n-slack algorithm described in [2]\n");
  printf("                        1: n-slack algorithm with shrinking heuristic\n");
  printf("                        2: 1-slack algorithm (primal) described in [5]\n");
  printf("                        3: 1-slack algorithm (dual) described in [5]\n");
  printf("                        4: 1-slack algorithm (dual) with constraint cache [5]\n");
  printf("                        9: custom algorithm in svm_struct_learn_custom.c\n");
  printf("         -e float    -> epsilon: allow that tolerance for termination\n");
  printf("                        criterion (default %f)\n",DEFAULT_EPS);
  printf("         -k [1..]    -> number of new constraints to accumulate before\n"); 
  printf("                        recomputing the QP solution (default 100)\n");
  printf("                        (-w 0 and 1 only)\n");
  printf("         -f [5..]    -> number of constraints to cache for each example\n");
  printf("                        (default 5) (used with -w 4)\n");
  printf("         -b [1..100] -> percentage of training set for which to refresh cache\n");
  printf("                        when no epsilon violated constraint can be constructed\n");
  printf("                        from current cache (default 100%%) (used with -w 4)\n");
  printf("SVM-light Options for Solving QP Subproblems (see [3]):\n");
  printf("         -n [2..q]   -> number of new variables entering the working set\n");
  printf("                        in each svm-light iteration (default n = q). \n");
  printf("                        Set n < q to prevent zig-zagging.\n");
  printf("         -m [5..]    -> size of svm-light cache for kernel evaluations in MB\n");
  printf("                        (default 40) (used only for -w 1 with kernels)\n");
  printf("         -h [5..]    -> number of svm-light iterations a variable needs to be\n"); 
  printf("                        optimal before considered for shrinking (default 100)\n");
  printf("         -# int      -> terminate svm-light QP subproblem optimization, if no\n");
  printf("                        progress after this number of iterations.\n");
  printf("                        (default 100000)\n");
  printf("Kernel Options:\n");
  printf("         -t int      -> type of kernel function:\n");
  printf("                        0: linear (default)\n");
  printf("                        1: polynomial (s a*b+c)^d\n");
  printf("                        2: radial basis function exp(-gamma ||a-b||^2)\n");
  printf("                        3: sigmoid tanh(s a*b + c)\n");
  printf("                        4: user defined kernel from kernel.h\n");
  printf("         -d int      -> parameter d in polynomial kernel\n");
  printf("         -g float    -> parameter gamma in rbf kernel\n");
  printf("         -s float    -> parameter s in sigmoid/poly kernel\n");
  printf("         -r float    -> parameter c in sigmoid/poly kernel\n");
  printf("         -u string   -> parameter of user defined kernel\n");
  printf("Output Options:\n");
  printf("         -a string   -> write all alphas to this file after learning\n");
  printf("                        (in the same order as in the training set)\n");
  printf("Application-Specific Options:\n");
  print_struct_help();
  wait_any_key();

  printf("\nMore details in:\n");
  printf("[1] T. Joachims, Learning to Align Sequences: A Maximum Margin Aproach.\n");
  printf("    Technical Report, September, 2003.\n");
  printf("[2] I. Tsochantaridis, T. Joachims, T. Hofmann, and Y. Altun, Large Margin\n");
  printf("    Methods for Structured and Interdependent Output Variables, Journal\n");
  printf("    of Machine Learning Research (JMLR), Vol. 6(Sep):1453-1484, 2005.\n");
  printf("[3] T. Joachims, Making Large-Scale SVM Learning Practical. Advances in\n");
	printf("    Kernel Methods - Support Vector Learning, B. Sch�lkopf and C. Burges and\n");
  printf("    A. Smola (ed.), MIT Press, 1999.\n");
  printf("[4] T. Joachims, Learning to Classify Text Using Support Vector\n");
  printf("    Machines: Methods, Theory, and Algorithms. Dissertation, Kluwer,\n");
  printf("    2002.\n");
  printf("[5] T. Joachims, T. Finley, Chun-Nam Yu, Cutting-Plane Training of Structural\n");
  printf("    SVMs, Machine Learning Journal, to appear.\n");
}



