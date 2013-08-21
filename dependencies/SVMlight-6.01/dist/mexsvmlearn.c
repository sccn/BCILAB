/*
* =============================================================
* mexsvmlearn - link MATLAB with SVM-Lite by Thorstenn Joachims.
* This is a MEX-file for MATLAB.  
* MEX-Interface:
*    Thomas Briggs (c) 2004
* SVM-Lite:
*    Thorsten Joachims
* =============================================================
*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "svm_common.h"
#include "svm_learn.h"
#include "mexsvmlearn.h"
#include "global.h"
#include "mexcommon.h"

char docfile[200];           /* file with training examples */
char modelfile[200];         /* file for resulting classifier */
char restartfile[200];       /* file with initial alphas */

/* call as  model = mexsvmlearn(data,labels,options) */
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	char **argv;
	int argc;
	DOC **docs;  /* training examples */
	long totwords,totdoc,i;
	double *target;
	double *alpha_in=NULL;
	KERNEL_CACHE *kernel_cache;
	LEARN_PARM learn_parm;
	KERNEL_PARM kernel_parm;
	MODEL model;

	/* check for valid calling format */
	if ((nrhs != 3)  || (nlhs != 1))
		mexErrMsgTxt(ERR001);

	if (mxGetM(prhs[0]) != mxGetM(prhs[1]))
		mexErrMsgTxt(ERR002);

	if (mxGetN(prhs[1]) != 1)
		mexErrMsgTxt(ERR003);

	/* reset static variables -- as a .DLL, static things are sticky  */
	global_init( );

	/* convert the parameters (given in prhs[2]) into an argv/argc combination */
	argv = make_argv((mxArray *)prhs[2],&argc); /* send the options */



	/* this was originally supposed to be argc, argv, re-written for MATLAB ...
	its cheesy - but it workss :: convert the options array into an argc, 
	argv pair  and let svm_lite handle it from there. */

	read_input_parameters(argc,argv,docfile,modelfile,restartfile,&verbosity, 
		&learn_parm,&kernel_parm);

	extract_user_opts((mxArray *)prhs[2], &kernel_parm);

	totdoc = mxGetM(prhs[0]);
	totwords = mxGetN(prhs[0]);

	/* prhs[0] = samples (mxn) array
	prhs[1] = labels (mx1) array */
	mexToDOC((mxArray *)prhs[0], (mxArray *)prhs[1], &docs, &target, NULL, NULL);

	/* TODO modify to accept this array 
	if(restartfile[0]) alpha_in=read_alphas(restartfile,totdoc); */

	if(kernel_parm.kernel_type == LINEAR) { /* don't need the cache */
		kernel_cache=NULL;
	}
	else {
		/* Always get a new kernel cache. It is not possible to use the
		same cache for two different training runs */
		kernel_cache=kernel_cache_init(totdoc,learn_parm.kernel_cache_size);
	}


	if(learn_parm.type == CLASSIFICATION) {
		svm_learn_classification(docs,target,totdoc,totwords,&learn_parm,
			&kernel_parm,kernel_cache,&model,alpha_in);

	}
	else if(learn_parm.type == REGRESSION) {
		svm_learn_regression(docs,target,totdoc,totwords,&learn_parm,
			&kernel_parm,&kernel_cache,&model);
	}
	else if(learn_parm.type == RANKING) {
		svm_learn_ranking(docs,target,totdoc,totwords,&learn_parm,
			&kernel_parm,&kernel_cache,&model);
	}
	else if(learn_parm.type == OPTIMIZATION) {
		svm_learn_optimization(docs,target,totdoc,totwords,&learn_parm,
			&kernel_parm,kernel_cache,&model,alpha_in);
	}
	else {
		mexErrMsgTxt(ERR004);
	}

	if(kernel_cache) {
		/* Free the memory used for the cache. */
		kernel_cache_cleanup(kernel_cache);
	}

	/* **********************************
	* After the training/learning portion has finished,
	* copy the model back to the output arrays for MATLAB 
	* ********************************** */
	store_model(&model, plhs);


	/* Warning: The model contains references to the original data 'docs'.
	If you want to free the original data, and only keep the model, you 
	have to make a deep copy of 'model'. */
	/* deep_copy_of_model=copy_model(model); */
	/*write_model(modelfile,model); */

	/* my_free(alpha_in);
	freeModelMem(&model);
	for(i=0;i<totdoc;i++) 
		free_example(docs[i],1);
	my_free(docs);
	my_free(target);   */

	global_destroy( );	


}





/*---------------------------------------------------------------------------*/

void read_input_parameters(int argc,char **argv,char *docfile,char *modelfile,
			   char *restartfile,long *verbosity,
			   LEARN_PARM *learn_parm,KERNEL_PARM *kernel_parm)
{
	long i;
	char type[100];
	char msg[512], msg1[255],msg2[255];


	/* set default */
	strcpy (modelfile, "svm_model");
	strcpy (learn_parm->predfile, "trans_predictions");
	strcpy (learn_parm->alphafile, "");
	strcpy (restartfile, "");
	(*verbosity)=1;

	learn_parm->biased_hyperplane=1;
	learn_parm->sharedslack=0;
	learn_parm->remove_inconsistent=0;
	learn_parm->skip_final_opt_check=0;
	learn_parm->svm_maxqpsize=10;
	learn_parm->svm_newvarsinqp=0;
	learn_parm->svm_iter_to_shrink=-9999;
	learn_parm->maxiter=100000;
	learn_parm->kernel_cache_size=40;
	learn_parm->svm_c=0.0;
	learn_parm->eps=0.1;
	learn_parm->transduction_posratio=-1.0;
	learn_parm->svm_costratio=1.0;
	learn_parm->svm_costratio_unlab=1.0;
	learn_parm->svm_unlabbound=1E-5;
	learn_parm->epsilon_crit=0.001;
	learn_parm->epsilon_a=1E-15;
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
		case 'z': i++; strcpy(type,argv[i]); break;
		case 'v': i++; (*verbosity)=atol(argv[i]); break;
		case 'b': i++; learn_parm->biased_hyperplane=atol(argv[i]); break;
		case 'i': i++; learn_parm->remove_inconsistent=atol(argv[i]); break;
		case 'f': i++; learn_parm->skip_final_opt_check=!atol(argv[i]); break;
		case 'q': i++; learn_parm->svm_maxqpsize=atol(argv[i]); break;
		case 'n': i++; learn_parm->svm_newvarsinqp=atol(argv[i]); break;
		case '#': i++; learn_parm->maxiter=atol(argv[i]); break;
		case 'h': i++; learn_parm->svm_iter_to_shrink=atol(argv[i]); break;
		case 'm': i++; learn_parm->kernel_cache_size=atol(argv[i]); break;
		case 'c': i++; learn_parm->svm_c=atof(argv[i]); break;
		case 'w': i++; learn_parm->eps=atof(argv[i]); break;
		case 'p': i++; learn_parm->transduction_posratio=atof(argv[i]); break;
		case 'j': i++; learn_parm->svm_costratio=atof(argv[i]); break;
		case 'e': i++; learn_parm->epsilon_crit=atof(argv[i]); break;
		case 'o': i++; learn_parm->rho=atof(argv[i]); break;
		case 'k': i++; learn_parm->xa_depth=atol(argv[i]); break;
		case 'x': i++; learn_parm->compute_loo=atol(argv[i]); break;
		case 't': i++; kernel_parm->kernel_type=atol(argv[i]); break;
		case 'd': i++; kernel_parm->poly_degree=atol(argv[i]); break;
		case 'g': i++; kernel_parm->rbf_gamma=atof(argv[i]); break;
		case 's': i++; kernel_parm->coef_lin=atof(argv[i]); break;
		case 'r': i++; kernel_parm->coef_const=atof(argv[i]); break;
		case 'l': i++; strcpy(learn_parm->predfile,argv[i]); break;
		case 'a': i++; strcpy(learn_parm->alphafile,argv[i]); break;
		case 'y': i++; strcpy(restartfile,argv[i]); break;
		default: 
			sprintf(msg,"\nUnrecognized option %s!\n\n",argv[i]);
			mexErrMsgTxt(msg);
		}
	}


	if ((i>1) && (i>argc)) {
		printf ("Errr?  i=%d, argc=%d\n", i, argc);
		mexErrMsgTxt("Not enough input parameters!\n\n");
	}

	/* don't read a docfile!!
	/*  strcpy (docfile, argv[i]);
	if((i+1)<argc) {
	strcpy (modelfile, argv[i+1]);
	}
	*/

	if(learn_parm->svm_iter_to_shrink == -9999) {
		if(kernel_parm->kernel_type == LINEAR) 
			learn_parm->svm_iter_to_shrink=2;
		else
			learn_parm->svm_iter_to_shrink=100;
	}
	if(strcmp(type,"c")==0) {
		learn_parm->type=CLASSIFICATION;
	}
	else if(strcmp(type,"r")==0) {
		learn_parm->type=REGRESSION;
	}
	else if(strcmp(type,"p")==0) {
		learn_parm->type=RANKING;
	}
	else if(strcmp(type,"o")==0) {
		learn_parm->type=OPTIMIZATION;
	}
	else if(strcmp(type,"s")==0) {
		learn_parm->type=OPTIMIZATION;
		learn_parm->sharedslack=1;
	}
	else {
		sprintf(msg,"\nUnknown type '%s': Valid types are 'c' (classification), 'r' regession, and 'p' preference ranking.\n",type);
		mexErrMsgTxt(msg);
	}    
	if((learn_parm->skip_final_opt_check) 
		&& (kernel_parm->kernel_type == LINEAR)) {
			printf("\nIt does not make sense to skip the final optimality check for linear kernels.\n\n");
			learn_parm->skip_final_opt_check=0;
		}    
		if((learn_parm->skip_final_opt_check) 
			&& (learn_parm->remove_inconsistent)) {
				mexErrMsgTxt("It is necessary to do the final optimality check when removing inconsistent \nexamples.\n");

			}    
			if((learn_parm->svm_maxqpsize<2)) {
				sprintf(msg,"\nMaximum size of QP-subproblems not in valid range: %ld [2..]\n",learn_parm->svm_maxqpsize); 
				mexErrMsgTxt(msg);
			}
			if((learn_parm->svm_maxqpsize<learn_parm->svm_newvarsinqp)) {
				sprintf(msg1,"Maximum size of QP-subproblems [%ld] must be larger than the number of\n",learn_parm->svm_maxqpsize); 
				sprintf(msg2,"new variables [%ld] entering the working set in each iteration.\n",learn_parm->svm_newvarsinqp); 
				sprintf(msg,"%s\n%s\n", msg1,msg2);
				mexErrMsgTxt(msg);
			}
			if(learn_parm->svm_iter_to_shrink<1) {
				sprintf(msg,"\nMaximum number of iterations for shrinking not in valid range: %ld [1,..]\n",learn_parm->svm_iter_to_shrink);
				mexErrMsgTxt(msg);
			}
			if(learn_parm->svm_c<0) {
				mexErrMsgTxt("\nThe C parameter must be greater than zero!\n\n");
			}
			if(learn_parm->transduction_posratio>1) {
				mexErrMsgTxt("\nThe fraction of unlabeled examples to classify as positives must be less than 1.0 !!!\n\n");
			}
			if(learn_parm->svm_costratio<=0) {
				mexErrMsgTxt("\nThe COSTRATIO parameter must be greater than zero!\n\n");
			}
			if(learn_parm->epsilon_crit<=0) {
				mexErrMsgTxt("\nThe epsilon parameter must be greater than zero!\n\n");
			}
			if(learn_parm->rho<0) {
				printf("\nThe parameter rho for xi/alpha-estimates and leave-one-out pruning must\n");
				printf("be greater than zero (typically 1.0 or 2.0, see T. Joachims, Estimating the\n");
				printf("Generalization Performance of an SVM Efficiently, ICML, 2000.)!\n\n");
				mexErrMsgTxt("ending.\n");
			}
			if((learn_parm->xa_depth<0) || (learn_parm->xa_depth>100)) {
				printf("\nThe parameter depth for ext. xi/alpha-estimates must be in [0..100] (zero\n");
				printf("for switching to the conventional xa/estimates described in T. Joachims,\n");
				printf("Estimating the Generalization Performance of an SVM Efficiently, ICML, 2000.)\n");
				mexErrMsgTxt("ending\n");
			}
}


void extract_user_opts(mxArray *options, KERNEL_PARM *kernel_parm)
{
	char *buf;
	int bufflen;
	char *user;

	if (! mxIsChar(options) )
		mexErrMsgTxt("Error: options not specified as a string.");


	bufflen = mxGetNumberOfElements(options)+1;
	buf = my_malloc(bufflen * sizeof(char));
	myBzero(buf, bufflen);
	mxGetString(options,buf, bufflen);

	/* -------  buf now contains the given user arguments ----------- */
	/* look for the -u parameter (if here) and strip it out */
	user = strstr(buf,"-u");
	if (user != NULL) {
		strcpy((*kernel_parm).custom,user+2);
	}

} 


char **make_argv(mxArray *options, int *argc)
{

	const int *dims;
	int i;
	char *buf,*ptr,*ptr2 ;
	int number_of_dimensions, bufflen;
	char **argv;
	char *user;


	argv = (char **)my_malloc(255 * sizeof(char *));

	if (! mxIsChar(options) )
		mexErrMsgTxt("Error: options not specified as a string.");



	bufflen = mxGetNumberOfElements(options)+1;
	buf = my_malloc(bufflen * sizeof(char));
	myBzero(buf, bufflen);
	mxGetString(options,buf, bufflen);

	/* -------  buf now contains the given user arguments ----------- */
	/* look for the -u parameter (if here) and strip it out */
	user = strstr(buf,"-u");
	if (user != NULL) {
		int userlen = strlen(user);
		int newlen = user - buf;
		char *newbuff = my_malloc(newlen * sizeof(char));

		strncpy(newbuff,buf,newlen);
		newbuff[newlen] = 0x00;
		buf = newbuff;
	}

	/* -------------------------------------------------------------- */


	argv[0] = my_malloc(13 * sizeof(char));
	myBzero(argv[0],13);

	i = 1;
	ptr = strtok(buf," ");
	while ((ptr != NULL) && (i<200))  {

		ptr2 = (char *)my_malloc((strlen(ptr)+1) * sizeof(char));
		strcpy (ptr2, ptr);

		argv[i] = ptr2;
		i++;

		ptr = strtok(NULL," ");
	}
	argv[i] = 0x00;
	*argc = i;

	if (user != NULL)
		my_free(buf);

	return argv;
}

void freeArgv(int argc, char **argv)
{
	int i;

	for (i = 0; i < argc; i++)
		my_free(argv[i]);

	my_free(argv);
}


void freeModelMem(MODEL *model)
{
	long i;

	if(model->supvec) {
		for(i=1;i<model->sv_num;i++) {
			if (model->supvec[i])
				free_example(model->supvec[i],1);
		}
	}
}
