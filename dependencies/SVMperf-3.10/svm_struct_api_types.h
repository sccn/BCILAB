/***********************************************************************/
/*                                                                     */
/*   svm_struct_api_types.h (instantiated for SVM-perform)             */
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

#ifndef svm_struct_api_types
#define svm_struct_api_types

# include "svm_common.h"
# include "svm_learn.h"

# define INST_NAME          "SVM-perf"
# define INST_VERSION       "V3.00"
# define INST_VERSION_DATE  "15.07.2009"

/* Identifiers for loss functions */
#define ZEROONE      0
#define FONE         1
#define ERRORRATE    2
#define PRBEP        3
#define PREC_K       4
#define REC_K        5
#define SWAPPEDPAIRS 10
#define AVGPREC      11

/* default precision for solving the optimization problem */
# define DEFAULT_EPS         0.1 
/* default loss rescaling method: 1=slack_rescaling, 2=margin_rescaling */
# define DEFAULT_RESCALING   2
/* default loss function:    */
# define DEFAULT_LOSS_FCT    ERRORRATE
/* default optimization algorithm to use: */
# define DEFAULT_ALG_TYPE    9
/* store Psi(x,y) once instead of recomputing it every time: */
# define USE_FYCACHE         1
/* decide whether to evaluate sum before storing vectors in constraint
   cache: 
   0 = NO, 
   1 = YES (best, if sparse vectors and long vector lists), 
   2 = YES (best, if short vector lists),
   3 = YES (best, if dense vectors and long vector lists) */
# define COMPACT_CACHED_VECTORS 1
/* minimum absolute value below which values in sparse vectors are
   rounded to zero. Values are stored in the FVAL type defined in svm_common.h 
   RECOMMENDATION: assuming you use FVAL=float, use 
     10E-15 if COMPACT_CACHED_VECTORS is 1 
     10E-10 if COMPACT_CACHED_VECTORS is 2 or 3 
*/
# define COMPACT_ROUNDING_THRESH 10E-15


typedef struct pattern {
	/* this defines the x-part of a training example, e.g. the structure
     for storing a natural language sentence in NLP parsing */
	DOC **doc;   /* set of example vectors */
	int totdoc;  /* size of set */
} PATTERN;

typedef struct label {
	/* this defines the y-part (the label) of a training example,
     e.g. the parse tree of the corresponding sentence. */
	double *class; /* vector of labels */
	int totdoc;    /* size of set */
} LABEL;

typedef struct structmodel {
	double *w;          /* pointer to the learned weights */
	MODEL  *svm_model;  /* the learned SVM model */
	long   sizePsi;     /* maximum number of weights in w */
  double walpha;
	/* other information that is needed for the stuctural model can be
     added here, e.g. the grammar rules for NLP parsing */
	long   sparse_kernel_type;   /* Selects the kernel type for sparse
	 kernel approximation. The values are
	 the same as for the -t
	 option. LINEAR (i.e. 0) means that
	 sparse kernel approximation is not
	 used. */
	long expansion_size;         /* Number of vectors in sparse kernel
	 expansion */
	DOC **expansion;             /* Vectors in sparse kernel expansion */
  long reducedset_size;        /* Number of vectors in reduced set
				  expansion */
  SVECTOR **reducedset;        /* Vectors in reduced set expansion */
  float **reducedset_kernel;   /* Kernel values between reduced set
				  expansion and training examples */
  MATRIX *reducedset_gram;     /* Gram matrix of reduced set expansion */
  MATRIX *reducedset_cholgram; /* Cholesky decomposition of Gram
				  matrix of reduced set expansion */
	MATRIX *invL;                /* Inverse of Cholesky decomposition of
	 Gram matrix over the vectors in the
	 sparse kernel expansion */
} STRUCTMODEL;

typedef struct struct_learn_parm {
	double epsilon;              /* precision for which to solve
	 quadratic program */
	double newconstretrain;      /* number of new constraints to
	 accumulate before recomputing the QP
	 solution */
	int    ccache_size;          /* maximum number of constraints to
	 cache for each example (used in w=4
	 algorithm) */
  double batch_size;           /* size of the mini batches in percent
				  of training set size (used in w=4
				  algorithm) */
	double C;                    /* trade-off between margin and loss */
	char   custom_argv[20][300]; /* string set with the -u command line option */
	int    custom_argc;          /* number of -u command line options */
	int    slack_norm;           /* norm to use in objective function
	 for slack variables; 1 -> L1-norm, 
	 2 -> L2-norm */
	int    loss_type;            /* selected loss rescaline type from -r
	 command line option. Select between
	 slack rescaling (1) and margin
	 rescaling (2) */
	int    loss_function;        /* select between different loss
	 functions via -l command line
	 option */
	/* further parameters that are passed to init_struct_model() */
	int num_features;
  int    truncate_fvec;        /* should test example vectors be
				  truncated to the number of features
				  seen in the training data? */
	double bias;                 /* value for artificial bias feature */
	long   bias_featurenum;      /* id number of bias feature */
	double prec_rec_k_frac;      /* fraction of training set size to use
	 as value of k for Prec@k and
	 Rec@k */
	long   sparse_kernel_type;   /* Same as -t option */
	long   sparse_kernel_size;   /* Number of basis functions to select
	 from the set of basis functions for
	 training with approximate kernel
	 expansion.  */
	char   sparse_kernel_file[999];/* File that contains set of basis
	 functions for training with
	 approximate kernel expansion.  */
  int    sparse_kernel_method; /* method for selecting sparse kernel
				  subspace (1 random sampling, 2
				  incomplete cholesky) */
	int    shrinking;            /* Selects whether shrinking heuristic
	 is used in the custom algorithm for
	 minimizing error rate. */
  double rset_precision;       /* (only used internally) minimum
				  improvement in euclidian distance so
				  that the preimage is added to the
				  reduced set expansion */
  int    preimage_method;      /* method to use for finding preimages */
  int    recompute_rset;       /* selects whether the reduced set is
				  recomputed and the method restarted
				  again, after it has converged for
				  the first time */
  int    classify_dense;       /* uses a dense vector representation
				  when classifying new examples in
				  svm_perf_classify. This uses more
				  memory, but is faster if the support
				  vectors in the model are dense. */
} STRUCT_LEARN_PARM;

typedef struct struct_test_stats {
	/* you can add variables for keeping statistics when evaluating the
     test predictions in svm_struct_classify. This can be used in the
     function eval_prediction and print_struct_testing_stats. */
	long   test_data_unlabeled;
	double errorrate;
	double precision;
	double recall;
	double fone;
	double prbep;
	double rocarea;
	double avgprec;
} STRUCT_TEST_STATS;

typedef struct struct_id_score {
	int id;
	double score;
	double tiebreak;
} STRUCT_ID_SCORE;

#endif
