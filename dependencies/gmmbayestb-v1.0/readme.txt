                         GMMBayes Toolbox
     Gaussian mixture model learning and Bayesian classification
                            Version 1.0


Authors: Pekka Paalanen <pekka.paalanen@lut.fi>
         Joni Kämäräinen <joni.kamarainen@lut.fi>
	 Jarmo Ilonen <jarmo.ilonen@lut.fi>

Website: http://www.it.lut.fi/project/gmmbayes/

A document describing the used algorithms and test results can be
found on the website. Note that the document is specific to the
GMMBayes Toolbox version 0.1.

Information about density quantile can be found in
   [1] Paalanen, P., Kamarainen, J.-K., Ilonen, J., Kälviäinen, H.,
    Feature Representation and Discrimination Based on Gaussian Mixture Model
    Probability Densities - Practices and Algorithms, Research Report 95,
    Lappeenranta University of Technology, Department of Information
    Technology, 2005.

For description of package contents see Contents.m file.
Descriptions of used data structures are at the end of this file.
All Matlab m-files should contain sufficient documentation.
For acknowledgements and copyright information see Contents.m file.


Contents
1. Major changes
2. The new classifier API
3. Changelog
4. Data structures



                        --- 1. MAJOR CHANGES ---

- The gmmb_create and gmmb_classify APIs have changed in version 0.2 and 0.3.

- A new API has replaced the gmmb_classify interface, which remains
  as an unsupported legacy interface, in version 0.3.

- A new fractile system in version 0.3. Fractiles are called density
  quantiles in version 1.0.


                    --- 2. THE NEW CLASSIFIER API ---

Let's assume training is already done and you have a bayesS structure.
How to proceed to classify new data?

First we compute PDF values for each class from each data point:
pdfmat = gmmb_pdf(data, bayesS);

If you are doing the usual Bayesian classification, you will
apply priors, normalize the sums and make the decisions:
postprob = gmmb_normalize( gmmb_weightprior(pdfmat, bayesS) );
labels = gmmb_decide(postprob);

The resulting labels vector will contain class indexes, as coming
from bayesS struct array. If you instead want to minimize the Bayes risk,
you can multiply the postprob matrix with a cost matrix (note that rows
correspond to data points) and get the labels:
[dummy, labels] = min(postprob*cost, [], 2);

The new API supports also density quantiles. For quantile related operations
you will need a histS structure matching the used PDFs; either generated
from bayesS structure (recommended), or computed from training data
with the gmmb_hist function.

histS structure can be generated like this:
histS = gmmb_generatehist(bayesS, N);

Since histS is a vector approximation of several functions, you need to
tell the amount of sampling: the number of generated random points N.
N should be fairly big, at least in the order of thousands. The bigger N,
the better approximation, but also bigger memory consumption.

You can convert PDF values to density quantile values:
fracmat = gmmb_lhood2frac(histS, pdfmat);

You can threshold the data points with respect to uniform quantile value,
i.e., you can identify outliers with respect to each distribution:
outlier_mask = gmmb_fracthresh(pdfmat, histS, thr);

If you want to do the Bayesian classification with thresholding,
you can now do:
postprob(outlier_mask) = 0;
labels = gmmb_decide(postprob);



                        --- 3. CHANGELOG ---


    ***   FINAL RELEASE GMMBayes Toolbox v1.0   ***

- Changed the documentation to correspond the terminology in [1]
	(in the beginning of this readme). Fractiles are now called density
	quantiles.

- Software cooldown period is over, this is the final version.


    ***   BUG FIXES in GMMBayes Toolbox v0.3   ***

- The complex Gaussian PDF was incorrect; it is now fixed. Also the real
	and complex cases need to be handled separately where applicable.
	NOTE: the complex range is detected from training data or the bayesS
	struct field 'mu'. Complex training data should have abs(Im)>0,
	otherwise it is truncated to real values.

- Missed a warning() call in gmmb_covfixer.m, changed to warning_wrap().

- gmmb_fj: fixed variable names related to covfixer and verbose output.

- gmmbvl_kmeans: fixed several function call names.

- gmmb_pdf() has been changed to be more friendly to use.

- The field histS has been removed from the bayesS structure.
	The boolean parameter has been removed from gmmb_create API.

- Some tweaks in gmmbvl_em_step_partial and gmmbvl_kmeans to prevent
	the Matlab C compiler from whining. Still, it won't work
	ascompiled C application because of problems with the version()
	function.


   ***   NEW FEATURES in GMMBayes Toolbox v0.3   ***

- The new classifier API, described in the previous section.
	gmmb_demo01 uses this new API.

- gmmb_cmvnpdf.m has been replaced with our own optimized code.

- The interface gmmb.m and newgmmb.m have been removed.

- gmmb_classify: fractile thresholding is done before classification.
	This way the "trash class" is the last choice, not second.
	As gmmb_classify no longer relies on bayesS.histS, it will call
	gmmb_generatehist() to create a histS structure.

- gmmb_em, gmmb_fj and gmmb_gem:
	Functions will now print warnings about possibly insufficient amount
	of training data, but it does not directly imply that the estimation
	could fail. The limit is quite loose.

- gmmb_frac2lhood, gmmb_lhood2frac, gmmb_fracthresh:
	The new fractile system, performing mapping between PDF-values
	(thresholds) and fractiles (probability mass of areas having
	PDF-value over a threshold).
	The fractile system is good for finding outliers.

- gmmb_generatehist:
	Generate a histS structure from bayesS structure with random
	sampling. The fractile system needs a histS structure to work.



    ***   BUG FIXES in GMMBayes Toolbox v0.2   ***

- gmmb_em: the convergence threshold is now relative log-likelihood change,
	like in the other estimation methods.

- infinite loop fixes in several places including: covfixer2, gmmb_fj

- finity checks added in several places including: covfixer2, gmmb_em, gmmb_fj

- covfixer2 is partly rewritten: the case when a value on the diagonal
	is less or equal to zero. The minimum correction amount should now
	be large enough to fix the covariance matrix.

- use of covfixer2 in gmmb_em and gmmb_fj:
	covfixer2 can modify the covariance matrix so that the change of
	log-likelihood is unpredictable; this could lead to premature
	stop of the algorithm.
	Now if a covariance matrix is fixed, the threshold stop condition
	is not evaluated on that iteration. To prevent from going into
	infinite estimate/fix loop, a maximum of 20 consecutive
	covariance fixes are allowed; otherwise the algorithm is stopped.

- gmmb_demo01: when done classification, show the classification result,
	not the true classes of the test data.



   ***   NEW FEATURES in GMMBayes Toolbox v0.2   ***

- gmmb_fj: added option 'broken' that is enabled by default (the old behaviour).
	When enabled, the number of free parameters is computed based on
	real-valued data (means and covariance matrices).
	If disabled and if data is complex valued, use a different
	computation for the number of free parameters.
	See the project report on the website for more information about
	handling complex data.

- extended algorithm logging capabilities added to the following functions:
	covfixer2, gmmb_em, gmmb_fj, gmmb_create

- tested with Matlab R12, should be compatible.
	NOTE: one change has to be made into getargs.m if using
	Matlab R12; comment a line and uncomment another.
	It is documented there.

- added field 'hist' to bayesS struct to allow "trash class" classification

- Optional trash class classification.
	When using Bayesian classification, you can give a 'fractile' value
	that tells how much of a class' training data is to be classified into
	the class in one-class classification scenario. Value of 1 says that
	all training data should be classified as belonging to the class, but
	no sample that has a lower likelihood than the worst sample in the
	class training data should be classified as belonging to the class.
	In multiclass scenarios, the one-class classification is performed 
	first for each of the known classes, and the classes with which the
	sample passed, are used in Bayesian classification.

- gmmb_em: added two new c-means based EM initialization methods,
	so that EM estimation no longer requires the Fuzzy Logic Toolbox.

- gmmb_fj: added 'maxloops' optional parameter to control the maximum
	iterations allowed in a CEM run.

- gmmb_demo01: added an example of trash class classification.



                    --- 4. DATA STRUCTURES ---

Specification for the so called bayesS struct in gmmb_* functions.

K number of classes
D number of dimensions
C number of GMM components (class specific, can vary)
N number of data points

bayesS - struct array 1 x K
  mu D x C		mean vectors
  sigma D x D x C	covariance matrices
  weight C x 1		component weights
  apriories 1 x 1	class a priori probability


When referring to "data" and "type" parameters:
data N x D matrix
type N x 1 vector (class label, positive integer)
type is assumed to have at least one of every integer 1, ..., K.

histS - cell array 1 x K of vectors
  Each vector contains an ordered list of PDF values of
  points from the distribution k described by the k'th PDF.


$Name:  $
$Id: readme.txt,v 1.4 2005/04/14 10:33:34 paalanen Exp $
