function [signal,state] = flt_ica(varargin)
% Annotate the Signal with a spatial decomposition into independent components (using ICA)
% [Signal,State] = flt_ica(Signal, Variant, CleaningLevel, TransformData, OutputCleanedData, ResumePrevious)
%
% The IC-decomposed signal [1] can be considered to have "better" channels than the raw sensor
% signal, for several reasons:
%  * The components are less correlated than the raw channels, which allows for easier-to-handle
%    statistical assumptions in higher levels (e.g. using diagonal covariance matrices).
%  * The relevant information can be assumed to be localized in only a small subset of the
%    components, where the majority of components will carry less relevant data; this allows for the
%    assumption of (channel-wise) sparsity in the derived features (e.g. using l1-regularized
%    classifiers and/or feature extractors).
%  * The components have a higher degree of semantic meaning than channels (which needs to be
%    computed, though), such as the presence of eye artifacts, muscle artifacts, brain activity,
%    etc.
%  * A fraction of components can be localized in the brain using dipole fitting, beamforming,
%    sparse bayesian learning, and other methods, which gives access to semantics that can be
%    derived/computed from the component's location in the brain (e.g., postcentral gyrus -> high
%    chance of touch-related brain activity, etc.).
%
% In:
%   Signal     : a continuous data set
%
%   Variant    : type of ICA model to run (default: 'amica'), possible values are:
%                'amica' : AMICA [2] is currently the best-of-breed ICA method and generally
%                          preferred. (default) It uses a flexible model of source signal densities
%                          (generalized gaussian scale mixtures) which allow it to obtain better
%                          solutions for EEG data (in terms of statistical independence) than the
%                          other available methods; furthermore, multiple models can be learned,
%                          which allow it to capture non-stationarities in a principled manner. If
%                          available, a binary implementation is used. Options: help runamica
%                'infomax' : (extended) Infomax [3] is the second-best approach to ICA for EEG data,
%                            and may be faster than the MATLAB implementation of AMICA, or possibly
%                            easier to handle. If available, a binary implementation is used.
%                            Options: help runica
%                'beamica' : essentially Infomax in its default setting, but offers the option of 
%                            constraining solutions with the help of beamforming. The fastest ICA 
%                            implementation in the toolbox (especially when run on a fast GPU).
%                'fastica' : FastICA [4] is a relatively simple ICA implementation, which is
%                            converges relatively quickly, though at the expense of solution
%                            quality. In many cases it can not attain results as good as extended
%                            Infomax or AMICA, but for repeated computations (e.g. inside a
%                            cross-validation) it can be a reasonable choice due to its speed.
%                            Options: help fastica
%                'rica' : (overcomplete) Reconstruction ICA [7] is a novel fast ICA approach that allows
%                         to learn arbitrarily over-complete decompositions.
%                'kernelica' : KernelICA [5] has a similarly flexible model of source densities as
%                              AMICA, using a kernel approach, but requires massive computation
%                              time, so KernelICA can realistically at best be used on small data
%                              sets. Options: help kernel_ica_options
%                'sphere' : just the spering matrix (second-order)
%                'robust_spere' : the robust sphering matrix (estimated under super-Gaussian noise)
%                others, if in path: 'jader','jadeop','jade_td_p','MatlabshibbsR','tica','erica',
%                                    'simbec','unica','amuse','fobi','evd','evd24','sons','sobi',
%                                    'ng_ol','acsobiro','acrsobibpf','pearson_ica','egld_ica','eeA',
%                                    'tfbss','icaML','icaMS'.
%
%                Note: each of these variants has its own set of sub-parameters, which can be
%                specified by passing Variant as a cell array, e.g., {'amica', 'num_models',3, 'max_iter',1000}
%
%                For a full list of these parameters, review either the argument specification below
%                or the review/edit approach panel (under ICA)
%
%   DataCleaning : Parameters for data cleaning prior to running an EEG; this is a cell array of 
%                  parameters for the function flt_clean_settings. In the simplest case, it is just
%                  a setting string (e.g. 'seated', 'noisy', 'walking', 'running')
%
%   TransformData : whether to place the decomposition result in the actual channel data instead of
%                   in an additional annotation field of the output data set (.icaact). (default: false)
%
%   OutputCleanedData : whether to output the cleaned version of the data instead of the original
%                       version of the data; by default, the cleaned data is only used to compute a
%                       better decomposition, while the (decomposed) original data is what is returned
%                       (default: false)
%
%   ResumePrevious : whether to try to resume from a previous computation (on the same data), if
%                    applicable (default: true)
%
%   State      : state, for online updates (default: [])
%
% Out:
%   Signal : continuous data set annotated with an ICA decomposition, and optionally with data
%            transformed into IC activations
%
% Notes:
%   Only the first two arguments can be specified without passing them by name.
%
% Examples:
%   % annotate the data set with an ICA decomposition, using default settings
%   eeg = flt_ica(eeg)
%
%   % do an ICA decomposition using Infomax
%   eeg = flt_ica(eeg, 'infomax')
%
%   % do an ICA decomposition using Infomax and transform the actual channel data (so that
%   % channel-based methods end up operating on components)
%   eeg = flt_ica(eeg, 'infomax', 'TransformData',true)
%
%   % do an ICA decomposition using Infomax and pass a specific cleaning level
%   eeg = flt_ica(eeg, 'infomax', 'CleaningLevel','hardcore')
%
%   % do an ICA decomposition using Infomax, pass a specific cleaning level, and override some
%   parameters to infomax
%   eeg = flt_ica(eeg, {'infomax' 'MaxIterations',300}, 'CleaningLevel','hardcore')
%
%   % as before, but pass all arguments by name (recommended to avoid confusion)
%   eeg = flt_ica('Signal',eeg, 'Variant',{'infomax' 'MaxIterations',300}, 'CleaningLevel','hardcore')
%
%   % do an ICA decomposition using amica and specify some of its parameters (and use a custom cleaning level)
%   eeg = flt_ica('Signal',eeg, 'Variant',{'amica', 'NumModels',3, 'MaxIterations',1000}, 'CleaningLevel','hardcore')
%
%   % run an AMICA and don't try to use the cluster (i.e., run locally)
%   eeg = flt_ica('Signal',eeg, 'Variant',{'amica', 'UseGridEngine','off'})
%
%   % do a multi-model amica decomposition using 3 models using 16 slots on the cluster, and do some
%   % moderate artifact handling
%   eeg = flt_ica('Signal',eeg, 'Variant',{'amica', 'NumModels',3, 'NumProcessors',16}, 'CleaningLevel','strong')
%
%   % do a quick-and-dirty FastICA
%   eeg = flt_ica(eeg, {'fastica', 'MaxIterations',100})
%
%   % do a super-slow Kernel ICA on highly cleaned data
%   eeg = flt_ica(eeg,'kernelica','CleaningLevel','ultrahardcore')
%
%
% References:
%   [1] Makeig S., Bell A.J., Jung T-P. and Sejnowski T.J. 1995. "Independent Component Analysis of Electroencephalographic Data"
%       in Mozer M. et al (eds) Advances in Neural Information Processing Systems 8, MIT press
%   [2] J. A. Palmer, S. Makeig, K. Kreutz-Delgado, and B. D. Rao, "Newton Method for the ICA Mixture Model".
%       In Proceedings of the 33rd IEEE International Conference on Acoustics and Signal Processing (ICASSP 2008), Las Vegas, NV, pp. 1805-1808, 2008.
%   [3] Bell, A. J., and Sejnowski, T. J. "An information-maximization approach to blind separation and blind deconvolution."
%       Neural Comput. 7, 6 (1995), 1129?1159.
%   [4] A. Hyvaerinen. "Fast and Robust Fixed-Point Algorithms for Independent Component Analysis."
%       IEEE Transactions on Neural Networks 10(3):626-634, 1999.
%   [5] Francis R. Bach, Michael I. Jordan. "Kernel Independent Component Analysis",
%       Journal of Machine Learning Research, 3, 1-48, 2002
%   [6] H. Shen, S. Jegelka and A. Gretton. "Fast Kernel ICA using an approximate Newton method."
%       AISTATS 2007.
%   [7] Q.V. Le, A. Karpenko, J. Ngiam, A.Y. Ng. "ICA with Reconstruction Cost for Efficient Overcomplete Feature Learning."
%       NIPS 2011
%
% See also:
%   runica, runamica11, fastica, kernel_ica_options, flt_clean_channels, flt_clean_windows, flt_clean_peaks
%
% TODO:
%   Add robust and CUDA ICA's.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-17

% flt_ica_version<1.2.11> -- for the cache


if ~exp_beginfun('filter') return; end

% has its own highpass filter, sometimes applied on re-referenced data
declare_properties('name','ICA', 'precedes',{'flt_fir','flt_iir'}, 'follows','flt_reref', 'independent_trials',false, 'independent_channels',false);

arg_define([0 2],varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg_subswitch({'variant','Variant'},'infomax',{ ...
    'noica', {}, ...  % noICA is a testing variant with no arguments
    'amica', { ...
        arg({'amica_version','AmicaVersion','version'},'stable11',{'devel','stable11','stable12','stable','matlab'},'Amica version to use. The stable12 is cross-platform, stable11 has been tested extensively on the SCCN cluster, and the matlab version should run everywhere.','cat','Miscellaneous'), ...
        arg({'max_iter','MaxIterations'},2500,uint32([1 10000]),'Maximum number of iterations to perform.','cat','Core Parameters'), ...
        arg({'num_models','NumModels'},1, uint32([1 20]),'Number of models to learn. Per model, approx. 100,000 data points are needed at 100 channels.','cat','Core Parameters'), ...
        arg({'num_mix_comps','NumMixtureComponents'},3, uint32([1 15]),'Number of mixture components to learn. This is per source.','cat','Core Parameters'), ...
        arg({'pdftype','TypeOfPDF'},'GeneralizedGaussian',{'GeneralizedGaussian','ExtendedInfomax','Gaussian','Logistic'},'Probability density type for source model.','cat','Core Parameters'), ...
        arg({'share_comps','ShareComponents'},false,[],'Flag to share components when num_models > 1.','cat','Core Parameters'), ...
        ...
        arg({'max_threads','MaxThreads'},[],[],'Number of threads per node. If this is too high, churn kicks in (especially between users on a shared cluster).','cat','Compute Resources','shape','scalar'), ...
        arg({'useqsub','UseGridEngine','qsub'},quickif(isdeployed,'off','on'),{'on','off'},'Use Sun Grid Engine cluster. Assumes qsub; the alternative is to run locally.','cat','Compute Resources'), ...
        arg({'use_queue','UseQueue'},{'q3:32','q4:32','q5:32','q6:32','q7:32','q8:32','q9:32','q10:32','qa1:64','qa2:64','qa3:64','qa4:64'},[],'Grid Engine queue to use. If this is a cell array of queue names (of the form name:xx, where xx is the number of processors required to be available), a free queue will be identified automatically (and the respective number of available processors will be used as numprocs). flt_ica will wait and print a notification if no queue is available.','cat','Compute Resources'), ...
        arg({'use_pe','ParallelEnvironment'},'autodetect',{'autodetect','mpich','orte','make'},'Parallel environment to use. Only recognized by the ''stable11'' and ''stable12'' versions.','cat','Compute Resources'), ...
        arg({'numprocs','NumProcessors','num_procs'},32, uint32([1 128]),'Number or processors to use. If use_queue was used with auto-detection, the # of slots in the queue will be used as numprocs.','cat','Compute Resources'), ...
        arg({'poll_interval','PollInterval'},10,[0 Inf],'Check Amica status every n secs. Used to monitor Amica''s progress, if running on a cluster.','cat','Compute Resources'), ...
        arg({'max_start_waiting','MaxStartWaiting'},600,[0 Inf],'Maximum wait time until Amica reschedule. If running on a cluster, maximum waiting time (if job not launched) until it is assumed that the AMICA job could not be scheduled, in seconds.','cat','Compute Resources'), ...
        arg({'max_init_waiting','MaxInitWaiting'},2000,[0 Inf],'Maximum wait time until Amica reschedule. If running on a cluster, maximum waiting time (if computation not initialized) until it is assumed that the AMICA job does not make progress (e.g., due to a straggler), in seconds.','cat','Compute Resources'), ...
        arg({'reduce_factor','ReductionFactor'},0.75,[0.1 1],'Processor reduction per restart. If running on a cluster, factor by which the number of processors is reduced when the maximum waiting time is exceeded.','cat','Compute Resources'), ...
        arg({'scheduler','Scheduler'},'computing',[],'Name of grid scheduler/master. If running on a cluster, name of the Sun Grid Engine master (the ICA runner may need to ssh into this).','cat','Compute Resources'), ...
        arg({'max_restarts','MaxRestarts'},20,uint32([0 1000]),'Maximum number of restarts. If amica has been restarted more than this many times (due to optimistic scheduling or NaN outputs), fall back to the MATLAB implementation.','cat','Compute Resources'), ...
        arg({'fallback_reduce','FallbackReductionFactor'},0.5,[0.1 1],'Iteration reduction after fallback. If falling back, reduce maximum number of iterations by this factor (for time reasons).','cat','Compute Resources'), ...
        arg({'measure_window','ThroughputMeasurementWindow'},5*60,[0 Inf],'Throughput Measurement Window. Amica computational throughput is measured within this past window (in seconds) to determine whether a re-schedule can significantly improve the expected running time.','cat','Compute Resources'), ...
        arg({'benefit_threshold','SpeedupRestartThreshold'},3,[1 100],'Minimum speedup to justify a restart. If the expected speedup due to a restart is at least this high, restart the computation.','cat','Compute Resources'), ...
        arg({'native_threshold','NativeFallbackThreshold'},10,[1 100],'Minimum speedup to justify a fallback. If the expected speedup due to a restart is at least this high, but we have already re-started MaxRestarts times, we fall back to a local MATLAB-based computation.','cat','Compute Resources'), ...
        arg({'flaky_cluster','HaveFlakyCluster'},true,[],'Cluster is flaky. If this is true, computation never falls back to local computation and just retries to schedule or jobs until they eventually get through.','cat','Compute Resources'), ...
        arg({'reuse_margin','DirectoryReuseSafetyMargin'},5,[0 Inf],'Safety margin for directory reuse. When selecting a directory to write results to, make sure that it has not been updated with (incomplete) results within the last n minutes.','cat','Compute Resources'), ...
        arg({'min_numprocs_fraction','MinimumNumProcsFraction'},0.33,[0 1],'Minimum fraction of NumProcessors. If amica gives NaN''s, this is sometimes due to the number of processors chosen - so we use a random # of processors between this fraction and your original number if it crashes, until it eventually works.','cat','Compute Resources'), ...
        ...
        arg({'lrate','LRate'},0.1,[0 1],'Initial learning rate for natural gradient. Extended Infomax only.','cat','Computation'), ...
        arg({'lratefact','LRateFactor'},0.5,[0 1],'Learning rate reduction. Multiplicative factor by which to decrease learning rate (extended Infomax only).','cat','Computation'), ...
        arg({'minlrate','MinLRate'},1e-8,[0 1],'Minimum learning rate. If lower, convergence is assumed; extended Infomax only.','cat','Computation'), ...
        arg({'rholrate','ShapeLRate'},0.05,[0 1],'Initial learning rate for shape parameters. Only for generalized Gaussian mode.','cat','Computation'), ...
        arg({'rholratefact','ShapeLRateFactor'},0.5,[0 1],'Shape learning rate reduction. Multiplicative factor by which to decrease the shape learning rate. Generalized Gaussian only.','cat','Computation'), ...
        arg({'rho0','InitialShape'},1.5,[0 Inf],'Initial shape parameter value. Generalized Gaussian only.','cat','Computation'), ...
        arg({'minrho','MinShape'},1.0,[0 Inf],'Minimum shape parameter value. Generalized Gaussian only.','cat','Computation'), ...
        arg({'maxrho','MaxShape'},2.0,[0 Inf],'Maximum shape parameter value. Generalized Gaussian only.','cat','Computation'), ...
        arg({'do_newton','DoNewton'},true,[],'Use the Newton method.','cat','Computation'), ...
        arg({'newt_start','NewtonStartIter'},50,uint32([1 10000]),'Iterations without Newton method. Only after this many iterations, the Newton method is used.','cat','Computation'), ...
        arg({'newtrate','NewtonRate'},1.0,[0 1],'Learning rate for Newton iterations.','cat','Computation'), ...
        arg({'newt_ramp','NewtonRampUp'},10,uint32([1 10000]),'Iterations to ramp up Newton method.','cat','Computation'), ...
        arg({'comp_thresh','ComponentThreshold'},0.98,[0 1],'Correlation Threshold to share components.','cat','Computation'), ...
        arg({'share_start','ShareStartIter'},100,uint32([1 10000]),'Iteration to start component sharing.','cat','Computation'), ...
        arg({'share_int','ShareInterval'},100,uint32([1 10000]),'Iterations between component sharing checks.','cat','Computation'), ...
        ...
        arg({'outdir','OutputDirectory'},[],[],'Output directory. Results of the computation are written there (empty: to BCILAB''s temp directory).','cat','File IO'), ...
        arg({'load_final','LoadFinalModels'},true,[],'Load final models. If the output directory contains a final model, try to load it.','cat','File IO'), ...
        arg({'writestep','WriteInterval'},10,uint32([1 10000]),'Iterations between output writes.','cat','File IO'), ...
        arg({'write_nd','WriteHistory'},true,[],'Write component update magnitudes. Per iteration.','cat','File IO'), ...
        arg({'write_LLt','WriteLoglikes'},true,[],'Write model log-likelihoods. Per time point.','cat','File IO'), ...
        arg({'indir','InputDirectory'},'',[],'Directory of initial model, if any. Optional input directory from which to load the initial model (during AMICA runtime, models are stored in directories).','cat','File IO','type','char','shape','row'), ...
        arg({'load_param','ReadParams'},false,[],'Read parameters from input directory.','cat','File IO'), ...
        arg({'load_rej','ReadRejections'},false,[],'Read rejection info from input directory. This is for the Amica-internal rejection.','cat','File IO'), ...
        arg({'load_comp_list','ReadCompList'},false,[],'Read component assigments.','cat','File IO'), ...
        ...
        arg({'do_reject','DoReject'},false,[],'Do online time-point rejection.','cat','Artifact Handling'), ...
        arg({'numrej','RejectCycles'},3,uint32([0 100]),'Number of rejections to perform in total.','cat','Artifact Handling'), ...
        arg({'rejsig','RejectSigma'},3,[0 Inf],'Likelihood threshold for rejection (stddev). This is the number of standard deviations of likelihood below which to reject data.','cat','Artifact Handling'), ...
        arg({'rejstart','RejectStart'},2,uint32([0 10000]),'Iteration at which to start rejection.','cat','Artifact Handling'), ...
        arg({'rejint','RejectInterval'},3,uint32([0 10000]),'Iterations between successive rejections.','cat','Artifact Handling'), ...
        arg({'kurt_start','KurtosisStart'},3,uint32([0 10000]),'Iterations without kurtosis calculations. This is for Extended Infomax only.','cat','Artifact Handling'), ...
        arg({'num_kurt','KurtosisCycles'},5,uint32([0 10000]),'Number of kurtosis calcs to perform. This is for Extended Infomax only.','cat','Artifact Handling'), ...
        arg({'kurt_int','KurtosisInterval'},1,uint32([0 10000]),'Iterations between successive kurtosis calculations. This is for Extended Infomax only.','cat','Artifact Handling'), ...
        ...
        arg({'update_A','UpdateWeights'},true,[],'Update the mixing matrix.','cat','Parameter Updates'), ...
        arg({'update_c','UpdateCenters'},true,[],'Update the model centers.','cat','Parameter Updates'), ...
        arg({'update_gamma','UpdateProbabilities'},true,[],'Update the model probabilities.','cat','Parameter Updates'), ...
        arg({'update_alpha','UpdateProportions'},true,[],'Update the source mixture proportions.','cat','Parameter Updates'), ...
        arg({'update_mu','UpdateLocations'},true,[],'Update the source mixture locations.','cat','Parameter Updates'), ...
        arg({'update_sbeta','UpdateScales'},true,[],'Update the source mixture scales.','cat','Parameter Updates'), ...
        arg({'do_rho','UpdateShapes'},true,[],'Update the shape parameters.','cat','Parameter Updates'), ...
        ...
        arg({'decwindow','DetectWindow'},1,[],'Moving average window for likelihood check. Likelihood decrease is detected in this window (presumably over iterations - Jason?).','cat','Miscellaneous'), ...
        arg({'invsigmax','MaxInverseSigma'},100,[],'Maximum value of inverse scale parameters.','cat','Miscellaneous'), ...
        arg({'invsigmin','MinInverseSigma'},1e-8,[],'Minimum value of inverse scale parameters.','cat','Miscellaneous'), ...
        arg({'do_mean','Centering'},true,[],'Remove the mean from the data.','cat','Miscellaneous'), ...
        arg({'do_sphere','Sphering'},true,[],'Sphere (whiten) the data.','cat','Miscellaneous'), ...
        arg({'doPCA','PCAReduction'},true,[],'Do a PCA dimensionality reduction.','cat','Miscellaneous'), ...
        arg({'pcakeep','RetainPCAs'},[],[],'Number of PCA components to keep. If PCA is enabled (if empty: #channels).','cat','Miscellaneous','shape','scalar','type','uint32'), ...
        arg({'doscaling','Rescale'},true,[],'Rescale unmixing matrix to unit norm.','cat','Miscellaneous'), ...
        arg({'scalestep','ScaleInterval'},1,uint32([1 10000]),'Iterations between unmixing matrix rescaling.','cat','Miscellaneous'), ...
        arg({'block_size','BlockSize'},128,uint32([10 10000]),'Matrix block size for block multiplication.','cat','Miscellaneous'), ...
        arg({'verbose','Verbose'},true,[],'Show progress updates.','cat','Miscellaneous')}, ...
    'infomax',{ ...
        arg({'maxsteps','MaxIterations'},512,uint32([1 10000]),'Maximum number of ICA training steps.','cat','Core Parameters'), ...
        arg({'extended','ExtendedInfomax'},int32(3),[],'Perform Extended-ICA sign estimations for N training blocks. If N > 0, automatically estimate the of sub-Gaussian sources. If N < 0, fix number of sub-Gaussian comps to -N [faster than N>0].','cat','Core Parameters'), ...
        arg({'pca','RetainPCAs'},0,uint32([0 10000]),'Do a PCA dimensionality reduction. Value is the number of PCs to retain (0=off)','cat','Core Parameters'), ...
        arg({'sphering','Sphering'},true,[],'Whether to sphere the data.','cat','Core Parameters'), ...
        ...
        arg({'lrate','LearningRate'},[],[],'Initial ICA learning rate (<< 1). If empty, use a heuristic.','cat','Computation'), ...
        arg({'anneal','Annealing'},[],[0 1],'Annealing constant, controls speed of convergence. If empty, use 0.90 for regular or 0.98 for extended Infomax.','cat','Computation'), ...
        arg({'annealdeg','AnnealChange'},60,[0 90],'Weight change threshold for annealing steps. In degrees.','cat','Computation'), ...
        arg({'pstop','MinWeightChange','stop'},0.000001,[0 Inf],'Minimum weight change to signal convergence.','cat','Computation'), ...
        ...
        arg({'pweights','InitialWeights','weights'},[],[],'Initial weight matrix, if any. If empty, use eye()/spher().','cat','Miscellaneous'), ...
        arg({'pblock','BlockSize','block'},[],[],'ICA block size (<< datalength). If empty, use a heuristic.','cat','Miscellaneous','shape','scalar','type','int32'), ...
        arg({'bias','BiasAdjust'},true,[],'Perform bias adjustment.','cat','Miscellaneous'), ...
        arg({'momentum','Momentum'},0,[0 1],'Training momentum. Technique to prevent getting stuck in local minima.','cat','Miscellaneous'), ...
        arg({'pposact','PosActivations','posact'},false,[],'Make all component activations net-positive. Requires time and memory; posact() may be applied separately.','cat','Miscellaneous'), ...
        arg({'verbose','Verbose'},true,[],'Show progress updates.','cat','Miscellaneous'), ...
        arg({'logfile','Logfile'},'',[],'Log file to save messages to. Messages will still be shown on screen.','cat','Miscellaneous','shape','row','type','char'), ...
        arg({'interupt','InterruptBtn'},false,[],'Show interupt button (slow).','cat','Miscellaneous')}, ...
    'beamica', { ...
        arg({'max_iter','MaxIterations'},1200,uint32([1 10000]),'Maximum number of iterations.'), ...
        arg({'anchorlabels','AnchorLabels'},{}, {'Precentral_L', 'Precentral_R', 'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', 'Rolandic_Oper_L', 'Rolandic_Oper_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R', 'Insula_L', 'Insula_R', 'Cingulum_Ant_L', 'Cingulum_Ant_R', 'Cingulum_Mid_L', 'Cingulum_Mid_R', 'Cingulum_Post_L', 'Cingulum_Post_R', 'Hippocampus_L', 'Hippocampus_R', 'ParaHippocampal_L', 'ParaHippocampal_R', 'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R', 'Lingual_L', 'Lingual_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'Occipital_Mid_L', 'Occipital_Mid_R', 'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Postcentral_L', 'Postcentral_R', 'Parietal_Sup_L', 'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Angular_L', 'Angular_R', 'Precuneus_L', 'Precuneus_R', 'Paracentral_Lobule_L', 'Paracentral_Lobule_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', 'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Inf_L', 'Temporal_Inf_R', 'Olfactory_L', 'Olfactory_R', 'Rectus_L', 'Rectus_R', 'Amygdala_L', 'Amygdala_R', 'Caudate_L', 'Caudate_R', 'Thalamus_L', 'Thalamus_R', 'Heschl_L', 'Heschl_R'}, ...
            'Cortical anchor locations. List of locations to which components shall be constrained. The first k components are encouraged to lie close to the given locations, in the order of appearance. This is experimental and currently requires a) 10-20 locations and b) Guido Nolte''s source analysis toolbox (not included).','experimental',true), ...
        arg({'reference','Reference'},'nasion',{'nasion','common_average'},'Referencing scheme. Only needed if anchor locations are selected.','experimental',true), ...
        arg({'tradeoff','BeamPenaltyStrength'},0,[0 1],'Beamformer penalty strength. Larger values emphasize the beamformer constraint over the Infomax cost function.','experimental',true), ...
        arg({'lrate','LearningRate'},0.5,[0 1],'Learning rate. The maximum is 1.0, but lower rates (down to 0.1 or below) can be used to ensure stability.'), ...
        arg({'verbose','VerboseOutput'},true,[],'Show verbose output.'), ...
        arg({'usegpu','TryUseGPU','UseGPU'},true,[],'Try to use the GPU.'), ...
        arg({'convergence_check','ConvergenceCheck'},false,[],'Force convergence check if on GPU. This is slow on the GPU, but can help ensuring that the solution converges.')}, ...
    'fastica', { ...
        arg({'maxNumIterations','MaxIterations'},1000,uint32([1 10000]),'Maximum number of iterations.','cat','Core Parameters'), ...
        arg({'approach','Approach'},'symm',{'symm','defl'},'The decorrelation approach used. Can be symmetric, i.e. estimate all the independent component in parallel, or deflation, i.e. estimate independent component one-by-one like in projection pursuit.','cat','Core Parameters'), ...
        arg({'numOfIC','NumICs'},[],[],' Number of estimated independent components. Default equals the dimension of data.','cat','Core Parameters','shape','scalar','type','uint32'), ...
        arg({'g','Nonlinearity'},'tanh',{'pow3','tanh','gauss','skew'},'Nonlinearity to use. Pow3 is g(u)=u^3, tanh is g(u)=tanh(a1*u), gauss is g(u)=u*exp(-a2*u^2/2), and skew is g(u)=u^2.','cat','Core Parameters'), ...
        arg({'finetune','FineTuning'},'tanh',{'off','pow3','tanh','gauss','skew'},'Nonlinearity for fine-tuning. Pow3 is g(u)=u^3, tanh is g(u)=tanh(a1*u), gauss is g(u)=u*exp(-a2*u^2/2), and skew is g(u)=u^2.','cat','Core Parameters'), ...
        arg({'a1','TanhParameter'},1,[],'Parameter a1 for tanh nonlinearity.','cat','Core Parameters'), ...
        arg({'a2','GaussParameter'},1,[],'Parameter a2 for gauss nonlinearity.','cat','Core Parameters'), ...
        ...
        arg({'stepsize','StepSize','mu'},1,[],'Step size. If the value is other than 1, then the program will use the stabilized version of the algorithm (see also parameter "Stabilization").','cat','Convergence'), ...
        arg({'stabilization','Stabilization'},false,[],'Use an adaptive step size. Serves to stabilize the convergence.','cat','Convergence'), ...
        arg({'epsilon','Epsilon'},0.001,[],'Stopping criterion.','cat','Convergence'), ...
        arg({'maxFinetune','MaxFinetune'},100,uint32([1,10000]),'Maximum number of iterations in fine-tuning.','cat','Convergence'), ...
        arg({'sampleSize','SampleSize'},1,[0 1],'Percentage of samples used per iteration. Samples are chosen randomly.','cat','Convergence'), ...
        arg({'initGuess','InitialGuess'},[],[],'Initial guess for A, if any. Default is random.','cat','Convergence'), ...
        ...
        arg({'verbose','Verbose'},true,[],'Report progress of algorithm.','cat','Miscellaneous'), ...
        arg({'displayMode','DisplayMode'},'off',{'signals','basis','filters','off'},'Plot running estimates of independent components.','cat','Miscellaneous'), ...
        arg({'displayInterval','DisplayInterval'},1,[],'Number of iterations between plots.','cat','Miscellaneous'), ...
        arg({'firstEig','FirstEigenvalue'},1,[],'Skip n largest components. This and "LastEigvalue" specify the range for eigenvalues that are retained, "FirstEigenvalue" is the index of largest eigenvalue to be retained.','cat','Miscellaneous'), ...
        arg({'lastEig','LastEigenvalue'},[],[],'Skip n smallest components. This is the index of the last (smallest) eigenvalue to be retained. Default equals the dimension of data.','cat','Miscellaneous')}, ...    
     'rica', { ...
        arg({'numFeatures','NumComponents'},[],[],'Number of components to learn. Can be larger than the number of channels. If empty, use # channels.'), ...
        arg({'dict_criterion','DictionaryCriterion'},'reconstruction',{'reconstruction','cortically_anchored'},'Dictionary learning criterion. The cortically anchored mode allows to specify a cell array of anchor locations.'), ...
        arg({'lambda','Lambda'},0.05,[0 Inf],'Sparsity/accuracy tradeoff. The lambda parameter in sparse dictionary learning.'), ...
        arg({'gamma','Gamma'},0.01,[0 Inf],'Anchoring tradeoff. This is the tradeoff parameter associated with the anchor constraint term.'), ...
        arg({'theta','Theta'},1,[0 Inf],'Reconstruction tradeoff. This is the tradeoff parameter associated with the reconstruction cost.'), ...
        arg({'anchors','AnchorLocations'},[],[],'Cortical anchor locations. Cell array of dipolar projection triplets, one cell per location that shall be anchored. The first k components are cortically constrained. The associated tuning parameter is Gamma.'), ...
        arg({'anchor_init','AnchorInitialization'},'from_subspace',{'from_subspace','bestmatch','perpendicular'},'Anchor component initialization. When using from_subspace a random linear combination from the subspace will be used. When using bestmatch, the best-matching component in the initialization will be used.'), ...
        arg({'max_restarts','MaxRestarts'},20,uint32([0,1000]),'Maximum # of restarts. When weights blow up or etc.'), ...
        arg({'cov_blocksize','CovarianceBlocksize'},10,uint32([1,100000]),'Robust blocksize. Blocksize for robust estimation (currently only covariance).'), ...
        arg({'initialization','Initialization'},'random_dipoles',{'random_dipoles','random','radial'},'Initialization. Either randomly or using radial (sphering) components.'), ...
        arg({'cov_rejection','CovarianceRejection'},5,[0 100],'Covariance-based rejection. Rejects samples that are beyond this many std-devs from a robustly estimated data distribution.'), ...
        arg({'temporal_normalization','TemporalNormalization'},false,[],'Temporal normalization. Whether the data should also be normalized in time.'), ...
        arg({'random_init_scale','RandomInitScale'},0.1,[],'Initial randomness level. This is the fraction of random noise used to determine the initial solution (e.g., eye+n*randn().'), ...
        arg({'randseed','RandomSeed'},10,[],'Random seed. Use the same seed for reproducible results.'), ...
        arg({'epsilon','Epsilon'},1e-5,[],'Epsilon parameter in sparsity cost.'), ...
        arg({'topoplot','ShowTopoplot'},false,[],'Show topoplots per update. For tracking of convergence, etc.'), ...
        arg_nogui({'chan_labels','ChannelLabels'}), ...
        arg_sub({'solverOptions','SolverOptions'},{}, { ...
            arg({'MaxIter','MaxIterations'},75,uint32([1 10000]),'Maximum number of iterations.'), ...
            arg({'Method','Optimizer'},'Scaled Non-Linear Conjugate Gradient',{'Steepest Descent','Cyclic Steepest Descent','Barzilai and Borwein Gradient','Non-Linear Conjugate Gradient','Scaled Non-Linear Conjugate Gradient','Preconditionined Non-Linear Conjugate Gradient','Quasi-Newton with Limited-Memory BFGS Updating','Hessian-Free Newton','Preconditioned Hessian-Free Newton ','Quasi-Newton Hessian approximation','Newton with Hessian update every k''th step'},'Optimization method to use. Later methods use fewer step sizes but espectially the Newton-type methods are more expensive, and some advanced methods might be too brittle for ICA.'), ...
            arg({'MaxFunEvals','MaxFunctionEvals'},1000,uint32([1 100000]),'Maximum function evaluations. Maximum number of function evaluations allowed (including during line searches).'), ...
            arg({'Display','VerbosityLevel'},'iter',{'off','final','iter','full','excessive'},'Verbosity level.'), ...
            arg({'optTol','OptimalityTolerance'},[],[],'First-order tolerance. Termination tolerance on the first-order optimality.'), ...
            arg({'progTol','ProgressTolerance'},[],[],'Progress tolerance. Termination tolerance on progress in terms of function/parameter changes.'), ...
            arg({'c1','ArmijoParam'},[],[],'Armijo parameter. Sufficient Decrease for Armijo condition.'), ...
            arg({'c2','WolfeParam'},[],[],'Wolfe parameter. Curvature Decrease for Wolfe conditions. If [], this defaults to 0.2 for CG methods and 0.9 otherwise.'), ...
            arg({'LS_init','LineSearchInit'},'Default',{'Default','Always1','AsPrevious','QuadraticInterp','TwicePrevious','ScaledConj'},'Line search initialization. The options are: Always try an initial step length of 1 (default for all except ''sd'' and ''cg''); Use a step similar to the previous step; Quadratic Initialization using previous function value and new; The minimum between 1 and twice the previous step length; The scaled conjugate gradient step length (may accelerate conjugate gradient methods, but requires a Hessian-vector product, default for ''scg'').'), ...
            arg({'LS_type','LineSearchType'},'Default',{'Default','BacktrackingArmijo','BracketingWolfe','MatlabOptTB'},'Line search type. The options are: A backtracking line-search based on the Armijo condition (default for ''bb''); A bracekting line-search based on the strong Wolfe conditions (default for all other methods); The line-search from the Matlab Optimization Toolbox (requires Matlab''s linesearch.m to be added to the path).'), ...
            arg({'LS_interp','LineSearchInterpolation'},'Default',{'Default','DoubleOrBisect','Cubic','MixedQuadraticCubic'},'Line search interpolation. For the Wolfe condition. The options are: Step Size Doubling and Bisection; Cubic interpolation/extrapolation using new function and gradient values; Mixed quadratic/cubic interpolation/extrapolation. MinFunc documents that Cubic is the default, while in practice it uses the Mixed mode.'), ...
            arg({'LS_multi','LineSearchExtraPoints'},'Default',{'Default','SameOrder','HigherOrder'},'Line search extra-points rule. How to handle extra points during the line search. The options are: Keep the same polynomial order regardless, or use a higher order (cubic,quadratic,quintic) if sufficient points are present.'), ...
            arg({'useMex','UseMex'},true,[],'Use mex functions. Where applicable, use compiled mex files to speed things up (may not be available on every system).')},'Control options for the optimizer (minFunc).')}, ...
     'dictica', @dictlearn, ...
     'kernelica', { ...
        arg({'contrastfun','ContrastFunc','contrast'},'kcca',{'kcca','kgv'},'Contrast function. Either Kernel Canonical Correlation Analysis or Kernel Generalized Variance.','cat','Core Parameters'), ...
        arg({'contrasttype','ContrastType'},'full',{'full','oneunit'},'Type of the contrast function.','cat','Core Parameters'), ...
        arg({'polish','Finetune'},true,[],'Double the precision. Finish with a half sigma value (gives better estimates).','cat','Core Parameters'), ...
        arg({'restarts','NumRestarts'},1,uint32([0 1000]),'Number of restarts.','cat','Core Parameters'), ...
        arg({'kernel','KernelFunc'},'gaussian',{'gaussian','poly','hermite'},'Type of kernel for contrast function.','cat','Core Parameters'), ...
        ...
        arg({'GaussianSigma'},1,[0 Inf],'Bandwidth parameter for the Gaussian kernel.','cat','Miscellaneous'), ...
        arg({'r','PolyR'},1,[],'Parameter r for the polynomial kernel.','cat','Miscellaneous'), ...
        arg({'s','PolyS'},1,[],'Parameter s for the polynomial kernel.','cat','Miscellaneous'), ...
        arg({'d','PolyDegree'},3,uint32([1 10]),'Degree of the polynomial kernel.','cat','Miscellaneous'), ...
        arg({'p','HermiteP'},3,[],'Parameter p of the Hermite Kernel.','cat','Miscellaneous'), ...
        arg({'HermiteSigma'},2,[],'Parameter sigma of the Hermite Kernel.','cat','Miscellaneous'), ...
        arg({'kap','Kappa'},0.01,[],'Regularization Parameter.','cat','Miscellaneous'), ...
        arg({'dodisp','Verbose','disp'},true,[],'Regularization Parameter.','cat','Miscellaneous')}, ...
    'fastkernelica', { ...
        arg({'maxiter','MaxIterations'},20,uint32([1,10000]),'Maximum number of iterations.','cat','Core Parameters'), ...
        arg({'psigma','KernelSize','sigma'},0.5,[],'Gaussian kernel size. One is a reasonable default; smaller values (e.g., 0.5) allow for higher precision, but can run into local minima.','cat','Core Parameters'), ...
        arg({'thresh','ConvergenceThreshold'},1e-7,[],'Convergence threshold. The algorithm terminates when the difference in subsequent values of the kernel independence measure is lower than this.','cat','Core Parameters'), ...
        arg({'restarts','NumRestarts'},10,uint32([0 1000]),'Number of restarts. Multiple restarts are necessary to find the global optimum, especially with small kernel sizes.','cat','Core Parameters')} ...
    'sphere' {}, ...
    'robust_sphere', {}, ...
    }, 'ICA variant. AMICA is the highest quality (but slowest, except if run on a cluster), Infomax is second-highest quality, FastICA is fastest (but can fail to converge and gives poorer results), KernelICA is experimental.'), ...
    arg_sub({'data_cleaning','DataCleaning','CleaningLevel','clean'},{}, @flt_clean_settings,'Optional data cleaning prior to running an ICA. The computed ICA solution will be applied to the original uncleaned data.'), ...
    arg({'do_transform','TransformData','transform'},false,[],'Transform the data rather than annotate. By default, ICA decompositions are added as annotations to the data set.'),...
    arg({'retain_labels','RetainLabels'},true,[],'Retain labels when transforming. If this is false the channel labels will be replaced by 1:k for k components, if TransformData is checked.'),...
    arg({'clear_after_trans','ClearAfterTransform'},true,[],'Clear .icaweights after transform. This is so that later functions do not attempt to transform the already transformed data.'),...
    arg({'do_calcact','CalculateActivation'},false,[],'Calculate component activations. If true, the .icaact field will be populated.'),...
    arg({'cleaned_data','OutputCleanedData'},false,[],'Emit cleaned data. Whether the cleaned data, instead of the original data should be output (note: this is not applicable for online use, since most cleaning filters cannot be run online).'),...
    arg({'doresume','ResumePrevious','resume'},true,[],'Try to resume previous computations if possible.'), ...
    arg({'doforce','ForceComputation','force'},false,[],'Force computation. Recompute ICA even if the input data set already has an attached ICA solution.'), ...
    arg({'dodebug','DebugMode','debug'},false,[],'Debug mode. Stops and waits for user input in case of an exception.','guru',true), ...
    arg({'normalize_weights','NormalizeWeights','NormalizeWights'},false,[],'Normalize weights.'), ...
    arg_nogui({'state','State'}));

if ~isempty(state)
    % online case: annotate the data
    signal.icasphere = state.icasphere;
    signal.icaweights = state.icaweights;
    signal.icachansind = state.icachansind;
    signal.icawinv = state.icawinv;
    if isfield(state,'amica')
        signal.etc.amica = state.amica; end
else
    logfile = env_translatepath('home:/.bcilab/logs/ica_datalog.log');

    % offline case: check if we actually need to compute an ICA solution (or if the signal already
    % has one)
    if doforce || ~(isfield(signal,'icaweights') && ~isempty(signal.icaweights))
        
        % first pre-process the data (using a sequence of non-causal data reductions)
        pre = signal;
        [chns,pnts,trials] = size(pre.data);
        
        if trials ~= 1
            % epoched dataset... reshape it
            pre.data = reshape(pre.data,chns,[],1);
            [pre.chns,pre.pnts,pre.trials] = size(pre.data);
            % delete event & epoch information (would not survive artifact rejection, anyway)
            pre.epochs = [];
            pre.event = [];
        end
        
        % clean the data
        pre = exp_eval(flt_clean_settings(data_cleaning,'signal',pre));
                
        % get the underlying chanlocs
        root_chanlocs = pre.chanlocs;
        global debug_chanlocs; debug_chanlocs = pre.chanlocs;
        
        % run ICA on preprocessed data
        switch variant.arg_selection %#ok<*NODEF>
            case 'noica'
                % this is just for testing - creates random weights
                pre.icaweights = randn(length(pre.chanlocs));
                pre.icasphere = eye(size(pre.data,1));
                
            case 'amica'
                % determine a unique identifier for this computation
                variant_core = hlp_struct2varargin(variant,'restrict',{ ...
                    'amica_version','max_iter','num_models','num_mix_comps','pdftype','share_comps','lrate','lratefact','minlrate','rholrate','rholratefact','rho0','minrho', ...
                    'maxrho','do_newton','newt_start','newtrate','newt_ramp','comp_thresh','share_start','share_int','do_reject','numrej','rejsig','rejstart','rejint','kurt_start', ...
                    'num_kurt','kurt_int','update_A','update_c','update_gamma','update_alpha','update_mu','update_sbeta','do_rho','decwindow','invsigmax','invsigmin','do_mean', ...
                    'do_sphere','doPCA','pcakeep','doscaling','scalestep','block_size'});
                tag = hlp_fingerprint({pre.tracking.expression,variant_core});
                
                % set the output directory
                if isempty(variant.outdir)
                    variant.outdir = env_translatepath(['temp:/amicaout-' num2str(tag)]);
                else
                    % sanitize output directory
                    variant.outdir = env_translatepath(variant.outdir);
                    if variant.outdir(end) ~= filesep
                        variant.outdir = [variant.outdir filesep]; end
                end
                if ~isempty(variant.indir)
                    % sanitize input directory
                    variant.indir = env_translatepath(variant.indir);
                    if variant.indir(end) ~= filesep
                        variant.indir = [variant.indir filesep]; end
                end
                
                % check if a final result for this problem is already there from a previous run
                if variant.load_final && exist([variant.outdir filesep 'out.txt'],'file')
                    try
                        disp('Checking for previous solutions for this computation...');
                        fileinfo = dir([variant.outdir filesep 'out.txt']);
                        time_difference = (now - fileinfo.datenum)*24*60;
                        if time_difference < variant.reuse_margin
                            % if the result file is too fresh, we might run into a
                            % conflict trying to write to a file that is still in use...
                            conflict_potential = true; end %#ok<NASGU>
                        
                        % cycle through all alternative storage locations for this tag
                        dirs = dir([variant.outdir '*']);
                        k = 1;
                        while ~exist('r','var') && k <= length(dirs)
                            curdir = env_translatepath(['temp:/' dirs(k).name]);
                            try
                                % an out.txt was written, check if it stems from a completed run
                                t = fopen([curdir filesep 'out.txt']);
                                while ~exist('r','var')
                                    line = fgetl(t);
                                    if ~ischar(line)
                                        break; end
                                    if ~isempty(strfind(line,'done.'))
                                        % found the 'done.' line - sanity-check the model
                                        switch variant.amica_version
                                            case 'stable12'
                                                tmp = loadmodout12(curdir);
                                            case 'stable11'
                                                tmp = loadmodout11(curdir);
                                            case 'devel'
                                                tmp = loadmodout10(curdir);
                                            case 'stable'
                                                tmp = loadmodout(curdir);
                                        end
                                        if isfield(tmp,'W') && size(tmp.W,3) == tmp.num_models
                                            disp('Found a previously successful Amica solution for the same problem.');
                                            r = tmp;
                                        end
                                    end
                                end
                                fclose(t);
                            catch
                                disp('Found an out.txt file of a previous Amica solution, but failed to load the model.');
                                try
                                    fclose(t);
                                catch,end
                            end
                            k = k+1;
                        end
                    catch
                        disp('Found directories with alternative solutions, but failed trying to read them.');
                    end
                end
                
                % if no previous solution is not already there, go ahead!
                if ~exist('r','var')
                    
                    if exist('conflict_potential','var')
                        % if there is conflict potential, we preferably choose a different output directory
                        fprintf('Target directory has results of a computation that was updated less than %i minutes ago; choosing a different output directory.\n',ceil(time_difference));
                        variant.outdir = [variant.outdir '-alt' num2str(mod(tic,100000))];
                    end
                    % translate pdftype into the form expected by amica
                    variant.pdftype = hlp_rewrite(variant.pdftype,'GeneralizedGaussian',0,'ExtendedInfomax',1,'Gaussian',2,'Logistic',3);
                    % translate the PE auto-detection
                    variant.use_pe = hlp_rewrite(variant.use_pe,'autodetect','');
                    % translate all booleans into doubles
                    for fn=fieldnames(variant)'
                        if islogical(variant.(fn{1}))
                            variant.(fn{1}) = double(variant.(fn{1})); end
                    end
                    % check pcakeep
                    if isempty(variant.pcakeep)
                        variant.pcakeep = size(pre.data,1); end
                    
                    % keep track of this to generate random system configurations if amica fails with NaNs...
                    original_numprocs = variant.numprocs;
                    
                    % adapt argument list for the different supported amica versions...
                    switch variant.amica_version
                        case {'stable11','stable12'}
                            suppress_args = {'scheduler','do_mean','do_sphere','doPCA','scalestep','kurt_start','num_kurt','kurt_int','load_comp_list','block_size'};
                            rewrite_args = {};
                            if isempty(variant.max_threads)
                                variant.max_threads = uint32(999); end
                        case 'stable'
                            suppress_args = {'share_comps','share_int','share_start','comp_thresh','doPCA','do_mean','do_sphere','kurt_int','num_kurt','kurt_start','load_comp_list','load_rej','scalestep','block_size','use_queue','use_pe'};
                            rewrite_args = {'update_A','update_W'};
                            if isempty(variant.max_threads)
                                variant.max_threads = uint32(4); end
                        case 'devel'
                            suppress_args = {'do_mean','do_sphere','doPCA','scalestep','kurt_start','num_kurt','kurt_int','load_comp_list','use_queue','use_pe'};
                            rewrite_args = {};
                            if isempty(variant.max_threads)
                                variant.max_threads = uint32(4); end
                        case 'matlab'
                            suppress_args = {'scheduler'};
                            rewrite_args = {};
                        otherwise
                            error('Unsupported Amica version.')
                    end
                    % ...and apply a final set of argument transformations that hold across all version
                    arglist = hlp_struct2varargin(variant,'suppress',[{'arg_selection','arg_direct','poll_interval','max_start_waiting','max_init_waiting','reduce_factor','fallback_reduce','verbose','max_restarts','amica_version','measure_window','benefit_threshold','native_threshold','load_final','flaky_cluster','reuse_margin','min_numprocs_fraction'},suppress_args],'rewrite',[{'useqsub','qsub'}, rewrite_args]);
                    
                    if strcmp(variant.amica_version,'matlab')
                        % run native MATLAB version
                        disp('Running native MATLAB implementation of AMICA');
                        r = amica_native(pre,variant);
                    else
                        % run fast binary version
                        disp(['AMICA output is in ' variant.outdir]);
                        try
                            
                            % try to delete all files in this directory
                            if ~(doresume && variant.load_param)
                                try
                                    files = dir(variant.outdir);
                                    for f=1:length(files)
                                        if ~isdir(files{f}.name)
                                            delete(files{f}.name); end
                                    end
                                catch, end
                            end
                            
                            % initiate...
                            if strcmp(variant.useqsub,'off')
                                % run locally (without using the cluster)
                                switch variant.amica_version
                                    case 'stable12'
                                        r = runamica12(pre.data,[],arglist{:});
                                    case 'stable11'
                                        r = runamica11(pre.data,[],arglist{:});
                                    case 'devel'
                                        r = runamica10(pre.data,[],size(pre.data,1),size(pre.data,2), arglist{:});
                                    case 'stable'
                                        r = runamica(pre.data,[],size(pre.data,1),size(pre.data,2), arglist{:});
                                end
                            else
                                % find out the median Amica running time from past history
                                try
                                    standard_normtime = median(load(env_translatepath('resources:/amica_runtimes.txt')));
                                catch
                                    disp('No past Amica performance history available.');
                                    standard_normtime = [];
                                end
                                
                                % schedule a run on the cluster
                                [job,r,t0,iter_times,maxiters,lastlines,job_started,early_finish] = schedule_amica(pre.data,variant.scheduler,variant.amica_version,logfile,arglist{:});
                                
                                % while not finished/terminated: monitor the computation...
                                num_restarts = 0;
                                while maxiters < variant.max_iter && ~early_finish
                                    pause(variant.poll_interval);
                                    
                                    % scan the output that amica has produced so far and derive a few
                                    % variables
                                    try
                                        % job is already started?
                                        if ~job_started && exist([variant.outdir filesep 'out.txt'],'file')
                                            job_started = true; end
                                        
                                        t = fopen([variant.outdir filesep 'out.txt']);
                                        curline = 0;
                                        got_nan = false;
                                        lastiters = maxiters;
                                        while 1
                                            % get next line; display
                                            line = fgetl(t);
                                            if ~ischar(line)
                                                break; end
                                            curline = curline + 1;
                                            if curline > lastlines
                                                disp(line); end
                                            
                                            % amica finished?
                                            if ~isempty(strfind(line,'done.'))
                                                if maxiters < variant.max_iter
                                                    disp('Amica finished early.'); end
                                                early_finish = true;
                                            end
                                            
                                            % we got NaN's?
                                            if ~isempty(strfind(line,'NaN'))
                                                got_nan = true; end
                                            
                                            try
                                                strs = hlp_split(line,' ');
                                                % we got an iteration output?
                                                if strcmp(strs{1},'iter')
                                                    % remember the iteration number
                                                    if 2 <= length(strs)
                                                        maxiters = str2num(strs{2}); end
                                                    try
                                                        % ... and collect samples of the iteration times
                                                        idx = find(strcmp(strs,'s,'),1)-1;
                                                        if ~isempty(idx) && idx <= length(strs)
                                                            iter_times(curline) = str2num(strs{idx}); end
                                                    catch,end
                                                end
                                            catch,end
                                        end
                                        
                                        if lastiters == 0
                                            % measure the time it took to start the computation
                                            % (for statistics & decision-making...)
                                            startup_time = toc(t0) - sum(iter_times); end
                                        
                                        lastlines = curline;
                                        fclose(t);
                                    catch
                                        try
                                            fclose(t);
                                        catch,end
                                    end
                                    
                                    % check if amica is running so slowly that we need to restart it
                                    if ~isempty(iter_times) && ~isempty(standard_normtime) && ~(any(strcmp(variant.amica_version,{'stable11','stable12'})) && ~isempty(variant.use_queue))
                                        % estimate computational complexity for the current configuration
                                        [C,S] = size(pre.data);
                                        M = variant.num_models;
                                        K = variant.num_mix_comps;
                                        P = variant.numprocs;
                                        config_complexity = ((M*(C*C*S + K*C*S + C*S*S))/P);
                                        
                                        % compute complexity-normalized runtime per iteration
                                        iter_normtime = iter_times / config_complexity;
                                        
                                        % measure the average computational throughput in the past measurement window
                                        cutoff = find(cumsum(iter_times(end:-1:1)) > variant.measure_window,1);
                                        if ~isempty(cutoff)
                                            mean_normtime = mean(iter_normtime(end-cutoff+1:end));
                                            % compute expected runtime with the current compute throughput
                                            expected_runtime = startup_time + sum(iter_times) + (mean_normtime * config_complexity * (variant.max_iter - maxiters));
                                            % compute expected runtime with the mean compute throughput
                                            standard_runtime = startup_time + (standard_normtime * config_complexity * variant.max_iter);
                                            
                                            % is it better to reschedule?
                                            if expected_runtime / standard_runtime > variant.benefit_threshold && num_restarts < variant.max_restarts
                                                disp('Computation on the cluster is very slow; re-starting to get a better node allocation');
                                                if length(job) >= 2
                                                    if isa(job{2},'onCleanup')
                                                        % clearing the job will implicitly lead to its deletion
                                                        clear job;
                                                    else isa(job{2},'function_handle')
                                                        % old MATLAB version: we invoke the function in job{2} to delete it explicitly
                                                        job{2}();
                                                    end
                                                else
                                                    disp('couldn''t delete the previous jobs; please delete them manually.');
                                                end
                                                variant.benefit_threshold = variant.benefit_threshold + 1;
                                                % restart computation...
                                                num_restarts = num_restarts + 1;
                                                [job,r,t0,iter_times,maxiters,lastlines,job_started,early_finish] = schedule_amica(pre.data,variant.scheduler,variant.amica_version,logfile,arglist{:},'numprocs',variant.numprocs);
                                                pause(startup_time);
                                            elseif expected_runtime / standard_runtime > variant.native_threshold && num_restarts >= variant.max_restarts
                                                if variant.flaky_cluster
                                                    disp('The cluster is apparently overloaded; not doing any further restarts from now on.');
                                                else
                                                    disp('Computation on the cluster is extremely slow; reverting to native MATLAB computatation on local machine.');
                                                    % fall back to native Amica implementation...
                                                    error('fall back');
                                                end
                                            end
                                        end
                                    end
                                    
                                    % check if amica is starting to produce NaN values
                                    if got_nan
                                        % sometimes AMICA gives NaN outputs; need to restart in this case
                                        disp('Got NaN outputs; trying to restart the job...');
                                        if length(job) >= 2
                                            if isa(job{2},'onCleanup')
                                                % clearing the job will implicitly lead to its deletion
                                                clear job;
                                            else isa(job{2},'function_handle')
                                                % old MATLAB version: we invoke the function in job{2} to delete it explicitly
                                                job{2}();
                                            end
                                        else
                                            disp('couldn''t delete the previous jobs; please delete them manually.');
                                        end
                                        if num_restarts > variant.max_restarts && ~variant.flaky_cluster
                                            disp('Amica consistently gives NaN outputs; falling back to MATLAB implementation.')
                                            error('fall back');
                                        else
                                            if strcmp(variant.amica_version,'stable11')
                                                disp('Choosing a random number of processors.');
                                                if ~exist('numprocs_tried','var')
                                                    numprocs_tried = []; end
                                                while true
                                                    variant.numprocs = floor(original_numprocs*variant.min_numprocs_fraction + original_numprocs*(1-variant.min_numprocs_fraction)*rand);
                                                    if ~any(variant.numprocs == numprocs_tried)
                                                        numprocs_tried(end+1) = variant.numprocs;
                                                        break;
                                                    elseif length(numprocs_tried) >= (original_numprocs*(1-variant.min_numprocs_fraction))-1
                                                        disp('Tried all possible system configurations; giving up and running locally...');
                                                        error('fall back');
                                                    end
                                                end
                                            end
                                            
                                            % re-start the job
                                            num_restarts = num_restarts + 1;
                                            [job,r,t0,iter_times,maxiters,lastlines,job_started,early_finish] = schedule_amica(pre.data,variant.scheduler,variant.amica_version,logfile,arglist{:},'numprocs',variant.numprocs);
                                            exist([variant.outdir filesep 'out.txt'],'file')
                                        end
                                    end
                                    
                                    % check if starting the job takes suspiciously long (indicating
                                    % that the job size exceeds the cluster quota)
                                    if (~job_started && toc(t0) > variant.max_start_waiting) || (maxiters == 0 && toc(t0) > variant.max_init_waiting)
                                        % in this case, we reduce that and delete the job
                                        if ~job_started
                                            disp('Apparently, the AMICA jobs didn''t get scheduled; restarting with ');
                                        else
                                            disp('Apparently, the AMICA computation didn''t init properly; restarting with ');
                                        end
                                        variant.numprocs = ceil(variant.numprocs*variant.reduce_factor);
                                        disp(['N = ' num2str(variant.numprocs)]);
                                        if length(job) >= 2
                                            if isa(job{2},'onCleanup')
                                                % clearing the job will implicitly lead to its deletion
                                                clear job;
                                            else isa(job{2},'function_handle')
                                                % old MATLAB version: we invoke the function in job{2} to delete it explicitly
                                                job{2}();
                                            end
                                        else
                                            disp('couldn''t delete the previous jobs; please delete them manually.');
                                        end
                                        
                                        if variant.numprocs <= 4 && ~variant.flaky_cluster
                                            disp('Could not acquire cluster resource; falling back to local computation. Note NaN checks cannot be performed in this mode.');
                                            
                                            % try it locally...
                                            switch variant.amica_version
                                                case 'devel'
                                                    r = runamica10(pre.data,[],size(pre.data,1),size(pre.data,2), arglist{:},'numprocs',variant.numprocs,'qsub','off');
                                                case 'stable'
                                                    r = runamica(pre.data,[],size(pre.data,1),size(pre.data,2), arglist{:},'numprocs',variant.numprocs,'qsub','off');
                                                case 'stable11'
                                                    r = runamica11(pre.data,[],arglist{:},'numprocs',variant.numprocs,'qsub','off');
                                                case 'stable12'
                                                    r = runamica12(pre.data,[],arglist{:},'numprocs',variant.numprocs,'qsub','off');
                                            end
                                            
                                            % ... and leave the while loop with a valid r
                                            break;
                                        else
                                            % re-start the job...
                                            [job,r,t0,iter_times,maxiters,lastlines,job_started,early_finish] = schedule_amica(pre.data,variant.scheduler,variant.amica_version,logfile,arglist{:},'numprocs',variant.numprocs);
                                        end
                                    end
                                end
                                
                                % done runnning on the cluster - load the results
                                if ~isfield(r,'W')
                                    switch variant.amica_version
                                        case 'stable12'
                                            r = loadmodout12(variant.outdir);
                                        case 'stable11'
                                            r = loadmodout11(variant.outdir);
                                        case 'devel'
                                            r = loadmodout10(variant.outdir);
                                        case 'stable'
                                            r = loadmodout(variant.outdir);
                                    end
                                end
                                
                                % append the iteration time samples to the current list of amica runtimes...
                                try
                                    f = fopen(env_translatepath('resources:/amica_runtimes.txt'),'a+');
                                    fprintf(f,' %d',iter_normtime);
                                    fclose(f);
                                catch
                                    try
                                        fclose(f);
                                    catch,end
                                end
                            end
                            
                            % check if the result data is broken, such that we need to fall back
                            % to a native MATLAB implementation
                            if isempty(r)
                                error('fall back'); end
                            if ~isfield(r,'S')
                                error('fall back'); end
                        catch e
                            % run the pure MATLAB variant
                            disp(['note: runamica failed to run; reason: ' e.message]);
                            disp('      falling back to MATLAB variant; available parameters: ');
                            disp('      num_models, num_mix_comps, do_mean/do_sphere, do_newton');
                            variant.max_iter = ceil(variant.max_iter * variant.fallback_reduce);
                            try
                                r = amica_native(pre,variant);
                            catch e
                                disp('Native MATLAB implementation of amica ran into an error. Giving up now.');
                                rethrow(e);
                            end
                        end
                    end
                    quicklog(logfile,'=== success; outdir=%s ===',variant.outdir);
                end
                
                % translate properties from r (the result) to pre (the data set)
                pre.icasphere  = r.S;
                if size(r.W,3) > 1
                    % we have multiple models...
                    % take the most likely model as weights and sphere
                    [x,bestmod] = max(r.mod_prob); %#ok<ASGLU>
                    pre.icaweights = r.W(:,:,bestmod);
                    % but retain all the rest in .etc.amica
                    if exist('mask','var')
                        r.sample_mask = mask; end
                else
                    pre.icaweights = r.W;
                end
                signal.etc.amica = r;                
            case 'infomax'
                % translate all booleans into 'on'/'off'
                variant = structfun(@(p) quickif(islogical(p),quickif(p,'on','off'),p),variant,'UniformOutput',false);
                arglist = hlp_struct2varargin(variant,'suppress',{'arg_selection','arg_direct'},'rewrite',{'pstop','stop','pweights','weights','pblock','block','pposact','posact'});
                [pre.icaweights,pre.icasphere] = hlp_diskcache('icaweights',@runica,pre.data,arglist{:});
            case 'beamica'
                sig = cov(pre.data');
                if ~isempty(variant.anchorlabels)
                    % calc beamformer penalty matrix for selected locations
                    [B,init_W,chanmask] = hlp_diskcache('filterdesign',@calc_beamformer_constraints,{pre.chanlocs.labels},variant.anchorlabels,sig,variant.reference); %#ok<ASGLU>
                    % remove missing channels
                    if ~any(chanmask)
                        error('None of your channels is in the head model; the AnchorLabels option can currently only be used for data with 10-20 labels.'); end
                    pre.data = pre.data(chanmask,:,:,:,:,:,:,:);
                    pre.chanlocs = pre.chanlocs(chanmask);
                    pre.nbchan = size(pre.data,1);
                else
                    B = {};
                end
                [pre.icaweights,pre.icasphere] = hlp_diskcache('icaweights',@beamica,pre.data,B,[],[],[],variant.max_iter,variant.lrate,variant.tradeoff,variant.verbose,variant.usegpu,variant.convergence_check);
            case 'fastica'
                % translate all booleans into 'on'/'off'
                variant = structfun(@(p) quickif(islogical(p),quickif(p,'on','off'),p),variant,'UniformOutput',false);
                arglist = hlp_struct2varargin(variant,'suppress',{'arg_selection','arg_direct'},'rewrite',{'stepsize','mu'});
                [dummy,pre.icaweights] = hlp_diskcache('icaweights',@fastica,pre.data,arglist{:}); %#ok<ASGLU>
                pre.icasphere = eye(size(pre.data,1));
            case 'rica'                
                variant.solverOptions.Method = hlp_rewrite(variant.solverOptions.Method, ...
                    'Steepest Descent','sd','Cyclic Steepest Descent','csd','Barzilai and Borwein Gradient','bb',...
                    'Non-Linear Conjugate Gradient','cg','Scaled Non-Linear Conjugate Gradient','scg','Preconditionined Non-Linear Conjugate Gradient','pcg',...
                    'Quasi-Newton with Limited-Memory BFGS Updating','lbfgs','Hessian-Free Newton','newton0',...
                    'Preconditioned Hessian-Free Newton','pnewton0','Quasi-Newton Hessian approximation','qnewton', ...
                    'Newton with Hessian update every k''th step','mnewton');
                variant.solverOptions.LS_init = hlp_rewrite(variant.solverOptions.LS_init,'Default',[],'Always1',0,'AsPrevious',1,'QuadraticInterp',2,'TwicePrevious',3,'ScaledConj',4);
                variant.solverOptions.LS_type = hlp_rewrite(variant.solverOptions.LS_type,'Default',[],'BacktrackingArmijo',0,'BracketingWolfe',1,'MatlabOptTB',2);
                variant.solverOptions.LS_interp = hlp_rewrite(variant.solverOptions.LS_interp,'Default',[],'DoubleOrBisect',0,'Cubic',1,'MixedQuadraticCubic',2);
                variant.solverOptions.LS_multi = hlp_rewrite(variant.solverOptions.LS_multi,'Default',[],'SameOrder',0,'HigherOrder',1);
                variant.chan_labels = {pre.chanlocs.labels};
                variant.chan_locs = pre.chanlocs;
                if doforce
                    [pre.icaweights,pre.icasphere,pre.icachansind] = rica(pre.data,variant,variant.solverOptions);
                else
                    [pre.icaweights,pre.icasphere,pre.icachansind] = hlp_diskcache('icaweights',@rica,pre.data,variant,variant.solverOptions);
                end
            case 'dictica'
                % use dictionary learning
                 [pre.icadict,pre.icasphere] = hlp_diskcache('icaweights',@dictlearn,variant,'X',pre.data,'chanlocs',pre.chanlocs);
            case 'kernelica'
                % translate all booleans into 0/1
                variant = structfun(@(p) quickif(islogical(p),double(p),p),variant,'UniformOutput',false);
                if strcmp(variant.kernel,'gaussian')
                    variant.sig = variant.GaussianSigma;
                else
                    variant.sig = variant.HermiteSigma;
                end
                arglist = hlp_struct2varargin(variant,'suppress',{'arg_selection','arg_direct','GaussianSigma','HermiteSigma'},'rewrite',{'contrastfun','contrast','dodisp','disp'});
                pre.icaweights = hlp_diskcache('icaweights',@kernel_ica_options,pre.data,arglist{:});
                pre.icasphere = eye(size(pre.data,1));
            case 'fastkernelica'
                % fast KernelICA is a euphemism
                C = size(pre.data,1);
                pre.icasphere = 2.0*inv(sqrtm(cov(pre.data'))); %#ok<MINV>
                guess = eye(C);
                for k=1:variant.restarts
                    tmpdata = (pre.icasphere * pre.data);
                    [Ws{k}, dummy, HSIC] = hlp_diskcache('icaweights',@fastkica,tmpdata, guess, variant.maxiter, variant.psigma, variant.thresh); %#ok<ASGLU>
                    guess = eye(C) + randn(C);
                    HSICs(k) = min(HSIC);
                end
                [dummy,bestidx] = min(HSICs); %#ok<ASGLU>
                pre.icaweights = Ws{bestidx};
            case 'sphere'
                 pre.icasphere = inv(sqrtm(cov(pre.data')));
                 pre.icaweights = eye(size(pre.data,1));
            case 'robust_sphere'
                 pre.icasphere = inv(sqrtm(hlp_diskcache('icaweights',@cov_blockgeom,pre.data')));
                 pre.icaweights = eye(size(pre.data,1));
            otherwise
                % let pop_runica handle all the rest
                arglist = hlp_struct2varargin(variant,'suppress',{'arg_selection','arg_direct'});
                pre = pop_runica(pre,'icatype',variant.label,arglist{:});
        end
        
        
        % add channel indices, if necessary
        try
            if isempty(pre.icachansind)
                pre.icachansind = 1:size(pre.data,1); end
            % add weight inverse (mixing matrix), if necessary
            if isempty(pre.icawinv)
                pre.icawinv = pinv(pre.icaweights*pre.icasphere); end;
            % normalize winv, for better plotting
            if normalize_weights
                normalizer = sqrt(1./sum(pre.icawinv.*pre.icawinv));
                pre.icawinv = single(pre.icawinv .* (ones(size(pre.icawinv,1),1)*normalizer));
                pre.icaweights = diag(1./normalizer)*pre.icaweights;
            end          
            
            if cleaned_data
                signal = pre;
                if trials ~= 1
                    warning('BCILAB:flt_ica:data_loss','The original epoch and event information has been erased from the cleaned data during artifact rejection; cannot restore.');
                    signal.data = reshape(signal.data,chans,pnts,trials);
                end
            else
                % carry over the fields to the output data set
                signal.icaweights = pre.icaweights;
                signal.icasphere = pre.icasphere;
                signal.icawinv = pre.icawinv;
                % if channel rejection was used, we remember from which channels the ica map was derived
                [dummy,idxeeg,idxpre] = intersect({signal.chanlocs.labels},{pre.chanlocs.labels}); %#ok<ASGLU>
                [dummy,idxpre_order] = sort(idxpre);  %#ok<ASGLU>
                signal.icachansind = idxeeg(idxpre_order);
            end
        catch e
            if dodebug
                disp('Error in flt_ica code; waiting for user input. You can disable this behavior by setting the DebugMode option to false.');
                keyboard;
            else
                rethrow(e);
            end
        end
    else
        disp('The dataset already contains an IC decomposition; skipping... (note: set parameter ''doforce'' to true to force re-computation).');
    end
    
    % generate root chanlocs, if not already there
    if ~exist('root_chanlocs','var')
        root_chanlocs = signal.chanlocs(signal.icachansind); end
    
    % store the ICA decomposition in the state
    state = struct('icasphere',{signal.icasphere}, 'icaweights',{signal.icaweights}, 'icawinv',{signal.icawinv}, 'icachansind',{signal.icachansind},'root_chanlocs',{root_chanlocs});
    if isfield(signal.etc,'amica')
        state.amica = signal.etc.amica; end
end


% keep track of the last ICA decomposition for inspection
global tracking; tracking.inspection.ica = state;

if do_calcact
    % just populate icaact
    signal.icaact = (signal.icaweights*signal.icasphere)*signal.data(signal.icachansind,:); end

if do_transform
    % transform the data itself, if necessary
    if isfield(signal.etc,'amica') && size(signal.etc.amica.W,3) > 1
        warn_once('Note: The signal will only be transformed according to the 1st amica model.'); end
    signal.data = (signal.icaweights*signal.icasphere)*signal.data(signal.icachansind,:);
    if retain_labels && nnz(signal.icachansind) == size(signal.data,1)
        signal.chanlocs = signal.chanlocs(signal.icachansind);
    else
        signal.chanlocs = struct('labels',cellfun(@num2str,num2cell(1:length(signal.icachansind),1),'UniformOutput',false)); 
    end
    signal.nbchan = size(signal.data,1);
    if clear_after_trans
        signal.icaweights = [];
        signal.icawinv = [];
        signal.icasphere = [];
    end
end

% remember which cleaning level was used
try
    signal.etc.clean_settings = data_cleaning;
catch
end

exp_endfun;


% schedule an AMICA run on an Sun Grid Engine cluster...
function [job,result,t0,iter_times,maxiters,lastlines,job_started,early_finish] = schedule_amica(X,scheduler,amica_version,logfile,varargin)

% schedule the job, get console output
switch amica_version
    case 'devel'
        [conout,dummy] = evalc('runamica10(X,[],size(X,1),size(X,2), varargin{:})'); %#ok<NASGU>
    case 'stable'
        [conout,dummy] = evalc('runamica(X,[],size(X,1),size(X,2), varargin{:})'); %#ok<NASGU>
    case {'stable11','stable12'}
        % resolve the queue auto-detection
        args = hlp_varargin2struct(varargin);
        if iscell(args.use_queue)
            startwait = now;
            while true
                fprintf('Scanning available queues...');
                [status,info] = system(['ssh ' scheduler ' -x "qstat -g c"']); %#ok<ASGLU>
                lines = hlp_split(info,sprintf('\n'));
                keys = hlp_split(lines{1},' ');
                availrow = find(strcmpi(keys,'avail'))-1;
                lines = lines(3:end);
                qavail = [];
                qname = {};
                qdata = {};
                for l=1:length(lines)
                    line = lines{l};
                    vals = hlp_split(line,' ');
                    qname{l} = vals{1};
                    qavail(l) = str2num(vals{availrow});
                    qdata{l} = [vals{1} ':' vals{availrow}];
                end
                % optionally backtrack to a 50% smaller queue
                if ~isempty(qavail) && ~any(qavail>=args.numprocs) && any(qavail>=args.numprocs/2)
                    args.numprocs = args.numprocs/2;
                    retain = qavail>=args.numprocs;
                    qdata = qdata(retain);
                    qname = qname(retain);
                    qavail = qavail(retain);
                end
                % our options
                options = intersect(args.use_queue,qdata);
                if isempty(options)
                    enough = find(qavail>=32);
                    if (now - startwait) > 6*3600 && ~isempty(enough)
                        % waited for more than a few hours
                        newargs = args;
                        newargs.use_queue = qname{enough(1)};
                        varargin = hlp_struct2varargin(newargs);
                        if newargs.numprocs == 0
                            newargs.numprocs = qavail; end
                        fprintf('patience exhausted. Using queue %s (%i processors).\n',newargs.use_queue,newargs.numprocs);
                        break;
                    else
                        fprintf('no available queue found; waiting...\n');
                        % wait for 10 minutes...
                        pause(10*60);
                    end
                else
                    % pick a random one
                    pick = options{1+floor(rand*length(options)-eps)};
                    parts = hlp_split(pick,':');
                    % rebuild varargin
                    newargs = args;
                    newargs.use_queue = parts{1};
                    if newargs.numprocs == 0
                        newargs.numprocs = qavail; end
                    varargin = hlp_struct2varargin(newargs);
                    fprintf('using queue %s (%i processors).\n',newargs.use_queue,newargs.numprocs);
                    break;
                end
            end
        end
        quicklog(logfile,'=== new attempt/%s: [%i x %i]; numprocs=%i; num_models=%i; outdir=%s',amica_version,size(X,1),size(X,2),newargs.numprocs,newargs.num_models, newargs.outdir);
        if strcmp(amica_version,'stable11')
            [conout,dummy] = evalc('runamica11(X,args.outdir,varargin{:})'); %#ok<NASGU>
        else
            [conout,dummy] = evalc('runamica12(X,args.outdir,varargin{:})'); %#ok<NASGU>
        end
end

try
    % take apart the console output of AMICA to find the qsub id....
    lines = hlp_split(conout,10);
    id = [];
    for l = 1:length(lines)
        if ~isempty(strfind(lines{l},'qsub id ='))
            parts = hlp_split(lines{l},'=');
            id = strtrim(parts{2});
            if isempty(id)
                error('Amica could not submit a job to the queue. Please check whether you are on a computer from which you can submit to that queue.'); end
            id = str2num(id);
            try
                disp(['current qsub id: ' num2str(id)]);
            catch,end
            break;
        end
    end
    if ~id
        error('didn''t find the qsub id...'); end
    % get the job's outdir
    args = hlp_varargin2struct(varargin);
    clean_function = @()delete_job(id,scheduler,args.outdir);
    job = {id, onCleanup(clean_function)};
catch
    disp('couldn''t identify job id, jobs will not be auto-deleted.');
    job = [];
end
result = [];
% also set a few other things that are associated with a new amica run
t0 = tic;
iter_times = [];
maxiters = 0;
lastlines = 0;
job_started = false;
early_finish = false;



% native MATLAB implementation of AMICA
function r = amica_native(pre,variant)
if ~isfield(variant,'do_mean')
    variant.do_mean = 1; end
if ~isfield(variant,'do_sphere')
    variant.do_sphere = 1; end
[r.W r.A r.c r.LL r.LLt r.mod_prob] = amica10(pre.data, variant.num_models, variant.num_mix_comps, variant.max_iter, variant.do_mean|variant.do_sphere, variant.do_newton);
% recompute some missing parameters
r.S = eye(size(pre.data,1));
r.v = zeros(size(r.LLt));
for m = 1:size(r.v,1)
    r.v(m,:) = 1./sum(exp(bsxfun(@minus,r.LLt,r.LLt(m,:))),1); end



% delete a running (or possibly completed) amica job
function delete_job(id,scheduler,outdir)
% first delete the job from qsub
system(['ssh ' scheduler ' qdel ' num2str(id)]);
% try to clean up the out.txt file, if the job did not get done...
job_done = false;
file_exists = false;
try
    t = fopen([outdir filesep 'out.txt']);
    file_exists = (t ~= -1);
    if file_exists
        while 1
            line = fgetl(t);
            if ~ischar(line)
                break; end
            if ~isempty(strfind(line,'done.'))
                job_done = true;
                break;
            end
        end
        fclose(t);
    end
catch
    try
        fclose(t);
    catch,end
end

if ~job_done && file_exists
    try
        if exist([outdir filesep 'out.txt'],'file')
            delete([outdir filesep 'out.txt']); end
        disp('Deleted out.txt file.');
    catch
        disp('Could not delete out.txt file.');
    end
end
