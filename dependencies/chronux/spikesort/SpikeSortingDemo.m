%% Demonstration code for the spike sorter derived from the Fee et al.
%% (J. Neurosci Meth (1996) 69: 175-88) algorithm.  Requires Matlab
%% version 6.1 or later and the Statistics Toolbox.
%%                                          Last Modified: 8/13/06
%%                                          samar mehta, sbmehta@ucsd.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTRO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The sorting procedure itself is fairly easy to perform (although it is
%%  currently mostly command line driven -- not much GUI yet -- and you 
%%  have to perform your own spike detection to create the inputs described
%%  below).  Automated spike sorting, however, always requires a human
%%  'sanity check' and can sometimes also need intervention.  To this end,
%%  several visualization tools are included.  This file steps through an
%%  extended example of their use.
%%
%% We start with an overview of how to sort if you want to get started
%%  quickly.  This is followed by example code and example data which
%%  will hopefully clarify the details.  The example does not, however,
%%  cover everything that the code can do . . . if there is something that
%%  you feel is not working or information that you feel would be useful,
%%  please feel free to look at the code itself (which is fully documented)
%%  and to use HELP.  (or email sbmehta@ucsd.edu)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OVERVIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The input to the sorting code is a Matlab structure with the following
%%  fields:    {WAVEFORMS, SPIKETIMES, FS, THRESHT, THRESHV}
%%     WAVEFORMS is an N x M matrix of spike data, where each row contains
%%      the voltage values from a single spike (this can be from a single
%%      electrode or e.g., the concatenated data from a tetrode).  More
%%      generally, any feature (not necessarily voltage) can appear in this
%%      matrix, although it is your responsibility to make sure that they
%%      are normalized such that differing units between features do not
%%      dominate the variance of the data.
%%     SPIKETIMES is an N x 1 vector of spike times (in seconds).  The nth
%%      entry should correspond to the nth row of 'waveforms'.
%%     FS is the sampling rate (in Hz) of the recording.
%%     THRESHT is the column index in WAVEFORMS where threshold crossing
%%      occurred (i.e., we are assuming the the spikes were extracted from a
%%      raw data trace via threshold detection).  For example, if the THRESHT
%%      is 10 and WAVEFORMS is N x 30, then the 10th column of WAVEFORMS
%%      contains the sample that crossed threshold, with 9 samples before and
%%      20 samples after.  The threshold detection must be done in such a way
%%      that no two detected spikes are ever less than M samples apart.
%%     THRESHV is a pair of numbers corresponding to the high and low voltage
%%      thresholds used to extract the spikes.  If only one threshold was used,
%%      the other can be given as +/- Inf.  For example, [-50, Inf] means that
%%      spikes were detected when the voltage dropped below -50 (arbitrary units)
%%      and no positive going threshold was defined.
%%  The WAVEFORMS, SPIKETIMES and FS are required.  The code will run
%%   if THRESHT and THRESHV are not defined, but the dejittering step will
%%   not work and the height/width graphical displays may produce errors.
%%  Sorting options can also be defined in the SPIKES structure to direct the
%%   algorithm.  The currently defined options are:
%%     SPIKES.OPTIONS.OUTLIER_CUTOFF     (see help for SS_OUTLIERS)
%%     SPIKES.OPTIONS.REFRACTORY_PERIOD  (see help for SS_AGGREGATE)
%%   These are not mandatory; if not defined, default values are used.
%%
%% The end product of the sorting is an N x 1 vector of assignments, with any
%%   outliers assigned to a cluster with label 0.  If you prefer to have the
%%   outliers segregated from the data, see HELP SS_AGGREGATE; in this case,
%%   the assignment vector will be of size P x 1, where P is the number of
%%   spikes not treated as outliers, and the WAVEFORMS and SPIKETIMES variables
%%   will also be modified to have P rows.  The remaining spikes (i.e., the
%%   outliers) can then be found in a new field called OUTLIERS.
%%   The assignments are numerical labels with no inherent meaning; their only
%%   interpretation is that spikes labeled with, e.g., 1 belong to the same cluster.
%%
%% If the input information is collected into a Matlab structure called 'spikes',
%%  the following code will produce the output assignments vector:
%          spikes = ss_dejitter(spikes);
%          spikes = ss_outliers(spikes);
%          spikes = ss_kmeans(spikes);
%          spikes = ss_energy(spikes);
%          spikes = ss_aggregate(spikes);
%          assignments = spikes.hierarchy.assigns;
%%
%% If this seems cryptic, read on ...
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP-BY-STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The above code will sort the spikes but provides little feedback along
%%  the way.  Until this is all wrapped into a GUI, effective sorting requires
%%  the use of the included command line tools to give you a sense for how
%%  the data look.  Note that this some of these steps might be easier to
%%  follow if you browse through the Fee et al. paper first.
%%
%% If you have put your own data into a 'spikes' data structure as described
%%  above, you can run the code on your own data.  This code, however, includes
%%  three demonstration data sets chosen to highlight the functionality.  i
%%  suggest that that you try one of these data sets first to make sure that
%%  things run the same on your computer as they do on mine.

%% *** The directory containing this file, the supporting code, the demo data
%%     and all of its subdirectories must be added to the Matlab path.  The
%%     following line accomplishes this if you haven't done it already.
% addpath(genpath(pwd), '-end');
%%     You can also just run the startup.m script in this directory.

%%%%% LOADING DATA
% If you are using your own data, load it into a variable named 'spikes' now.
% If you want to use the demo data, you can just uncomment and run one of the
% following lines:

% load spikes1;    % Real data (barrel cortex).  Good sig to noise but electrode drifts 3/4 way through.
% load spikes2;    % Simulated data (Bionics simulator).  High sig to noise (ideal case).
% load spikes3;    % Real data (barrel cortex).  Mediocre sig to noise; tough to sort manually.


%%%%% LOOKING AT THE RAW DATA
% The following commands plot the raw data in two ways (there are, of course, 
% many others).  The top plot is simply all of the raw waveforms overlaid on top
% of one another.  If you have more than 100-200 spikes, it becomes impossible to
% determine the density of the spikes, however.  So the bottom plot displays the
% same information but as a density after binning in time/voltage (see help for
% 'histxt').  The demo file 'spikes2' is a good example of why this helps.
figure(1); colormap jet;
subplot(2,1,1); plot(spikes.waveforms'); axis tight; title('Raw Data');
subplot(2,1,2); histxt(spikes.waveforms); h = gca;

% You can also use the following two functions at any point throughout the 
% following code.  Try running them at various points to visualize the data
% cluster in 2-D or 3-D as it is sorted.  Also try:
%    Click on the axis labels to change the features used for the projection.
%    Double-click outside the plot to draw the current data as a density.
%    For SSG_DATABROWSE2D, clicking on a data point will make that point's
%        cluster easier to see.
%    For SSG_DATABROWSE3D, clicking and dragging on the axis allows you to
%        rotate the data in 3 dimensions.
ssg_databrowse2d(spikes);
ssg_databrowse3d(spikes);

%%%%% ALIGNING SPIKES
% Noise on the electrode can jitter the exact time at which threshold crossing
% occurs.  If you are using this method to extract your spikes, this can be
% a significant source of variability for spikes from the same neuron.  The
% next step 'de-jitters' the data by aligning all spikes on the same sample.
% (Note: this requires THRESHT and THRESHV to be defined as described above).
% The data are then replotted just as they were in figure 1.  Put both plots
% side by side.  The density plot for the dejittered data should look tighter.
spikes = ss_dejitter(spikes);
figure(2); colormap jet;
subplot(2,1,1); plot(spikes.waveforms'); axis tight; title('Centered Data');
subplot(2,1,2); histxt(spikes.waveforms); clim = get(gca, 'CLim'); if(ishandle(h)), set(h, 'CLim', clim); end;

%%%%% REMOVING OUTLIERS
% The density estimation techniques used in spike sorting routines are typically
% not robust -- that is, they are sensitive to outliers -- so we need to remove
% outliers.  'outliers' are spikes that look so much unlike the other events in
% the data that they cannot be sorted reliably.  This does not mean that they are
% not interesting (they often include overlapping spikes/doublets in addition to
% electrical artifact), it just means that you'll have to look at them by hand.
% The following code removes the worst offenders.  See the help for 'ss_outliers'
% if you want to be more/less conservative.  (We don't replot density because
% this won't change much since outliers are only a small percent of the data).
spikes = ss_outliers(spikes);
figure(3); 
plot(spikes.waveforms'); axis tight; title('Centered Data w/ Outliers Removed');
% (note: if you want to look at the outliers themselves, you can do this as follows:)
% figure(9); plot(spikes.outliers.waveforms'); axis tight;
% In this plot, some of the spikes might not look like outliers.  This is often
% a visual artifact arising from plotting overlapping waveforms -- try looking at
% them one by one before resorting to the following:
% NOTE: If the above function is removing spikes that you do not consider
%         outliers, you have two choices.  You can manually assign these waveforms
%         to clusters after the automatic clustering is done or you can lower the
%         sensitivity of the outliers function.  Lowering the sensitivity may
%         cause problems, though, because if severe enough outliers are left in
%         the data, they will confuse the quality of the automatic clustering.
%         With that warning, to change sensitivity, set the following:
%                   spikes.options.outlier_cutoff = cutoff;
%         before the ss_outliers call above.  Cutoff is normally (1 - 1/M),
%         where M is the total number of spikes; this means that _using the
%         outliers statistical heuristics_ on average one spike in the data set
%         will be rejected incorrectly.  In practice, the number can be anything
%         close to 1.0 -- e.g., to make the algorithm throw away fewer waveforms,
%         try a cutoff of (1 - 0.0001/M).

%%%%% INITIAL CLUSTERING
% The Fee algorithm deals with possibly non-Gaussian data (e.g., bursting neurons)
% by doing the sorting in two steps.  The first step fits many local Gaussian
% spheres to the data to identify groups of spikes with similar shapes; these will
% later be combined into spike assignments.  This two step procedure is a good place
% to do a sanity check; do the results of the local clustering look like the
% algorithm is capturing local density?  The following code plots the waveforms
% (now colored according to local similarity) and the type of height-width that
% is often used for manual sorting (colored similarly).
spikes = ss_kmeans(spikes);
figure(4);  set(gcf, 'Renderer', 'OpenGL');
clusterXT(spikes, spikes.overcluster.assigns);  title('Local Clusters');
% (note: this is a good time to check stationarity of your data.  For example,
%         the data set in 'spikes1' contains electrode drift towards the end of the
%         record.  You can see this be looking at the local clusters vs. time:)
% figure(9); plot(spikes.spiketimes, spikes.overcluster.assigns, '.'); xlabel('Time (sec)'); ylabel('Cluster #');
% (note: you can get much of the information from the above two plots by re-running
%         the SSG_DATABROWSE functions.  e.g.,)
% ssg_databrowse2d(spikes);

%%%%% FINAL CLUSTERING
% The local clusters are now combined into larger groups corresponding (with luck)
% to neurons.  This step is typically the most time consuming (especially for
% larger data sets).  After the aggregation, the final results are summarized in
% a series of plots.
% Figure 5 contains one row of graphs for each final cluster.  The left plot is
% a density image and the two right plots are interspike interval histograms
% (at two time scales).  If the clustering is good, there should be few events
% with interspike intervals less than 2 msec; the ISI score displayed in the
% middle column reflects this (smaller numbers are better).  
% Figure 6 looks like Figure 4 from above, but recolored after aggregation.
% The bottom plot also contains a legend matching the colors to the cluster #
% labels in Figure 5.
% Figure 7 shows how the local clusters were combined into the final clusters.
% The color and length of the lines in this tree are described in the help
% text for SS_AGGREGATE.
% Figure 8 shows cross-correlations in spike timing between the final clusters.
spikes = ss_energy(spikes); spikes = ss_aggregate(spikes);
figure(5); colormap jet;
showclust(spikes, spikes.hierarchy.assigns);
figure(6); set(gcf, 'Renderer', 'OpenGL');
clusterXT(spikes, spikes.hierarchy.assigns); title('Final Clusters');
figure(7); 
aggtree(spikes); title('Aggregation Tree');
figure(8);
correlations(spikes);  title('Auto- and Cross- Correlations');

% This is, once again, a good time to run the SSG_DATABROWSE function to look at
% projections of the clustered data and to check if the results seem reasonable.
% ssg_databrowse3d(spikes);


% (note: For neurons with relatively high firing rates, the algorithm can watch
%         for refractory periods to determine when to stop the aggregation.  Since
%         this does not work well for neurons with very low firing rates, it can
%         sometimes make mistakes towards the end of the aggregation, joining
%         clusters that should not be joined or vice versa.  A person looking at
%         figures 5-7 can usually see this fairly quickly.  The problem can be
%         corrected using the following types of commands:)
% spikes = merge_clusters(spikes, clusternumber1, clusternumber2);
% spikes = split_cluster(spikes, clusternumber);
%         After you do this, replot figures 5-7 above and see if things look better).

%%%%% Getting the cluster assignments
% The vector of cluster assignments can be found in the following location:
%                          spikes.hierarchy.assigns
