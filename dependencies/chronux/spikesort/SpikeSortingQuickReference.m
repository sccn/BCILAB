%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SORTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load spikes2

%% The following commands run the Fee algorithm
spikes = ss_dejitter(spikes);
spikes = ss_outliers(spikes);
spikes = ss_kmeans(spikes);
spikes = ss_energy(spikes);
spikes = ss_aggregate(spikes);

%% The following commands merge/split clusters
% spikes = merge_clusters(spikes, clusternumber1, clusternumber2);
% spikes = split_cluster(spikes, clusternumber);


%% The following commands run Ken Harris' KlustKwik
projection = pcasvd(spikes.waveforms);
[assigns,log] = kkm(spikes.waveforms(:,1:3),0);

%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General visualization tools in 2 and 3 dimensions -- results depend on
%% whether you have run ss_kmeans or ss_aggregate.  Can also specify the
%% assignments list as in:  ssg_databrowse2d(spikes,assigns);
ssg_databrowse2d(spikes);
ssg_databrowse3d(spikes);


%% You can always look at:
figure;  plot(spikes.waveforms');  axis tight;  title('Spike Waveforms');
figure;  histxt(spikes.waveforms);  colormap jet(256);

%% After running the ss_outliers step, try
figure; plot(spikes.outliers.waveforms'); axis tight;

%% After running the ss_kmeans step, try
figure;  set(gcf, 'Renderer', 'OpenGL');  clusterXT(spikes, spikes.overcluster.assigns);  title('Local Clusters');

%% After running the ss_energy step, try
figure; set(gcf, 'Renderer', 'OpenGL');  clusterXT(spikes, spikes.hierarchy.assigns); title('Final Clusters');
figure; showclust(spikes, spikes.hierarchy.assigns);
figure; correlations(spikes, spikes.hierarchy.assigns);
figure; aggtree(spikes); title('Aggregation Tree');
