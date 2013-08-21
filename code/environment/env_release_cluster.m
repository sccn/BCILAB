function env_release_cluster
% Release any currently acquired cluster resources.
%
% This stops the heartbeat signal; it is the workers' responsibility to clean themselves up when
% they are no longer needed.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-04-12

disp('Releasing cluster resources...');
global tracking;
% this will lead to the deletion of any involved heartbeat timer.
tracking.cluster_requested = [];
