function env_clear_memcaches
% Clear the in-memory caches of BCILAB.
%
% BCILAB caches intermediate results in a memory cache. Under certain circumstances, especially when
% helper functions have been edited, BCILAB may fail to re-execute the just changed code and instead
% pull its outputs from the cache. This function clears that cache.
%
% Note that there is also a disk cache, which is being considered when computations take a relatively
% long time to finish. When debugging signal processing functions in BCILAB, it may be useful to 
% disable these caches in the startup configuration, and/or delete the cache folder/files.
%
% This function also clears various other small caches, and can be used to purge stray instances
% of Paradigm classes from memory, so that changes can be made to them without having to restart 
% MATLAB.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-04-15

global tracking;

tracking.cache.data = struct();
tracking.cache.times = struct();
tracking.cache.sizes = struct();

% clear all micro-caches
hlp_microcache('clear');

% and get rid of class instances, too...
env_clear_classes;
