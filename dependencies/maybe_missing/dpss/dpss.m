function dps_seq = dpss(seq_length,time_halfbandwidth,k)

% DPSS approximate the discrete prolate spheroidal (DPSS), or Slepian sequences
% using a precomputed version on disk

persistent w
if isempty(w)
  % load the precomputed tapers only once
  load precompute_dpss
end

n = 1000;
s = 1:0.5:50;

% find the nearest match
i = nearest(s, time_halfbandwidth, true, true);
if k>size(w{i},2)
    error('More tapers were requested than have been precomputed.'); end

% interpolate onto the desired number of samples
dps_seq = interp1(linspace(0,1,n), w{i}(:,1:k), linspace(0,1,seq_length));
