function [data,age,ndiscarded] = onl_lslreceive(inlet,timeout,maxage)
if nargin<2
    timeout = -1;
end
if nargin<3
    maxage = Inf;
end
ndiscarded = 0;

[tmpsample,timestamp] = inlet.pull_sample(timeout,true);

if isempty(timestamp) || timestamp == 0
    data = [];
    age = NaN;
    return;
end

% tic;
% calculate the age of the sample
age = lsl_local_clock(inlet.get_libhandle()) - (timestamp+inlet.time_correction());
% fprintf('age: %0.5g\n',age);
while age > maxage
    % discard samples that have expired
    [tmpsample,timestamp] = inlet.pull_sample(timeout,true);
    age = lsl_local_clock(inlet.get_libhandle()) - (timestamp+inlet.time_correction());
    ndiscarded = ndiscarded + 1;
end
% toc

if isempty(tmpsample)
    data = [];
    age = NaN;
    return;
end

% now we have a sufficiently current sample
% ...decode the data
% tic;
data = hlp_deserialize(tmpsample{1});
% toc