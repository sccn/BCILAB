function isOpen = hlp_pll_openPool(numWorkers,profileName)
% open a parallel toolbox pool

if nargin<2
    profileName = 'local';
end

sz = matlabpool('size');
isOpen = (sz>0);
% check if matlabpool is open...
if ~isOpen
    % ... if not, open one
    matlabpool(profileName,numWorkers);
elseif sz~=numWorkers
    % wrong number of workers
    % ...restart pool
    matlabpool('close');
    matlabpool(profileName,numWorkers);
end