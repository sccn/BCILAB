function [ci,bootstat]  = stat_bootci(nboot,bootfun,varargin)
%BOOTCI Bootstrap Confidence Interval
%   CI = BOOTCI(NBOOT,BOOTFUN,...) computes the 95 percent BCa bootstrap
%   confidence interval of the statistic defined by the function BOOTFUN.
%   NBOOT is a positive integer indicating the number of bootstrap data
%   samples used in the computation. BOOTFUN is a function handle specified
%   with @. The third and later input arguments to BOOTCI are data
%   (scalars, column vectors, or matrices) that are used to create inputs
%   to BOOTFUN. BOOTCI creates each bootstrap sample by sampling with
%   replacement from the rows of the non-scalar data arguments (these must
%   have the same number of rows). Scalar data are passed to BOOTFUN
%   unchanged. 
%
%   If BOOTFUN returns a scalar, CI is a vector containing the lower and
%   upper bounds of the confidence interval. If BOOTFUN returns a vector of
%   length M, CI is an array of size 2-by-M, where CI(1,:) are lower bounds
%   and CI(2,:) are upper bounds. If BOOTFUN returns an array of size
%   M-by-N-by-P-by-..., CI is an array of size 2-by-M-by-N-by-P-by-...,
%   where CI(1,:,:,:,...) is an array of lower bounds and CI(2,:,:,:,...)
%   is an array of upper bounds.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},'alpha',ALPHA) computes the 100*(1-ALPHA)
%   percent BCa bootstrap confidence interval of the statistic defined by the
%   function BOOTFUN. ALPHA is a scalar between 0 and 1. The default value of
%   ALPHA is 0.05. A cell array groups BOOTFUN and the arguments used to create
%   inputs to it. ALPHA and any other arguments to BOOTCI appear outside the
%   cell array.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},...,'type',TYPE) computes the bootstrap
%   confidence interval of the statistic defined by the function BOOTFUN.
%   TYPE is the confidence interval type, specifying different methods of
%   computing the confidence interval. TYPE is a string chosen from
%       'norm' or 'normal':               normal approximated interval with
%                                         bootstrapped bias and standard
%                                         error;                                        
%       'per' or 'percentile':            basic percentile method; 
%       'cper' or 'corrected percentile': bias corrected percentile method;
%       'bca' :                           bias corrected and accelerated 
%                                         percentile method;
%       'stud' or 'student':              studentized confidence interval.
%   The default value of TYPE is 'bca'.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},...,'type','stud','nbootstd',NBOOTSTD)
%   computes the studentized bootstrap confidence interval of the statistic
%   defined by the function BOOTFUN. The standard error of the bootstrap
%   statistics is estimated using bootstrap with NBOOTSTD bootstrap data
%   samples. NBOOTSTD is a positive integer value. The default value of
%   NBOOTSTD is 100.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},...,'type','stud','stderr',STDERR)
%   computes the studentized bootstrap confidence interval of statistics
%   defined by the function BOOTFUN. The standard error of the bootstrap
%   statistics is evaluated by the function STDERR. STDERR is a function
%   handle created using @. STDERR should take the same arguments as
%   BOOTFUN and return the standard error of the statistic computed by
%   BOOTFUN.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},...,'Weights',WEIGHTS) specifies
%   observation weights. WEIGHTS must be a vector of non-negative numbers
%   with at least one positive element. The number of elements in WEIGHTS
%   must be equal to the number of rows in non-scalar input arguments to
%   BOOTFUN. To obtain one bootstrap replicate, BOOTSTRP samples N out of N
%   with replacement using these weights as multinomial sampling
%   probabilities.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},...,'Options',OPTIONS) contains options
%   that specify whether to compute bootstrap iterations in parallel, and how
%   to use random numbers during the bootstrap sampling. OPTIONS is a struct
%   that can be created by a call to STATSET. BOOTCI uses the following fields:
%       'UseParallel'
%       'UseSubstreams'
%       'Streams'
%   For information on these fields see PARALLELSTATS.
%   NOTE: if 'UseParallel' is 'always' and 'UseSubstreams' 
%   is 'never', then the length of Streams must equal the number 
%   of processors used by BOOTCI. There are two possibilities. 
%   If a MATLAB pool is open, then Streams is the same length as
%   the size of the MATLAB pool. If a MATLAB pool is not open,
%   then Streams must supply a single random number stream.
%   
%   [CI,BOOTSTAT] = BOOTCI(...) also returns the bootstrapped statistic
%   computed for each of the NBOOT bootstrap replicate samples.  Each row of
%   BOOTSTAT contains the results of applying BOOTFUN to one bootstrap sample.
%   If BOOTFUN returns a matrix or array, then this output is converted to a
%   row vector for storage in BOOTSTAT.
%
%   Example:
%     Compute the confidence interval for the capability index in
%     statistical process control:
%          y = normrnd(1,1,30,1);                  % simulated process data
%          LSL = -3;  USL = 3;                     % process specifications
%          capable = @(x) (USL-LSL)./(6* std(x));  % process capability
%          bootci(2000,capable, y)                 % Bca confidence interval
%          bootci(2000,{capable, y},'type','per')  % basic percentile method
%
%   See also: BOOTSTRP, JACKKNIFE, STATSET, STATGET, 
%   RANDSAMPLE, PARFOR, PARALLELSTATS.

% The BCa method is described in the following references:
%
%     T.J. DiCicio and B. Efron (1996), "Bootstrap confidence intervals,"
%     Statistical Science, v. 11, n. 3, pp. 189-228.
% 
%     B. Efron and R.J. Tibshirani (1993), An Introduction to the
%     Bootstrap, Chapman & Hall, New York.
%
% Their formula involves a Z0 factor that is computed using the proportion
% of bootstrap values less than the original sample value.  In order to get
% reasonable results when the sample is lumpy, we include half of the
% bootstrap values that are tied with the original sample value when we
% compute Z0.

% Copyright 2005-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.8 $  $Date: 2011/07/20 00:08:12 $

if nargin<2
    error(message('stats:bootci:TooFewInputs'));
end;
if nboot<=0 || nboot~=round(nboot)
    error(message('stats:bootci:BadNboot'))
end; 

if ~iscell(bootfun) % default syntax
    type = 'bca';
    alpha = .05;
    fun = bootfun;
    data = varargin;
    weights = [];
    bootstrpOptions = statset('bootstrp');
else % syntax with optional type, alpha, nbootstd, and stderrfun name/value pairs
    fun = bootfun{1};
    data = bootfun(2:end);
    pnames = {'type', 'alpha', 'stderr', 'nbootstd', 'weights', 'options'};
    dflts =  {'bca', .05, [], 100, [], statset('bootstrp')};
    [type,alpha,stderrfun,nbootstd,weights,bootstrpOptions] = ...
                           internal.stats.parseArgs(pnames, dflts, varargin{:});
end

% error check for the bootfun
try 
    obsstat = fun(data{:});    
catch ME
    m = message('stats:bootci:BadBootfun',func2str(fun));
    throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
end
if any(~isfinite(obsstat))
        error(message('stats:bootci:NonfiniteBootfun'));
end; 

% Save the original size of stat and reshape it
sz = [];
if ~isvector(obsstat)
    sz = size(obsstat);
end
obsstat = obsstat(:)'; % turn into a row-vector

% call subfunctions to compute the intervals 
switch (lower(type))
    case {'norm','normal'}
        [ci,bootstat] = bootnorm(obsstat,nboot,fun,alpha,weights,bootstrpOptions,data{:});
    case {'per','percentile'}
        [ci,bootstat] = bootper(nboot,fun,alpha,weights,bootstrpOptions,data{:});
    case {'cper', 'corrected percentile'}
        [ci,bootstat] = bootcper(obsstat,nboot,fun,alpha,weights,bootstrpOptions,data{:});
    case 'bca'
        [ci,bootstat] = bootbca(obsstat,nboot,fun,alpha,weights,bootstrpOptions,data{:});
    case {'stud','student'}
        [ci,bootstat] = bootstud(obsstat,nboot,fun,alpha,nbootstd,stderrfun,...
            weights,bootstrpOptions,data{:});
    otherwise
        error(message('stats:bootci:BadType'))
end;

% Reshape
if ~isempty(sz)
    ci = reshape(ci,[2 sz]);
end

end   % bootci()
 
%-------------------------------------------------------------------------    
function [ci,bstat] = bootnorm(stat,nboot,bootfun,alpha,weights,bootstrpOptions,varargin)
% normal approximation interval
% A.C. Davison and D.V. Hinkley (1996), p198-200
 
bstat = bootstrp(nboot,bootfun,varargin{:},'weights',weights,'Options',bootstrpOptions);

se = std(bstat,0,1);   % standard deviation estimate
bias = mean(bsxfun(@minus,bstat,stat),1);
za = norminv(alpha/2);   % normal confidence point
lower = stat - bias + se*za; % lower bound
upper = stat - bias - se*za;  % upper bound

% return
ci = [lower;upper];        
end   % bootnorm() 
 
%-------------------------------------------------------------------------
function [ci,bstat] = bootper(nboot,bootfun,alpha,weights,bootstrpOptions,varargin)
% percentile bootstrap CI
 
bstat = bootstrp(nboot,bootfun,varargin{:},'weights',weights,'Options',bootstrpOptions);

pct1 = 100*alpha/2;
pct2 = 100-pct1;
lower = prctile(bstat,pct1,1); 
upper = prctile(bstat,pct2,1);

% return
ci =[lower;upper];
end % bootper() 

%-------------------------------------------------------------------------
function [ci,bstat] = bootcper(stat,nboot,bootfun,alpha,weights,bootstrpOptions,varargin)
% corrected percentile bootstrap CI
% B. Efron (1982), "The jackknife, the bootstrap and other resampling
% plans", SIAM.
 
bstat = bootstrp(nboot,bootfun,varargin{:},'weights',weights,'Options',bootstrpOptions);

% stat is transformed to a normal random variable z0.
% z0 = invnormCDF[ECDF(stat)]
z_0 = fz0(bstat,stat);
z_alpha = norminv(alpha/2); % normal confidence point
 
% transform z0 back using the invECDF[normCDF(2z0-za)] and
% invECDF[normCDF(2z0+za)] 
pct1 = 100*normcdf(2*z_0-z_alpha); 
pct2 = 100*normcdf(2*z_0+z_alpha);

% inverse ECDF
m = numel(stat);
lower = zeros(1,m);
upper = zeros(1,m);
for i=1:m
    lower(i) = prctile(bstat(:,i),pct2(i),1);
    upper(i) = prctile(bstat(:,i),pct1(i),1);
end

% return
ci = [lower;upper];
end % bootcper() 
 
%-------------------------------------------------------------------------
function [ci,bstat] = bootbca(stat,nboot,bootfun,alpha,weights,bootstrpOptions,varargin)
% corrected and accelerated percentile bootstrap CI
% T.J. DiCiccio and B. Efron (1996), "Bootstrap Confidence Intervals",
% statistical science, 11(3)
 
bstat = bootstrp(nboot,bootfun,varargin{:},'weights',weights,'Options',bootstrpOptions);

% same as bootcper, this is the bias correction
z_0 = fz0(bstat,stat);

% apply jackknife
try
    jstat = jackknife(bootfun,varargin{:},'Options',bootstrpOptions);
catch ME
    m = message('stats:bootci:JackknifeFailed',func2str(bootfun));
    throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
end
N = size(jstat,1);
if isempty(weights)
    weights = repmat(1/N,N,1);
else
    weights = weights(:);
    weights = weights/sum(weights);
end

% acceleration finding, see DiCiccio and Efron (1996)
mjstat = sum(bsxfun(@times,jstat,weights),1); % mean along 1st dim.
score = bsxfun(@minus,mjstat,jstat); % score function at stat; ignore (N-1) factor because it cancels out in the skew
iszer = all(score==0,1);
skew = sum(bsxfun(@times,score.^3,weights),1) ./ ...
    (sum(bsxfun(@times,score.^2,weights),1).^1.5) /sqrt(N); % skewness of the score function
skew(iszer) = 0;
acc = skew/6;  % acceleration

% transform back with bias corrected and acceleration
z_alpha1 = norminv(alpha/2);
z_alpha2 = -z_alpha1;
pct1 = 100*normcdf(z_0 +(z_0+z_alpha1)./(1-acc.*(z_0+z_alpha1)));
pct1(z_0==Inf) = 100;
pct1(z_0==-Inf) = 0;
pct2 = 100*normcdf(z_0 +(z_0+z_alpha2)./(1-acc.*(z_0+z_alpha2)));
pct2(z_0==Inf) = 100;
pct2(z_0==-Inf) = 0;

% inverse of ECDF
m = numel(stat);
lower = zeros(1,m);
upper = zeros(1,m);
for i=1:m
    lower(i) = prctile(bstat(:,i),pct2(i),1);
    upper(i) = prctile(bstat(:,i),pct1(i),1);
end

% return
ci = sort([lower;upper],1);
end % bootbca()
 
%-------------------------------------------------------------------------
function [ci,bstat] = bootstud(stat,nboot,bootfun,alpha,nbootstd,stderrfun,weights,bootstrpOptions,varargin)
% studentized bootstrap CI with bootstrp to estimate the se
% T.J. DiCiccio and B. Efron (1996), "Bootstrap Confidence Intervals",
% statistical science, 11(3)

% Should we use stderrfun to compute st. dev. or an inner bootstrap loop?
if isempty(stderrfun) % studentized with bootstrap error
    if nbootstd<=0 || nbootstd~=round(nbootstd)
        error(message('stats:bootci:BadNbootstd'))
    end
else  % studentized with stderrfun fun
    % error check for stderrfun
    try
        out=stderrfun(varargin{:});
    catch ME
        m = message('stats:bootci:BadStderr',func2str(stderrfun));
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    end
    if any(~isfinite(out))
        error(message('stats:bootci:NonfiniteStderr'));
    end
    if isvector(out)
        out = out(:)';
    end
    if size(stat)~=size(out)
        error(message('stats:bootci:BadStderrSize'));
    end
end;

% bootstrap
[bstat,bootsam] = bootstrp(nboot,bootfun,varargin{:}, ...
    'weights',weights,'Options',bootstrpOptions); % bootstrap statistics
             
% find non-scalar data in varargin
la = length(varargin);
scalard = zeros(la,1);
for k = 1:la
   [row,col] = size(varargin{k});
   if max(row,col) == 1
      scalard(k) = 1;
   end
   if row == 1 && col ~= 1
      varargin{k} = varargin{k}(:);
   end
end

% bootstrap generated bootstrap replica to get student errors
N = size(bootsam,1);
sd_t = zeros(nboot,numel(stat));
for b=1:nboot
    db = cell(la,1);
    for k = 1:la
        % store the bootstrap data samples in a cell array
        if scalard(k) == 0
            db{k} = varargin{k}(bootsam(:,b),:);
        else
            db{k} = varargin{k};
        end
    end
    if isempty(weights)
        w = [];
    else
        w = weights(bootsam(:,b));
    end
    if isempty(stderrfun)
        bstatstd = bootstrp(nbootstd,bootfun,db{:},'weights',w,...
            'Options',bootstrpOptions);
        sd_t(b,:) = std(bstatstd,0,1);
    else
        bstatstd = stderrfun(db{:});
        sd_t(b,:) = reshape(bstatstd,1,numel(stat));
    end
end

% studentized statistics
above0 = sd_t>0;
mbstat = bsxfun(@minus,bstat,stat);
tstat = zeros(size(mbstat));
tstat(above0) = mbstat(above0)./sd_t(above0);

% percentiles for the studentized stats are computed.
lower = prctile(tstat,100*alpha/2,1);
upper = prctile(tstat,100*(1-alpha/2),1);

% back to the original stats from the studentized stats
lower = lower.*std(bstat,0,1) + stat;
upper = upper.*std(bstat,0,1) + stat;

% return
ci = [lower; upper];
end 
 
% -------------------------
function z0=fz0(bstat,stat)
% Compute bias-correction constant z0
z0 = norminv(mean(bsxfun(@lt,bstat,stat),1) + mean(bsxfun(@eq,bstat,stat),1)/2);
end   % fz0()
