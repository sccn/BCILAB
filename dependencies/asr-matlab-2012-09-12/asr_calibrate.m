function state = asr_calibrate(X,srate,cutoff,blocksize,B,A,window_len,window_overlap,max_dropout_fraction,min_clean_fraction)
% Calibration function for the Artifact Subspace Reconstruction (ASR) method.
% State = asr_calibrate(Data,SamplingRate,Cutoff,BlockSize,FilterB,FilterA,WindowLength,WindowOverlap,MaxDropoutFraction,MinCleanFraction)
%
% The input to this data is a multi-channel time series of calibration data. In typical uses the
% calibration data is clean resting EEG data of ca. 1 minute duration (can also be longer). One can
% also use on-task data if the fraction of artifact content is below the breakdown point of the
% robust statistics used for estimation (50% theoretical, ~30% practical). If the data has a
% proportion of more than 30-50% artifacts then bad time windows should be removed beforehand. This
% data is used to estimate the thresholds that are used by the ASR processing function to identify
% and remove artifact components.
%
% The calibration data must have been recorded for the same cap design from which data for cleanup
% will be recorded, and ideally should be from the same session and same subject, but it is possible
% to reuse the calibration data from a previous session and montage to the extent that the cap is
% placed in the same location (where loss in accuracy is more or less proportional to the mismatch
% in cap placement).
%
% The calibration data should have been high-pass filtered (for example at 0.5Hz or 1Hz using a
% Butterworth IIR filter).
%
% In:
%   Data : Calibration data [#channels x #samples]; *zero-mean* (e.g., high-pass filtered) and
%          reasonably clean EEG of not much less than 30 seconds length (this method is typically
%          used with 1 minute or more).
%
%   SamplingRate : Sampling rate of the data, in Hz.
%
%
%   The following are optional parameters (the key parameter of the method is the RejectionCutoff):
%
%   RejectionCutoff: Standard deviation cutoff for rejection. Data portions whose variance is larger
%                    than this threshold relative to the calibration data are considered missing
%                    data and will be removed. The most aggressive value that can be used without
%                    losing too much EEG is 2.5. A quite conservative value would be 5. Default: 5.
%
%   Blocksize : Block size for calculating the robust data covariance and thresholds, in samples;
%               allows to reduce the memory and time requirements of the robust estimators by this 
%               factor (down to Channels x Channels x Samples x 16 / Blocksize bytes). Default: 10
%
%   FilterB, FilterA : Coefficients of an IIR filter that is used to shape the spectrum of the signal
%                      when calculating artifact statistics. The output signal does not go through
%                      this filter. This is an optional way to tune the sensitivity of the algorithm
%                      to each frequency component of the signal. The default filter is less
%                      sensitive at alpha and beta frequencies and more sensitive at delta (blinks)
%                      and gamma (muscle) frequencies. Default: 
%                      [b,a] = yulewalk(8,[[0 2 3 13 16 40 min(80,srate/2-1)]*2/srate 1],[3 0.75 0.33 0.33 1 1 3 3]);
%
%   WindowLength : Window length that is used to check the data for artifact content. This is 
%                  ideally as long as the expected time scale of the artifacts but short enough to 
%				   allow for several 1000 windows to compute statistics over. Default: 0.5.
%
%   WindowOverlap : Window overlap fraction. The fraction of two successive windows that overlaps.
%                   Higher overlap ensures that fewer artifact portions are going to be missed (but
%                   is slower). Default: 0.66
%
%   MaxDropoutFraction : Maximum fraction of windows that can be subject to signal dropouts 
%                        (e.g., sensor unplugged), used for threshold estimation. Default: 0.1
%
%   MinCleanFraction : Minimum fraction of windows that need to be clean, used for threshold
%                      estimation. Default: 0.25
%
%
% Out:
%   State : initial state struct for asr_process
%
% Notes:
%   This can run on a GPU with large memory and good double-precision performance for faster processing 
%   (e.g., on an NVIDIA GTX Titan or K20), but requires that the Parallel Computing toolbox is
%   installed.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-08-31

% asr_calibrate_version<1.03> -- for the cache

% UC Copyright Notice
% This software is Copyright (C) 2013 The Regents of the University of California. All Rights Reserved.
% 
% Permission to copy, modify, and distribute this software and its documentation for educational,
% research and non-profit purposes, without fee, and without a written agreement is hereby granted,
% provided that the above copyright notice, this paragraph and the following three paragraphs appear
% in all copies.
% 
% Permission to make commercial use of this software may be obtained by contacting:
% Technology Transfer Office
% 9500 Gilman Drive, Mail Code 0910
% University of California
% La Jolla, CA 92093-0910
% (858) 534-5815
% invent@ucsd.edu 
% 
% This software program and documentation are copyrighted by The Regents of the University of
% California. The software program and documentation are supplied "as is", without any accompanying
% services from The Regents. The Regents does not warrant that the operation of the program will be
% uninterrupted or error-free. The end-user understands that the program was developed for research
% purposes and is advised not to rely exclusively on the program for any reason.
% 
% IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
% THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
% CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
% MODIFICATIONS.

[C,S] = size(X);

if nargin < 3 || isempty(cutoff)
    cutoff = 5; end
if nargin < 4 || isempty(blocksize)
    blocksize = 10; end
blocksize = max(blocksize,ceil((C*C*S*8*3*2)/hlp_memfree));
if nargin < 6 || isempty(A) || isempty(B)
    try
        % try to use yulewalk to design the filter (Signal Processing toolbox required)
        [B,A] = yulewalk(8,[[0 2 3 13 16 40 min(80,srate/2-1)]*2/srate 1],[3 0.75 0.33 0.33 1 1 3 3]);
    catch e %#ok<NASGU>
        % yulewalk not available (maybe no toolbox installed) -- use precomputed filter
        % coefficients depending on sampling rate
        switch srate
            case 100
                [B,A] = deal([0.9314233528641650 -1.0023683814963549 -0.4125359862018213  0.7631567476327510  0.4160430392910331 -0.6549131038692215 -0.0372583518046807  0.1916268458752655  0.0462411971592346],[1.0000000000000000 -0.4544220180303844 -1.0007038682936749  0.5374925521337940  0.4905013360991340 -0.4861062879351137 -0.1995986490699414  0.1830048420730026  0.0457678549234644]);
            case 128
                [B,A] = deal([1.1027301639165037 -2.0025621813611867  0.8942119516481342  0.1549979524226999  0.0192366904488084  0.1782897770278735 -0.5280306696498717  0.2913540603407520 -0.0262209802526358],[1.0000000000000000 -1.1042042046423233 -0.3319558528606542  0.5802946221107337 -0.0010360013915635  0.0382167091925086 -0.2609928034425362  0.0298719057761086  0.0935044692959187]);
            case 200
                [B,A] = deal([1.4489483325802353 -2.6692514764802775  2.0813970620731115 -0.9736678877049534  0.1054605060352928 -0.1889101692314626  0.6111331636592364 -0.3616483013075088  0.1834313060776763],[1.0000000000000000 -0.9913236099393967  0.3159563145469344 -0.0708347481677557 -0.0558793822071149 -0.2539619026478943  0.2473056615251193 -0.0420478437473110  0.0077455718334464]);
            case 256
                [B,A] = deal([1.7587013141770287 -4.3267624394458641  5.7999880031015953 -6.2396625463547508  5.3768079046882207 -3.7938218893374835  2.1649108095226470 -0.8591392569863763  0.2569361125627988],[1.0000000000000000 -1.7008039639301735  1.9232830391058724 -2.0826929726929797  1.5982638742557307 -1.0735854183930011  0.5679719225652651 -0.1886181499768189  0.0572954115997261]);
            case 300
                [B,A] = deal([1.9153920676433143  -5.7748421104926795   9.1864764859103936 -10.7350356619363630   9.6423672437729007  -6.6181939699544277   3.4219421494177711  -1.2622976569994351   0.2968423019363821],[1.0000000000000000 -2.3143703322055491  3.2222567327379434 -3.6030527704320621  2.9645154844073698 -1.8842615840684735  0.9222455868758080 -0.3103251703648485  0.0634586449896364]);
            case 500
                [B,A] = deal([2.3133520086975823 -11.9471223009159130  29.1067166493384340 -43.7550171007238190  44.3385767452216370 -30.9965523846388000  14.6209883020737190  -4.2743412400311449   0.5982553583777899],[1.0000000000000000  -4.6893329084452580  10.5989986701080210 -14.9691518101365230  14.3320358399731820  -9.4924317069169977   4.2425899618982656  -1.1715600975178280   0.1538048427717476]);
            case 512
                [B,A] = deal([2.3275475636130865 -12.2166478485960430  30.1632789058248850 -45.8009842020820410  46.7261263011068880 -32.7796858196767220  15.4623349612560630  -4.5019779685307473   0.6242733481676324],[1.0000000000000000  -4.7827378944258703  10.9780696236622980 -15.6795187888195360  15.1281978667576310 -10.0632079834518220   4.5014690636505614  -1.2394100873286753   0.1614727510688058]);
            otherwise
                error('repair_bursts:NoYulewalk','The yulewalk() function was not found and there is no pre-computed spectral filter for your sampling rate. If you would like to use the default spectral filter please try to resample to one of the supported rates (100,128,200,256,300,500,512) or get the appropriate toobox license (you can also disable the spectral weighting feature or supply your own precalculated IIR filter coefficients).');
        end
    end
end
if nargin < 8 || isempty(window_len)
    window_len = 0.5; end
if nargin < 9 || isempty(window_overlap)
    window_overlap = 0.66; end
if nargin < 10 || isempty(max_dropout_fraction)
    max_dropout_fraction = 0.1; end
if nargin < 11 || isempty(min_clean_fraction)
    min_clean_fraction = 0.25; end

X(~isfinite(X(:))) = 0;

% apply the signal shaping filter and initialize the IIR filter state
[X,iirstate] = filter(B,A,double(X),[],2); X = X';
if any(~isfinite(X(:)))
    error('The IIR filter diverged on your data. Please try using either a more conservative filter or removing some bad sections/channels from the calibration data.'); end

% calculate the sample covariance matrices U (averaged in blocks of blocksize successive samples)
U = zeros(length(1:blocksize:S),C*C);
for k=1:blocksize
    range = min(S,k:blocksize:(S+k-1));
    U = U + reshape(bsxfun(@times,reshape(X(range,:),[],1,C),reshape(X(range,:),[],C,1)),size(U));
end

% get the mixing matrix M
M = sqrtm(real(reshape(block_geometric_median(U/blocksize),C,C)));

% window length for calculating thresholds
N = round(window_len*srate);

% get the threshold matrix T
fprintf('Determining per-component thresholds...');
[V,D] = eig(M); %#ok<NASGU>
X = abs(X*V);
for c = C:-1:1
    % compute RMS amplitude for each window...
    rms = X(:,c).^2;
    rms = sqrt(sum(rms(bsxfun(@plus,round(1:N*(1-window_overlap):S-N),(0:N-1)')))/N);
    % fit a distribution to the clean part
    [mu(c),sig(c)] = fit_eeg_distribution(rms,min_clean_fraction,max_dropout_fraction);
end
T = diag(mu + cutoff*sig)*V';
disp('done.');

% initialize the remaining filter state
state = struct('M',M,'T',T,'B',B,'A',A,'cov',[],'carry',[],'iir',iirstate,'last_R',[],'last_trivial',true);



function y = block_geometric_median(X,blocksize,varargin)
% Calculate a blockwise geometric median (faster and less memory-intensive 
% than the regular geom_median function).
%
% This statistic is not robust to artifacts that persist over a duration that
% is significantly shorter than the blocksize.
%
% In:
%   X : the data (#observations x #variables)
%   blocksize : the number of successive samples over which a regular mean 
%               should be taken
%   tol : tolerance (default: 1.e-5)
%   y : initial value (default: median(X))
%   max_iter : max number of iterations (default: 500)
%
% Out:
%   g : geometric median over X
%
% Notes:
%   This function is noticably faster if the length of the data is divisible by the block size.
%   Uses the GPU if available.
% 

if nargin < 2 || isempty(blocksize)
    blocksize = 1; end

if blocksize > 1
    [o,v] = size(X);       % #observations & #variables
    r = mod(o,blocksize);  % #rest in last block
    b = (o-r)/blocksize;   % #blocks
    if r > 0
        X = [reshape(sum(reshape(X(1:(o-r),:),blocksize,b*v)),b,v); sum(X((o-r+1):end,:))*(blocksize/r)];
    else
        X = reshape(sum(reshape(X,blocksize,b*v)),b,v);
    end
end

try
    y = gather(geometric_median(gpuArray(X),varargin{:}))/blocksize;
catch
    y = geometric_median(X,varargin{:})/blocksize;
end



function y = geometric_median(X,tol,y,max_iter)
% Calculate the geometric median for a set of observations (mean under a Laplacian noise distribution)
% This is using Weiszfeld's algorithm.
%
% In:
%   X : the data, as in mean
%   tol : tolerance (default: 1.e-5)
%   y : initial value (default: median(X))
%   max_iter : max number of iterations (default: 500)
%
% Out:
%   g : geometric median over X

if ~exist('tol','var') || isempty(tol)
    tol = 1.e-5; end
if ~exist('y','var') || isempty(y)
    y = median(X); end
if ~exist('max_iter','var') || isempty(max_iter)
    max_iter = 500; end

for i=1:max_iter
    invnorms = 1./sqrt(sum(bsxfun(@minus,X,y).^2,2));
    [y,oldy] = deal(sum(bsxfun(@times,X,invnorms)) / sum(invnorms),y);
    if norm(y-oldy)/norm(y) < tol
        break; end
end



function result = hlp_memfree
% Get the amount of free physical memory, in bytes
result = java.lang.management.ManagementFactory.getOperatingSystemMXBean().getFreePhysicalMemorySize();
