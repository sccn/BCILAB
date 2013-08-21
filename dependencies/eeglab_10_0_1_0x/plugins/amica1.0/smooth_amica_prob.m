% Applies a Hanning window of specified size in sec, to the log-likelihood (LLt) resulting in llt_smooth and then
% converts it to probability of models v_smooth.
% Inputs:
%   sRate          - Sampling rate of data (EEG.srate)
%   LLt            - Log-likelihood of data
%   smoothlength   - Length of smoothing window in sec. 

function [v_smooth llt_smooth] = smooth_amica_prob(sRate,LLt, smoothlength)

if nargin<2
    smoothlength = 2; % default smoothing length (in sec.);
end

% Hanning smoothing window generation
smoothwnd = hann(round(smoothlength*sRate)); smoothwnd = smoothwnd/sum(smoothwnd);

if size(LLt,3) == 1
    
    llt_smooth = filtfilt_fast(fastif(smoothlength == 0, 1,smoothwnd),1,LLt')';
    v_smooth = LLt2v(llt_smooth);
    
else
    llt = reshape(LLt,size(LLt,1),size(LLt,2)*size(LLt,3));
    llt_smooth = filtfilt_fast(fastif(smoothlength == 0, 1,smoothwnd),1,llt')';
    v_smooth = LLt2v(llt_smooth);
    llt_smooth = reshape(llt_smooth,size(LLt,1),size(LLt,2),size(LLt,3));
    v_smooth = reshape(v_smooth,size(LLt,1),size(LLt,2),size(LLt,3));
end