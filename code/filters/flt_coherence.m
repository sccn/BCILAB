function signal = flt_coherence(varargin)
% Calculate between-channel / component coherence.
% Signal = flt_coherence(Signal, TimeBandwidth, Tapers, Padding, IncludePhase, Normalized, LogTransform)
%
% The fourier representation is practical if there is a highly complex (relevant) spectral structure
% in the data. In any representation except for 'complex', it is necessary that an appropriate
% spatial filter has been applied to the data beforehand (e.g., ICA with the 'transform' option set
% to 1, or the surface laplacian). The complex representation can be spatially transformed after
% running the Fourier filter, since it is a linear operator. If classifiers operate directly on the
% fourier representation, it is usually necessary to chose a representation that can be linearly
% mapped to a prediction, such as 'amplitude', 'polar' or 'phase'; 'complex' requires a non-linear
% classification, which is almost certain to overfit any other random non-linear dependencies in the
% data. Some interesting assumptions that can be imposed on Fourier data include sparsity (l1) and
% group sparsity (if independent component activity was transformed).
%
% Note: The computation time of coherence scales quadratically with the number of components under
%       consideration; 16 components is still fast enough for real-time use. Also note that coherence
%       applied to channels is not nearly as useful as applied to spatially filtered source signals.
%
% In:
%   Signal :   Epoched data set to be processed
%
%   TimeBandwidth : Spectral smoothing. Controls the bias vs. variance of the spectral estimation.
%                   (default: 3)
%
%   Tapers : Number of tapers. Should be an integer smaller than 2*TimeBandwith (default: 2*TimeBandwidth-1)
%
%   Padding : FFT padding factor. Controls the oversampling of the spectrum; 0 is the next largest
%             power of two, 1 is 2x as much, etc. (default: 0)
%
%   IncludePhase : Include phase information. Whether to include the phase of coherence in the result.
%                  (default: true)
%
%   Normalized  : Normalize the spectrum by 1/f. Doing this has benefits for classifiers that work
%                 best with naturally normalized features (e.g. some regularized classifiers).
%                 (default: true)
%
%   LogTransform : Log-Transform. Whether to transform the resulting spectral data into the log
%                  domain; can facilitates the use of simple (linear) classifiers. (default: false)
%
% Out: 
%   Signal  :   processed data set
%
% Examples:
%   % use default settings (basically contains all information, but not log-transformed PSDs)
%   eeg = flt_coherence(eeg)
%
%   % use a custom bandwith
%   eeg = flt_coherence(eeg,4)
%
%   % turn off the phase information
%   eeg = flt_coherence('Signal',eeg,'IncludePhase',false)
%
%   % turn on log-transformation of the PSD's (diagonals)
%   eeg = flt_coherence('Signal',eeg,'LogTransform',true)
%
%   % turn off the 1/f normalization (giving raw spectra)
%   eeg = flt_coherence('Signal',eeg,'Normalized',false)
%
% See also:
%   cohmatrixc
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% flt_coherence_version<0.95> -- for the cache

if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','CoherenceTransform', 'depends','set_makepos', 'follows',{'flt_reconstruct','flt_project','flt_window'}, 'independent_channels',false, 'independent_trials',true);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'mtbandwidth','TimeBandwidth'},3,[0 Inf],'Spectral smoothing. Controls the bias vs. variance of the spectral estimation.'), ...
    arg({'mttapers','Tapers'},[],uint32([1 1000]),'Number of tapers. Should be an integer smaller than 2*TimeBandwith; default 2*TimeBandwidth-1','guru',true), ...
    arg({'mtpadding','Padding'},0,[],'FFT padding factor. Controls the oversampling of the spectrum; 0 is the next largest power of two, 1 is 2x as much, etc.','guru',true), ...
    arg({'includephase','IncludePhase'}, true, [], 'Include phase information. Whether to include the phase of coherence in the result.'), ...
    arg({'normalized','Normalized'}, true, [], 'Normalize the spectrum by 1/f. Doing this has benefits for classifiers that work best with naturally normalized features (e.g. some regularized classifiers).'), ...
    arg({'logtransform','LogTransform'}, false, [], 'Log-Transform. Whether to transform the resulting spectral data into the log domain; can facilitates the use of simple (linear) classifiers.'));

    
[C,S,T] = size(signal.data);

% compute the Coherence & PSD using Slepian tapers
if isempty(mttapers)
    mttapers = 2*mtbandwidth-1; end
tmp = cell(1,size(signal.data,3));
for t = 1:length(tmp)
    % calc coherence magnitude and phase
    [coh,phi] = cohmatrixc(signal.data(:,:,t)',struct('tapers',[2*mtbandwidth mttapers],'pad',mtpadding,'Fs',signal.srate));
    % also calc the spectrum
    [spec,F] = mtspectrumc(signal.data(:,:,t)',struct('tapers',[2*mtbandwidth mttapers],'pad',mtpadding,'Fs',signal.srate));
    transp_coh = cell(1,size(signal.data,1));
    transp_phi = cell(1,size(signal.data,1));
    for c=1:length(transp_coh)
        if normalized
            coh(:,c,c) = spec(:,c) .* ([1 1:size(spec,1)-1]')/size(spec,1);
        else
            coh(:,c,c) = spec(:,c);
        end
        transp_coh{c} = coh(:,:,c)';
        transp_phi{c} = phi(:,:,c)';
    end
    if includephase
        tmp{t} = [vertcat(transp_coh{:}) ; vertcat(transp_phi{:})];
    else
        tmp{t} = vertcat(transp_coh{:});
    end
end
signal.data = cat(3,tmp{:});
signal.pnts = size(signal.data,2);
signal.freqs = F;

if includephase
    % replicate chanlocs
    idx = [1:signal.nbchan; 1:signal.nbchan];
    signal.chanlocs = signal.chanlocs(idx);
    for k=1:length(signal.chanlocs)
        if mod(k,2)
            signal.chanlocs(k).labels = [signal.chanlocs(k).labels '_mag'];
        else
            signal.chanlocs(k).labels = [signal.chanlocs(k).labels '_ang'];
        end
    end
    % recalc channel number
    signal.nbchan = size(signal.data,1);
end

if logtransform
    if includephase
        signal.data(1:2:end,:,:) = log(signal.data(1:2:end,:,:)); 
    else
        signal.data = log(signal.data); 
    end
    signal.data(~isfinite(signal.data(:))) = 0;
end

exp_endfun;
