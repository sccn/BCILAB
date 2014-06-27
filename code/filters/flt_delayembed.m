function signal = flt_delayembed(varargin)
% Apply delay embedding to epoched data.
% Signal = flt_delayembed(Signal,NumLags)
%
% Delay embedding is essentially appending to each multi-channel samples the subsequent k
% multi-channel samples, and thereby multiplies the number of channels by k. Delay embedding is a
% practical tool to extend linear spatial models (e.g., spatial filters or independent components)
% to linear spatio-temporal (and therefore implicitly spatio-spectral) models, just by applying
% those models to delay-embedded data. As a result, approaches that can learn optimal spatial
% filters (finding sources of interest) can be repurposed to learning optimal spatio-spectral
% filters (jointly finding sources and frequencies of interest).
%
% The tradeoff associated with delay-embedding is that the complexity of the models increases and
% they become harder to estimate, which might require more data or better
% constraints, priors, or regularization.
%
% In:
%   Signal : Epoched data set to be processed
%   
%   NumLags : the number of lags that shall be used for delay-embedding (default: 1)
%
%   IncludeIntermediates : Include intermediate lags. If this is set to false, only the 0''th and
%                          the N''th lag will be embedded. (default: true)
%
% Out:
%   Signal : the processed signal; will have more channels
%
% Notes:
%   The temporal filters that can be designed for a small number of lags are often limited to 
%   high-frequency responses; the frequency range can often be extended without increasing model
%   complexity by first resampling the data to the lowest acceptable sampling rate (e.g., 60 Hz).
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-11-17

% flt_delayembed_version<1.00> -- for the cache

if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','DelayEmbedding', 'depends','set_makepos', 'follows',{'flt_project','flt_window'}, 'independent_channels',false, 'independent_trials',true);

% declare arguments
arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'numlags','NumLags'}, 1, uint32([1 1 20 1000]), 'Number of lags. For delay-embedding.'), ...
    arg({'includeIntermediates','IncludeIntermediates'}, true, [], 'Include intermediate lags. If this is set to false, only the 0''th and the N''th lag will be embedded.'));

for k=quickif(includeIntermediates,0:numlags,[0 numlags])
    tmp{k+1} = signal.data(:,k+(1:end-numlags),:); end
signal.data = cat(1,tmp{:});
[signal.nbchan,signal.pnts,signal.trials] = size(signal.data);
signal.chanlocs = struct('labels', cellfun(@num2str,num2cell(1:signal.nbchan),'UniformOutput',false));

exp_endfun;
