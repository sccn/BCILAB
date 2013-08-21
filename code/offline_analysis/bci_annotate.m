function data = bci_annotate(varargin)
% Annotate a raw dataset with a bci-derived channel
% Signal = bci_annotate(Model,Data,OutputFormat,SamplingRate,Interpolation,SmoothingWindow,StoreOutputs)
% 
% This function allows to add channels (or streams) to a continuous data set which
% contain the  estimates of a predictive model applied to the data. The resulting set can then be
% analyzed and visualized using a variety of tools (e.g. some of EEGLAB's plot tools).
%
% In:
%   Model : A predictive model (previously learned with bci_train).
%
%   Data : A continuous and not further filtered EEGLAB data set, or a stream bundle.
%
%   OutputFormat : Output format. The format in which each individual prediction is represented -- 
%                  can be a discrete/continuous probability distribution, the most likely value
%                  (mode), or the expected value (expectation). (default: 'distribution')
%
%   SamplingRate : Output sampling rate. The rate at which predictions are made, in Hz.
%                  (default: 10)
%
%   Interpolation : Output interpolation. Determines how the predictions are upsampled to the sampling 
%                   rate of the data set (or the first stream in a bundle). The most conservative
%                   approaches are ''constant'' and ''nans''. Resample is the only smooth output
%                   that is causal, but but can delay the BCI output significantly (esp. at low
%                   output sampling rates). The non-causal methods give a smooth output that is
%                   aligned with the time course of the predictions but exhibit "false" pre-ringing
%                   (and other effects) that cannot be physically realized for online execution.
%                   Among those, noncausallinear is simplest to understand and noncausalpchip
%                   results in a good-looking shape-preserving smoothing. (default: 'constant')
%
%   StoreOutputs : Output behavior. Determines whether the BCI outputs are attached to the data set or 
%                  stream as additional channels, or as an additional stream, or whether an EEGLAB
%                  set containing just the outputs should be returned. (default: 'newchannels')
%
% Out:
%   Data: the new data, with BCI-derived predictions (or a data set that contains only the BCI 
%         predictions)
%
% Notes:
%   if this function is applied to epoched data, there will be a period at the beginning of each
%   epoch for which there exist no meaningful predictions, because each model usually requires a 
%   minumum amount of samples to operate.
%
% Examples:
%   % given a predictive model and a continuous data set, append a channel which encodes the BCI 
%   % model's output
%   eeg = bci_annotate(lastmodel,eeg)
%
%   % as before, but use a different output format (namely: expected value instead of probability 
%   % distribution)
%   eeg = bci_annotate(lastmodel,eeg,'OutputFormat','expectation')
%
%   % as before, but this time use a higher sampling rate (slower, but more precise)
%   eeg = bci_annotate(lastmodel,eeg,'SamplingRate',60)
%
%   % as before, but use a smooth (but lagged) interpolation method
%   eeg = bci_annotate(lastmodel,eeg,'Interpolation','resample')
% 
% See also:
%   bci_predict, onl_simulate, bci_train
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% read arguments
opts = arg_define([0 2],varargin, ...
    arg({'model','Model'},mandatory,[],'Predictive model. This is a model as previously computed via bci_train.'), ...
    arg({'data','Data'},mandatory,[],'Data set. EEGLAB data set, or stream bundle, or cell array of data sets / stream bundles to use for prediction.'), ...
    arg({'format','OutputFormat'},'distribution',{'distribution','expectation','mode'},'Output format. The format in which each individual prediction is represented -- can be a discrete/continuous probability distribution, the most likely value (mode), or the expected value (expectation).'), ...
    arg({'srate','SamplingRate'},10,[],'Output sampling rate. The rate at which predictions are made, in Hz.'), ...
    arg({'interpolation','Interpolation'},'constant',{'nans','constant','laggedresample','noncausalpchip','noncausallinear','noncausalresample'},'Output interpolation. Determines how the predictions are upsampled to the sampling rate of the data set (or the first stream in a bundle). The most conservative approaches are ''constant'' and ''nans''. Resample produces a smoothed output, but delays the BCI output significantly (esp. at low output sampling rates). The non-causal methods give a smooth output that is aligned with the time course of the predictions but exhibit "false" pre-ringing (etc) that cannot be physically realized for online execution. Among those, linear is simplest to understand and pchip results in a good shape-preserving smoothing.'), ...
    arg({'storeoutputs','StoreOutputs'},'newchannels',{'newchannels','newstream','newset'},'Output behavior. Determines whether the BCI outputs are attached to the data set / stream as additional channels, or as an additional stream, or whether an EEGLAB set containing just the outputs should be returned.'), ...
    arg({'name','ChannelName'},'BCI',[],'Name of new channels. If there are multiple channels, they will be consecutively numbered.'));

% uniformize the data
data = opts.data;
was_simple = all(isfield(data,{'data','srate'})) || all(isfield(data,{'head','parts'}));  % whether the data was a simple EEGLAB data set
if iscell(data)
    data = struct('streams',{data}); end
if was_simple
    data = struct('streams',{{data}}); end
data = utl_check_bundle(data);
for s=1:length(data.streams)
    data.streams{s} = exp_eval_optimized(data.streams{s}); end

% use onl_simulate to get outputs, for each epoch
[preds,lats] = onl_simulate(data,opts.model,'sampling_rate',opts.srate,'format',opts.format);
preds(~isfinite(preds(:))) = 0;
lats = 1 + round(lats*data.streams{1}.srate);

[C,S,E] = size(data.streams{1}.data);   % #channels & #samples of the first stream
NC = size(preds,2);                     % # of newly-added channels

% interpolate the data at the signal's sampling rate
switch opts.interpolation
    case 'nans'
        % generate NaN dataset
        tmp = exp_eval(set_new('data',nan(NC,S),'srate',data.streams{1}.srate));
        % fill in the predictions at the appropriate locations
        tmp.data(:,lats) = preds';
    case 'constant'
            % generate NaN dataset
            tmp = exp_eval(set_new('data',nan(NC,S),'srate',data.streams{1}.srate));
            % fill in the predictions at the appropriate locations
            tmp.data(:,lats) = preds';
            % go replace NaNs by the last valid output
            last = tmp.data(:,1);
            for k=2:size(tmp.data,2)
                cur = tmp.data(:,k);
                if all(isnan(cur))
                    tmp.data(:,k) = last;
                else
                    last = cur;
                end
            end
    case {'laggedresample','noncausalresample'}        
        if strcmp(opts.interpolation,'laggedresample')
            % causal version
            tmp = exp_eval(flt_resample(exp_eval(set_new('data',preds','srate',opts.srate)),data.streams{1}.srate));
        else
            % non-causal version
            tmp = pop_resample(exp_eval(set_new('data',preds','srate',opts.srate)),data.streams{1}.srate);
        end
        % fix exact sample range to align BCI outputs with the rest of the data
        if lats(1) > 1
            % prepend constant data, if necessary
            tmp.data = [repmat(tmp.data(:,1),1,lats(1)-1) tmp.data];
            tmp.pnts = size(tmp.data,2);
        end
        if tmp.pnts < S
            % pad with constant data, if necessary
            tmp.data = [tmp.data repmat(tmp.data(:,end),1,S - tmp.pnts)];
            tmp.pnts = size(tmp.data,2);
        elseif tmp.pnts > S
            % or cut samples from the end, if necessary
            tmp.data = tmp.data(:,1:S);
            tmp.pnts = size(tmp.data,2);
        end
    otherwise
        % one of the other non-causal modes
        tmp = zeros(NC,S);
        for d=1:NC
            % interpolate using interp1
            tmp(d,:) = interp1(lats,preds(:,d),1:S,opts.interpolation(10:end))';
            tmp(d,[1:min(lats) max(lats):end]) = 0;
        end
        % turn into a dataset
        tmp = exp_eval(set_new('data',tmp,'srate',data.streams{1}.srate));
end

% add chanlocs to tmp
if NC > 1
    for k=1:NC
        tmp.chanlocs(k).labels = [opts.name num2str(k)]; 
        tmp.chanlocs(k).type = 'Latent';
    end
else
    tmp.chanlocs.labels = opts.name;
    tmp.chanlocs.type = 'Latent';
end

switch opts.storeoutputs
    case 'newset'
        data = tmp;
    case 'newchannels'
        % append the results as new channel(s) to the data set
        data.streams{1}.data(C+(1:NC),:) = tmp.data;
        % update chanlocs, nbchan, ...
        data.streams{1}.nbchan = size(data.streams{1}.data,1);
        for k=1:NC
            data.streams{1}.chanlocs(k+C).labels = tmp.chanlocs(k).labels; 
            data.streams{1}.chanlocs(k+C).type = tmp.chanlocs(k).type;             
        end
        if was_simple
            data = data.streams{1}; end
    case 'newstream'
        data.streams{end+1} = tmp;
    otherwise
        error('Unsupported output mode selected.');
end
