function [signal,state] = flt_bandpower(varargin)
% Compute logarithmic bandpower features.
% [Signal,State] = flt_bandpower(Signal, Bands, Smoothing, State)
%
% TODO: detailed description
%
% In:
%   Signal       :   continuous data set to be filtered
%
%   Bands        :   bands specification:
%                    * if all channels have the same bands, use a cell array containing the
%                      frequency bands, e.g. {[8, 10], [12, 16], [22, 30]} (in Hz)
%                      this example creates 3 identical bands for each channel
%                    * if you want individual bands for each channel, use a
%                      two-dimensional cell array containing the frequency bands
%                      (second dimension) of each channel (first dimension),
%                      e.g. {{[7, 11], [13, 18]}, {[6, 35]}, {[12, 15], [20, 22], [24, 35]}} (in Hz)
%                      this example creates 2 bands for the first channel, 1
%                      band for the second channel, and 3 bands for the third
%                      channel
%                    * default: {[10, 12], [16, 24]}
%
%   Smoothing    :   length of smoothing windows (in s) (default: 1)
%
%   Order        :   filter order of IIR Butterworth filter (default: 'auto')
%
%   State        :   previous filter state, as obtained by a previous execution of flt_bandpower on an
%                    immediately preceding data set (default: [])
%
% Out:
%   Signal       :  filtered, continuous EEGLAB data set
%
%   State        :  state of the filter, can be used to continue on a subsequent portion of the data
%
% Examples:
%   % calculate bandpower features with two identical bands per channel
%   signal = flt_bandpower(eeg,{[10, 12], [14, 22]})
%
%   % calculate bandpower features with individual bands per channel
%   signal = flt_bandpower(eeg,{{[10, 12], [14, 22]}, {[7, 35]}, {[6, 10], [12, 16], [20, 30]}})
%
% References:
%  [1] TODO
%
%                                Clemens Brunner, Swartz Center for Computational Neuroscience, UCSD
%                                2011-06-29

% flt_bandpower_version<0.95> -- for the cache

if ~exp_beginfun('filter') return; end

% makes no sense on epoched data
declare_properties('name', 'BandPower', 'experimental',true, 'cannot_follow', 'set_makepos', 'follows', 'flt_ica', 'precedes', 'flt_iir', 'independent_channels', false, 'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'bands','Bands'}, {[10,12],[16,24]}, [], 'Frequency bands for each channel (in Hz).'), ...  % TODO: Modify description
    arg({'smoothing','SmoothingWindow'}, 1, [], 'Length of smoothing window (in s).'), ...
    arg({'logpower','LogPower'},true,[],'Return logarithm of power (dB)'), ...
    arg({'order', 'FilterOrder'}, 'auto', [], 'Filter order of IIR Butterworth filter.'), ...
    arg_norep({'state','State'},unassigned));

if size(signal.data,3) > 1
    error('flt_bandpower is supposed to be applied to continuous (non-epoched) data.');  end

% parse frequency bands parameter
% the goal is to reshape the bands parameter to the most generic form
% bands{n}{m}
%   n is the channel index (i = 1...signal.nbchan)
%   m is the band index for channel i (m = 1...n_bands_i)
%
% however, a simpler version exists when bands are identical for all channels.
% then, this parameter reduces to bands{m}, but the following code converts the
% format to the generic one described above.
if ~iscell(bands)
    error('bands must be a cell array.');
elseif ~iscell(bands{1})  % these frequency bands apply to all channels
    for k = 1:signal.nbchan
        temp{k} = bands;  end
    bands = temp;
    all_channels = true;  % filters applied to all channels
    clear('temp');
else  % there are different frequency bands for each channel in a 2D cell array
    all_channels = false;  % individual filters for each channel
end

for k = 1:length(bands)
    n_bands(k) = length(bands{k});  end

% loop over filters
if all_channels  % apply filters to all channels
    signal_filtered = cell(1, n_bands(1));  % will contain filtered signals, the bands are applied to all channels
    for n_filter = 1:n_bands(1)
        
        % need to create dfilt state object?
        if ~exist('state','var') || length(state) < n_filter
            if ~exist('dfilt','file')
                error('You need the Signal Processing toolbox to make use of IIR filters in BCILAB.'); end
            
            % create 1 Hz transition bands for filter order calculation
            [Wp,Ws,label] = deal(2 * bands{1}{n_filter}/signal.srate, 2 * [bands{1}{n_filter}(1)-1, bands{1}{n_filter}(2)+1]/signal.srate, {});
            
            if ischar(order) && strcmp(order, 'auto')
                [n,Wn] = buttord(Wp,Ws,0.5,50);
            else
                n = order;
                Wn = bands{1}{n_filter} ./ (signal.srate / 2);
            end
            
            % compute filter coefficients (in Zero-Pole-Gain form, to prevent instability)
            [z,p,k] = butter(n,Wn,label{:});
            
            [sos,g] = zp2sos(z,p,k);
            state{n_filter} = dfilt.df2sos(sos,g);
            state{n_filter} = dfilt.df2sos(sos,g/max(abs(freqz(state{n_filter}))));
            set(state{n_filter},'PersistentMemory',true);
        else
            % make a deep copy of the state
            state{n_filter} = copy(state{n_filter});
        end

        % apply the filter, smooth with moving average, and take the logarithm
        signal_filtered{n_filter} = filter(ones(round(smoothing*signal.srate),1),round(smoothing*signal.srate),transpose(filter(state{n_filter},double(signal.data),2).^2));
        if logpower
            signal_filtered{n_filter} = transpose(log(signal_filtered{n_filter}));
        else
            signal_filtered{n_filter} = transpose(signal_filtered{n_filter});
        end;
    end
    
    % Merge filtered signals into one EEGLAB signal structure
    temp = signal;
    temp.data = zeros(sum(n_bands), signal.pnts);
    temp.data = vertcat(signal_filtered{:});
    temp.nbchan = sum(n_bands);
    temp.chanlocs = repmat(signal.chanlocs, 1, n_bands(1));
    % make channel names unique by appending '_X'
    counter = 1;
    for k = 1:n_bands(1)
        for l = 1:length(bands)
            temp.chanlocs(counter).labels = [temp.chanlocs(counter).labels, '_', num2str(k)];
            counter = counter + 1;
        end
    end
    signal = temp;
    clear('temp');
    
else  % apply filters to each channel individually
    signal_filtered = cell(1, length(bands));  % will contain filtered signals, each channel has separate bands
    for k = 1:length(signal_filtered)
        signal_filtered{k} = cell(1, n_bands(k));  end
    
    % create empty state variable
    state = signal_filtered;
    
    for n_channel = 1:length(bands)
        for n_filter = 1:n_bands(n_channel)
            
            % need to create dfilt state object?
            if isempty(state{n_channel}{n_filter})
                if ~exist('dfilt','file')
                    error('You need the Signal Processing toolbox to make use of the bandpower filter in BCILAB.'); end
                
                % create 1 Hz transition bands for filter order calculation
                [Wp,Ws,label] = deal(2 * bands{n_channel}{n_filter}/signal.srate, 2 * [bands{n_channel}{n_filter}(1)-1, bands{n_channel}{n_filter}(2)+1]/signal.srate, {});
                
                if ischar(order) && strcmp(order, 'auto')
                    [n,Wn] = buttord(Wp,Ws,0.5,50);
                else
                    n = order;
                    Wn = bands{n_channel}{n_filter} ./ (signal.srate / 2);
                end
                
                % compute filter coefficients (in Zero-Pole-Gain form, to prevent instability)
                [z,p,k] = butter(n,Wn,label{:});
                
                [sos,g] = zp2sos(z,p,k);
                state{n_channel}{n_filter} = dfilt.df2sos(sos,g);
                state{n_channel}{n_filter} = dfilt.df2sos(sos,g/max(abs(freqz(state{n_channel}{n_filter}))));
                set(state{n_channel}{n_filter},'PersistentMemory',true);
            else
                % make a deep copy of the state
                state{n_channel}{n_filter} = copy(state{n_channel}{n_filter});
            end;         
            
            % apply the filter, smooth with moving average, and take the logarithm
            signal_filtered{n_channel}{n_filter} = filter(ones(round(smoothing*signal.srate),1),round(smoothing*signal.srate),transpose(filter(state{n_channel}{n_filter},double(signal.data),2).^2));
            if logpower
                signal_filtered{n_channel}{n_filter} = transpose(log(signal_filtered{n_channel}{n_filter}));
            else
                signal_filtered{n_channel}{n_filter} = transpose(signal_filtered{n_channel}{n_filter});
            end;
        end;   
        
    end;
    
    % Merge filtered signals into one EEGLAB signal structure
    temp = signal;
    for k = 1:length(signal_filtered)
        temp_data{k} = vertcat(signal_filtered{k}{:}); end;
    temp.data = vertcat(temp_data{:});
    clear('temp_data');
    temp.nbchan = sum(n_bands);
    for k = 1:length(signal_filtered)
        temp_chanlocs{k} = repmat(signal.chanlocs(k), n_bands(k), 1); end;
    temp.chanlocs = vertcat(temp_chanlocs{:});
    % make channel names unique by appending '_X'
    counter = 1;
    for k = 1:length(bands)  % loop over channels
        for l = 1:n_bands(k)  % loop over bands per channel
            temp.chanlocs(counter).labels = [temp.chanlocs(counter).labels, '_', num2str(l)];
            counter = counter + 1;
        end
    end
    temp.chanlocs = temp.chanlocs';
    clear('temp_chanlocs');
    signal = temp;
    clear('temp');
    
end;


exp_endfun;
