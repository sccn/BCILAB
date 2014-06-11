function signal = flt_fourier_bandpower(varargin)
% Estimate bandpower for multiple frequencies using a filter bank.

% Compute logarithmic bandpower features.
% [Signal,State] = flt_fourier_bandpower(Signal, Bands)
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
% Out:
%   Signal       :  EEGLAB data set with annotated/modified fields:
%                   EEG.data         contains the bandpower for each band
%                                    and channel
%                   signal.spectrum: contains the complete power spectrum
%                   signal.freqs:    contains all frequencies used to
%                                    compute the spectrum
%%
% Examples:
%   % calculate bandpower features with two identical bands per channel
%   signal = flt_bandpower(eeg,{[10, 12], [14, 22]})
%
%   % calculate bandpower features with individual bands per channel
%   signal = flt_bandpower(eeg,{{[10, 12], [14, 22]}, {[7, 35]}, {[6, 10], [12, 16], [20, 30]}})
%
% References:
%
% Based on flt_bandpower() by Clemens Brunner
%
%                                Tim Mullen, Swartz Center for Computational Neuroscience, UCSD
%                                2013-04-28

if ~exp_beginfun('filter') return; end

% makes no sense on epoched data
declare_properties('name', 'FourierBandPower', 'cannot_follow', 'set_makepos', 'follows', {'flt_ica', 'flt_iir'},'independent_channels', false, 'independent_trials',true);


arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg_sub({'filterSpec','Filtering'},{'rep','multitaper'},@flt_fourier,'Filtering approach'), ...
    arg({'bands','Bands'}, {[10, 12], [14, 22]}, [], 'Frequency-band selection (Hz). Can be specified for all channels as {[low high], [low high], ...} or for each channel individually as {{[low high], [low high], ...} {[low high], [low high], ...}}. Frequencies in Hz.'), ...
    arg({'freqcol','FreqCollapse'},'sum',{'sum','integrate','mean','max'},'Method for collapsing over frequency within each band'), ...
    arg({'avgChans','AverageChannels'},false,[],'Average power across channels'), ...
    arg({'preserveSpec','PreserveSpectrum'},true,[],'Keep original spectrum. If true, spectrum and frequencies are stored in signal.spectrum and signal.freqs'), ...
    arg({'anotdata','AnnotateData'},false,[],'Annotate or replace signal.data') ...
    );

if size(signal.data,3) > 1
    error('flt_fourier_bandpower must be applied to continuous (non-epoched) data.');  end


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

% obtain spectrum for each channel
signal_fourier = exp_eval(flt_fourier('signal',signal,filterSpec,'arg_direct',true));

% loop over filters
if all_channels  % apply filters to all channels
    signal_filtered = cell(1, n_bands(1));  % will contain filtered signals, the bands are applied to all channels
    for n_filter = 1:n_bands(1)
        % extract power over frequency band
        fidx = hlp_getindex(signal_fourier.freqs,bands{1}{n_filter});
        signal_filtered{n_filter} = hlp_collapsefreqs(signal_fourier.data,freqcol,fidx(1):fidx(2));
    end
    
    % Merge filtered signals into one EEGLAB signal structure
    if avgChans
        % average bandpower over channels
        for n_filter=1:n_bands(1)
            signal_filtered{n_filter} = mean(signal_filtered{n_filter},1);
        end
        
        if ~anotdata
            temp = signal;
            temp.chanlocs = repmat(temp.chanlocs(1),1,n_bands(1));
            temp.pnts = 1;
            temp.data = zeros(n_bands(1), temp.pnts);
            temp.data = vertcat(signal_filtered{:});
            temp.nbchan = size(temp.data,1);
            % rename channel labels to bands
            temp.chanlocs = temp.chanlocs(1:n_bands(1));
            for k = 1:n_bands(1)
                temp.chanlocs(k).labels = ['Hz ' num2str(bands{1}{k}(1)) '-' num2str(bands{1}{k}(2))];
            end
            signal = temp;
        else
            signal.fourier_bandpower = vertcat(signal_filtered{:});
        end
        clear('temp');
    else
        if ~anotdata
            temp = signal;
            temp.pnts = 1;
            temp.data = zeros(sum(n_bands), temp.pnts);
            temp.data = vertcat(signal_filtered{:});
            temp.nbchan = sum(n_bands);
            temp.chanlocs = repmat(signal.chanlocs, 1, n_bands(1));
            % make channel names unique by appending '_X'
            temp.chanlocs = hlp_microcache('chanlocs',@make_chanlocs1,signal.chanlocs,n_bands,bands);
            signal = temp;
        else
            signal.fourier_bandpower = vertcat(signal_filtered{:});
        end
    end    
else  % apply filters to each channel individually
    signal_filtered = cell(1, length(bands));  % will contain filtered signals, each channel has separate bands
    for k = 1:length(signal_filtered)
        signal_filtered{k} = cell(1, n_bands(k));  end    
    for n_channel = 1:length(bands)
        for n_filter = 1:n_bands(n_channel)
            % extract power over frequency band
            fidx = hlp_getindex(signal_fourier.freqs,bands{n_channel}{n_filter});
            signal_filtered{n_channel}{n_filter} = hlp_collapsefreqs(signal_fourier.data,freqcol,fidx(1):fidx(2));
        end
    end;    
    % Merge filtered signals into one EEGLAB signal structure
    
    for k = 1:length(signal_filtered)
        temp_data{k} = vertcat(signal_filtered{k}{:}); end;
    
    if ~anotdata
        temp = signal;
        temp.pnts = 1;
        temp.data = vertcat(temp_data{:});
        temp.nbchan = sum(n_bands);
        temp.chanlocs = hlp_microcache('chanlocs',@make_chanlocs2,signal.chanlocs,length(signal_filtered),n_bands,bands);
        signal = temp;
    else
        signal.fourier_bandpower = vertcat(temp_data{:});
    end
end;

if preserveSpec
    signal.freqs    = signal_fourier.freqs;
    signal.spectrum = signal_fourier.data;
else
    signal.freqs = [];
    signal.spectrum = [];
end
exp_endfun;

function chanlocs = make_chanlocs1(chanlocs,n_bands,bands)
chanlocs = repmat(chanlocs, 1, n_bands(1));
counter = 1;
for k = 1:n_bands(1)
    for l = 1:length(bands)
        chanlocs(counter).labels = sprintf('%s_%i',chanlocs(counter).labels,k);
        counter = counter + 1;
    end
end

        
function chanlocs = make_chanlocs2(chanlocs,lsf,n_bands,bands)
for k = 1:lsf
    chanlocs{k} = repmat(chanlocs(k), n_bands(k), 1); end;
chanlocs = vertcat(chanlocs{:});
counter = 1;
for k = 1:length(bands)  % loop over channels
    for l = 1:n_bands(k)  % loop over bands per channel
        chanlocs(counter).labels = [chanlocs(counter).labels, '_', num2str(l)];
        counter = counter + 1;
    end
end
chanlocs = chanlocs';
