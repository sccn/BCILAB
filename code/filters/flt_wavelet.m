function signal = flt_wavelet(varargin)
% Transform an epoched data set into a per-channel (multi-level) wavelet representation.
% Signal = flt_wavelet(Signal, Family, Levels)
%
% The wavelet representation allows to extract local variations in a signal (both in time and 
% frequency). The captured frequency range depends on the number of decomposition levels, and the
% shape of captured local variations (e.g. symmetric, smooth, ...) depends on the chosen wavelet 
% family. Some signals are sparse in a wavelet representation, which allows to use sparse classifiers
% on them.
%
% In:
%   Signal : epoched data set to be processed
%
%   Family : wavelet familty (default: 'db10')
%
%   Levels : decomposition levels (default: 10)
%
%   PacketDecomposition : do a wavelet packet decomposition instead of a multi-level wavelet 
%                         decomposition (default: false)
%
% Out: 
%   Signal : processed data set
%
% TODO/Notes:
%   This filter is currently too slow for online use (except on very few channels); needs to be optimized.
%
% Examples:
%   % use default settings (using Daubechies-10 wavelets)
%   eeg = flt_wavelet(eeg)
%
%   % use a different wavelet family
%   eeg = flt_wavelet(eeg,'sym4')
%
%   % do a multi-level decomposition (5 levels)
%   eeg = flt_wavelet(eeg,'sym4',5)
%
%   % as before, but pass all arguments by name
%   eeg = flt_wavelet('Signal',eeg,'NumLevels',5,'Family','sym4')
%
% See also:
%   wfilters, waveinfo
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% flt_wavelet_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','WaveletTransform', 'depends','set_makepos', 'follows',{'flt_fourier'}, 'independent_channels',true, 'independent_trials',true);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'wfamily','Family','family'}, 'db10', {'db2','db3','db4','db5','db6','db7','db8','db9','db10','db11','db12','db13','db14','db15','db16','db17','db18','db19','db20','coif2','coif3','coif4','coif5','sym2','sym3','sym4','sym5','sym6','sym7','sym8','sym9','sym10','sym11','sym12','sym13','sym14','sym15','sym16','sym17','sym18','sym19','sym20'}, 'Wavelet Family. The wavelet family to use; Typical examples are db2 - db20, sym2 - sym45, and coif2 - coif5.'), ...
    arg({'wlevels','NumLevels','numlevels'}, 10, uint32([0 1 20 100]), 'Number of decomposition levels. If this is > 1, a multi-level wavelet decomposition will be done.'));

% generate coefficients if necessary
if ischar(wfamily) %#ok<USENS>
	[lo,hi] = hlp_diskcache('filterdesign',@wfilters,wfamily,'d');
else
    [lo,hi] = deal(wfamily{:});
end
lf = length(lo);

for f = utl_timeseries_fields(signal) %#ok<NODEF>
    if ~isempty(signal.(f{1}))
        % flip dimensions so we can filter along the first dimension
        [X,dims] = spatialize_transpose(signal.(f{1}));
        C = [];
        % for each decomposition level...
        for i = 1:wlevels
            % figure out the sub-range to retain
            lx = size(X,1); first = 2; last = lx+lf-1;
            % pad the signal
            Y = X([lf-1:-1:1, 1:lx, lx:-1:lx-(lf-1)+1],:);
            % compute approximation coefficients
            Z = conv2(Y,lo(:),'valid');
            X = Z(first:2:last,:);
            % compute & append detail coefficients
            Z = conv2(Y,hi(:),'valid');
            C = [Z(first:2:last,:); C]; %#ok<AGROW>
        end
        C = [X; C]; %#ok<AGROW>
        % unflip dimensions and write back the result
        signal.(f{1}) = unspatialize_transpose(C,dims);
    end
end
signal.pnts = size(signal.data,2);

exp_endfun('append_online',{'wfamily',{lo,hi}});
