function signal = flt_rej_impedance(varargin)
% Reject channels based on an impedance threshold; requires hardware support.
% this assumes that impedence is encoded in the peak-to-peak amplitdue  
% difference between every EncodingFactor samples
%
% Author: Tim Mullen, SCCN/INC/UCSD 2013

if ~exp_beginfun('filter') return; end

% used as a tool to select channel subsets before these ops are applied
declare_properties('name',{'RejectHighImpedanceChans'}, 'experimental',false, 'precedes',{'flt_laplace','flt_ica','flt_reref','flt_movavg'}, 'independent_trials',true, 'independent_channels',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'},[],[],'Signal structure. Must contain .data field with channel data'), ...
    arg({'impedanceThreshold','ImpedanceThreshold'},2,[0 Inf],'Impedance threshold (Mohms). Channels with impedence higher than this will be rejected.'), ...
    arg({'imp_period','ImpedancePeriod'},4,[],'Period (samples) for impedance calculation'), ...
    arg({'verb','Verbosity'},false,[],'Verbose output'), ...
    arg_norep('removed_channels',unassigned) ...
    );

% flag channels
if ~exist('removed_channels','var')
    
    if signal.pnts < imp_period*2
        return; end

    ISTIM = 0.000000024;
    TO_Z  = 1.4/(ISTIM*2.0);
    TO_V  = (1/1e6);

    npts =  signal.pnts-mod(signal.pnts,imp_period);
    nchs =  signal.nbchan;
    % get the peak and trough for every cycle of data
    datmax = max(reshape(signal.data(:,1:npts)',imp_period,nchs*(npts/imp_period)),[],1)';
    datmin = min(reshape(signal.data(:,1:npts)',imp_period,nchs*(npts/imp_period)),[],1)';
    datmax = reshape(datmax,(npts/imp_period),nchs);
    datmin = reshape(datmin,(npts/imp_period),nchs);
    % convert average peak-to-peak microvoltage to impedence
    imped  = mean(datmax-datmin)'*TO_V*TO_Z;
    removed_channels = find(imped > impedanceThreshold*1e6);
    
    if verb && any(removed_channels)
        labels = {signal.chanlocs.labels};
        if isempty(labels)
            labels = 1:signal.nbchan; 
        end
        removed_channels_lbl = labels(removed_channels);
        
        fprintf('Impedance check:\n');
        fprintf('Impedances:\t%s\n',hlp_tostring(imped));
        fprintf('Bad channels (%d):\t%s\n',length(removed_channels_lbl),hlp_tostring(removed_channels_lbl));
    end
end

% execute
if ~isempty(removed_channels)
    retain_channels = true(1,size(signal.data,1)); 
    retain_channels(removed_channels) = false;
    signal.data = signal.data(retain_channels,:,:,:,:,:,:,:);
    signal.chanlocs = signal.chanlocs(retain_channels);
    signal.nbchan = size(signal.data,1);
end

exp_endfun('append_online',{'removed_channels',removed_channels});