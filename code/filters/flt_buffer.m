function [signal,state] = flt_buffer(varargin)
% Buffer samples and emit chunks.
% function [Signal,State] = flt_buffer(Signal,ChunkLength)
%
% In:   
%   Signal      : EEGLAB data set, either continuous or epoched
%
%   ChunkLength : chunk length to emit, in samples (default: 128)
%
%   State : input state
%
% Out:
%   Signal : chunked signal
%
%   State : output state
%
% Examples:
%   % use default settings
%   eeg = flt_buffer(eeg)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-02-12

if ~exp_beginfun('filter') return; end

declare_properties('name','Buffer', 'follows','flt_selchans', 'independent_channels',true, 'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'chunk_len','ChunkLength'}, 128, uint32([1 1000000]), 'Chunk length to buffer.'),...
    arg_nogui({'state','State'}));

if isempty(state)
    state.buffer = []; end

% prepend buffer to data
state.buffer = [state.buffer signal.data];
% emit the part that's ready
signal.data = state.buffer(:,1:end-mod(end,chunk_len));
% buffer the rest
state.buffer = state.buffer(:,end-mod(end,chunk_len)+1:end);

exp_endfun;
