function varargout = hlp_progBar(verbLevel,varargin)
% simple switch function to determine whether to generate/update 
% text-based progress bar or graphical progress bar.
% 
% Inputs:
%   verbLevel:  verbosity level:
%               1: text-based progress bar (progress.m)
%               2: graphical progress bar (multiWaitbar.m)
%   varargin:   arguments to either progress() or multiWaitbar()

switch verbLevel
    case 1
        % generate text-based progress bar
        N = progress(varargin{:});
    case 2
        % generate graphical progress bar
        multiWaitbar(varargin{:});
end

if exist('N','var')
    varargout{1} = N;
end
