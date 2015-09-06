function out_opt = setopt(varargin)

%SETOPT Set options.
%   out_opt = SETOPT(optname1, optvalue1, optname2, optvalue2, ...)
%   options:
%       usec - set to true to use algorithms in c code.
%       inv - inverse method used for Kalman filter - 'ge', 'sy', or 'po'.
%       disp - 'off', 'notify', 'final' or 'iter'.
%       tol - tolerance.
%       fmin - minimization algorithm to use ('simplex' or 'bfgs').
%       maxiter - maximum number of iterations.
%       nsamp - number of samples to use in simulation.
%       antithetic - number of antithetic variables to use in simulation.
%       hideinit - hide values returned for exact initialization.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

persistent opt;

if isempty(opt)
    %% Set default values %%
    opt.debug       = true;
    opt.usec        = true;
    opt.inv         = 0; % 'ge'
    opt.disp        = 1; % 'notify'
    opt.tol         = 10^-7;
    opt.fmin        = 'simplex';
    opt.maxiter     = 1000;
    opt.nsamp       = 1000;
    opt.antithetic  = 1;
    opt.hideinit    = true;
end

out_opt = opt;

if nargin > 0
    change  = ~isempty(varargin{1});
    if ~change, varargin(1) = []; end
    for i = 1 : 2 : length(varargin)
        optname  = varargin{i};
        optvalue = varargin{i+1};
        if ~ischar(optname) || ~isfield(out_opt, optname), error('ssm:setopt:InvalidOptName', 'Invalid option name.'); end
        if strcmp(optname, 'disp') && ischar(optvalue)
            switch optvalue
                case 'off',     out_opt.disp = 0;
                case 'notify',  out_opt.disp = 1;
                case 'final',   out_opt.disp = 2;
                case 'iter',    out_opt.disp = 3;
                otherwise,      out_opt.disp = 1;
            end
        elseif strcmp(optname, 'inv') && ischar(optvalue)
            switch optvalue
                case 'ge', out_opt.inv = 0;
                case 'sy', out_opt.inv = 1;
                case 'po', out_opt.inv = 2;
                otherwise, out_opt.inv = 0;
            end
        else out_opt.(optname) = optvalue;
        end
    end
    if change, opt = out_opt; end
end


