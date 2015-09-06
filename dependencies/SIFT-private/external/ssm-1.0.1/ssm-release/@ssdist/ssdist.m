function distribution = ssdist(type, matf, logpf, n, p)

%@SSDIST/SSDIST State space non-Gaussian distribution class constructor.
%   distribution = SSDIST(type[, matf, logpf]) creates a general
%       non-Gaussian distribution.
%   distribution = SSDIST(type, matf, logpf[, n, p]) creates a constant
%       general non-Gaussian distribution.
%       type is 0 for exponential family distributions and 1 for additive
%           noise distributions.
%       matf is the function that generates the approximated Gaussian variance
%           matrix given observations and signal or disturbances.
%       logpf is the function that calculates the log probability of
%           observations given observation and signal or disturbances.
%       n and p is the number of time points and number of variables of the
%           non-Gaussian distribution specified by matf and logpf.
%       (if not both matf and logpf is provided, they're assumed to be
%           variable.)
%       (Only one non-Gaussian distribution can be specified, combine multiple
%           SSDIST objects for multiple distributions.)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin == 0
    %% Default constructor %%
    parent                  = ssmat;
    distribution.type       = [];
    distribution.matf       = {};
    distribution.logpf      = {};
    distribution.diagmask   = false(0);
    distribution.dmask      = false(1, 0);
elseif isa(type, 'ssdist')
    %% Copy constructor %%
    distribution            = type;
    return;
else
    if isa(type, 'ssmat')
        if size(type, 1) ~= size(type, 2), error('ssm:ssdist:ssdist:UnableToConvert', 'Non-square SSMATs cannot be converted to SSDIST.'); end
        %% Convert from SSMAT %%
        parent                  = type;
        distribution.type       = [];
        distribution.matf       = {};
        distribution.logpf      = {};
        distribution.diagmask   = false(size(type, 1), 0);
        distribution.dmask      = false(1, 0);
    else
        if nargin < 1, error('ssm:ssdist:ssdist:NotEnoughInputs', 'Insufficient input arguments.'); end
        if ~isnumeric(type) || ~isscalar(type) || all(type ~= [0 1]), error('ssm:ssdist:ssdist:InputError', 'type must be 0 or 1.'); end

        if nargin < 3
            %% Variable general non-Gaussian distribution %%
            parent                  = ssmat(0, [], true, 0);
            distribution.type       = type;
            distribution.matf       = {0}; % not set until updated
            distribution.logpf      = {0}; % not set until updated
            distribution.diagmask   = true;
            distribution.dmask      = true;
        else
            %% Constant general non-Gaussian distribution %%
            if ~isa(matf, 'function_handle') && ~isequal(matf, 0), error('ssm:ssdist:ssdist:InputError', 'matf must be a function handle.'); end
            if ~isa(logpf, 'function_handle') && ~isequal(logpf, 0), error('ssm:ssdist:ssdist:InputError', 'logpf must be a function handle.'); end
            if nargin < 4, n = 1; elseif ~isnumeric(n) || ~isscalar(n), error('ssm:ssdist:ssdist:InputError', 'n must be a scalar.'); end
            if nargin < 5, p = 1; elseif ~isnumeric(p) || ~isscalar(p), error('ssm:ssdist:ssdist:InputError', 'p must be a scalar.'); end
            parent                  = ssmat(zeros(p), [], true(p), zeros(p^2, n));
            distribution.type       = type;
            distribution.matf       = {matf};
            distribution.logpf      = {logpf};
            distribution.diagmask   = true(p, 1);
            distribution.dmask      = false;
        end
    end
end

superiorto('ssmat');

%% Register object instance %%
distribution    = class(distribution, 'ssdist', parent);

