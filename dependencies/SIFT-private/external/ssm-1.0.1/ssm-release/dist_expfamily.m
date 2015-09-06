function dist = dist_expfamily(b, d2b, id2bdb, c)

%DIST_EXPFAMILY Create SSDIST object for general exponential family distribution.
%   dist = DIST_EXPFAMILY(b, d2b, id2bdb, c)
%       b is a function of theta.
%       d2b is the second derivative of b.
%       id2bdb is the inverse of d2b times db.
%       c is a function of y.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4, error('ssm:dist_expfamily:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isa(b, 'function_handle'), error('ssm:dist_expfamily:InputError', 'b must be a function handle.'); end
if ~isa(d2b, 'function_handle'), error('ssm:dist_expfamily:InputError', 'd2b must be a function handle.'); end
if ~isa(id2bdb, 'function_handle'), error('ssm:dist_expfamily:InputError', 'id2bdb must be a function handle.'); end
if ~isa(c, 'function_handle'), error('ssm:dist_expfamily:InputError', 'c must be a function handle.'); end

    function [H y] = matf_expfamily(y, theta)
        H   = 1./d2b(theta);
        y   = theta - id2bdb(theta) + H.*y;
    end

    function logp = logpf_expfamily(y, theta)
        logp    = diag(y'*theta)' - b(theta) + c(y);
    end

dist    = ssdist(0, @matf_expfamily, @logpf_expfamily);
end

