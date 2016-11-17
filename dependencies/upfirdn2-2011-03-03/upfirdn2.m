function [Y,Zf] = upfirdn2(X,H,P,Q,Zi)
% Upsample by P, appy given FIR filter H, and downsample by Q.
% [Y,Zf] = upfirdn2mex(X,H,P,Q,Zi)
%
% Like upfirdn, except that final filter conditions Zf are being returned, which can
% be used as initial conditions of a subsequent run of upfirdn2.
%
% In:
%   X : the signal to be resampled; operating on columns
%
%   H : the filter kernel to use; column vector
%
%   P : upsampling factor, must be a positive integer
%
%   Q : downsampling factor, must be a positive integer
%
%   Zi : optional initial filter conditions
%
% Out:
%   Y : resampled singnal; column vector (note: Y may be delayed depending on the
%       delay of the FIR filter H
%
%   Zf : final conditions of the filter; this is a 3-element cell array
%
% See also:
%   upfirdn
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-03-03

% make sure that the binary is compiled...
persistent is_compiled;
if isempty(is_compiled)
    is_compiled = exist('upfirdn2mex','file'); end

% check inputs / outputs
warning off MATLAB:nargchk:deprecated 
error(nargchk(4,5,nargin,'struct'));
error(nargoutchk(0,2,nargout,'struct'));
v = size(X,2);
if ~isreal(X)
    error('X cannot be complex-valued.'); end
if ~isvector(H)
    error('H must be a vector.'); end
H = H(:);

if is_compiled
    try
        % run mex file on each column of X
        Y = cell(1,v);
        if exist('Zi','var') && ~isempty(Zi)
            [Zi,Ti,Oi] = deal(Zi{:});
            if nargout > 1
                Zf = cell(1,v);
                for c = 1:v
                    [Y{c},Zf{c},Tf,Of] = upfirdn2mex(X(:,c),H,P,Q,Zi(:,c),Ti,Oi); end
                Zf = {[Zf{:}],Tf,Of};
            else
                for c = 1:v
                    Y{c} = upfirdn2mex(X(:,c),H,P,Q,Zi(:,c),Ti,Oi); end
            end
        else
            if nargout > 1
                Zf = cell(1,v);
                for c = 1:v
                    [Y{c},Zf{c},Tf,Of] = upfirdn2mex(X(:,c),H,P,Q); end
                Zf = {[Zf{:}],Tf,Of};
            else
                for c = 1:v
                    Y{c} = upfirdn2mex(X(:,c),H,P,Q); end
            end
        end
        Y = [Y{:}];
    catch
        is_compiled = false;
        warn_once('The MATLAB compiler on your platform produces incorrect code for the function "upfirdn2mex", which is required for online resampling. You can work around this by making sure that the sampling rate of your training data matches the rate produced by your online equipment.');
        if size(X,1) == 1
            % workaround for upfirdn shortcoming if there's 1 sample
            Y = X;
        else 
            Y = upfirdn(double(X),double(H),P,Q);
        end            
        Zf = {[],[],[]};
    end
else
    warn_once('Some code needs to be compiled to allow correct online resampling. Please consider setting up a supported compiler using the command "mex -setup".');
    if size(X,1) == 1
        % workaround for upfirdn shortcoming if there's 1 sample
        Y = X;
    else
        Y = upfirdn(double(X),double(H),P,Q);
    end
    Zf = {[],[],[]};
end
