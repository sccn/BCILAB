function [AR, PE, argsout]=mvar_arfit(varargin)
% Algorithm:  ARfit
%
% Description:  
% 
% Stepwise least squares estimation of the parameters 
% [A1, ... Ap] and noise covariance matrix C of a 
% multivariate autoregressive (VAR) model of order p:
%
% v(t) = mu + A1*v(t-1) +...+ Ap*v(t-p) + noise(C)
%
% Author Credits:
% 
% This function is a wrapper for arfit() from the 
% ARFIT toolbox [1-3]
%
% References and Code:
%
% [1] A. Neumaier and T. Schneider. ACM TOMS, 27(1):27-57, 2001.
% [2] T. Schneider and A. Neumaier. ACM TOMS, 27(1):58-65, 2001.
% [3] http://www.gps.caltech.edu/~tapio/arfit/
%
% Dependencies: arfit()
% ------------------------------------------------------------------------

g = arg_define(varargin, ...
            arg_norep({'data','Data'},mandatory,[],'[chans x time x trials] data matrix'), ...
            arg({'morder','ModelOrder','p','pmax'},10,[],'Maximum model order'), ...
            arg({'selector','OrderSelector'},'sbc',{'sbc','fpe'},'Model order selection criterion'), ...
            arg({'no_const','NoConstantTerm'},true,[],'Exclude constant term in model. Setting this implicitly assumes the process mean is zero.') ...
            );

arg_toworkspace(g);
        
optargs = fastif(no_const,{'zero'},{});

% evaluate the subfunction
[argsout.mu, AR, PE, sbc, fpe, argsout.th]=arfit(permute(data,[2 1 3]), morder, morder, selector, optargs{:});
