function out=locfit_all(varargin)

% Smoothing noisy data using Local Regression and Likelihood.
%
% This is a combination of the locfit and predict functions
%

% Minimal input validation    
if nargin < 1
   error( 'At least one input argument required' );
end

predict_args = {};

locfit_args = varargin{1};

if nargin==2
predict_args = varargin{2};
end;

fit = locfit( locfit_args{:} );

predict_out = predict( fit, predict_args{:} );

out = {fit predict_out};
