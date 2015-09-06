function [data sigma mu] = sim_genTVARdata(varargin)
%  Simulate n time steps of the VAR(p) process
%
%     v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + eta(k,:)',
%
%  where A=[A1 ... Ap] is the p-lag coefficient matrix, and w is a vector
%  of intercept terms that is included to allow for a nonzero mean of the
%  process. The vectors eta(k,:) are noise vectors with mean zero and
%  covariance matrix C. The noise can have either a hyperbolic secant
%  distribution or a generalized Gaussian distribution with shape
%  parameter beta (def=2) and scale parameter alpha (def=1).
%  If beta=2 (default), the noise is normally distributed.
%  If beta < 2, the noise is super-gaussian (beta=1 is the laplacian).
%  If beta > 2, the noise is super-gaussian
%  (beta=Inf is the uniform distribution).
%
%  Input A may be a matrix or a cell array of size equal to n (the length
%  of the process). In this case, A{t} is interpreted as the time-varying
%   VAR coefficient matrix (same dimension as above) at time t.
%
%  The p vectors of initial values for the simulation are taken to
%  be equal to the mean value of the process. (The process mean is
%  calculated from the parameters A and w.) To avoid spin-up effects,
%  the first 10^3 time steps are discarded. Alternatively,
%  the user can specify how many initial samples to discard.
%
%
% -------------------------------------------------------------------------
% Input                Information
% -------------------------------------------------------------------------
% ModelCoeffs:         Model coefficients
%                      Input Range  : Unrestricted
%                      Default value: MANDATORY INPUT
%                      Input Data Type: string
%
% NoiseCovMat:         Noise covariance matrix
%                      If scalar, assumes uniform noise variance (cov
%                      matrix is a diagonal matrix with this value on
%                      main diagonal)
%                      Input Range  : Unrestricted
%                      Default value: 1
%                      Input Data Type: real number (double)
%
% ProcessMean:         Process mean
%                      If empty, zero-mean process is assumed. If
%                      scalar, all processes have same mean. Otherwise
%                      can be a [1 x M] row vector specifying the mean
%                      for each of M processes
%                      Input Range  : Unrestricted
%                      Default value: n/a
%                      Input Data Type: real number (double)
%
% NumSamplesToDiscard: Number of 'burn-in' samples to discard
%                      Input Range  : Unrestricted
%                      Default value: 1000
%                      Input Data Type: real number (double)
%
% TrialLength:         Trial length in samples
%                      If empty, determined automatically from the
%                      length of cell array A
%                      Input Range  : [1  Inf]
%                      Default value: n/a
%                      Input Data Type: real number (double)
%
% NumTrials:           Possible values: Number of trials (realizations)
%                      Default value  : 'd'
%                      Input Data Type: real number (double)
%
% NoiseDistribution:   Noise distribution
%                      hsec: Hyperbolic secant distribution.
%                        gengauss: Generalized gaussian distribution
%                      Possible values: 'gengauss','hsec'
%                      Default value  : 'gengauss'
%                      Input Data Type: string
%
%     | ScaleParam:    Scale parameter
%                      This determines the variance.
%                        var =
%                      alpha^2*(gamfunc(3/beta))/gamfunc(1/beta).
%                      Input Range  : [0  Inf]
%                      Default value: 1
%                      Input Data Type: real number (double)
%
%     | ShapeParam:    Shape parameter
%                      This determines the kurtosis.
%                        If beta=2 (default), the distribution is
%                      Gaussian.
%                        If beta < 2, the distribution is super-gaussian
%                      (beta=1 is the laplacian).
%                        If beta > 2, the distribution is sub-gaussian.
%                        beta=Inf is the uniform distribution.
%                      Input Range  : [0  Inf]
%                      Default value: 2
%                      Input Data Type: real number (double)
%
% -------------------------------------------------------------------------
% Output                Information
% -------------------------------------------------------------------------
% Data:                Data matrix             [nchs x npnts x ntr]
% sigma:               Noise covariance matrix [nchs x nchs]
% mu:                  Process mean vector     [1 x nchs]
%
%
% Author: Tim Mullen 2011-2013, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

g = arg_define([0 Inf],varargin,...
    arg_norep({'A','ModelCoeffs'},mandatory,[],'Model coefficients'), ...
    arg({'sigma','NoiseCovMat'},1,[],'Noise covariance matrix. If scalar, assumes uniform noise variance (cov matrix is a diagonal matrix with this value on main diagonal)','shape','matrix'), ...
    arg({'mu','ProcessMean'},[],[],'Process mean. If empty, zero-mean process is assumed. If scalar, all processes have same mean. Otherwise can be a [1 x M] row vector specifying the mean for each of M processes','shape','row'), ...
    arg({'ndisc','NumSamplesToDiscard'},1000,[],'Number of ''burn-in'' samples to discard'), ...
    arg({'Nl','TrialLength'},[],[1 Inf],'Trial length in samples. If empty, determined automatically from the length of cell array A'), ...
    arg({'Nr','NumTrials'},100,[1 Inf],'Number of trials (realizations)'), ...
    arg_subswitch({'noiseDistrib','NoiseDistribution'},'gengauss',...
        {'gengauss', ...
        {...
        arg({'alpha','ScaleParam'},1,[0 Inf],{'Scale parameter.', sprintf('\nThis determines the variance.\n  var = alpha^2*(gamfunc(3/beta))/gamfunc(1/beta).')}), ...
        arg({'beta','ShapeParam'},2,[0 Inf],{'Shape parameter.', sprintf('\nThis determines the kurtosis.\n  If beta=2 (default), the distribution is Gaussian.\n  If beta < 2, the distribution is super-gaussian (beta=1 is the laplacian).\n  If beta > 2, the distribution is sub-gaussian.\n  beta=Inf is the uniform distribution.')}), ...
        }, ...
    'hsec' {}, ...
    'user_defined', {...
        arg({'pdf','PDF'},@(x,a,b) a + b*log(tan(pi*x/2)),[],'Function handle for pdf. Must be of form f(x,param1,param2,...) accepting uniform random input x') ...
        arg({'params'},{1,pi/2},[],'Cell array of parameters for pdf') ...
        }}, ...
    {'Noise distribution.',sprintf('\n  hsec: Hyperbolic secant distribution.\n  gengauss: Generalized gaussian distribution')}) ...
    );

% get number of channels
if iscell(g.A)
    M = size(g.A{1},1);
else
    M = size(g.A,1);
end
% determine trial length
if isempty(g.Nl)
    if iscell(g.A)
        g.Nl = length(g.A);
    else
        error('SIFT:sim_genTVARdata:badInput', ...
            'Cannot automatically determine the trial length. Please specify the TrialLength manually.');
    end
end
% check noise covariance
if isempty(g.sigma)
    g.sigma = 1;
end
if isscalar(g.sigma)
    g.sigma = g.sigma*eye(M);
elseif any(size(g.sigma) ~= M)
    error('SIFT:sim_genTVARdata:wrongArraySize',...
        'Noise covariance matrix must be scalar or an [M x M] square matrix, where M is the number of channels in the model');
end
% check mean
if isempty(g.mu)
    g.mu = zeros(1,M);
elseif isscalar(g.mu)
    g.mu = g.mu(ones(1,M));
elseif size(g.mu,2) ~= M
    error('SIFT:sim_genTVARdata:wrongArraySize',...
        'Process mean must be scalar or a [1 x M] vector, where M is the number of channels in the model');
end

% simulate the data
switch g.noiseDistrib.arg_selection
    case 'hsec'
        % hyperbolic secant noise
        data = tvarsim(g.mu, g.A, g.sigma, [g.Nl g.Nr], g.ndisc, 2/pi, 0, 'hsec');
    case 'gengauss'
        % generalized gaussian noise
        data = tvarsim(g.mu, g.A, g.sigma, [g.Nl g.Nr], g.ndisc, g.noiseDistrib.beta, g.noiseDistrib.alpha, 'gengauss');
    case 'user_defined'
        % user-defined distribution
        data = tvarsim(g.mu, g.A, g.sigma, [g.Nl g.Nr], g.ndisc, [], [], [{g.noiseDistrib.pdf} g.noiseDistrib.params]);
end
% return data as [nchs x npnts x ntr]
data = permute(data,[2 1 3]);

if nargout > 2
    sigma  = g.sigma;
    mu = g.mu;
end