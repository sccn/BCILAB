function res = hlp_findscaling(X, scaling)
% Obtain information necessary to scale the given data. Works with hlp_applyscaling.
% Scale-Info = hlp_findscaling(Data, Scale-Mode)
%
% This is just a convenience tool to implement simple data (e.g. feature) scaling operations.
%
% In:
%   Data        : data matrix of [Observations x Variables]
%   Scale-Mode  : scaling mode, with the following options:
%                 'center' : shift the data to make them zero-mean
%                 'std' : center and standardize the data
%                 'std-nocenter': just standardize but do not center the data
%                 'minmax' : scale the data to the range 0-1
%                 'whiten' : whiten/sphere the data in a multivariate manner (rotates into PCA space)
%                 'rescale' : rescale the data by 1 / median standard deviation
%                 'decorrelate' : sphere/whiten/decorrelate the data without rotating
%                 'decorrelate_shrinkage' : sphere/whiten/decorrelate the data without rotating, using shrinkage regularization
%                 'decorrelate_robust' : sphere/whiten/decorrelate the data without rotating, using outlier-robust estimation
%
% Out:
%   Scale-Info  : scaling structure that can be used with hlp_applyscaling, to scale data
%
% Examples:
%   scaleinfo = hlp_findscaling(data,'whiten')
%   hlp_applyscaling(data,scaleinfo)
%
% See also:
%   hlp_applyscaling
%
%               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%               2010-03-28

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if nargin < 2
    scaling = 'minmax'; end
if ~ischar(scaling)
    error('The given scaling mode must be a string.'); end

switch scaling
    case 'center'
        res = struct('add',{-mean(X)});
    case 'std'
        res = struct('add',{-mean(X)}, 'mul',{1./ std(X)});
        res.mul(~isfinite(res.mul(:))) = 1;
    case 'std-nocenter'
        res = struct('add',{zeros(1,size(X,2))}, 'mul',{1./ std(X)});
        res.mul(~isfinite(res.mul(:))) = 1;
    case 'minmax'
        res = struct('add',{-min(X)}, 'mul',{1./ (max(X) - min(X))});
        res.mul(~isfinite(res.mul(:))) = 1;
    case 'whiten'
        [Uc,Lc] = eig(cov(X));
        res = struct('add',{-mean(X)},'project',{Uc * sqrt(inv(Lc))'});
        res.project(~isfinite(res.project(:))) = 1;
    case 'decorrelate'
        res = struct('add',{-mean(X)},'project',{cov(X)^(-1/2)});
        res.project(~isfinite(res.project(:))) = 1;
    case 'decorrelate_shrinkage'
        res = struct('add',{-mean(X)},'project',{cov_shrink(X)^(-1/2)});
        res.project(~isfinite(res.project(:))) = 1;
    case 'decorrelate_robust'
        res = struct('add',{-median(X)},'project',{cov_blockgeom(X)^(-1/2)});
        res.project(~isfinite(res.project(:))) = 1;
    case 'rescale'
        res = struct('add',{zeros(1,size(X,2))}, 'mul',{ones(1,size(X,2)) ./ median(std(X))});
        res.mul(~isfinite(res.mul(:))) = 1;
    otherwise
        if ~isempty(scaling) && ~strcmp(scaling,'none')
            error('hlp_findscaling: unknown scaling mode specified: %s',scaling); end
        res = struct();
end
