function [featureshape,trials,vectorize_trials] = utl_determine_featureshape(trials,shape,multitask)
% Uniformize the given trials and shape information.
% [FeatureShape,Trials,VectorizeTrials] = utl_determine_featureshape(Trials,Shape,Multitask)
%
% This function deals with the fact that trials can be represented either as [NxF] matrix of
% vectorized features (N=#observations, F=#features), or as [AxBxCx...xN] array of tensor-shaped
% features (A,B,C, ... = #elements in the respective dimension), which is a bit inconsistent for
% historical reasons.
%
% Takes trials in a variety of shapes, and an (optionally non-empty) desired shape vector and
% uniformizes them into a [#trials x #features] matrix of vectorized trials, the deduced feature
% shape, and a boolean flag of whether trials had to be vectorized. By default also handles a cell
% array of multiple trial arrays, which are assumed to stem from multiple compatibly shaped tasks.
%
% In:
%   Trials : A NxF matrix of observations with pre-vectorized features, or a [AxBxCx...xN] array
%            of tensor-shaped features. Can also be a cell array of multiple such arrays where all
%            dimensions except N must be identical.
%
%   Shape : A vector that is the desired output of size() applied to a single observation, that is,
%           a vector of sizes. If [], shape information will be deduced from Trials, otherwise
%           overrides the shape of Trials.
%
%   Multitask : If true, the output will be a cell array of trial matrices. Otherwise it will be 
%               a matrix (and an error will be thrown if there was more than one task in the
%               input).
%
% Out:
%   FeatureShape : Final 1xD vector that describes the feature shape (same as size(one_trial)).
%
%   Trials : Vectorized version of input Trials (#trials x #features).
%
%   VectorizeTrials : whether trials had to be vectorized.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-02-04

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
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


if ~iscell(trials)
    trials = {trials}; end
featureshape = cell(1,length(trials));
vectorize_trials = cell(1,length(trials));
% for each task (each of which has multiple trials)...
for t=1:length(trials)
    if ndims(trials{t}) > 2 %#ok<ISMAT>
        featureshape{t} = size(trials{t}); featureshape{t} = featureshape{t}(1:end-1);
        if ~isempty(shape) && ~isequal(shape,featureshape{t})
            if prod(featureshape{t}) == prod(shape)
                warning('You are specifying a shape property but also multidimensional features of a different shape; using the explicit shape parameter.');
                featureshape{t} = shape;
            else
                error('You are specifying a shape property but also features with an incompatible number of elements. Please correct.');
            end
        end
        trials{t} = double(reshape(trials{t},[],size(trials{t},ndims(trials{t})))');
        vectorize_trials{t} = true;
    else
        if ~isempty(shape)
            if size(shape,1) > 1
                if all(all(bsxfun(@eq,shape(1,:),shape)))
                    featureshape{t} = [shape(1,1),shape(1,2),size(shape,1)];
                    if prod(featureshape{t}) == size(trials{t},2)
                        % the reason is that it is much more efficient to operate on a dense 3d array than a very sparse 2d array
                        warn_once('This method will by convention reshape block-diagonalized feature matrices with identical blocks into a 3d tensor. This warning will only come up once.');
                    else
                        error('Your shape parameter has a different number of features than your data.');
                    end
                else
                    % we don't implement block-diagonalization in here
                    error('This method does not handle implicitly block-diagonal features; please either reformulate in tensor form or pass a large sparse data matrix (pre-blockdiagonalized). Note that the tensor form is likely several times faster.');
                end
            elseif prod(shape) ~= size(trials{t},2)
                error('Your shape parameter has a different number of features than data.');
            else
                featureshape{t} = shape;
            end
        else
            featureshape{t} = [size(trials{t},2),1];
        end
        vectorize_trials{t} = false;
    end
end
vectorize_trials = unique([vectorize_trials{:}]);
if length(vectorize_trials)>1 || ~all(cellfun(@(x)isequal(x,featureshape{1}),featureshape))
    error('The number or shape of features must be the same for each task.'); end    
featureshape = featureshape{1};
if nargin >=3 && ~multitask
    if length(trials) > 1
        error('Multi-task learning is disabled by multiple tasks were given.');
    else
        trials = trials{1};
    end
end
