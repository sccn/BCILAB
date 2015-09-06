function varargout = findpeaksND(x,minDist, exculdeEqualPoints, numpeaks)
% function varargout = localMaximum(x,minDist, exculdeEqualPoints)
%
% pkval = findpeaksND(x,3); returns the height of all peaks with a minimum
%           between-peak distance of 3
% [pkval idx] = findpeaksND(x,3);  additionally returns the linear
%               indices of the peaks
% [pkval idx1 idx2 ... idxN] = findpeaksND(x,3);  returns the 
%               [idx1 idx2 ... idxN] matrix indices of the peaks
% 
% 
% This function returns the indexes\subscripts of local maximum in the data x.
% x can be a vector or a matrix of any dimension
%
% minDist is the minimum distance between two peaks (local maxima)
% minDist should be a vector in which each argument corresponds to it's
% relevant dimension OR a number which is the minimum distance for all
% dimensions
%
% exculdeEqualPoints - is a boolean definning either to recognize points with the same value as peaks or not
% x = [1     2     3     4     4     4     4     4     4     3     3     3     2     1];  
%  will the program return all the '4' as peaks or not -  defined by the 'exculdeEqualPoints'
% localMaximum(x,3)
% ans = 
%      4     5     6     7     8     9    11    12
%
%  localMaximum(x,3,true)
% ans =
%      4     7    12
%      
%
% Example:
% a = randn(100,30,10);
% minDist = [10 3 5];
% peaks = localMaximum(a,minDist);
% 
% To recieve the subscript instead of the index use:
% [xIn yIn zIn] = localMaximum(a,minDist);
%
% To find local minimum call the function with minus the variable:
% valleys = localMaximum(-a,minDist);
%
% DOWNLOADED FROM: http://www.mathworks.com/matlabcentral/fileexchange/14498
%
    if nargin < 3
        exculdeEqualPoints = false;
        if nargin < 2
            minDist = size(x)/10;
        end       
    end
    if nargin<4
        numpeaks = [];
    end
    % if all values are identical, just return the index of the first one
    if all(diff(x(:))==0)
        varargout{1} = x(1);
        if nargout>1
            varargout(2:nargout)=num2cell(ones(1,nargout-1));
        end
        return;
    end
    
    if isempty(minDist)
        minDist = size(x)/10;
    end
    
    dimX = length ( size(x) );
    if length(minDist) ~= dimX
        % In case minimum distance isn't defined for all of x dimensions
        % I use the first value as the default for all of the dimensions
        minDist = minDist( ones(dimX,1) );
    end
    
    % validity checks
    minDist = ceil(minDist);
    minDist = max( [minDist(:)' ; ones(1,length(minDist))] );
    minDist = min( [minDist ; size(x)] );

    % ---------------------------------------------------------------------
    if exculdeEqualPoints
        % this section comes to solve the problem of a plato
        % without this code, points with the same hight will be recognized as peaks
        y = sort(x(:));
        dY = diff(y);
        % finding the minimum step in the data
        minimumDiff = min( dY(dY ~= 0) );   
        %adding noise which won't affect the peaks
        x = x + rand(size(x))*minimumDiff;
    end
    % ---------------------------------------------------------------------
    
    
    se = ones(minDist);
    X = imdilate(x,se);
    f = find(x == X);
    peakvals = x(f);
    
    if ~isempty(numpeaks)
        % select only the largest numpeaks 
        [peakvals idxx] = sort(peakvals,'descend');
        peakvals = peakvals(1:min(numpeaks,length(peakvals)));
        f = f(idxx(1:min(numpeaks,length(peakvals))));
    end
    
    varargout{1} = peakvals;
    if nargout > 2
        [varargout{2:nargout}] = ind2sub( size(x), f );
    elseif nargin > 1
        varargout{2} = f;
    end
