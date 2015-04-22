% [] = XXXX()
%
% -----------------------------Definition---------------------------------%
%
% description
%
% usage : 
%
% -----------------------------Input--------------------------------------%
%
% XX : description + dimension
%
% -----------------------------Output-------------------------------------%
%
% XX : description + dimension
%
% -----------------------------References---------------------------------%
%
% [1] : XXX
%
%
%   Project : BCI-EEG
%
%   author : A. Barachant
%   date : 2011-XXXX
%   version : 1.0 
%   status : a terminer, termin�   
%   CEA/GRENOBLE-LETI/DTBS
%
%   See also distance_riemann, distance_kullback, distance_logeuclid, norm.

% [EOF: XXX.m]

function d = distance(C1,C2,method_dist,arg_dist)

if (nargin<3)||(isempty(method_dist))
    method_dist = 'euclid';
end
if (nargin<4)
    arg_dist = {};
end

switch method_dist
    case 'riemann'
        d = distance_riemann(C1,C2);
    case 'kullback'
        d = distance_kullback(C1,C2);
    case 'logeuclid'
        d = distance_logeuclid(C1,C2);
    case 'opttransp'
        d = distance_opttransp(C1,C2);
    case 'ld'
        d = distance_ld(C1,C2);
    otherwise
        d = sqrt(norm(C1-C2,'fro'));
end
