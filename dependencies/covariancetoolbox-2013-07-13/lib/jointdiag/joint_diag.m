% V = joint_diag(C,method_diag,arg_diag)
%
% -----------------------------Definition---------------------------------%
%
% description : performs a joint diagonalization of a set of covariance
% matrices.
%
% usage : joint_diag(C,'uwedge',{})
%
% -----------------------------Input--------------------------------------%
%
% C : a cell array of Nc-by-Nc covariance matrices.
% method_diag : the joint diagonalization method. 'uwedge', 'jade' methods
% are available. default is 'uwedge'.
% arg_diag : parameters for the method, in a cell array. Default is empty
%
% -----------------------------Output-------------------------------------%
%
% V : the matrix which best diagonalize the set of covariances matrices. Is
% a Nc-by-Nc matrix. The diagonalisation is givent by V'*C{1}*V.
%
% -----------------------------References---------------------------------%
%
%
%   Project : BCI-EEG
%
%   author : A. Barachant
%   date : 2011-10-26
%   version : 1.0 
%   status : a terminer
%   CEA/GRENOBLE-LETI/DTBS
%
%   See also uwedge, jade.

function V = joint_diag(C,method_diag,arg_diag)

if (nargin<2)||(isempty(method_diag))
    method_diag = 'uwedge';
end
if (nargin<3)
    arg_diag = {};
end

switch method_diag
    case 'jade'
        V = jade(C,arg_diag);
    case 'riemann'
        C = reshape(cell2mat(C)',size(C{1},1),size(C{1},1),[]);
        V = riemann_diag(C,arg_diag)';
    otherwise
        V = uwedge(cell2mat(C)');
end
% [EOF: joint_diag.m]
