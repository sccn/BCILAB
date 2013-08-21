function Y = CORE_leadingedges(X)
%CORE_LEADINGEDGES Core computational routine for LEADINGEDGES.
%   Y = CORE_LEADINGEDGES(X) takes an M x N matrix X and returns a M x N
%   matrix Y, containing a logical 1 in each location Y(j,k) such that
%   X(j,k) AND ~X(j,k-1) for k>1 and X(j,k) for k==1.
%
%   The matrix Y is of type logical.
%
%   CONDITIONS
%   ----------
%   X must be a matrix of type: DOUBLE, UINT8, LOGICAL.
%   X can not be sparse.
%   X should not contain NaN values; the results (i.e., whether NaN==0)
%     are compiler dependent.
%     

Y = (X(1,:) ~= 0);
Y = [Y; X(2:end,:) & ~X(1:end-1,:)];

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = uint8((rand(1e6,3)*255) > 250);
% tic;  trig = CORE_leadingedges(data);  t(1) = toc;
% tic;  trig2 = [(data(1,:)~=0) ; [(data(2:end,:)~=0) & (data(1:end-1,:)==0)]];  t(2) = toc;
% printf('\nCORE_leadingedges took %5.3f sec and equivalent native code took %5.3f sec.', t(1), t(2));
% if (~isequal(trig,trig2))
% 	printf('The two calls did not produce the same results.');
% end
