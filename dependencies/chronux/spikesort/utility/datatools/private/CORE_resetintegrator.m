function Y = CORE_resetintegrator(X,R)
%CORE_RESETINTEGRATOR Core computational routine for RESETINTEGRATOR.
%   Y = CORE_RESETINTEGRATOR(X,R) takes length N vectors X and R, where X
%   is of type double and R is of type LOGICAL and returns a length N
%   vector Y such that Y(j) = / Y(j-1) + X(j), if R(j) == 1
%                             \ 0            , if R(j) == 0
%
%   The matrix Y is of type double.
%
%   CONDITIONS
%   ----------
%   X must be a matrix of type DOUBLE.
%   X can not be sparse.
%   X should not contain NaN or Inf values; the results will be compiler
%     dependent.    

Y = zeros(size(X));
if (R(1)),  Y(1) = X(1);
else        Y(1) = 0;
end
for t = 2:length(X)
	if (R(t)),   Y(t) = Y(t-1) + X(t);
    else         Y(t) = 0;
	end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 1e6;  data = rand(N,1);   reset = (rand(N,1) > 0.01);
% tic;  int = CORE_resetintegrator(data,reset);  t = toc;
% err = [diff([0; int]) - data].^2;   err(~reset) = NaN;
% printf('\nCORE_resetintegrator took %5.3f sec with MSE %g.', t, nanmean(err));
% if (~all(int(~reset)==0)),  printf('Resetting failed.');  end;
