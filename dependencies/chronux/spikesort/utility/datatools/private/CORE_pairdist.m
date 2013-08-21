function D = CORE_pairdist(X,Y,takesqrt,reuse,useSafe)
%CORE_PAIRDIST     Core computational routine for PAIRDIST.
%   D = CORE_PAIRDIST(X,Y,TAKESQRT,REUSEMEM),
%   given an M x P matrix X and an N x P matrix Y, returns the M x N 
%   matrix D such that D(i,j) = || X(i,:) - Y(j,:) ||, where ||.|| means
%   the Euclidean 2-norm.  The next 2 values are required (0 or 1) flags
%   such that:
%       TAKESQRT - If 0, D(i,j) = ||X(i,:)-Y(j,:)||^2 rather than the
%                  magnitude of the distance.
%       REUSE    - If 1 and if M and N are the same as (or smaller than)
%                  their values in the previous call to this function, 
%                  attempt to reuse memory rather than reallocating space
%                  for D.  **See WARNING below**
%
%   WARNING:  Setting the REUSE flag is not guaranteed to have an effect.
%   If it does work, setting REUSE to 1 can result in unexpected behavior
%   for the returned Matlab variable D.  First, clearing D will not free
%   the memory associated with it (i.e., the memory will not become
%   available for other Matlab variables).  To fully free this memory
%   without restarting Matlab, type "clear CORE_pairdist".  Second, 
%   CORE_PAIRDIST retains a handle to D and can alter its contents on
%   later calls, potentially unexpectedly.  See the help for PAIRDIST for
%   an example.
%
%   CONDITIONS
%   ----------
%   X and Y must be REAL 2-D arrays of type DOUBLE.
%   X and Y must have the same number of columns.
%   TAKESQRT and REUSEMEM must each be either 0 or 1 and of type DOUBLE.

%%%%%%%%%%%%%%%%%%%%%%%%%%% Prep inputs %%%%%%%%%%%%%%%%%%%%%%%%%%
normsqrX = sum(X.^2,2);
normsqrY = sum(Y.^2,2);

[N,P1] = size(X);   % not going to check P1==P2 out of stubbornness --
[M,P2] = size(Y);   %   CORE_ functions are described as not doing error checking


%%%%%%%%%%%%%%%%%%%%%%%%% Distance Calculation %%%%%%%%%%%%%%%%%%%%%%%
%         dist(x,y)^2 = (x-y)'(x-y) = x'x + y'y - 2 x'y,
% Note that we do this in two steps for memory efficiency (i.e., computing
% z = x'x + y'y - 2x'y all at once requires 4 matrices of size z to be in 
% memory at the same time, while breaking it up only requires 2 matrices at
% any given time.  Still inefficient, i know . . . thats why there's a MEX 
% alternative).
D = (normsqrX * ones(1,M));
D = D + (ones(N,1) * normsqrY');
D = D - (2 * X * Y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Postprocess %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (takesqrt),  D = sqrt(D);  end;


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% (requires the Statistics toolbox)
% X = randn(1000,10);
% tic;  D1 = CORE_pairdist(X,X,1,0);            t(1) = toc;
% tic;  D2 = squareform(pdist(X,'euclidean'));  t(2) = toc; 
% printf('\nCORE_pairdist took %5.3f sec and equivalent native code took %5.3f sec.', t(1), t(2));
% printf('The RMS error between the two results was %6.4g.\n', sqrt(mean((D1(:)-D2(:)).^2)));



