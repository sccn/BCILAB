function C = CORE_histxt(X,d)
%CORE_HISTXT       Core computational routine for HISTXT. 
%   C = CORE_HISTXT(X,D), where X is an M x N matrix, returns a D x N
%   matrix C such that C(i,j) = #[X(:,j)=i].
%
%   CONDITIONS
%   ----------
%   X must only contain integer values between 1 and D.
%   X must be of type double.

C = hist(X,1:d);


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = ceil(256*rand(1e5,10));
% tic; counts = CORE_histxt(data,256);  t(1) = toc;
% tic; counts2 = hist(data,1:256); t(2) = toc;
% printf('\nCORE_histxt took %5.3f sec and equivalent native code took %5.3f sec.', t(1), t(2));
% if (~isequal(counts,counts2))
%     printf('The two calls did not produce the same results.');
% end
