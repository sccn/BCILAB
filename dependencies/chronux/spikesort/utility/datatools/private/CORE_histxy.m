function Z = CORE_histxy(c,r,cmax,rmax)
%CORE_HISTXY       Core computational routine for HISTXY.
%   Z = CORE_HISTXY(C,R,CMAX,RMAX) returns an RMAX x CMAX matrix Z such
%   that Z(i,j) = #[R(a)=i,C(b)=j].
%
%   CONDITIONS
%   ----------
%   R and C must be column vectors of the same length.
%   R and C must only contain integer values between 1 and RMAX or CMAX
%       respectively.
%   R and C must be of type DOUBLE.
%   RMAX and CMAX must be integer-valued and of type DOUBLE.

Z = full(sparse(r,c,1,rmax,cmax));


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rows = ceil(256*rand(1e6,1));  cols = ceil(128*rand(1e6,1));
% tic; counts = CORE_histxy(cols,rows,128,256);  t(1) = toc;
% tic; counts2 = full(sparse(rows,cols,1,256,128));  t(2) = toc;
% printf('\nCORE_histxy took %5.3f sec and equivalent native code took %5.3f sec.', t(1), t(2));
% if (~isequal(counts,counts2))
%     printf('The two calls did not produce the same results.');
% end
