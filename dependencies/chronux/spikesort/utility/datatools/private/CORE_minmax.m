function [xmn,xmx,mni,mxi] = CORE_minmax(x)
%CORE_MINMAX       Core computational routine for MINMAX.
%   [XMN,XMX,MNI,MXI] = CORE_MINMAX(X) returns scalars such that
%   XMN = X(MNI) = min(X) and XMX = X(MXI) = max(X).  Ties in indices are
%   broken in favor of the lowest magnitude index.
%   
%   NaN values are ignored unless the input is all NaN.  In this case, XMN
%   and XMX are set to NaN, while MNI and MXI are set to 1.  This mimics
%   the behavior of the Matlab native MIN/MAX functions. 
%
%   CONDITIONS
%   ----------
%   X must be a real vector of type DOUBLE.  An N-D array X is treated as X(:).
%   Infinite values are allowed.   NaN's are ignored (see above).

[xmn,mni] = min(x(:));
[xmx,mxi] = max(x(:));

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 5e6;
% array = randn(N,1);  array(1:1e4:end) = NaN;
% % array = [1:N];     array(1:10:end) = NaN;  % much slower
% tic; [mminval,mmaxval,mminind,mmaxind] = CORE_minmax(array);  t(1) = toc;
% tic; [minval,minind] = min(array);   [maxval,maxind] = max(array); t(2) = toc;
% printf('\nCORE_minmax took %5.3f sec and Matlab code took %5.3f sec.', t(1), t(2));
% if (~isequal([mminval,mmaxval,mminind,mmaxind], [minval,maxval,minind,maxind]))
%     printf('The two calls did not produce the same results.');
% end

