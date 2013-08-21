function Z = CORE_pxcorr(x,y,maxlag)
%CORE_PXCORR       Core computational routine for PXCORR.
%   Z = CORE_PXCORR(X,Y,MAXLAG), where X and Y are length M and N
%   vectors respectively, returns the length (2*MAXLAG+1) vector
%   Z such that Z(i) = #(X(j)-Y(k) == i), |i| <= MAXLAG.
%
%   CONDITIONS (including requirements for corresponding MEX file)
%   ----------
%   X and Y must contain sorted integer values.
%   X and Y must be row vectors of type DOUBLE.

Z = zeros(1,2*maxlag+1);

% We do this efficiently by stepping along X and keeping track of those
% indices in Y that are near by (i.e., within a distance of MAXLAG).
limit = length(x);
a = 1;  c = 1;   

for b = 1:length(y)
    while((y(b)-x(a)) > maxlag),        % move left bracket until less than MAXLAG before x(b)
            a=a+1;   if (a > limit), return; end;
    end
    if (c < a), c = a; end;             % catch up the right bracket if it was passed by
    if (c <= limit)
        while((x(c)-y(b)) <= maxlag),   % move right bracket until more than MAXLAG after x(b)
            c=c+1;   if (c > limit), break; end;
        end
    end

    offset = -y(b)+maxlag+1;            % add 'em up
    for bb = a:(c-1)
        ind = x(bb)+offset;
        Z(ind) = Z(ind) + 1;
    end
end


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fs = 1000; N = 5000;  mu = 0.05;  maxlag = 1;
% intervals = exprnd(mu,2,N) + 2/Fs;     lags = [-maxlag:1/Fs:maxlag];
% events1 = cumsum(intervals(1,:),2);  events1(events1>(mu*N)) = [];  events1 = round(events1*Fs);
% events2 = cumsum(intervals(2,:),2);  events2(events2>(mu*N)) = [];  events2 = round(events2*Fs);
% tic;  cross  = CORE_pxcorr(events1,events2,ceil(maxlag*Fs));  t(1) = toc;
% tic;  cross2 = xcorr(double(rasterize(events1)), double(rasterize(events2)),ceil(maxlag*Fs)); t(2) = toc;
% printf('\nCORE_pxcorr took %5.3f sec and equivalent native code took %5.3f sec.', t(1), t(2));
% printf('The MSE between the two results was %6.4f.', mean((cross-cross2).^2));
