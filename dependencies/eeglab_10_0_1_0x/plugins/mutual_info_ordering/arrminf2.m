function [mo,ord] = arrminf(mi0,maxpass,fignum,ord0)
global n mi ord
mi = mi0;
[m,n] = size(mi);
if nargin < 4
    ord = 1:n;
else
    ord = ord0;
end
minimprov = 1;
still_changing = 1;
pass = 0;
if nargin < 2
    maxpass = 3;
end
if nargin < 3
    f = figure;
else
    f = fignum;
end
figure(f), imagesc(mi);%, colormap hot;

if 0 && nargin < 4
sorted = 0;
while (1 & sorted == 0)
    sorted = 1;
    for i = 1:(n-1)
       if sum(mi(:,i)) < sum(mi(:,i+1))
           doswap(i,i+1);
           %tmp = ord(i);
           %ord(i) = ord(i+1);
           %ord(i+1) = tmp;
           sorted = 0;
       end     
    end
    %numpasses = numpasses + 1;
    %disp(['finished pass ' int2str(numpasses)]);
end
end
%odg = inf*ones(n,n);
%for i = 1,n
%    odg(i,i) = 0;
%end
while still_changing & pass < maxpass
    pass = pass+1;
    %disp(['pass ' int2str(pass)]);
    still_changing = 0;
    for k = 1:n
        % find the swap that minimizes new offdiagonality
        ok = offdiag(k,k);
        odg(k,k) = 0;
        for t = 1:n
            ot = offdiag(t,t);
            for s = 1:n
                %if s ~= t
                    % get the change in offdiagonality of these swaps
                    if ~(s == k & t == k)
                        odg(t,s) = offdiag(k,s) + offdiag(t,k) + offdiag(s,t) - ok - ot - offdiag(s,s);
                    end
                %end
            end
        end
        [mn,indi] = min(odg);
        [mn2,indj] = min(mn);
        if ~(indi(indj) == k & indj == k)
            doswap(k,indi(indj));
            doswap(indi(indj),indj);
            still_changing = 1;
        end
        figure(f), imagesc(mi);
        pause(0.3);
    end
    cost(pass) = 0;
    for kk = 1:n
       cost(pass) = cost(pass) + offdiag(kk,kk);
    end 
    disp(['cost = ' num2str(cost(pass))]);

    %figure(f), imagesc(mi);
    %pause(1);
end

mo = mi;


function od = offdiag(i,r)
global n mi
% off diagonality of ith row/col in rth positiion
od = 0;
for j = 1:n
    if j ~= i
        od = od + mi(i,j) * abs(j-r)^0.9  / r^0;
    else
        od = od + mi(i,r) * abs(j-r)^0.9  / r^0;
    end
end


function doswap(i,j)
global mi ord
swaprow(i,j);
swapcol(i,j);
tmp = ord(i);
ord(i) = ord(j);
ord(j) = tmp;

function swaprow(i,j)
global mi
tmp = mi(i,:);
mi(i,:) = mi(j,:);
mi(j,:) = tmp;

function swapcol(i,j)
global mi
tmp = mi(:,i);
mi(:,i) = mi(:,j);
mi(:,j) = tmp;