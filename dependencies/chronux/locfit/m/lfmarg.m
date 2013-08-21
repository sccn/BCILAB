function xfit = lfmarg(fit)

% computes grid margins from a locfit object, used for plotting.
% this function is called from lfplot and lfband, will not normally
% be called directly by users.

box = fit.fit_points.fit_limits;
d = size(box,1);
xfit = cell(1,d);
mg = 10*ones(1,d);
if (d==1) mg = 200; end;
if (d==2) mg = [51 50]; end;
for i=1:d
    xfit{i} = box(i,1) + (box(i,2)-box(i,1))*(0:mg(i))'/mg(i);
end;

return;
