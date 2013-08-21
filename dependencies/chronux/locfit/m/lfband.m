function lfband(fit,varargin)

% adds confidence bands around the plot of a locfit() fit.
%
% for 2-d fits, produces separate surface plots of upper and
% lower confidence limits.
%
% Bands are based on 95% pointwise coverage, using a single
% (i.e. global) estimate of sigma^2.

xfit = lfmarg(fit);
% placing 'band','g' before varargin{:} ensures that
% user-provided 'band' has precedence.
ypp = predict(fit,xfit,'band','g',varargin{:});
yfit = ypp{1};
se = ypp{2};
bands = ypp{3};

data = fit.data;
xdata = data.x;
p = size(xdata,2);
cv = 1.96;
fali = fit.fit_points.family_link;
cl = invlink(bands(:,1),fali);
cu = invlink(bands(:,2),fali);

if (p==1)
  hold on;
  plot(xfit{1},cu,':');
  plot(xfit{1},cl,':');
  hold off;
end;

if (p==2)
  x1 = xfit{1};
  x2 = xfit{2};
  figure(1);
  surf(x1,x2,reshape(cl,length(x1),length(x2))');
  figure(2);
  surf(x1,x2,reshape(cu,length(x1),length(x2))');
end;

return;
