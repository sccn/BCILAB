function lfplot(varargin)
 
% Plot (for one or two dimensions) a locfit() fit.
%
% Usage:
%   fit = locfit(x,y);
%   lfplot(fit);
%
% Plot the fitted smooth curve, and add a scatterplot of the data.
%
% Required argument:
%   fit  (produced by locfit()).
%
% Optional arguments:
%   'nodata'  - don't add data to plot.
%   'contour' - for 2-d predictors, use contour instead of surf.
%   'direct'  - fit directly, instead of using interpolation
%               (see the predict() function).
%   'what'    - locfit what argument ('coef', 'infl', 'vari', 'band' etc).
%   Any additional arguments are passed to Matlab's plot(), contour()
%     or surf() function as appropriate.
%
% To add confidence bands, use the lfband() function.
%
% Author: Catherine Loader.

fit = varargin{1};
data = fit.data;
xdata = data.x;
n = size(xdata,1);
d = size(xdata,2);
fali = fit.fit_points.family_link;
ydata = data.y;
wdata = data.weights;
cdata = data.censor;
if (length(cdata)==1) cdata = zeros(n,1); end;
showdata = (fit.evaluation_structure.derivative==0);
ppargs = {};
plotargs = {};

type = 's';
na = 2;
while na <= length(varargin)
  inc = 0;
  if (strcmp(varargin{na},'contour'))
    type = 'c';
    inc = 1;
  end;
  if (strcmp(varargin{na},'what'))
    ppargs = {ppargs{:}, 'what', varargin{na+1}};
    showdata = 0;
    inc = 2;
  end;
  if (strcmp(varargin{na},'nodata'))
    showdata = 0;
    inc = 1;
  end;
  if (strcmp(varargin{na},'direct'))
    ppargs = {ppargs{:} 'direct'};
    inc = 1;
  end;
  if (inc==0)
    plotargs = {plotargs{:} varargin{na}};
    inc = 1;
  end;
  na = na+inc;
end;

xfit = lfmarg(fit);
yfit = predict(fit,xfit,ppargs{:});
yfit = invlink(yfit,fali);
fam = mod(fali(1),64);
if (fam>4)
  ydata = ydata ./ wdata;
end;

if (d==1)
  plot(xfit{1},yfit,plotargs{:});
  if (showdata)
    hold on;
    if (length(ydata)==1) ydata = zeros(n,1); end;
    plotbyfactor(xdata,ydata,cdata);
    hold off;
  end;
end;

if (d==2)
  x1 = xfit{1};
  x2 = xfit{2};
  yfit = reshape(yfit,length(x1),length(x2));
  if (type=='c')
    [C h] = contour(x1,x2,yfit',plotargs{:});
    clabel(C,h);
    if (showdata)
      hold on;
      plotbyfactor(xdata(:,1),xdata(:,2),cdata);
      hold off;
    end;
  else
    surf(x1,x2,yfit',plotargs{:});
  end;
end;

return;
