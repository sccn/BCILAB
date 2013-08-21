function z=smooth_lf(x,y,varargin)

% must (unlike R smooth.lf() function) give x and y.
% also R's direct=T is automatic.
%

xev = x;
if (k>1)
  if (strcmp(varargin{1},'xev'))
    xev = varargin{2};
    varargin(1:2) = [];
  end;
end;
fit = locfit(x,y,varargin{:},'ev',xev,'module','simple');
z = lfknots(fit);
fv = invlink(z(:,1),fit.fit_points.family_link);

z = { xev, fv };

return;
