function z=backtr(y,fit)
%
% Inverse-link transformation of y.

fali = fit.fit_points.family_link;
link = fali(2);
switch link
  case 3
    z=y;
  case 4
    z=exp(y);
  case 5
    z = y;
    i = find(y<=0);
    if (length(i)>0)
      z(i) = exp(y(i))./(1+exp(y(i)));
    end;
    i = find(y>0);
    if (length(i)>0)
      z(i) = 1./(1+exp(-y(i)));
    end;
  case 6
    z=1./y;
  case 7
    z=y.*abs(y);
  case 8
    z=sin(y).^2;
  otherwise
     disp('Backtr: invalid link');
     z=y;
end;

return;
