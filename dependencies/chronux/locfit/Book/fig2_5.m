% Local Regression and Likelihood, Figure 2.5.
%
% Using different smoothing parameters.

load ethanol;

alp = [0.2 0.4 0.6 0.8];
f=['a','b','c','d'];
for (i=1:4)
  figure('Name', sprintf( 'figure2_5%c: Smoothing parameters', f(i) ) );
  fit = locfit(E,NOx,'alpha',alp(i));
  res = residuals(fit);
  fit2 = locfit(E,res,'alpha',0.2);
  lfplot(fit2);
  title(strcat('a = ',num2str(alp(i))));
  ylabel('Residual');
  hold on;
  plot([min(E) max(E)], [0 0], ':');
  hold off;
end;
