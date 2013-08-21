% Local Regression and Likelihood, Figure 2.3.
%
% Using different smoothing parameters.

load ethanol;

alp = [0.2 0.4 0.6 0.8];
f=['a','b','c','d'];
for (i=1:4)
  figure('Name', sprintf( 'figure2_3%c: Changing smoothing parameters', f(i) ) );
  fit = locfit(E,NOx,'alpha',alp(i));
  lfplot(fit);
  title(strcat('a = ',num2str(alp(i))));
end;
