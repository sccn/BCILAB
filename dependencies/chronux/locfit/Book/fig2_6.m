% Local Regression and Likelihood, Figure 2.6.
%
% Influence and variance functions.

load ethanol;
d = [2 3];
a = [0.49 0.59];
main = {'Local Quadratic','Local Cubic'};
f=['a','b'];
for (i=1:2)
  figure('Name', sprintf( 'figure2_6%c: Influence and variance functions', f(i) ) );
  fit = locfit(E,NOx,'nn',a(i),'deg',d(i),'ev','grid','mg',100);
  vecr(fit)
  plot(E,zeros(88,1),'o');
  title(main{i});
  hold on;
  lfplot(fit,'what','infl');
  lfplot(fit,'what','vari','--');
  xlabel('Equivalence Ratio');
  ylabel('Influence');
  legend('Data','Influence Function','Variance Function');
  hold off;
end;
