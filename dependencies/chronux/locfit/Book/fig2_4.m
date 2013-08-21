% Local Regression and Likelihood, Figure 2.4.
%
% Bivariate Local Regression.

load ethanol;

alp = [0.25 0.3 0.49 0.59];
lab = {'Local Constant' 'Local Linear' 'Local Quadratic' 'Local Cubic'};
f=['a','b','c','d'];
for (i=1:4)
  figure('Name', sprintf( 'figure2_4%c: Bivariate local regression', f(i) ) );
    fit = locfit(E,NOx,'deg',i-1,'alpha',alp(i));
    lfplot(fit);
    title(lab(i));
end;
