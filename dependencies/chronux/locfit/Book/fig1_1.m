% Local Regression and Likelihood, Figure 1.1.
%
% A linear regression and plot for Spencer's mortality
% data. Unlike the book, I use locfit for this!
%
% Setting 'ev','none' means that locfit computes only
% the parametric fit. In this case, 'deg',1; i.e. linear.

load spencer;
fit = locfit(age,mortality,'deg',1,'ev','none');
figure('Name','fig1_1: locfit linear regression');
lfplot(fit);
