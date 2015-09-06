function plot_lc(rho,eta,marker,ps,reg_param)
%PLOT_LC Plot the L-curve.
%
% plot_lc(rho,eta,marker,ps,reg_param)
%
% Plots the L-shaped curve of the solution norm
%    eta = || x ||      if   ps = 1
%    eta = || L x ||    if   ps = 2
% as a function of the residual norm rho = || A x - b ||.  If ps is
% not specified, the value ps = 1 is assumed.
%
% The text string marker is used as marker.  If marker is not
% specified, the marker '-' is used.
%
% If a fifth argument reg_param is present, holding the regularization
% parameters corresponding to rho and eta, then some points on the
% L-curve are identified by their corresponding parameter.

% Per Christian Hansen, IMM, 12/29/97.

% Set defaults.
if (nargin==2), marker = '-'; end  % Default marker.
if (nargin < 4), ps = 1; end       % Std. form is default.
np = 10;                           % Number of identified points.

% Initialization.
if (ps < 1 | ps > 2), error('Illegal value of ps'), end
n = length(rho); ni = round(n/np);

% Make plot.
loglog(rho(2:end-1),eta(2:end-1)), ax = axis;
if (max(eta)/min(eta) > 10 | max(rho)/min(rho) > 10)
  if (nargin < 5)
    loglog(rho,eta,marker), axis(ax)
  else
    loglog(rho,eta,marker,rho(ni:ni:n),eta(ni:ni:n),'x'), axis(ax)
    HoldState = ishold; hold on;
    for k = ni:ni:n
      text(rho(k),eta(k),num2str(reg_param(k)));
    end
    if (~HoldState), hold off; end
  end
else
  if (nargin < 5)
    plot(rho,eta,marker), axis(ax)
  else
    plot(rho,eta,marker,rho(ni:ni:n),eta(ni:ni:n),'x'), axis(ax)
    HoldState = ishold; hold on;
    for k = ni:ni:n
      text(rho(k),eta(k),num2str(reg_param(k)));
    end
    if (~HoldState), hold off; end
  end
end
xlabel('residual norm || A x - b ||_2')
if (ps==1)
  ylabel('solution norm || x ||_2')
else
  ylabel('solution semi-norm || L x ||_2')
end
title('L-curve')