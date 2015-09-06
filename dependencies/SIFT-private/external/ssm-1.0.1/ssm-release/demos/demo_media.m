y           = load('data/M3_21.dat')';

%% Airline seasonal adjustment %%
air         = estimate(y, ssm_airline, 0.1);
aircom      = ssmhtd(air);
alphahat    = statesmo(y, aircom);
ycom        = signal(alphahat, aircom);
airseas     = ycom(2, :);
% figure('Name', 'Airline model seasonal adjustment');
% plot(y, 'r:', 'DisplayName', 'data'), hold all, plot(y - airseas, 'b', 'DisplayName', 'sa: airline'), hold off, ylim([640 1180]), legend('show');
% if ispc, set(gcf, 'WindowStyle', 'docked'); end

%% Generalized airline 3-5-1(3) seasonal adjustment %%
param0      = air.param([1 2 2 3]);
param0(1:3) = -param0(1:3);
param0(2:3) = param0(2:3).^(1/12);
ga          = estimate(y, ssm_genair(3, 5, 3), param0);
gacom       = ssmhtd(ga);
alphahat    = statesmo(y, gacom);
ycom        = signal(alphahat, gacom);
gaseas      = ycom(2, :);
% figure('Name', 'Generalized airline 3-5-1(3) model seasonal adjustment');
% plot(y, 'r:', 'DisplayName', 'data'), hold all, plot(y - gaseas, 'b', 'DisplayName', 'sa: ga 3-5-1(3)'), hold off, ylim([640 1180]), legend('show');
% if ispc, set(gcf, 'WindowStyle', 'docked'); end
figure('Name', 'Seasonal adjustment comparison');
plot(y, 'r:', 'DisplayName', 'data'), hold all, plot(y - airseas, 'g-.', 'DisplayName', 'sa: airline'), plot(y - gaseas, 'b', 'DisplayName', 'sa: ga 3-5-1(3)'), hold off, ylim([640 1180]), legend('show');
if ispc, set(gcf, 'WindowStyle', 'docked'); end

N       = length(y);
lambda  = (0:0.01:6)/12;
[airwtfil airwtsmo] = weights(aircom, N, floor(N/2), true);
[gawtfil gawtsmo]   = weights(gacom, N, floor(N/2), true);
clear('pi');
figure('Name', 'Concurrent Squared Gain');
plot(1 - abs(sum(diag(airwtsmo(2, :))*exp(2i*pi*(floor(N/2)-(1:N))'*lambda))).^2); hold all;
plot(1 - abs(sum(diag(gawtsmo(2, :))*exp(2i*pi*(floor(N/2)-(1:N))'*lambda))).^2); hold off;
if ispc, set(gcf, 'WindowStyle', 'docked'); end
drawnow;

%% Out-of-sample forecasts %%
% [yfair1 errair1 SSair1]     = oosforecast(y, air, 24, 1);
% [yfair12 errair12 SSair12]  = oosforecast(y, air, 24, 12);
% [yfga1 errga1 SSga1]        = oosforecast(y, ga, 24, 1);
% [yfga12 errga12 SSga12]     = oosforecast(y, ga, 24, 12);
% forcerrstat1    = 24*(SSga1 - SSair1)/SSair1(end);
% forcerrstat12   = 13*(SSga12 - SSair12)/SSair12(end);
[yfair errair SSair]    = oosforecast(y, air, 24, [1 12]);
[yfga errga SSga]       = oosforecast(y, ga, 24, [1 12]);
forcerrstat1    = 24*(SSga(1, :) - SSair(1, :))/SSair(1, end);
forcerrstat12   = 13*(SSga(2, :) - SSair(2, :))/SSair(2, end);
figure('Name', 'Out-of-sample forecast comparison');
plot(1:24, forcerrstat1, 'DisplayName', 'One step ahead'), hold all, plot(1:24, forcerrstat12, 'DisplayName', '12 steps ahead'), set(gca, 'YGrid', 'on'), legend('show'), title('Out-of-sample square forecast error diagnostic\newlineComparing 3-5-1(3) model with the airline model'), ylabel('Cumulative sums of squared errors');
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% %% Airline components with t-distribution %%
% airt        = ssm_sarimahtd(0, 1, 1, 0, 1, 1, 12, false);
% airt        = estimate(y, airt, [air.param 10], [], 'fmin', 'bfgs', 'disp', 'iter');
% alphahat    = statesmo(y, airt);
% ycom        = signal(alphahat, airt);
% airtseas    = ycom(2, :);
% figure('Name', 'Airline components model with t-distributed irregular seasonal adjustment');
% plot(y, 'r:', 'DisplayName', 'data'), hold all, plot(y - airtseas, 'b', 'DisplayName', 'sa: t airline'), hold off, ylim([640 1180]), legend('show');
% if ispc, set(gcf, 'WindowStyle', 'docked'); end
