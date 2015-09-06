% Load data
y       = load('data/nile.dat')';
time    = 1871:1970;
fprintf(1, '\n');

% Maximum loglikelihood estimation
llm         = estimate(y, ssm_llm, [10000 5000]);
[logL fvar] = loglik(y, llm);
fprintf(1, 'Loglikelihood = %g, variance = %g.\n', logL, fvar);
fprintf(1, 'epsilon variance = %g, eta variance = %g.\n', llm.param(1), llm.param(2));

% Kalman filtering
[a P v F]   = kalman(y, llm);
figure('Name', 'Filtered state');
subplot(2, 2, 1), plot(time, y, 'r:', 'DisplayName', 'nile'), hold all, plot([time 1971], a, 'b', 'DisplayName', 'filt. state'), plot([time 1971], [a+1.645*sqrt(P); a-1.645*sqrt(P)], 'g:', 'DisplayName', {'90% conf. +' '90% conf. -{}'}), hold off, title('Filtered state'), ylim([450 1400]), legend('show');
subplot(2, 2, 2), plot([time 1971], P), title('Filtered state variance'), ylim([5000 17500]);
subplot(2, 2, 3), plot(time, v), title('Prediction errors'), ylim([-450 400]);
subplot(2, 2, 4), plot(time, F), title('Prediction error variance'), ylim([20000 32500]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% State smoothing
[alphahat V r N]    = statesmo(y, llm);
figure('Name', 'Smoothed state');
subplot(2, 2, 1), plot(time, y, 'r:', 'DisplayName', 'nile'), hold all, plot(time, alphahat, 'DisplayName', 'smo. state'), plot(time, [alphahat+1.645*sqrt(V); alphahat-1.645*sqrt(V)], 'g:', 'DisplayName', {'90% conf. +' '90% conf. -'}), hold off, title('Smoothed state'), ylim([450 1400]), legend('show');
subplot(2, 2, 2), plot(time, V), title('Smoothed state variance'), ylim([2300 4100]);
subplot(2, 2, 3), plot(time, r), title('Smoothing cumulant'), ylim([-0.036 0.024]);
subplot(2, 2, 4), plot(time, N), title('Smoothing variance cumulant'), ylim([0 0.000105]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% Disturbance smoothing
[epshat etahat epsvarhat etavarhat] = disturbsmo(y, llm);
figure('Name', 'Smoothed disturbances');
subplot(2, 2, 1), plot(time, epshat), title('Observation error'), ylim([-360 280]);
subplot(2, 2, 2), plot(time, epsvarhat), title('Observation error variance'), ylim([2300 4100]);
subplot(2, 2, 3), plot(time, etahat), title('State error'), ylim([-50 35]);
subplot(2, 2, 4), plot(time, etavarhat), title('State error variance'), ylim([1225 1475]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% Simulation smoothing
[alphatilde epstilde etatilde alphaplus] = simsmo(y, llm, 1);
figure('Name', 'Simulation');
subplot(2, 2, 1), plot(time, alphahat, 'DisplayName', 'samp. state'), hold all, scatter(time, alphaplus+alphahat(1)-alphaplus(1), 8, 'r', 's', 'filled', 'DisplayName', 'nile'), hold off, title('Unconditioned sampled state'), legend('show');
subplot(2, 2, 2), plot(time, alphahat, 'DisplayName', 'disp. samp. state'), hold all, scatter(time, alphatilde, 8, 'r', 's', 'filled', 'DisplayName', 'nile'), hold off, title('Conditioned sampled state'), ylim([740 1160]), legend('show');
subplot(2, 2, 3), plot(time, epshat, 'DisplayName', 'smo. obs. disturb.'), hold all, scatter(time, epstilde, 8, 'r', 's', 'filled', 'DisplayName', 'samp. obs. disturb.'), hold off, title('Conditioned sampled observation error'), ylim([-440 280]), legend('show');
subplot(2, 2, 4), plot(time, etahat, 'DisplayName', 'smo. state disturb.'), hold all, scatter(time, etatilde, 8, 'r', 's', 'filled', 'DisplayName', 'samp. state disturb.'), hold off, title('Conditioned sampled state error'), ylim([-440 280]), legend('show');
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% Missing Observations
ymis                = y;
ymis([21:40 61:80]) = NaN;
[amis Pmis]         = kalman(ymis, llm);
[alphahatmis Vmis]  = statesmo(ymis, llm);
figure('Name', 'Filtering and smoothing of data with missing observations');
subplot(2, 2, 1), plot(time, ymis, 'r:', 'DisplayName', 'nile w/ miss. values'), hold all, plot([time 1971], amis, 'DisplayName', 'filt. state'), hold off, title('Filtered state (extrapolation)'), ylim([450 1400]), legend('show');
subplot(2, 2, 2), plot([time 1971], Pmis), title('Filtered state variance'), ylim([4000 36000]);
subplot(2, 2, 3), plot(time, ymis, 'r:', 'DisplayName', 'nile w/ miss. values'), hold all, plot(time, alphahatmis, 'DisplayName', 'smo. state'), hold off, title('Smoothed state (interpolation)'), ylim([450 1400]), legend('show');
subplot(2, 2, 4), plot(time, Vmis), title('Filtered state (extrapolation)'), ylim([2000 10000]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% Forecasting (equivalent to future missing values)
yforc                       = [y repmat(NaN, 1, 50)];
[aforc Pforc vforc Fforc]   = kalman(yforc, llm);
figure('Name', 'Forecasting');
subplot(2, 2, 1), plot([time 1972:2021], yforc, 'r:', 'DisplayName', 'nile'), hold all, plot([time 1972:2022], aforc, 'DisplayName', 'forecast'), plot([time 1972:2022], [repmat(NaN, 2, length(time)) [aforc(end-50:end)+0.675*sqrt(Pforc(end-50:end)); aforc(end-50:end)-0.675*sqrt(Pforc(end-50:end))]], 'g:', 'DisplayName', {'50% conf. +' '50% conf. -'}), title('State forecast'), xlim([1868 2026]), ylim([450 1400]), legend('show');
subplot(2, 2, 2), plot([time 1972:2022], Pforc), title('State variance'), xlim([1868 2026]), ylim([4000 80000]);
subplot(2, 2, 3), plot([time 1972:2022], aforc), title('Observation forecast'), xlim([1868 2026]), ylim([700 1200]);
subplot(2, 2, 4), plot([time 1972:2021], Fforc), title('Observation forecast variance'), xlim([1868 2026]), ylim([20000 96000]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

fprintf(1, '\n');
