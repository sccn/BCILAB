load('data/fMRItimeseries');
y               = y';
ft              = frametimes;
X               = X';
regar           = [ssm_arma(1, 0) ssm_reg(X)];
[regar logL]    = estimate(y, regar, [0.2 20], [], 'tol', 10^-4);
[alphahat V]    = statesmo(y, regar);
ycom            = signal(alphahat, regar);
fprintf(1, '[RegAR fit of fMRI time series]\n');
fprintf(1, 'Loglikelihood: %g\n', logL);
fprintf(1, 'AR phi: %g\n', regar.param(1));
fprintf(1, 'AR variance: %g\n', regar.param(2));
fline   = '  %-16.4g%-16.4g\n';
fprintf(1, fline, [alphahat(2:end, end) alphahat(2:end, end)./realsqrt(diag(V(2:end, 2:end, end)))]');

figure('Name', 'fMRI fit');
subplot(1, 2, 1), plot(ft, y, 'DisplayName', 'fMRI data'), hold all, plot(ft, ycom(2, :), 'DisplayName', 'Regression effect'), title('original fMRI series'), legend('show');
subplot(1, 2, 2), plot(ft, ycom(1, :)), title('fMRI series w/ regression effects removed');
