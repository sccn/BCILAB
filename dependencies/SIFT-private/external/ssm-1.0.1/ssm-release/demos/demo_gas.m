y       = load('data/gas.dat')';
time    = 1960:1/4:1986+3/4;
fprintf(1, '\n');

%% Analysis based on Gaussian model %%
stsm        = ssm_stsm('trend', 'dummy', 4);
[stsm logL] = estimate(y, stsm, 0.01, [], 'fmin', 'bfgs', 'disp', 'off');
fprintf(1, '[Gaussian STSM]\n');
fprintf(1, 'loglikelihood:      %g\n', logL);
fprintf(1, 'Irregular variance: %g\n', stsm.param(1));
fprintf(1, 'Trend variance:     %g\n', stsm.param(2));
fprintf(1, 'Level variance:     %g\n', stsm.param(3));
fprintf(1, 'Seasonal variance:  %g\n\n', stsm.param(4));

[alpha irr] = fastsmo(y, stsm);
ycom        = signal(alpha, stsm);
seas        = ycom(2, :);
figure('Name', 'Gaussian model analysis');
subplot(2, 1, 1), plot(time, seas), title('Seasonal'), xlim([1959 1988]), ylim([-0.9 0.8]);
subplot(2, 1, 2), plot(time, irr), title('Gaussian Irregular'), xlim([1959 1988]), ylim([-0.4 0.4]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end
drawnow;

%% Analysis based on t-model %%
stsmt           = [ssm_t ssm_llt ssm_seasonal('dummy', 4)];
randn('state', [3765023265; 2369472656]);
[stsmt logLt]   = estimate(y, stsmt, [stsm.param(1) 4 stsm.param(2:end)], alpha, 'fmin', 'bfgs', 'disp', 'off');
[alphat irrt]   = fastsmo(y, stsmt);
ycomt           = signal(alphat, stsmt);
lvlt            = ycomt(1, :);
seast           = ycomt(2, :);

fprintf(1, '[t-dist error STSM]\n');
fprintf(1, 'loglikelihood:      %g\n', logLt);
fprintf(1, 't-dist variance:    %g\n', stsmt.param(1));
fprintf(1, 't-dist df:          %g\n', stsmt.param(2));
fprintf(1, 'Trend variance:     %g\n', stsmt.param(3));
fprintf(1, 'Level variance:     %g\n', stsmt.param(4));
fprintf(1, 'Seasonal variance:  %g\n', stsmt.param(5));

figure('Name', 'Component comparison between Gaussian and t-dist. models');
subplot(2, 1, 1), plot(time, seas, 'r', 'DisplayName', 'Gaussian'), hold all, plot(time, seast, 'b', 'DisplayName', 't-distribution'), hold off, title('Seasonal component'), xlim([1959 1988]), ylim([-0.9 0.8]);
subplot(2, 1, 2), plot(time, irr, 'r', 'DisplayName', 'Gaussian'), hold all, plot(time, irrt, 'b', 'DisplayName', 't-distribution'), hold off, title('Irregular component'), xlim([1959 1988]), ylim([-0.4 0.4]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

fprintf(1, '\n');
