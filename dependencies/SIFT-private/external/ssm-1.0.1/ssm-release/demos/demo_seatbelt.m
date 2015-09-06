% Load data
seatbelt    = load('data/seatbelt.dat')';
time        = 69:1/12:85;
fprintf(1, '\n');

%% Analysis of drivers series %%
y           = seatbelt(1, :);

% Estimation of basic structural time series model
bstsm       = ssm_stsm('level', 'trig1', 12);
bstsm       = estimate(y, bstsm, [0.003 0.0009 5e-7], [], 'fmin', 'bfgs', 'disp', 'off');
[a P]       = kalman(y, bstsm);
alphahat    = statesmo(y, bstsm);

% Retrieve components
ycom        = signal(a, bstsm);
lvl         = ycom(1, :);
seas        = ycom(2, :);
ycomhat     = signal(alphahat, bstsm);
lvlhat      = ycomhat(1, :);
seashat     = ycomhat(2, :);
[irr etahat epsvarhat etavarhat] = disturbsmo(y, bstsm);

figure('Name', 'Estimated Components');
subplot(3, 1, 1), plot(time(1:end-1), y, 'r:', 'DisplayName', 'drivers'), hold all, plot(time(1:end-1), lvlhat, 'DisplayName', 'est. level'), hold off, title('Level'), ylim([6.875 8]), legend('show');
subplot(3, 1, 2), plot(time(1:end-1), seashat), title('Seasonal'), ylim([-0.16 0.28]);
subplot(3, 1, 3), plot(time(1:end-1), irr), title('Irregular'), ylim([-0.15 0.15]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

figure('Name', 'Data and level');
plot(time, lvl, 'DisplayName', 'filtered level'); hold all; plot(time(1:end-1), lvlhat, ':', 'DisplayName', 'smoothed level');
scatter(time(1:end-1), y, '+', 'DisplayName', 'drivers'); hold off; ylim([6.95 7.9]), legend('show');
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% Calculate standardized residuals
u       = irr./realsqrt(epsvarhat);
r       = zeros(12, size(y, 2));
for t = 1:size(y, 2), r(:, t) = inv(sqrtm(etavarhat(:,:, t)))*etahat(:, t); end
comres  = signal(r, bstsm);
lvlres  = comres(1, :);
figure('Name', 'Residuals');
subplot(3, 1, 1), plot(time(1:end-1), y - lvl(1:end-1) - seas(1:end-1)), title('One-step ahead prediction residuals'), xlim([68 86]), ylim([-0.35 0.25]);
subplot(3, 1, 2), plot(time(1:end-1), u), title('Auxiliary irregular residuals'), xlim([68 86]), ylim([-4.5 4.5]);
subplot(3, 1, 3), plot(time(1:end-1), lvlres), title('Auxiliary level residuals'), xlim([68 86]), ylim([-2.5 1.5]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% Adding explanatory variables and intervention to the model
petrol              = seatbelt(5, :);
bstsmir             = [bstsm ssm_intv(size(y, 2), 'step', 170) ssm_reg(petrol, 'petrol')];
[bstsmir logL]      = estimate(y, bstsmir, [0.004 0.00027 1e-6]);
[alphahatir Vir]    = statesmo(y, bstsmir);
irrir               = disturbsmo(y, bstsmir);
ycomir              = signal(alphahatir, bstsmir);
lvlir               = sum(ycomir([1 3 4], :), 1);
seasir              = ycomir(2, :);

fprintf(1, '[Analysis on drivers series]\n');
fprintf(1, 'Loglikelihood: %g\n', logL);
fprintf(1, 'Irregular variance: %g\n', bstsmir.param(1));
fprintf(1, 'Level variance: %g\n', bstsmir.param(2));
fprintf(1, 'Seasonal variance: %g\n', bstsmir.param(3));
fprintf(1, 'Variable             Coefficient     R. m. s. e.     t-value\n');
fprintf(1, 'petrol coefficient   %-14.5g  %-14.5g  %g\n', alphahatir(end, 1), realsqrt(Vir(end, end, end)), alphahatir(end, 1)/realsqrt(Vir(end, end, end)));
fprintf(1, 'level shift at 83.2  %-14.5g  %-14.5g  %g\n\n', alphahatir(end-1, 1), realsqrt(Vir(end-1, end-1, end)), alphahatir(end-1, 1)/realsqrt(Vir(end-1, end-1, end)));

figure('Name', 'Estimated Components w/ intervention and regression');
subplot(3, 1, 1), plot(time(1:end-1), y, 'r:', 'DisplayName', 'drivers'), hold all, plot(time(1:end-1), lvlir, 'DisplayName', 'est. level'), hold off, title('Level'), ylim([6.875 8]), legend('show');
subplot(3, 1, 2), plot(time(1:end-1), seasir), title('Seasonal'), ylim([-0.16 0.28]);
subplot(3, 1, 3), plot(time(1:end-1), irrir), title('Irregular'), ylim([-0.15 0.15]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end
drawnow;

% [irr etahat epsvarhat etavarhat] = disturbsmo(y, bstsmir);
% r       = zeros(12, size(y, 2));
% for t = 1:size(y, 2), r(:, t) = inv(sqrtm(etavarhat(:,:, t)))*etahat(:, t); end
% comres  = signal(r, bstsm);
% lvlres  = comres(1, :);
% figure('Name', 'Estimated Components w/ intervention and regression');
% subplot(2, 1, 1), plot(time(1:end-1), y, 'r-.', 'DisplayName', 'drivers'), hold all, plot(time(1:end-1), lvlir, 'DisplayName', 'est. level'), hold off, title('Level'), xlim([68 86]), ylim([6.875 8]);
% subplot(2, 1, 2), plot(time(1:end-1), seasir), title('Seasonal'), xlim([68 86]), ylim([-0.16 0.28]);
% figure('Name', 'Estimated Components w/ intervention and regression');
% subplot(2, 1, 1), plot(time(1:end-1), irrir), title('Irregular'), xlim([68 86]), ylim([-0.15 0.15]);
% subplot(2, 1, 2), plot(time(1:end-1), lvlres), title('Normalized level residuals'), xlim([68 86]), ylim([-1.5 1]);

%% Analysis of both front and rear seat passengers bivariate series %%
y2          = seatbelt(2:3, :);

% Bivariate basic structural time series model with regression variables
% petrol and kilometer travelled, before intervention
bibstsm         = ssm_mvstsm(2, [true true false], 'level', 'trig fixed', 12, false, seatbelt(4:5, :));
[bibstsm logL]  = estimate(y2(:, 1:169), bibstsm, [0.00531 0.0083 0.00441 0.000247 0.000229 0.000218], [], 'fmin', 'bfgs', 'disp', 'off');
fprintf(1, '[Parameters estimated w/o intervention on front and rear seat bivariate series]\n');
fprintf(1, 'Loglikelihood: %g.\n', logL);
fprintf(1, 'Irregular disturbance   Level disturbance\n');
fline   = '%-10.5g  %-10.5g  %-10.5g  %-10.5g\n';
fprintf(1, fline, bibstsm.param([1 3 4 6]));
fprintf(1, fline, bibstsm.param([3 2 6 5]));
fprintf(1, '\n');

alphahat    = statesmo(y2(:, 1:169), bibstsm);
comhat      = signal(alphahat, bibstsm);
lvlhat      = comhat(:, :, 1);
seashat     = comhat(:, :, 2);
reghat      = comhat(:, :, 3);
figure('Name', 'Estimated components w/o intervention on front and rear seat bivariate series');
subplot(2, 2, 1), plot(time(1:169), lvlhat(1, :)+reghat(1, :)), hold all, scatter(time(1:169), y2(1, 1:169), 8, 'r', 's', 'filled'), hold off, title('Front seat passenger level (w/o seasonal)'), xlim([68 85]), ylim([6 7.25]);
subplot(2, 2, 2), plot(time(1:169), lvlhat(1, :)), title('Front seat passenger level'), xlim([68 85]),% ylim([3.84 4.56]);
subplot(2, 2, 3), plot(time(1:169), lvlhat(2, :)+reghat(2, :)), hold all, scatter(time(1:169), y2(2, 1:169), 8, 'r', 's', 'filled'), hold off, title('Rear seat passenger level (w/o seasonal)'), xlim([68 85]), ylim([5.375 6.5]);
subplot(2, 2, 4), plot(time(1:169), lvlhat(2, :)), title('Rear seat passenger level'), xlim([68 85]),% ylim([1.64 1.96]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% Add intervention to both series
bibstsm2i           = [bibstsm ssm_mvintv(2, size(y2, 2), 'step', 170)];
[bibstsm2i logL2i]  = estimate(y2, bibstsm2i, [0.0054 0.00857 0.00445 0.000256 0.000232 0.000225], [], 'fmin', 'bfgs', 'disp', 'off');
[alphahat2i V2i]    = statesmo(y2, bibstsm2i);
fprintf(1, '[Parameters estimated w/ intervention on both series]\n');
fprintf(1, 'Loglikelihood: %g.\n', logL2i);
fprintf(1, 'Irregular disturbance   Level disturbance\n');
fprintf(1, fline, bibstsm2i.param([1 3 4 6]));
fprintf(1, fline, bibstsm2i.param([3 2 6 5]));
fprintf(1, 'Level shift intervention:\n');
fprintf(1, '        Coefficient     R. m. s. e.     t-value\n');
fprintf(1, 'front   %-14.5g  %-14.5g  %g\n', alphahat2i(end-1, end), realsqrt(V2i(end-1, end-1, end)), alphahat2i(end-1, end)/realsqrt(V2i(end-1, end-1, end)));
fprintf(1, 'rear    %-14.5g  %-14.5g  %g\n\n', alphahat2i(end, end), realsqrt(V2i(end, end, end)), alphahat2i(end, end)/realsqrt(V2i(end, end, end)));

% Add intervention only to front seat passenger series
bibstsmi            = [bibstsm ssm_mvintv(2, size(y2, 2), {'step' 'null'}, 170)];
[bibstsmi logLi]    = estimate(y2, bibstsmi, [0.00539 0.00856 0.00445 0.000266 0.000235 0.000232]);
[alphahati Vi]      = statesmo(y2, bibstsmi);
fprintf(1, '[Parameters estimated w/ intervention only on front seat series]\n');
fprintf(1, 'Loglikelihood: %g.\n', logLi);
fprintf(1, 'Irregular disturbance   Level disturbance\n');
fprintf(1, fline, bibstsmi.param([1 3 4 6]));
fprintf(1, fline, bibstsmi.param([3 2 6 5]));
fprintf(1, 'Level shift intervention:\n');
fprintf(1, '        Coefficient     R. m. s. e.     t-value\n');
fprintf(1, 'front   %-14.5g  %-14.5g  %g\n\n', alphahati(end, end), realsqrt(Vi(end, end, end)), alphahati(end, end)/realsqrt(Vi(end, end, end)));

comhati     = signal(alphahati, bibstsmi);
lvlhati     = comhati(:, :, 1);
seashati    = comhati(:, :, 2);
reghati     = comhati(:, :, 3);
intvhati    = comhati(:, :, 4);
figure('Name', 'Estimated components w/ intervention only on front seat series');
subplot(2, 2, 1), plot(time(1:end-1), lvlhati(1, :)+reghati(1, :)+intvhati(1, :), 'DisplayName', 'est. level'), hold all, scatter(time(1:end-1), y2(1, :), 8, 'r', 's', 'filled', 'DisplayName', 'front seat'), hold off, title('Front seat passenger level (w/o seasonal)'), xlim([68 85]), ylim([6 7.25]), legend('show');
subplot(2, 2, 2), plot(time(1:end-1), lvlhati(1, :)+intvhati(1, :)), title('Front seat passenger level'), xlim([68 85]),% ylim([3.84 4.56]);
subplot(2, 2, 3), plot(time(1:end-1), lvlhati(2, :)+reghati(2, :), 'DisplayName', 'est. level'), hold all, scatter(time(1:end-1), y2(2, :), 8, 'r', 's', 'filled', 'DisplayName', 'rear seat'), hold off, title('Rear seat passenger level (w/o seasonal)'), xlim([68 85]), ylim([5.375 6.5]), legend('show');
subplot(2, 2, 4), plot(time(1:end-1), lvlhati(2, :)), title('Rear seat passenger level'), xlim([68 85]),% ylim([1.64 1.96]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end

% figure('Name', 'Estimated components w/ intervention only on front seat series');
% subplot(2, 1, 1), plot(time(1:end-1), lvlhati(1, :)+reghati(1, :)+intvhati(1, :), 'DisplayName', 'est. level'), hold all, scatter(time(1:end-1), y2(1, :), 8, 'r', 's', 'filled', 'DisplayName', 'front seat'), hold off, title('Front seat passenger level (w/o seasonal)'), xlim([68 85]), ylim([6 7.25]);
% subplot(2, 1, 2), plot(time(1:end-1), lvlhati(2, :)+reghati(2, :), 'DisplayName', 'est. level'), hold all, scatter(time(1:end-1), y2(2, :), 8, 'r', 's', 'filled', 'DisplayName', 'rear seat'), hold off, title('Rear seat passenger level (w/o seasonal)'), xlim([68 85]), ylim([5.375 6.5]);
% if ispc, set(gcf, 'WindowStyle', 'docked'); end
