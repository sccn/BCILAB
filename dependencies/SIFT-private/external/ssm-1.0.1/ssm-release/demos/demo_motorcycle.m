y       = load('data/motorcycle.dat')';
time    = y(1, :);
delta   = y(2, [2:end 1]);
y       = y(3, :);

%% Spline smoothing model %%
spline              = ssm_spline(delta);
spline              = estimate(y, spline, [1 0.1]);
[alphahat V]        = statesmo(y, spline);
conf                = squeeze(1.96*realsqrt(V(1, 1, :)))';
[eps eta epsvar]    = disturbsmo(y, spline);

figure('Name', 'Motorcycle acceleration data analyzed by a cubic spline');
subplot(2, 1, 1), plot(time, alphahat(1, :), 'b'), hold all, plot(time, [alphahat(1, :)+conf; alphahat(1, :)-conf], 'b:'), scatter(time, y, 10, 'r', 's', 'filled'), hold off, title('Spline and 95% confidence intervals'), ylim([-140 80]), set(gca,'YGrid','on');
subplot(2, 1, 2), scatter(time, eps./realsqrt(epsvar), 10, 'r', 's', 'filled'), title('Standardized irregular'), set(gca,'YGrid','on');
if ispc, set(gcf, 'WindowStyle', 'docked'); end

%% Continuous local level model %%
abseps          = abs(eps);
contllm         = [ssm_gaussian ssmodel('continuous local level', 0, 1, 1, 1, ssmat(0, [], true, zeros(size(delta)), true), 'Qd', {@(X) exp(2*X)*delta}, {[]}, ssparam({'zeta var'}, '1/2 log'))];
[contllm logL]  = estimate(abseps, contllm, [1 0.1]);
alphahat        = statesmo(abseps, contllm);

figure('Name', 'Correction for heteroscedasticity');
subplot(3, 1, 1), plot(time, alphahat, 'b'), hold all, scatter(time, abseps, 10, 'r', 's', 'filled'), hold off, title('Absolute smoothed irregular and h^\ast_t'), ylim([0 87.5]);

%% Correction for heteroscedasticity %%
h2                  = (alphahat/alphahat(1)).^2;
splineh             = [ssmodel('Heteroscedastic noise', ssmat(0, [], true, zeros(size(h2)), true), zeros(1, 0), [], [], [], 'Hd', {@(X) exp(2*X)*h2}, {[]}, ssparam({'epsilon var'}, '1/2 log')) spline];
splineh             = estimate(y, splineh, [1 0.1]);
[alphahath Vh]      = statesmo(y, splineh);
confh               = squeeze(1.96*realsqrt(Vh(1, 1, :)))';
[epsh eta epsvarh]  = disturbsmo(y, splineh);

subplot(3, 1, 2), plot(time, alphahath(1, :), 'b'), hold all, plot(time, [alphahath(1, :)+confh; alphahath(1, :)-confh], 'b:'), scatter(time, y, 10, 'r', 's', 'filled'), hold off, title('Spline and 95% confidence intervals'), ylim([-140 80]), set(gca,'YGrid','on');
subplot(3, 1, 3), scatter(time, epsh./realsqrt(epsvarh), 10, 'r', 's', 'filled'), title('Standardized irregular'), set(gca,'YGrid','on');
if ispc, set(gcf, 'WindowStyle', 'docked'); end
