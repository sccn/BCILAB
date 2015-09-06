y           = load('data/sv.dat')';
fprintf(1, '\n');

%% Analysis with linear Gaussian approximation %%
ylg         = reallog(y.^2);
asv         = [ssm_gaussian ssm_const ssm_arma(1, 0)];
asv.name    = 'Approximate stochastic volatility model';
[asv logL]  = estimate(ylg, asv, [var(ylg) 0.5 0.1*var(ylg)], [], 'fmin', 'bfgs', 'disp', 'off');
alphahat    = statesmo(ylg, asv);
ylgcom      = signal(alphahat, asv);
hhat        = ylgcom(2, :);

fprintf(1, '[Approximate stochastic volatility model]\n');
fprintf(1, 'loglikelihood: %g\n', logL);
fprintf(1, 'Irr variance:  %g\n', asv.param(1));
fprintf(1, 'AR variance:   %g\n', asv.param(3));
fprintf(1, 'AR(1) phi:     %g\n\n', asv.param(2));

figure('Name', 'Estimated volatility for the Pound series'), plot(exp(hhat/2)), ylim([0.52 1.84]), title('Estimated volatility for the Pound series');
drawnow;

%% Analysis with zero-mean stochastic volatility model %%
zmsv                = [ssm_zmsv ssm_arma(1, 0)];
zmsv.name           = 'Zero-mean stochastic volatility model';
randn('state', [1252786190; 3322582042]);
[zmsv logL output]  = estimate(y, zmsv, [0.6392 0.9774 0.02358]);

fprintf(1, '[Zero-mean stochastic volatility model]\n');
fprintf(1, 'loglikelihood: %g\n', logL);
fprintf(1, 'Sigma:         %g\n', zmsv.param(1));
fprintf(1, 'AR variance:   %g\n', zmsv.param(3));
fprintf(1, 'AR(1) phi:     %g\n', zmsv.param(2));

fprintf(1, '\n');
