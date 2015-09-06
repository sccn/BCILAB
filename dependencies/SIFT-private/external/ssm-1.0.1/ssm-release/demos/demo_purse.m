y               = load('data/purse.dat')';
model           = [ssm_poisson ssm_llm];
randn('state', [1105946959; 3715058465]);
[model logL]    = estimate(y, model, exp(-6), [], 'fmin', 'bfgs', 'disp', 'iter');
fprintf(1, 'Loglikelihood: %g\nEta variance: %g\n', logL, model.param);

[alpha irr]     = fastsmo(y, model);
figure('Name', 'Purse data w/ poisson analysis');
plot(y, 'r:', 'DisplayName', 'Purse data'), hold all, plot(alpha, 'b', 'DisplayName', 'Estimated signal'), hold off; legend('show');
