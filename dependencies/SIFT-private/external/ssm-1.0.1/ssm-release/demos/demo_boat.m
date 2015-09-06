y       = load('data/boat.dat')';
time    = 1829:2000;
fprintf(1, '\n');

%% Binary distribution with random walk probability %%
llmb                = [ssm_binary ssm_llm];
randn('state', [1105946959; 3715058465]);
[llmb logL output]  = estimate(y, llmb, 1, [], 'fmin', 'bfgs', 'disp', 'off');
fprintf(1, '[Binary distribution random walk model]\n');
fprintf(1, 'Loglikelihood: %g\n', logL);
fprintf(1, 'Eta variance: %g\n', llmb.param);
llmb.param  = 0.521;
[theta V]   = statesmo(output.ytilde, llmb);
pi          = exp(theta)./(1+exp(theta));
pi0         = exp(theta+0.675*realsqrt(V))./(1+exp(theta+0.675*realsqrt(V)));
pi1         = exp(theta-0.675*realsqrt(V))./(1+exp(theta-0.675*realsqrt(V)));

figure('Name', 'Probability of Cambridge win w/ 50% confidence');
scatter(time, y, 's', 'filled'), hold all, plot(time, pi, 'b'), plot(time, [pi0; pi1], 'g:'), hold off;

a           = kalman(output.ytilde, llmb);
fprintf(1, 'Prob. of Cambridge win in 2001: %.2f.\n', a(length(y)+1));
fprintf(1, '\n');
