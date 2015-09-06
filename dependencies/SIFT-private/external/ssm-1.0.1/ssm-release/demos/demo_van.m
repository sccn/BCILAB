van     = load('data/van.dat')';
time    = 1969:1/12:1984+11/12;
y       = van(1, :);
fprintf(1, '\n');

%% Structural time series model with poisson density %%
bstsmp                  = [ssm_poisson ssm_llm ssm_seasonal('dummy fixed', 12) ssm_reg(van(2, :))];
bstsmp.name             = 'Poisson STSM';
randn('state', [3765023265; 2369472656]);
[bstsmp logL output]    = estimate(y, bstsmp, 0.0006, [], 'fmin', 'bfgs', 'disp', 'off');
[alphahat V]            = statesmo(output.ytilde, bstsmp);
ycom                    = signal(alphahat, bstsmp);
lvl                     = ycom(1, :) + ycom(3, :);
seas                    = ycom(2, :);

fprintf(1, '[Poisson STSM]\n');
fprintf(1, 'loglikelihood:     %g\n', logL);
fprintf(1, 'Level variance:    %g\n', bstsmp.param);
fprintf(1, 'Intervention coeff: %g\n', alphahat(end, end));
fprintf(1, '             st.d.: %g\n', realsqrt(V(end, end, end)));

figure('Name', 'Number of van drivers killed and estimated level');
subplot(2, 1, 1), plot(time, exp(lvl), 'DisplayName', 'est. level'), hold all, plot(time, y, 'r:', 'DisplayName', 'data'), hold off, ylim([1 18]), legend('show');
subplot(2, 1, 2), plot(time, exp(lvl), 'DisplayName', 'est. level'), hold all, plot(time, exp(reallog(y)-seas), 'r:', 'DisplayName', 'seas. adj. data'), hold off, ylim([0 19]), legend('show');

fprintf(1, '\n');
