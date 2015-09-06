internet    = load('data/internet.dat')';
y           = internet(2, 2:end);
fprintf(1, '\n');

P       = 5;
Q       = 5;

%% Model selection for complete series %%
fprintf(1, 'AIC for ARMA(p, q) models on complete data:\n');
fprintf(1, '    '); for q = 0:Q, fprintf(1, '%-12d', q); end; fprintf(1, '\n');
logL    = zeros(P+1, Q+1);
AIC     = zeros(P+1, Q+1);
arma    = cell(P+1, Q+1);
for p = 0 : P
    fprintf(1, '%-4d', p);
    for q = 0 : Q
        model                                   = ssm_arma(p, q);
        [arma{p+1, q+1} logL(p+1, q+1) output]  = estimate(y, model, [repmat(0.3, 1, model.w-1) 10], [], 'fmin', 'bfgs', 'disp', 'off');
        AIC(p+1, q+1)                           = output.AIC;
        fprintf(1, '%-12g', AIC(p+1, q+1));
    end
    fprintf(1, '\n');
end
[m i] = min(AIC(:)); temp = AIC(i); AIC(i) = Inf;
fprintf(1, 'ARMA(%d, %d) found to be the best model, ', mod(i-1, P+1), floor((i-1)/(P+1)));
[m j] = min(AIC(:)); AIC(i) = temp;
fprintf(1, 'ARMA(%d, %d) is the second best.\n\n', mod(j-1, P+1), floor((j-1)/(P+1)));

%% Model selection for data with missing values %%
ymis        = y; ymis([6 16 26 36 46 56 66 72 73 74 75 76 86 96]) = NaN;
fprintf(1, 'AIC for ARMA(p, q) models on data w/ missing observations:\n');
fprintf(1, '    '); for q = 0:Q, fprintf(1, '%-12d', q); end; fprintf(1, '\n');
logLmis     = zeros(P+1, Q+1);
AICmis      = zeros(P+1, Q+1);
armamis     = cell(P+1, Q+1);
for p = 0 : P
    fprintf(1, '%-4d', p);
    for q = 0 : Q
        model                                           = ssm_arma(p, q);
        [armamis{p+1, q+1} logLmis(p+1, q+1) output]    = estimate(ymis, model, [repmat(-0.1, 1, model.w-1) 1], [], 'fmin', 'bfgs', 'disp', 'off');
        AICmis(p+1, q+1)                                = output.AIC;
        fprintf(1, '%-12g', AICmis(p+1, q+1));
    end
    fprintf(1, '\n');
end
[m i] = min(AICmis(:)); temp = AICmis(i); AICmis(i) = Inf;
fprintf(1, 'ARMA(%d, %d) found to be the best model, ', mod(i-1, P+1), floor((i-1)/(P+1)));
[m j] = min(AICmis(:)); AICmis(i) = temp;
fprintf(1, 'ARMA(%d, %d) is the second best.\n\n', mod(j-1, P+1), floor((j-1)/(P+1)));

%% Forecast with ARMA(1, 1) on the complete data %%
[armafore logL] = estimate(y, ssm_arma(1, 1), 0.1);
yf              = signal(kalman([y repmat(NaN, 1, 20)], armafore), armafore);
figure('Name', 'Internet series forecast');
plot(yf, 'DisplayName', 'forecast'), hold all, scatter(1:length(y), y, 10, 'r', 's', 'filled', 'DisplayName', 'data'), hold off, ylim([-15 15]), legend('show');
if ispc, set(gcf, 'WindowStyle', 'docked'); end

%% Forecast with ARMA(1, 1) on the data w/ missing values %%
[armafore logLmis]  = estimate(ymis, ssm_arma(1, 1), 0.1);
ymisf               = signal(kalman([ymis repmat(NaN, 1, 20)], armafore), armafore);
figure('Name', 'Internet series in-sample one-step and out-of-sample forecasts');
plot(ymisf), hold all, scatter(1:length(ymis), ymis, 10, 'r', 's', 'filled'), hold off, ylim([-15 15]);
if ispc, set(gcf, 'WindowStyle', 'docked'); end
