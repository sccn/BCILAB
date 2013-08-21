function x = check_1d_distributed(x)
standard_distributions = {'bino','chi2','exp','ev','f','gam','gev','gp','geo','hyge','logn','nbin','ncf','nct','ncx2','norm', 'poiss','rayl','t','unif','unid','wbl'};
x = is_1d_distributed(x) && any(strcmp(x{1},standard_distributions)) && isnumeric(x{2});
