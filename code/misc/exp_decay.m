% exponental decay function
% t: time
% gamma: decay rate
% lambda_0: initial value
% t_0: initial time offset (usually 0)
% lambda_ss: asymptotic value of the function
% Tim Mullen, SCCN/INC UCSD 2013
function lambda = exp_decay(t,gamma,lambda_0,t_0,lambda_ss)
    lambda = lambda_0 ./ ((t + t_0) .^ gamma) + lambda_ss;   % forgetting rate