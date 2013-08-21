function [results] = vdpgm(given_data, opts)
%
% function [results] = vdpgm(given_data, opts)
%

start_time = clock;
if nargin == 1
  opts = struct();
end
if issparse(given_data)
  given_data = full(given_data);
end
if ~ isfield(opts, 'algorithm')
  % algorithm can be one of 'vdp', 'bj', 'cdp', 'csb' and 'non_dp'
  % vdp : variational DP
  % bj : Blei and Jordan
  % cdp : collapsed Dirichlet prior
  % csb : collapsed stick-breaking
  % non_dp : variational Bayes for Gaussian mixtures
  opts.algorithm = 'vdp';
end
if ~ isfield(opts, 'collapsed_means')
  opts.collapsed_means = 0;
end
if ~ isfield(opts, 'do_sort')
  opts.do_sort = '0';
end
if ~ isfield(opts, 'get_q_of_z')
  opts.get_q_of_z = 0;
end
if ~ isfield(opts, 'weight')
  opts.weight = 1;
end
if ~ isfield(opts, 'get_log_likelihood')
  opts.get_log_likelihood = 0;
end
if ~ isfield(opts, 'use_kd_tree')
  opts.use_kd_tree = 1;
end
if ~ isfield(opts, 'threshold')
  opts.threshold = 1.0e-5;
end
if ~ isfield(opts, 'sis')
  opts.sis = 0;
end
if ~ isfield(opts, 'initial_depth')
  opts.initial_depth = 3;
end
if ~ isfield(opts, 'initial_K')
  opts.initial_K = 1;
end
if ~ isfield(opts, 'ite')
  opts.ite = inf;
end
if ~ isfield(opts, 'do_split')
  opts.do_split = 0;
end
if ~ isfield(opts, 'do_merge')
  opts.do_merge = 0;
end
if ~ isfield(opts, 'do_greedy')
  opts.do_greedy = 1;
end
if ~ isfield(opts, 'max_target_ratio')
  opts.max_target_ratio = 0.5;
end
if ~ isfield(opts, 'init_of_split')
  % 'pc', 'rnd', 'rnd_close' or 'close_f'
  opts.init_of_split = 'pc';
end
if ~ isfield(opts, 'recursive_expanding_depth')
  opts.recursive_expanding_depth = 2;
end
if ~ isfield(opts, 'recursive_expanding_threshold')
  opts.recursive_expanding_threshold = 1.0e-1;
end
if ~ isfield(opts, 'recursive_expanding_frequency')
  opts.recursive_expanding_frequency = 3;
end
if isfield(opts, 'seed')
  rand('state', opts.seed);
else
  seed = rand('state');
  results.seed = seed;
end

data.given_data = given_data;
if opts.use_kd_tree
  partitions = init_kdtree_partitions(given_data, opts.initial_depth);
  data.kdtree = partitions;
end

% the hyperparameters of priors
hp_prior = mk_hp_prior(data, opts);

if isfield(opts, 'hp_posterior')
  opts.use_kd_tree = 0;
  if opts.get_q_of_z
    results.q_of_z = mk_q_of_z(data, opts.hp_posterior, opts.hp_prior, opts);
  end
  if opts.get_log_likelihood
    results.log_likelihood = mk_log_likelihood(data, opts.hp_posterior, opts.hp_prior, opts);
  end
  if isfield(opts, 'test_data')
    results.predictive_posterior = log_predictive_dist(data, opts.q_of_z, opts.hp_posterior, ...
                                                      opts.hp_prior, opts);
  end
  return
end

if opts.sis > 0
  q_of_z = sequential_importance_sampling(data, hp_prior, opts);
elseif isfield(opts, 'q_of_z')
  q_of_z = opts.q_of_z;
else
  q_of_z = rand_q_of_z(data, opts.initial_K, opts);
end

hp_posterior = mk_hp_posterior(data, q_of_z, hp_prior, opts);
if opts.do_greedy
  [free_energy, hp_posterior, data] = greedy(data, hp_posterior, hp_prior, opts);
else
  [free_energy, hp_posterior, data] = split_merge(data, hp_posterior, hp_prior, opts);
end

results.algorithm = opts.algorithm;
results.elapsed_time = etime(clock, start_time);
results.free_energy = free_energy;
results.hp_prior = hp_prior;
results.hp_posterior = hp_posterior;
results.K = length(hp_posterior.eta);
results.opts = opts;
if opts.get_q_of_z
  results.q_of_z = mk_q_of_z(data, hp_posterior, hp_prior, opts);
  if opts.use_kd_tree
    results.q_of_z = mk_non_kdtree_q_of_z(data, results.q_of_z);
  end
end
if opts.get_log_likelihood
  results.log_likelihood = mk_log_likelihood(data, hp_posterior, hp_prior, opts);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function log_likelihood = mk_log_likelihood(data, hp_posterior, hp_prior, opts);
[D,N] = size(data.given_data);
K = size(hp_posterior.m, 2);
log_likelihood = zeros(K,N);
E_pi = mk_E_pi(hp_posterior, hp_prior, opts);
for c=1:K
  mu = hp_posterior.m(:,c);
  f = hp_posterior.eta(c) + 1 - D;
  Sigma = (hp_posterior.xi(c)+1) / hp_posterior.xi(c) / f * hp_posterior.B{c};
  log_likelihood(c,:) = log_no_w(E_pi(c)) + logmvtpdf(data.given_data, mu, f, Sigma);
end
log_likelihood = log_sum_exp(log_likelihood, 1); % 1 by N
log_likelihood = sum(log_likelihood, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function log_prob = log_predictive_dist(data, q_of_z, hp_posterior, hp_prior, opts);
if opts.use_kd_tree
  N = length(data.kdtree);
  D = size(data.kdtree(1).sum_x, 1);
  Na = [data.kdtree(:).N];
  if isequal(opts.algorithm, 'vdp')
    true_Nc = Na*q_of_z; % 1*K
    q_of_z(:,end) = 0;
  end
  Nc = Na*q_of_z; % 1*K
  sum_x = [data.kdtree(:).sum_x] * q_of_z;
else
  [D,N] = size(data.given_data);
  if isequal(opts.algorithm, 'vdp')
    true_Nc = sum(q_of_z, 1); % 1*K
    q_of_z(:,end) = 0;
  end
  Nc = sum(q_of_z, 1); % 1*K
  sum_x = data.given_data * q_of_z;
end
E_pi = mk_E_pi(hp_posterior, hp_prior, opts);
m = (sum_x + repmat(hp_prior.xi0*hp_prior.m0, 1, size(sum_x,2))) ./ repmat(Nc, D, 1);
E_pi(find(Nc==0)) = 0;
n = size(opts.test_data, 2);
log_prob = zeros(size(m,2), n);
for c=1:size(m,2)
  f = hp_posterior.eta(c) + 1 - D;
  sigma = hp_posterior.B{c} * (hp_posterior.xi(c)+1) / hp_posterior.xi(c) / f;
  if E_pi(c) > 0
    log_prob(c,:) = log(E_pi(c)) + log_T_dist(opts.test_data, m(:,c), ...
                                              sigma, f);
  else
    log_prob(c,:) = -inf;
  end
end
log_prob = log_sum_exp(log_prob,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function log_prob = log_T_dist(X, m, sigma, f);
% X is a D by n matrix
[D,n] = size(X);
diff = X-repmat(m,1,n);
log_prob = gammaln((f+D)/2) ...
    - D/2*log(f*pi) ...
    - gammaln(f/2) ...
    - 0.5*detln(sigma) ...
    - (f+D)/2*log(1 + sum( diff .* ((f*sigma) \ diff), 1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quad_term = mk_quad_term(data, q_of_z, hp_prior, opts)
% calculate the quad_term in marginalizing out means
% quad_term : D by D by K
%
[N, K] = size(q_of_z);
D = size(data.given_data, 1);
Nc = sum(q_of_z, 1); % 1 by K
xi = Nc + hp_prior.xi0;
f0 = 1 ./ (xi.^2);
f1 = - 2 ./ (xi.^3);
f2 = 3 ./ (xi.^4);
v = q_of_z.*(1 - q_of_z); % N by K
c = v.*((1-q_of_z).^2 + q_of_z.^2);
q = v.*((1-q_of_z).^3 + q_of_z.^3);
p_x = data.given_data * q_of_z; % D by K
v_x = data.given_data * v;      % D by K
c_x = data.given_data * c;      % D by K
v_p_x = data.given_data * (q_of_z + v); % D by K
c_p_x = data.given_data * (q_of_z + c); % D by K
for t=1:K
  first_term = f0(t) ...
      *((repmat(v(:,t)', D, 1).*data.given_data)*data.given_data' ...
        + p_x(:,t)*p_x(:,t)');
  second_term = f1(t) ...
      *((repmat(c(:,t)', D, 1).*data.given_data)*data.given_data' ...
        + v_p_x(:,t)*v_p_x(:,t)' - v_x(:,t)*v_x(:,t)' - p_x(:,t)*p_x(:,t)');
  tmp = q(:,t) - 3*v(:,t).^2 + v(:,t)*sum(v(:,t));
  third_term = f2(t) ...
      *((repmat(tmp', D, 1).*data.given_data)*data.given_data' ...
        + c_p_x(:,t)*c_p_x(:,t)' - c_x(:,t)*c_x(:,t)' - p_x(:,t)*p_x(:,t)' ...
        + 2*v_x(:,t)*v_x(:,t)' ...
        + sum(v(:,t))*p_x(:,t)*p_x(:,t)');
  quad_term(:,:,t) = first_term + second_term + third_term;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_log_p_of_z_given_other_z = mk_E_log_p_of_z_given_other_z(hp_posterior, hp_prior, opts);
% returns E[log p(z_n|Z^-n)]_q(Z^-n)
% E_log_p_of_z_given_other_z : N by K
q_of_z = hp_posterior.q_of_z;
[N,K] = size(q_of_z);
[E_Nc_minus_n, V_Nc_minus_n] = mk_E_Nc_minus_n(q_of_z); % N by K
use_variance = 1;
if isequal(opts.algorithm, 'cdp')
  if use_variance
    E_log_p_of_z_given_other_z = ...
        log(E_Nc_minus_n + hp_prior.alpha/K) ...
        - 0.5*V_Nc_minus_n./((E_Nc_minus_n+hp_prior.alpha/K).^2) ...
        - log(N - 1 + hp_prior.alpha);
  else
    E_log_p_of_z_given_other_z = ...
        log(E_Nc_minus_n + hp_prior.alpha/K) ...
        - log(N - 1 + hp_prior.alpha);
  end
elseif isequal(opts.algorithm, 'csb')
  E_Nc_minus_n_cumsum_geq = fliplr(cumsum(fliplr(E_Nc_minus_n), 2));
  q_of_z_cumsum_geq = fliplr(cumsum(fliplr(q_of_z), 2));
  dummy = q_of_z_cumsum_geq.*(1-q_of_z_cumsum_geq);
  v_Nc_cumsum_geq = repmat(sum(dummy, 1), N, 1) - dummy;

  E_Nc_minus_n_cumsum = E_Nc_minus_n_cumsum_geq - E_Nc_minus_n;
  q_of_z_cumsum = q_of_z_cumsum_geq - q_of_z;
  dummy = q_of_z_cumsum.*(1-q_of_z_cumsum);
  v_Nc_cumsum = repmat(sum(dummy, 1), N, 1) - dummy;
  
  if use_variance
    first_term = ...
        (log(1+E_Nc_minus_n) ...
         - 0.5*V_Nc_minus_n./((1+E_Nc_minus_n).^2) ...
         - log(1+hp_prior.alpha+E_Nc_minus_n_cumsum_geq) ...
         + 0.5*v_Nc_cumsum_geq./((1+hp_prior.alpha+E_Nc_minus_n_cumsum_geq).^2));
    first_term(:,end) = 0;
    dummy = log(hp_prior.alpha+E_Nc_minus_n_cumsum) ...
            - 0.5*v_Nc_cumsum./((hp_prior.alpha+E_Nc_minus_n_cumsum).^2) ...
            - log(1+hp_prior.alpha+E_Nc_minus_n_cumsum_geq) ...
            + 0.5*v_Nc_cumsum_geq./((1+hp_prior.alpha+E_Nc_minus_n_cumsum_geq).^2);
    second_term = cumsum(dummy, 2) - dummy;
  else
    first_term = log(1+E_Nc_minus_n) - log(1+hp_prior.alpha+E_Nc_minus_n_cumsum_geq);
    first_term(:,end) = 0;
    dummy = log(hp_prior.alpha+E_Nc_minus_n_cumsum) ...
            - log(1+hp_prior.alpha+E_Nc_minus_n_cumsum_geq);
    second_term = cumsum(dummy, 2) - dummy;
  end
  E_log_p_of_z_given_other_z = first_term + second_term;
else
  error('unsupported algorithm');
end
% Gaussian approximation may not give a proper distribution.
% note. E[log p] is not need to be proper


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_Nc_minus_n, V_Nc_minus_n] = mk_E_Nc_minus_n(q_of_z)
% returns E[Nc^-n]; the expected value of Nc-n; N by K
%         V[Nc^-n]; the variance of Nc-n;       N by K
% q_of_z : N by K
N = size(q_of_z, 1);
E_Nc_minus_n = repmat(sum(q_of_z, 1), N, 1) - q_of_z;
tmp = q_of_z.*(1-q_of_z);
V_Nc = sum(tmp, 1);
V_Nc_minus_n = repmat(V_Nc, N, 1) - tmp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_pi = mk_E_pi(hp_posterior, hp_prior, opts);
if isequal(opts.algorithm, 'cdp') | isequal(opts.algorithm, 'non_dp')
  E_pi = hp_posterior.tilde_alpha / sum(hp_posterior.tilde_alpha);
elseif isequal(opts.algorithm, 'bj')
  second_term = [0 cumsum(log(hp_posterior.gamma(2,:)) - log(sum(hp_posterior.gamma)),2)]; % 1 by K
  E_pi = exp([(log(hp_posterior.gamma(1,:)) ...
               - log(sum(hp_posterior.gamma,1)) ...
               + second_term(1:end-1)) ...
              , second_term(end)]); % 1 by K
elseif isequal(opts.algorithm, 'vdp')
  opts_tmp = opts;
  opts_tmp.algorithm = 'bj';
  E_pi = mk_E_pi(hp_posterior, hp_prior, opts_tmp);
  E_pi(end) = max(1 - sum(E_pi(1:end-1)), 0);
elseif isequal(opts.algorithm, 'csb')
  [N,K] = size(hp_posterior.q_of_z);
  a = 1 + hp_posterior.Nc; % 1 by K
  b = hp_prior.alpha + N - cumsum(hp_posterior.Nc, 2); % 1 by K
  first_term = a ./ (a+b);
  first_term(end) = 1;
  second_term = cumprod(b./(a+b)) ./ (b./(a+b));
  E_pi = first_term .* second_term;
else
  error('not supported algorithm')
end

%
% p_TSB(z_test|Z)
%
% \begin{align*}
% p_{\text{TSB}}(Z) =& 
% \alpha^{T-1}
% \prod_{i<T}
% \frac
% {\Gamma(1+N_i) \Gamma(\alpha + N_{>i})}
% {\Gamma(1+\alpha+N_{\geq i})}
% \\
% p_{\text{TSB}}(z_{test}=t,Z) = &
% \left(
% \frac{1+N_t}{1+\alpha+N_{\geq t}}
% \right)^{\mathbb{I}(t<T)}
% \left\{
% \prod_{j<t}\frac
% {\alpha + N_{>j}}
% {1+\alpha+N_{\geq j}}
% \right\}
% \alpha^{T-1}
% \prod_{i<T}
% \frac
% {\Gamma(1+N_i) \Gamma(\alpha + N_{>i})}
% {\Gamma(1+\alpha+N_{\geq i})}
% \\
% p_{\text{TSB}}(z_{test}=t|Z) = &
% \left(
% \frac{1+N_t}{1+\alpha+N_{\geq t}}
% \right)^{\mathbb{I}(t<T)}
% \prod_{j<t}\frac
% {\alpha + N_{>j}}
% {1+\alpha+N_{\geq j}}
% \end{align*}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q_of_z = sequential_importance_sampling(data, hp_prior, opts);
disp('start SIS.')
N = size(data.given_data, 2);
for r = 1:opts.sis
  I = randperm(N);
%   I = I(1:ceil(N*0.1));
  data_r.given_data = data.given_data(:,I);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   data_r_sub.given_data = data_r.given_data(:,1);
%   q_of_z = ones(1,opts.initial_K) / opts.initial_K;
%   hp_posterior = mk_hp_posterior(data_r_sub, q_of_z, hp_prior, opts);
%   for n = 1:N
%     data_r_sub.given_data = data_r.given_data(:,1:n);
%     q_of_z = mk_q_of_z(data_r_sub, hp_posterior, hp_prior, opts);
%     q_of_z = sort_q_of_z(data, q_of_z, opts);  
%     hp_posterior = mk_hp_posterior(data_r_sub, q_of_z, hp_prior, opts);
%   end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  q_of_z = [];
  for n = 1:size(data_r.given_data,2)
    data_r_sub.given_data = data_r.given_data(:,1:n);
    q_of_z(n,:) = ones(1, opts.initial_K) / opts.initial_K;
    hp_posterior = mk_hp_posterior(data_r_sub, q_of_z, hp_prior, opts);
    q_of_z = mk_q_of_z(data_r_sub, hp_posterior, hp_prior, opts);
%     q_of_z = sort_q_of_z(data, q_of_z, opts);  
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  hist(r).hp_posterior = hp_posterior;
  hist(r).free_energy = mk_free_energy(data, hp_posterior, hp_prior, opts);
  [min_f, best_r] = min([hist(:).free_energy]);
  hp_posterior = hist(best_r).hp_posterior;
  disp(['SIS: ' num2str(r) ';  best f = ' num2str(min_f) ';  best Nc = ' num2str(hp_posterior.Nc)])
end
disp('SIS done.')
q_of_z = mk_q_of_z(data, hp_posterior, hp_prior, opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nonkdtree_q_of_z = mk_non_kdtree_q_of_z(data, q_of_z);
nonkdtree_q_of_z = zeros(size(data.given_data,2), size(q_of_z, 2)); % N*K
for a=1:length(data.kdtree);
  nonkdtree_q_of_z(data.kdtree(a).indices,:) = repmat(q_of_z(a,:), length(data.kdtree(a).indices), 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [free_energy, hp_posterior, data] = greedy(data, hp_posterior, hp_prior, opts);
free_energy = mk_free_energy(data, hp_posterior, hp_prior, opts);
disp_status(free_energy, hp_posterior, opts);
while 1
  disp('finding the best one....')
  [new_free_energy, new_hp_posterior, new_data, c] = find_best_splitting(data, ...
                                                    hp_posterior, ...
                                                    hp_prior, opts);
  if c == -1
    break
  end
  disp(['finding the best one.... done.  component ' num2str(c) ' was split.'])
  disp_status(new_free_energy, new_hp_posterior, opts);
  [new_free_energy, new_hp_posterior, new_data] = update_posterior2(new_data, ...
                                                    new_hp_posterior, ...
                                                    hp_prior, opts, 1, opts.ite);
  if free_energy_improved(free_energy, new_free_energy, 0, opts) == 0
    break
  end
  free_energy = new_free_energy;
  hp_posterior = new_hp_posterior;
  data = new_data;
end
disp_status(free_energy, hp_posterior, opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [free_energy, hp_posterior, data, c] = find_best_splitting(data, ...
                                                  hp_posterior, ...
                                                  hp_prior, opts);
c_max = 10;
K = size(hp_posterior.m, 2);
candidates = find(hp_posterior.Nc>2);
if isempty(candidates)
  free_energy = 0;
  c = -1;
  return
end
q_of_z = mk_q_of_z(data, hp_posterior, hp_prior, opts);
new_free_energy = ones(1,max(candidates))*inf;
%%%%%%%%%%%%%%%%%%%%%
fc = mk_E_log_q_p_eta(data, hp_posterior, hp_prior, opts);
log_lambda = mk_log_lambda(data, hp_posterior, hp_prior, opts);
%%%%%%%%%%%%%%%%%%%%%
for c = candidates(1:min(c_max, length(candidates)))
  disp(['splitting ' num2str(c) '...'])
  [new_data(c), new_q_of_z, info] = split(c, data, q_of_z, hp_posterior, hp_prior, opts, 1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  new_c = info.new_c;
  relating_n = find(sum(new_q_of_z(:,[c new_c]),2) > 0.5);
  if isempty(relating_n)
    continue
  end
  new_K = size(new_q_of_z, 2);
  sub_q_of_z = new_q_of_z(relating_n, [c new_c new_K]);
  if opts.use_kd_tree
    sub_data.kdtree = new_data(c).kdtree(relating_n);
  else
    sub_data.given_data = new_data(c).given_data(:,relating_n);
  end
  sub_hp_posteior = mk_hp_posterior(sub_data, sub_q_of_z, hp_prior, opts);
  [sub_f, sub_hp_posteior, dummy, sub_q_of_z] = update_posterior2(sub_data, ...
                                                    sub_hp_posteior, ...
                                                    hp_prior, opts, 0, 10, 0);
  if size(sub_q_of_z,2) < 3
    continue
  end
  if length(find(sum(sub_q_of_z,1)<1.0e-10)) > 1
    continue
  end
  new_log_lambda = log_lambda;
  if opts.use_kd_tree
    updated_data.kdtree = new_data(c).kdtree(info.updated_I);
    new_log_lambda(info.updated_I,:) = mk_log_lambda(updated_data, hp_posterior, hp_prior, opts);
  end
  sub_log_lambda = mk_log_lambda(new_data(c), sub_hp_posteior, hp_prior, opts);
  insert_indices = [c new_c new_K:(new_K+size(sub_q_of_z,2)-3)];
  new_log_lambda(:,insert_indices) = sub_log_lambda;
  new_fc = fc;
  new_fc(insert_indices) = mk_E_log_q_p_eta(sub_data, sub_hp_posteior, hp_prior, opts);
  new_free_energy(c) = mk_free_energy(new_data(c), sub_hp_posteior, hp_prior, opts, new_fc, new_log_lambda);
  new_q_of_z(relating_n,:) = 0;
  new_q_of_z(relating_n,insert_indices) = sub_q_of_z;
  new_q_of_z_cell{c} = new_q_of_z;
end
[free_energy, c] = min(new_free_energy);
if isinf(free_energy)
  c = -1;
  return
end
data = new_data(c);
hp_posterior = mk_hp_posterior(data, new_q_of_z_cell{c}, hp_prior, opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [free_energy, hp_posterior, data] = split_merge(data, hp_posterior, hp_prior, opts);
%
% split: 
% tdp: if an empty cluster exists, split is available.
%      otherwise, make a new empty cluster. =>  K<-K+1.
%
% non-dp: if 
%
c_max = 3;
[free_energy, hp_posterior] = update_posterior2(data, hp_posterior, hp_prior, opts, 0, opts.ite);
while 1
  K = size(hp_posterior.m, 2);
  new_free_energy = inf;
  %%% split
  if ~ opts.do_split
    break
  end
  disp(['### start splitting'])
  [dummy, candidates_for_splliting] = sort(hp_posterior.Nc, 2, 'descend');
  candidates_for_splliting(find(hp_posterior.Nc(candidates_for_splliting)<1)) = [];
  q_of_z = mk_q_of_z(data, hp_posterior, hp_prior, opts); % N*K
  for c=candidates_for_splliting(1:min(c_max, length(candidates_for_splliting)))
    disp(['### start splitting k=' num2str(c) '...'])
    [new_data, new_q_of_z] = split(c, data, q_of_z, hp_posterior, hp_prior, opts);
    new_hp_posterior = mk_hp_posterior(new_data, new_q_of_z, hp_prior, opts);
    [new_free_energy, new_hp_posterior, new_data] = update_posterior2(new_data, ...
                                                      new_hp_posterior, ...
                                                      hp_prior, ...
                                                      opts, 1);
    if free_energy_improved(free_energy, new_free_energy, 0, opts)
      disp(['### splitting k=' num2str(c) ' improved the free energy.  ' ...
           num2str(free_energy) ' -> ' num2str(new_free_energy)])
      break
    else
      disp(['### splitting k=' num2str(c) ' did not improve the free energy.'])
    end
  end
  disp('### end splitting')
  if free_energy_improved(free_energy, new_free_energy, 0, opts)
    free_energy = new_free_energy;
    hp_posterior = new_hp_posterior;
    data = new_data;
    continue
  end
  if ~ opts.do_merge
    break
  end
  %%% merge
  disp(['### start merging'])
  q_of_z = mk_q_of_z(data, hp_posterior, hp_prior, opts); % N*K
  candidates_for_merging = mkcandidates_to_merge(q_of_z'); % 2*#pairs
  for i=1:min(c_max, size(candidates_for_merging,2))
    pair = candidates_for_merging(:,i);
    c1 = pair(1); c2 = pair(2);
    disp(sprintf('### start merging c1=%d, c2=%d ...', c1, c2))
    new_q_of_z = q_of_z;
    new_q_of_z(:,c1) = new_q_of_z(:,c1) + new_q_of_z(:,c2);
    new_q_of_z(:,c2) = zeros(size(new_q_of_z,1),1);
    new_hp_posterior = mk_hp_posterior(data, new_q_of_z, hp_prior, opts);
    [new_free_energy, new_hp_posterior] = update_posterior2(data, ...
                                                      new_hp_posterior, ...
                                                      hp_prior, opts);
    if free_energy_improved(free_energy, new_free_energy, 0, opts)
      disp(sprintf('### merging c1=%d, c2=%d improved the free energy.  %0.5g -> %0.5g', ...
                   c1, c2, free_energy, new_free_energy))
      break
    else
      disp(sprintf('### merging c1=%d, c2=%d did not improve the free energy.', c1, c2))
    end
  end
  disp(['### end merging'])
  if free_energy_improved(free_energy, new_free_energy, 0, opts)
    free_energy = new_free_energy;
    hp_posterior = new_hp_posterior;
    continue
  end
  break
end % while 1
disp_status(free_energy, hp_posterior, opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function log_p_of_x_given_c = mk_map_log_p_of_x_given_c(data, clusters, hp_posterior, opts);
% clusters: e.g [1:K]
[D,N] = size(data);
K = length(clusters);
log_p_of_x_given_c = zeros(K,N);
for i = 1:K
  c = clusters(i);
  m = hp_posterior.m(:,c);
  precision = hp_posterior.inv_B(:,:,c)*hp_posterior.eta(c);
  d = data - repmat(m, 1, N);
  log_p_of_x_given_c(i,:) = (-D*0.5)*log(2*pi)+0.5*detln(precision)-0.5*sum(d.*(precision*d),1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_data, new_q_of_z, info] = split(c, data, q_of_z, ...
                                              hp_posterior, hp_prior, opts, update_kdtree)
% q_of_z: N*K
if nargin < 7
  update_kdtree = 1;
end
new_data = data;
if opts.use_kd_tree & update_kdtree
  [dummy, indices] = max(q_of_z,[],2);
  target_partitions = find(indices==c)';
  if ~ isempty(target_partitions)
    m = [data.kdtree(target_partitions).mean];
    log_p_of_m_given_c = mk_map_log_p_of_x_given_c(m, c, hp_posterior, opts); % 1*|target_partitions|
    [dummy, I] = sort(log_p_of_m_given_c, 2, 'descend');
    target_partitions = target_partitions(I(1:ceil(length(I)*opts.max_target_ratio)));
    %   target_partitions = find(q_of_z(:,c)>0.01)';
    [new_data, q_of_z, updated_I] = expand_all_nodes(new_data, q_of_z, hp_posterior, ...
                                                     hp_prior, target_partitions, opts);
    info.updated_I = updated_I;
  else
    info.updated_I = [];
  end
end
if isequal(opts.init_of_split, 'pc')  % principal eigenvector
  if opts.use_kd_tree
    arg1_data = [new_data.kdtree(:).mean];
  else
    arg1_data = new_data.given_data;
  end
  dir = divide_by_principal_component(arg1_data, ...
                                      hp_posterior.B{c}/hp_posterior.eta(c), ...
                                      hp_posterior.m(:,c));
  q_of_z_c1 = zeros(size(q_of_z,1),1);
  q_of_z_c2 = q_of_z(:,c);
  I = find(dir>=0);
  q_of_z_c1(I) = q_of_z(I,c);
  q_of_z_c2(I) = 0;
else
  q_of_z_c = q_of_z(:,c);
  if isequal(opts.init_of_split, 'rnd')  % random
    r = rand(size(q_of_z,1),1);
  elseif isequal(opts.init_of_split, 'rnd_close')  % make close clusters
    r = 0.5 + (rand(size(q_of_z,1),1)-0.5)*0.01;
  elseif isequal(opts.init_of_split, 'close_f')  % one is almost zero.
    r = 0.98 + rand(size(q_of_z,1),1)*0.01;
  else
    init_of_split = opts.init_of_split
    error('unknown option')
  end
  q_of_z_c1 = q_of_z_c.*r;
  q_of_z_c2 = q_of_z_c.*(1-r);
end
new_q_of_z = zeros(size(q_of_z,1), size(q_of_z,2)+1);
new_q_of_z(:,[1:end-2 end]) = q_of_z;
new_q_of_z(:,c) = q_of_z_c1;
new_c = size(new_q_of_z, 2) - 1;
new_q_of_z(:,new_c) = q_of_z_c2;
info.new_c = new_c;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_data, new_q_of_z, updated_I] = expand_all_nodes(data, ...
                                                  q_of_z, hp_posterior, ...
                                                  hp_prior, ...
                                                  target_partitions, opts);
if size(target_partitions, 1) ~= 1
  error('target_partitions must be a row vector.')
end
new_data = data;
for a=target_partitions
  if isempty(data.kdtree(a).children)
    children = mk_children_of_partition(data.kdtree(a));
    data.kdtree(a).children = children;
  else
    children = data.kdtree(a).children;
  end
  if length(children) == 1
    continue
  end
  new_data.kdtree(a) = children(1);
  new_data.kdtree(end+1) = children(2);
end % while a <= length(data.kdtree)
updated_I = [target_partitions length(data.kdtree)+1:length(new_data.kdtree)];
sub_data.given_data = data.given_data;
sub_data.kdtree = new_data.kdtree(updated_I);
sub_q_of_z = mk_q_of_z(sub_data, hp_posterior, hp_prior, opts);
new_q_of_z = zeros(length(new_data.kdtree),size(q_of_z,2));
new_q_of_z(1:size(q_of_z,1),:) = q_of_z;
new_q_of_z(updated_I,:) = sub_q_of_z;
disp(['### building kd-tree done; #partition ' num2str(length(data.kdtree)) ...
      ' -> ' num2str(length(new_data.kdtree))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, q_of_z] = expand_recursively_until_convergence(data, ...
                                                  q_of_z, hp_posterior, ...
                                                  hp_prior, opts);
kdtree_size = length(data.kdtree);
for a=1:length(data.kdtree)
  [data, q_of_z] = expand_recursively_until_convergence2(data, a, ...
                                                    q_of_z, hp_posterior, ...
                                                    hp_prior, opts, ...
                                                    opts.recursive_expanding_depth);
end
[prob, best_c] = max(q_of_z, [], 2);
% for c=1:size(q_of_z,2);
%   I_n = 
% end
disp(['### building kd-tree done; #partition ' num2str(kdtree_size) ...
      ' -> ' num2str(length(data.kdtree))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, q_of_z] = expand_recursively_until_convergence2(data, a, ...
                                                  q_of_z, hp_posterior, ...
                                                  hp_prior, opts, depth);
% q_of_z : N*K
if depth == 0
  return
end
if isempty(data.kdtree(a).children)
  children = mk_children_of_partition(data.kdtree(a));
  data.kdtree(a).children = children;
else
  children = data.kdtree(a).children;
end
if length(children) == 1
  return
end
sub_data.given_data = data.given_data;
sub_data.kdtree = children;
sub_q_of_z = mk_q_of_z(sub_data, hp_posterior, hp_prior, opts); % 2*K
diff = sub_q_of_z - repmat(q_of_z(a,:), 2, 1);
if sum(sum(diff.*diff, 1), 2)/prod(size(sub_q_of_z)) < opts.recursive_expanding_threshold
  return
end
b = length(data.kdtree)+1;
data.kdtree(a) = children(1);
data.kdtree(b) = children(2);
q_of_z([a b],:) = sub_q_of_z;
[data, q_of_z] =  expand_recursively_until_convergence2(data, a, ...
                                                  q_of_z, hp_posterior, ...
                                                  hp_prior, opts, depth-1);
[data, q_of_z] =  expand_recursively_until_convergence2(data, b, ...
                                                  q_of_z, hp_posterior, ...
                                                  hp_prior, opts, depth-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [free_energy, hp_posterior, data, q_of_z] = update_posterior2(data, hp_posterior, hp_prior, opts, upkdtree, ite, do_sort);
% update q_of_z: N*K
disp(['### updating posterior ...'])
free_energy = inf;
if nargin == 4
  upkdtree = 0;
end
if nargin < 6
  ite = inf;
end
if nargin < 7
  do_sort = 1;
end
i = 0;
last_Nc = 0;
start_sort = 0;
while 1
  i = i+1;
  [new_free_energy, log_lambda] = mk_free_energy(data, hp_posterior, hp_prior, opts);
  disp_status(new_free_energy, hp_posterior, opts);
  if (~isinf(ite) && i>=ite) || ...
        (isinf(ite) && free_energy_improved(free_energy, new_free_energy, 0, opts) == 0)
    free_energy = new_free_energy;
    if do_sort && opts.do_sort && ~ start_sort
      start_sort = 1;
    else
      break
    end
  end
  last_Nc = hp_posterior.Nc;
  free_energy = new_free_energy;
  [q_of_z, data] = mk_q_of_z(data, hp_posterior, hp_prior, opts, log_lambda);
  freq = opts.recursive_expanding_frequency;
  if opts.use_kd_tree & upkdtree & (freq==1 | mod(i,freq)==1)
    [data, q_of_z] = expand_recursively_until_convergence(data, ...
                                                      q_of_z, ...
                                                      hp_posterior, ...
                                                      hp_prior, opts);
  end
  % check if the last component is small enough
  if isequal(opts.algorithm, 'vdp') & sum(q_of_z(:,end)) > 1.0e-20
    q_of_z(:,end+1) = 0;
  end
  if start_sort
    q_of_z = sort_q_of_z(data, q_of_z, opts);
  end
  if isequal(opts.algorithm, 'vdp') & sum(q_of_z(:,end-1)) < 1.0e-10
    q_of_z(:,end-1) = [];
  end
  hp_posterior = mk_hp_posterior(data, q_of_z, hp_prior, opts);
end
% disp_status(free_energy, hp_posterior, opts);
disp(['### updating posterior ... done.'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q_of_z, I] = sort_q_of_z(data, q_of_z, opts);
disp('sorting...')
if opts.use_kd_tree
  Nc = [data.kdtree(:).N]*q_of_z; % 1*K
else
  Nc = sum(q_of_z, 1); % 1*K
end
if isequal(opts.algorithm, 'vdp')
  [dummy,I] = sort(Nc(1:end-1), 2, 'descend');
  I(end+1) = length(Nc);
else
  [dummy,I] = sort(Nc, 2, 'descend');
end
q_of_z = q_of_z(:,I);
disp('sorting... done.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disp_status(free_energy, hp_posterior, opts);
if isequal(opts.algorithm, 'vdp')
  Nc = hp_posterior.true_Nc;
else
  Nc = hp_posterior.Nc;
end
disp(['F=' num2str(free_energy) ...
      ';    Nc=[' num2str(Nc, ' %0.5g ') ...
      '];'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = free_energy_improved(free_energy, new_free_energy, warn_when_increasing, opts);
diff = new_free_energy - free_energy;
if abs(diff/free_energy) < opts.threshold
  bool = 0;
elseif diff > 0
  if warn_when_increasing
    if abs(diff/free_energy) > 1.0e-3
      error(['the free energy increased.  the diff is ' num2str(new_free_energy-free_energy)])
    else
      warning(['the free energy increased.  the diff is ' num2str(new_free_energy-free_energy)])
    end
  end
  bool = 0;
elseif diff == 0
  bool = 0
else
  bool = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partitions = init_kdtree_partitions(given_data, depth);
[D,N] = size(given_data);
root = mk_kdtree_partition(given_data, [1:N]);
partitions = expand_recursively(root, depth);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partitions = expand_recursively(partition, depth);
children = mk_children_of_partition(partition);
if depth == 1 | length(children) == 1
  partitions = children;
else
  partitions1 = expand_recursively(children(1), depth-1);
  partitions2 = expand_recursively(children(2), depth-1);
  partitions = [partitions1 partitions2];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function children = mk_children_of_partition(partition);
if partition.N == 1
  children = partition;
  return
end
dir = divide_by_principal_component(partition.given_data(:,partition.indices), ...
                                    partition.sigma, ...
                                    partition.mean);
positive_I = find(dir>=0);
negative_I = find(dir<0);
if length(positive_I) == 0 | length(negative_I) == 0
  children = partition;
  return
end
children(1) = mk_kdtree_partition(partition.given_data, partition.indices(positive_I));
children(2) = mk_kdtree_partition(partition.given_data, partition.indices(negative_I));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function direction = divide_by_principal_component(data, covariance, mean);
N = size(data, 2);
if size(data,1) <= 16
  [V,D] = eig(covariance);
  [eig_val, principal_component_i] = max(diag(D));
  principal_component = V(:,principal_component_i);
else
  [principal_component,eig_val] = power_method(covariance);
end
direction = sum((data - repmat(mean, 1, N)).*repmat(principal_component, 1, N), 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partition = mk_kdtree_partition(given_data, indices);
partition.N = length(indices);
if partition.N == 0
  error('indices must have at least one index.')
end
data_a = given_data(:,indices);
partition.given_data = given_data;
partition.indices = indices;
partition.sum_x = sum(data_a, 2); % D*1
mean = partition.sum_x / partition.N; % D*1
partition.mean = mean;
partition.sum_xx = data_a*data_a'; % D*D
partition.sigma = partition.sum_xx / partition.N - mean*mean';
partition.children = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q_of_z = rand_q_of_z(data, K, opts);
% q_of_z: N*K
if opts.use_kd_tree
  N = length(data.kdtree);
else
  N = size(data.given_data, 2);
end
if isequal(opts.algorithm, 'vdp')
  q_of_z = zeros(N, K+1);
else
  q_of_z = zeros(N, K);
end
q_of_z(:,1:K) = rand(N, K);
q_of_z = normalize(q_of_z, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hp_posterior = mk_hp_posterior(data, q_of_z, hp_prior, opts);
% the last component of q_of_z represents the infinite rest of components
% the last component is the prior.
% q_of_z: N*K
% q_of_z(:,end) is the rest of responsibilities.
threshold_for_N = 1.0e-200;
K = size(q_of_z, 2);
if opts.use_kd_tree
  N = length(data.kdtree);
  D = size(data.kdtree(1).sum_x, 1);
  Na = [data.kdtree(:).N];
  if isequal(opts.algorithm, 'vdp')
    true_Nc = Na*q_of_z; % 1*K
    q_of_z(:,end) = 0;
  end
  Nc = Na*q_of_z; % 1*K
  sum_x = [data.kdtree(:).sum_x] * q_of_z;
else
  [D,N] = size(data.given_data);
  if isequal(opts.algorithm, 'vdp')
    true_Nc = sum(q_of_z, 1); % 1*K
    q_of_z(:,end) = 0;
  end
  Nc = sum(q_of_z, 1); % 1*K
  sum_x = data.given_data * q_of_z;
end
I = find(Nc>threshold_for_N);
inv_Nc = zeros(1,K);
inv_Nc(I) = 1./Nc(I);
hp_posterior.eta = hp_prior.eta0 + Nc;
hp_posterior.xi = hp_prior.xi0 + Nc;
means = sum_x .* repmat(inv_Nc, D, 1); % D*K
hp_posterior.inv_B = zeros(D,D,K);
if opts.use_kd_tree
  t1 = reshape([data.kdtree(:).sum_xx],D,D,N);
  for c=1:K
    v0 = means(:,c) - hp_prior.m0;
    q_of_z_c = reshape(q_of_z(:,c), 1, 1, N);
    S = sum(repmat(q_of_z_c,[D,D,1]).*t1, 3) - Nc(c)*means(:,c)*means(:,c)';
    hp_posterior.B{c} = hp_prior.B0 ...
        + S ...
        + Nc(c)*hp_prior.xi0*v0*v0'/(hp_posterior.xi(c)); % D*D
  end
else
  if opts.collapsed_means
    quad_term = mk_quad_term(data, q_of_z, hp_prior, opts);
  end
  for c=1:K
    if opts.collapsed_means
      hp_posterior.B{c} = hp_prior.B0 ...
          + (repmat(q_of_z(:,c),1,D)'.*data.given_data)*data.given_data' ...
          + quad_term(:,:,c);
    else
      v = data.given_data - repmat(means(:,c),1,N); % D*N
      v0 = means(:,c) - hp_prior.m0;
      hp_posterior.B{c} = hp_prior.B0 ...
          + (repmat(q_of_z(:,c),1,D)'.*v)*v' ...
          + Nc(c)*hp_prior.xi0*v0*v0'/(hp_posterior.xi(c)); % D*D
    end
  end
end
for c=1:K
  hp_posterior.inv_B(:,:,c) = inv(hp_posterior.B{c});
end
hp_posterior.m = (sum_x + repmat(hp_prior.xi0*hp_prior.m0,1,K)) ...
    ./ repmat(Nc+hp_prior.xi0,D,1);
if isequal(opts.algorithm, 'vdp')
  % gamma: 2*K
  hp_posterior.gamma = zeros(2,K);
  hp_posterior.gamma(1,:) = 1 + true_Nc;
  hp_posterior.gamma(2,:) = hp_prior.alpha + sum(true_Nc) - cumsum(true_Nc,2);
elseif isequal(opts.algorithm, 'bj')
  hp_posterior.gamma = zeros(2,K-1);
  hp_posterior.gamma(1,:) = 1 + Nc(1:K-1);
  hp_posterior.gamma(2,:) = hp_prior.alpha + sum(Nc) - cumsum(Nc(1:K-1),2);
elseif isequal(opts.algorithm, 'non_dp') | isequal(opts.algorithm, 'cdp')
  hp_posterior.tilde_alpha = hp_prior.alpha/K + Nc;
elseif isequal(opts.algorithm, 'csb')
  1;
else
  error('unknown algorithm')
end

hp_posterior.Nc = Nc; 
if isequal(opts.algorithm, 'vdp')
  hp_posterior.true_Nc = true_Nc;
end
hp_posterior.q_of_z = q_of_z; % q_of_z is a N by K matrix where N is
                              % #given_data or #nodes in a kd-tree.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc = mk_E_log_q_p_eta(data, hp_posterior, hp_prior, opts);
% returns E[ log q(eta)/p(eta) ]_q
% fc : 1 by K
D = size(hp_posterior.m, 1);
K = size(hp_posterior.eta, 2);
log_det_B = zeros(1,K);
for c=1:K
  log_det_B(c) = detln(hp_posterior.B{c});
  d = hp_posterior.m(:,c)-hp_prior.m0; % D*1
  term_eta(1,c) = sum(sum(hp_posterior.inv_B(:,:,c).*(hp_prior.xi0*d*d'),1),2);
  term_eta(2,c) = sum(sum(hp_posterior.inv_B(:,:,c).*hp_prior.B0,1),2) - D;
end
E_log_q_p_mean = ...
    + 0.5*D*(hp_prior.xi0./hp_posterior.xi ...
             - log(hp_prior.xi0./hp_posterior.xi) ...
             - 1) ...
    + 0.5*(hp_posterior.eta).* term_eta(1,:);     

psi_sum = sum(psi( (repmat(hp_posterior.eta+1,D,1) - repmat([1:D]',1,K))*0.5 ), 1); % 1*K
E_log_q_p_cov = ...
    0.5*hp_prior.eta0*(log_det_B-detln(hp_prior.B0)) ...
    + 0.5*hp_posterior.Nc.*psi_sum ...
    + 0.5*(hp_posterior.eta).* term_eta(2,:) ...
    + gamma_multivariate_ln(hp_prior.eta0*0.5,D) ...
    - gamma_multivariate_ln(hp_posterior.eta*0.5,D);

%debug
if length(find(E_log_q_p_mean<-1.0e-8,1)) > 0
  E_log_q_p_mean
  error('E_log_q_p_mean is negative.')
end
if length(find(E_log_q_p_cov<-1.0e-8,1)) > 0
  E_log_q_p_cov
  error('E_log_q_p_mean is negative.')
end

fc = E_log_q_p_mean + E_log_q_p_cov;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [free_energy, log_lambda] = mk_free_energy(data, hp_posterior, ...
                                                  hp_prior, opts, ...
                                                  fc, log_lambda);
if nargin == 4
  fc = mk_E_log_q_p_eta(data, hp_posterior, hp_prior, opts); % 1*K
  log_lambda = mk_log_lambda(data, hp_posterior, hp_prior, opts); % N*K
end
[N,K] = size(log_lambda);
if isequal(opts.algorithm, 'vdp') || isequal(opts.algorithm, 'bj')
  % note when bj,  if full hp_posterior is given, len_gamma = K - 1.
  len_gamma = size(hp_posterior.gamma, 2);
  if isequal(opts.algorithm, 'bj') && len_gamma ~= K - 1
    error('invalid length')
  end
  E_log_p_of_V = ...
      gammaln(sum(hp_posterior.gamma, 1)) ...
      - gammaln(1+hp_prior.alpha) ...
      - sum(gammaln(hp_posterior.gamma), 1) ...
      + gammaln(hp_prior.alpha) ...
      + ((hp_posterior.gamma(1,:)-1) ...
         .*(psi(hp_posterior.gamma(1,:))-psi(sum(hp_posterior.gamma,1)))) ...
      + ((hp_posterior.gamma(2,:)-hp_prior.alpha) ...
         .*(psi(hp_posterior.gamma(2,:))-psi(sum(hp_posterior.gamma,1))));
  extra_term = sum(E_log_p_of_V);
elseif isequal(opts.algorithm, 'non_dp')
  E_log_p_of_pi = ...
      sum(gammaln(hp_prior.alpha/K) ...
          - gammaln(hp_posterior.tilde_alpha) ...
          + hp_posterior.Nc.*(psi(hp_posterior.tilde_alpha) ...
                              - psi(sum(hp_posterior.tilde_alpha)))) ...
      + gammaln(sum(hp_prior.alpha+size(data,2))) - gammaln(hp_prior.alpha);
  extra_term = E_log_p_of_pi;
elseif isequal(opts.algorithm, 'cdp') || isequal(opts.algorithm, 'csb')
  E_log_p_of_z_given_other_z = ...
      mk_E_log_p_of_z_given_other_z(hp_posterior, hp_prior, opts); % N by K
  q_of_z = mk_q_of_z(data,hp_posterior,hp_prior,opts,log_lambda);
%   q_of_z = hp_posterior.q_of_z;
  E_Nc = hp_posterior.Nc;
  V_Nc = sum(q_of_z.*(1-q_of_z), 1); % 1 by K
  if isequal(opts.algorithm, 'cdp')
    E_log_p_of_z = gammaln(hp_prior.alpha) - gammaln(N+hp_prior.alpha) ...
        + sum(gammaln(hp_prior.alpha/K + E_Nc) ...
              - 0.5*psi(1, hp_prior.alpha/K + E_Nc).*V_Nc - gammaln(hp_prior.alpha/K));
  else % csb
    E_Nc_geq_to_i = cumsum(E_Nc, 2);
    q_of_z_geq_to_i = cumsum(q_of_z, 2);
    V_Nc_geq_to_i = sum(q_of_z_geq_to_i.*(1-q_of_z_geq_to_i), 1);

    E_Nc_greater_than_i = cumsum(E_Nc, 2) - E_Nc;
    q_of_z_greater_than_i = q_of_z_geq_to_i - q_of_z;
    V_Nc_greater_than_i = sum(q_of_z_greater_than_i.*(1-q_of_z_greater_than_i), 1);
    
    tmp = gammaln(1 + E_Nc) - 0.5*psi(1, 1 + E_Nc).*V_Nc ...  % E[log p(1+Nc)]
          + (gammaln(1 + E_Nc_greater_than_i) ...             % E[log p(1+N_{>c})]
             - 0.5*psi(1, 1 + E_Nc_greater_than_i).*V_Nc_greater_than_i) ...
          - (gammaln(1 + hp_prior.alpha + E_Nc_geq_to_i) ...  % E[log p(1+alpha+N_{>=c})]
             - 0.5*psi(1, 1 + hp_prior.alpha + E_Nc_geq_to_i).*V_Nc_geq_to_i);
    E_log_p_of_z = sum(tmp(1:end-1));
  end
  extra_term = sum(sum(E_log_p_of_z_given_other_z.*q_of_z, 2), 1) - E_log_p_of_z;
else
  error('unknown algorithm')
end
if opts.use_kd_tree
  free_energy = extra_term + sum(fc) - [data.kdtree(:).N]*log_sum_exp(log_lambda, 2);
else
  free_energy = extra_term + sum(fc) - sum(log_sum_exp(log_lambda, 2), 1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function log_lambda = mk_log_lambda(data, hp_posterior, hp_prior, opts);
% log_lambda: N*K
% q(z_n=c|x_n) = lambda_n_c / sum_c lambda_n_c

if isequal(opts.algorithm, 'vdp')
  if abs(hp_posterior.gamma(2,end) - hp_prior.alpha) > 1.0e-5
    hp_posterior.gamma(2,end)
    hp_prior.alpha
    diff = hp_prior.alpha - hp_posterior.gamma(2,end)
    error('must be alpha')
  end
end

if opts.use_kd_tree
  N = length(data.kdtree);
  D = size(data.kdtree(1).sum_x, 1);
else
  [D,N] = size(data.given_data);
end
K = size(hp_posterior.eta, 2);

psi_sum = sum(psi( (repmat(hp_posterior.eta+1,D,1) - repmat([1:D]',1,K))*0.5 ), 1); % 1*K
log_lambda = zeros(N,K);
if opts.use_kd_tree
  t1 = reshape([data.kdtree(:).sum_xx],D,D,N);
  sum_x = repmat(reshape([data.kdtree(:).sum_x],D,1,N),[1,D,1]);
  Na = repmat(reshape([data.kdtree(:).N],1,1,N),[D,D,1]);
end
if isequal(opts.algorithm, 'csb') | isequal(opts.algorithm, 'cdp')
  E_log_p_of_z_given_other_z = mk_E_log_p_of_z_given_other_z(hp_posterior, hp_prior, opts);
end
for c=1:K
%   Precision = 0.5*hp_posterior.inv_B(:,:,c)*hp_posterior.eta(c);
%   if isequal(opts.algorithm, 'vdp')
%     constant = psi(hp_posterior.gamma(1,c)) - psi(sum(hp_posterior.gamma(:,c),1)) ...
%         + sum(psi(hp_posterior.gamma(2,[1:c-1])) - psi(sum(hp_posterior.gamma(:,[1:c-1]),1)), 2);
%   elseif isequal(opts.algorithm, 'bj')
%     if c < K
%       constant = psi(hp_posterior.gamma(1,c)) - psi(sum(hp_posterior.gamma(:,c),1)) ...
%           + sum(psi(hp_posterior.gamma(2,[1:c-1])) - psi(sum(hp_posterior.gamma(:,[1:c-1]),1)), 2);
%     else
%       constant = sum(psi(hp_posterior.gamma(2,[1:c-1])) - psi(sum(hp_posterior.gamma(:,[1:c-1]),1)), 2);
%     end
%   elseif isequal(opts.algorithm, 'non_dp')
%     % E[log pi] ; pi is the weight of mixtures.
%     constant = psi(hp_posterior.tilde_alpha(c)) - psi(sum(hp_posterior.tilde_alpha));
%   elseif isequal(opts.algorithm, 'cdp')
%     constant = 0;
%   elseif isequal(opts.algorithm, 'csb')
%     constant = 0;
%   else
%     error('unknown algorithm')
%   end
%   constant = constant - 0.5*D*log(pi) - 0.5*detln(hp_posterior.B{c}) ...
%       + 0.5*psi_sum(c) - 0.5*D/(hp_posterior.xi(c));
%   if opts.use_kd_tree
%     t2 = sum_x.*repmat(hp_posterior.m(:,c)',[D,1,N]);
%     term_dependent_on_n = (t1 - t2 - permute(t2,[2,1,3]))./Na + ...
%         repmat(hp_posterior.m(:,c)*hp_posterior.m(:,c)',[1,1,N]);
%     column_c = - sum(sum(repmat(Precision,[1,1,N]).*term_dependent_on_n,2),1) + constant;
%   else
%     d = data.given_data - repmat(hp_posterior.m(:,c),1,N);
%     column_c = - sum(d.*(Precision*d),1) + constant; % 1*N
%   end
%   if isequal(opts.algorithm, 'csb') | isequal(opts.algorithm, 'cdp')
%     column_c = column_c + E_log_p_of_z_given_other_z(:,c)';
%   end
%   log_lambda(:,c) = column_c;
  if isequal(opts.algorithm, 'vdp')
    E_log_p_of_z_given_other_z_c = ...
        psi(hp_posterior.gamma(1,c)) ...
        - psi(sum(hp_posterior.gamma(:,c),1)) ...
        + sum(psi(hp_posterior.gamma(2,[1:c-1])) - psi(sum(hp_posterior.gamma(:,[1:c-1]),1)), 2);
  elseif isequal(opts.algorithm, 'bj')
    if c < K
      E_log_p_of_z_given_other_z_c = ...
          psi(hp_posterior.gamma(1,c)) ...
          - psi(sum(hp_posterior.gamma(:,c),1)) ...
          + sum(psi(hp_posterior.gamma(2,[1:c-1])) - psi(sum(hp_posterior.gamma(:,[1:c-1]),1)), 2);
    else
      E_log_p_of_z_given_other_z_c = sum(psi(hp_posterior.gamma(2,[1:c-1])) ...
                                         - psi(sum(hp_posterior.gamma(:,[1:c-1]),1)), 2);
    end
  elseif isequal(opts.algorithm, 'non_dp')
    % E[log pi] ; pi is the weight of mixtures.
    E_log_p_of_z_given_other_z_c = psi(hp_posterior.tilde_alpha(c)) ...
        - psi(sum(hp_posterior.tilde_alpha));
  elseif isequal(opts.algorithm, 'csb') | isequal(opts.algorithm, 'cdp')
    E_log_p_of_z_given_other_z_c = E_log_p_of_z_given_other_z(:,c)';
  else
    error('unknown algorithm')
  end

  Precision = 0.5*hp_posterior.inv_B(:,:,c)*hp_posterior.eta(c);
  E_log_p_of_x = - 0.5*D*log(pi) - 0.5*detln(hp_posterior.B{c}) ...
      + 0.5*psi_sum(c) - 0.5*D/(hp_posterior.xi(c));
  if opts.use_kd_tree
    t2 = sum_x.*repmat(hp_posterior.m(:,c)',[D,1,N]);
    term_dependent_on_n = (t1 - t2 - permute(t2,[2,1,3]))./Na + ...
        repmat(hp_posterior.m(:,c)*hp_posterior.m(:,c)',[1,1,N]);
    E_log_p_of_x = - sum(sum(repmat(Precision,[1,1,N]).*term_dependent_on_n,2),1) + E_log_p_of_x;
  else
    d = data.given_data - repmat(hp_posterior.m(:,c),1,N);
    E_log_p_of_x = - sum(d.*(Precision*d),1) + E_log_p_of_x; % 1*N
  end
  log_lambda(:,c) = E_log_p_of_x + E_log_p_of_z_given_other_z_c;
end
if isequal(opts.algorithm, 'vdp')
  log_lambda(:,end) = log_lambda(:,end) - log(1- exp(psi(hp_prior.alpha) - psi(1+hp_prior.alpha)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q_of_z, data, log_lambda] = mk_q_of_z(data, hp_posterior, hp_prior, opts, log_lambda);
% q_of_z: N*K
if nargin == 4
  log_lambda = mk_log_lambda(data, hp_posterior, hp_prior, opts);
end
q_of_z = exp(normalizeln(log_lambda, 1));
if opts.weight ~= 1
  q_of_z = q_of_z .* repmat(opts.weight, 1, size(q_of_z,2));
else
  q_of_z = q_of_z;
end

function M =normalizeln(M ,dimension);
M = lpt2lcpt(M, dimension);

function lcpt=lpt2lcpt(lpt,dimension);
% function lcpt=lpt2lcpt(lpt,dimension);
%
% make a log conditional probability table with log probability table
%
% lpt is a matrix.
% dimension must be either 1 or 2.
%
% lcpt = log( exp(lpt) / sum(exp(lpt)) )
%      = lpt - log(sum(exp(lpt)))
%

the_other_dimension=-(dimension-1.5)+1.5;
lpt=permute(lpt,[dimension,the_other_dimension]);
% now we can calculate as if dimension=1.

log_sum_exp_lpt = log_sum_exp(lpt,2); % Mx1
lcpt = lpt - repmat(log_sum_exp_lpt,1,size(lpt,2));

lcpt=permute(lcpt,[dimension,the_other_dimension]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hp_prior = mk_hp_prior(data, opts)
if isfield(opts, 'xi0')
  hp_prior.xi0 = opts.xi0;
else
  hp_prior.xi0 = 0.01;
end
if isfield(opts, 'eta_p')
  eta_p = opts.eta_p;
else
  eta_p = 1;
end
if isfield(opts, 'alpha')
  hp_prior.alpha = opts.alpha;
else
  hp_prior.alpha = 1;
end
[D,N] = size(data.given_data);
if opts.use_kd_tree
  sum_xx = sum(reshape([data.kdtree(:).sum_xx],D,D,length(data.kdtree)),3);
  sum_x = sum([data.kdtree(:).sum_x],2);
  covariance = sum_xx/N - sum_x*sum_x'/(N*N);
else
  covariance = cov(data.given_data');
end
hp_prior.m0 = mean(data.given_data, 2);
if D > 16
  [dummy, max_eig] = power_method(covariance);
else
  max_eig = max(eig(covariance));
end
hp_prior.eta0 = eta_p * D + 1;
hp_prior.B0 = hp_prior.eta0 * max_eig * eye(D) * hp_prior.xi0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec,value]=power_method(A, start, precision)

switch nargin
 case 1
  start = ones(length(A),1);
  precision = 1e-10;
 case 2
  precision = 1e-10;
 case 3
 otherwise
  error 'invalid number of arguments'
end

%
% xp = A^p * x
% x(p+1) = lambda * xp   when p is big enough
% lambda = x(p+1)'*xp / xp'*xp
% 
diff=precision+1;
x=start;
n=norm(x)+diff;
i = 0;
while diff> precision
  i = i + 1;
  y=A*x;
  n2 = norm(x);
  diff=abs(n2-n);
  n=n2;
  if n < 1.0e-200
    x = zeros(length(A), 1);
    break
  else
    x=y/n;
  end
  if i > 100
    break
  end
end
n = norm(x);
if n < 1.0e-200
  vec = zeros(length(A), 1);
else
  vec = x/n;
end
value=n;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m=normalize(m,dim)
% return m normalized with 'dimension'
%
% e.g.
% m : i by j by k by ...
% m = sum(normalize(m,2), 2);
% m(i, :, k, ...) = ones(1, J, 1, ...)
%

dims = ones(1, ndims(m));
dims(dim) = size(m, dim);
m = m ./ repmat(sum(m, dim), dims);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y] = detln( X )
% Y = logdet( X )
% return a log determinant of X
[d err] = chol(X);
if err
  error('error in Choleski disposition for detln');
end
Y = sum(log(diag(d))) *2;

% Local Variables: ***
% mode: matlab ***
% End: ***
