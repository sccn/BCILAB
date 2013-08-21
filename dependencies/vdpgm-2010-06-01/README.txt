Variational Dirichlet Process Gaussian Mixture Models
written by Kenichi Kurihara
distributed under the modified BSD license.



ALGORITHMS (see Refrences, the end of this document)

1. accelerated variational Dirichlet process Gaussian mixture model
2. collapsed variational stick-breaking Dirichlet process Gaussian mixture model
3. variational Gaussian mixture model with a collapsed Dirichlet prior.
4. variational Dirichlet process Gaussian mixture model by Blei and Jordan.


USAGE

>> result = vdpgm(X, opts);

The first argument is data. Each data point is a column vector of X.

The second argument opts is the option of this program which
determines an algorithm and hyperparameters. You can set opts as you
want, or basic option generators are also available.

>> opts = mkopts_avdp;     % for the algorithm 1
>> opts = mkopts_csb(10);  % for the algorithm 2 with T=10 
>> opts = mkopts_cdp(10);  % for the algorithm 3 with K=10 
>> opts = mkopts_bj(10);   % for the algorithm 4 with T=10 
  
Although opts accepts many options, some options are exclusive.

The output result is a structure containing parameters for posteriors.
Maybe, the most useful result is result.q_of_z which is the posterior
probability of assignments. q_of_z is a N by K (or T) matrix
s.t. sum_c q_of_z(i,c) = 1 for any c. q_of_z is available only when
opts.get_q_of_z is set to 1.

Other useful stats:
- The expected value of the covariance of component c,
>> result.hp_posterior.B{c} / results.hp_posterior.eta(c)

- The expected value of the centroid of component c,
>> result.hp_posterior.m(:,c)

One may want to know the number of discovered clusters.
If opts.algorithms is 'vdp', it is 
>> results.K - 1

Otherwise, results.K is initialized by opts, and does not change.  In
these cases, the number of clusters is K s.t.
result.hp_posterior.Nc(K+1) is close enough to zero.



REFERENCES

    * Kenichi Kurihara, Max Welling and Yee Whye Teh,
      Collapsed Variational Dirichlet Process Mixture Models,
      the Twentieth International Joint Conference on Artificial Intelligence (IJCAI 2007). 

    * Kenichi Kurihara, Max Welling and Nikos Vlassis,
      Accelerated Variational Dirichlet Mixture Models,
      Advances in Neural Information Processing Systems 19 (NIPS 2006). 

    * David M. Blei and Michael I. Jordan,
      Variational Inference for Dirichlet Process Mixtures,
      Bayesian Analysis, Vol.1, No.1, 2005.
