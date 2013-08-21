%   [m,ga,b,z,nlZ] = DLI(X,y,s2,B,pot,tau,opts,G)
%
% DLI - Double Loop Inference - Bayesian Inference using a Double Loop algorithm
%       The variational criterion whose stationary point is found has the form
%         phi(ga,b) = ln|A| + h(ga,b) + \min_u R(u,ga,b), where
%         A         = X'*X/s2 + B'*diag(1./ga)*B, s=B*u,
%         h(ga,b)   = \sum_{j=1}^q h_j(ga_j,b_j)
%         R(u,ga,b) = (1/s2) * ||X*u-y||^2 + s'*diag(1./ga)*s - 2*b'*s.
%
% The probabilistic model consists of
%   (i)  Gaussian potentials y = X*u + e, e~N(0,s2*I)   and 
%   (ii) non-Gaussian potentials t_i(s_i), pot(s) = \prod_i t_i(tau_i*s_i)
% leading to a posterior P(u) of the form
%    P(u) = (1/Z) N(u|Xu,s2*I) * pot(s) where s = B*u and the normalisation
%    constant or partition function is Z = \int  N(u|Xu,s2*I) * pot(s) du.
%
% OUTER LOOP: z \approx diag(B'*inv(A)*B) = var(s)
%
% opts.innerType == 'VB'
% INNER LOOP: min_u  (1/s2)*||X*u-y||^2 + 2*pen(s) where r = sign(s)*sqrt(s^2+z)
%              pen(s) = tau*b*(r-s) - ln pot(tau*r), see pen/penVB.m.
%   Our variational inference relaxation uses Gaussian individual lower bounds 
%   pot(s) \ge exp( b*s -s^2/ga +h(ga)/2 ) to obtain a joint a lower bound
%     Z \ge exp(h(ga)/2) * \int  N(u|Xu,s2*I) * exp( b*s -s^2/ga ) du = Z_{VB}
%   having the form of a Gaussian integral that can be written as
%     Z_{VB} \ge C*exp( -phi(ga)/2 ), C = (2*pi)^(n/2) * (2*pi*s2)^(-m/2).
%   Here, the variational criterion phi(ga) is to be minimised w.r.t. the
%   variational parameters ga. For log-concave models, this constitutes a convex
%   variational optimisation problem (1).
%   Our approach to solve the problem uses Fenchel duality to decouple ln|A| by
%   upper bounding ln|A| \min_z z'*(1./ga) g*(z); a function that is just a sum
%   over individual components of ga. The algorithm consists of an outer loop 
%   where the decoupling bound with coefficients z is refit and an inner loop
%   where the variational criterion is jointly minimised in both u and ga 
%   (wich reduces to a PLS problem, pls/pls*.m).
%
% opts.innerType == 'EP'
% INNER LOOP:
%   We interleave parallel EP updating steps and posterior mean recomputations
%   to find a stationary point of the inner loop criterion.
%
% INPUT
% =====
%  X   [mxn]  measurement matrix or operator
%  y   [mx1]  measurement vector
%  s2  [1x1]  measurement variance
%  B   [qxn]  matrix or operator
%  pot        potential function handle or function name string from pot/pot*.m
%  tau [qx1]  scale parameters of the potentials
%  opts.      optimisation parameters
%       outerZinit   initial value for upper bound             [default 0.05]
%       outerGainit  initial value for variational parameter   [default 1]
%       outerNiter   number of outer loop iterations           [default 10]
%       outerMVM     number of MVMs/Lanczos vectors            [default 50]
%       outerExact   Do we try to compute exact outer loops?   [default false]
%                    If B is the identity and X numeric, we use the Woodbury 
%                    formula, otherwise A is expanded into a full matrix.
%                    This flag is also used for the posterior mean computation
%                    at the end and in the EP inner loop.
%       outerOutput  flag saying whether some output is shown  [default false]
%       innerMVM     number of inner loop MVMs or CG steps     [default 50]
%       innerIt      number of inner Newton steps for plsTN    [default 15]
%       innerOutput  flag saying whether some output is shown  [default false]
%       innerType    which inference method EP or VB           [default 'VB']
%       innerVBpls   inner loop PLS solver                  [default 'plsLBFGS']
%       => the opts struct is passed to the PLS solver which
%          allows to pass additional arguments to it
%       innerEPeta   scalar parameter for fractional EP        [default 1]
%                    only possible for potGauss and potLaplace
%  G   [qtxq] grouping matrix (for VB only)                     [default eye(q)]
%             this changes the inner loop penaliser to
%             pen(s) = -log( pot(tau*r) ) where r = sqrt( G*(s^2 + z) ).
%
% OUTPUT
% ======
%  m   [nx1]  posterior mean estimate                m = mean(u)
%  ga  [qx1]  optimal value of the variational width parameters
%  b   [qx1]  optimal value of the variational position parameters
%  z   [qx1]  posterior marginal variance estimates  z = var(s)
%  nlZ [1x1]  approxmation to the negative log marginal likelihood -log(Z)
%
%   See also PENFUNCTIONS.M, POTFUNCTIONS.M, INF/DLI.M.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 23