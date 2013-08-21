% A potential pot(s) is a scalar positive function IR -> IR_+.
% Applied to a vector valued input s, a potential function is understood as 
% pointwise evaluation of every component s(i).
%
% We currently support two classes of inference algorithms beyond Laplace's
% approximation aka MAP estimation: Variational bounding and expectation 
% propagation. They are specified by the type parameter in the general call
%   P = pot(s,theta,type,z)
% where theta contains additional (hyper)parameter(s). The input s is only
% used as s(:), that is, matrix structure is ignored.
%
% 1) Variational bounding (VB)  [for ExpPow, Gauss, Laplace, Logistic, Sech2, T]
%   P = pot(s,theta) = pot(s,theta,'VB')
%   - pot(s) Needs to be strongly super-Gaussian.
%   - pot(s) Has to be symmetric or symmetrisable i.e. there is a constant b
%            such that f(s) = pot(s)*exp(-b*s) is even or symmetric f(s)=f(-s).
%   The return argument is a matrix P of size [qx4] whose columns 
%   [lp,dlp,d2lp,b] are as long as there are elements in the input s.
%     lp(s)   = log( pot(s) ),                    [simple log evaluation]
%     dlp(s)  = d   lp(s) / ds,                   [first  derivative of the log]
%     d2lp(s) = d^2 lp(s) / ds^2, and             [second derivative of the log]
%     b s.t. f(s)=f(-s), f(s) = pot(s)*exp(-b*s). [linear symmetry parameter]
%
% 2) Expectation propagation (EP)          [for Gauss, Laplace, Logistic, Sech2]
%   P = pot(s,theta,'EP',z)
%   - pot(s) Does not need to be strongly super-Gaussian or symmetric; we need
%            Gaussian expectations w.r.t. pot(s) in order to run EP.
%   The return argument is a matrix P of size [qx3] whose columns 
%   [lZ,dlZ,d2lZ] are as long as there are elements in the input s.
%     lZ(s,z)   = log( \int N(t|s,z) pot(t) dt ),       [log partition function]
%     dlZ(s,z)  = d   lZ(s,z) / ds, and           [first  derivative of the log]
%     d2lZ(s,z) = d^2 lZ(s,z) / ds^2.             [second derivative of the log]
%
% We use a single matrix-valued output argument to allow very compact 
% specifications of mixed potentials, see the last line of the examples below. 
% This ease of specification, however, comes at the expense of some computional
% overhead since we always evaluate the full matrix P even if only the first
% column is needed.
%
% Currently, the following potentials are implemented in pot/pot<NAME>.m:
%   <NAME>         Analytic Expression     Param  Description
%   ----------------------------------------------------------------------------
%   ExpPow(s,al) = exp(-|s|^al)            al>0   Exponential power distribution
%   Gauss(s)     = exp(-s^2/2)                    Gaussian          distribution
%   Laplace(s)   = exp(-abs(s))                   Laplace           distribution
%   Logistic(s)  = 1/(1+exp(-s))                  Logistic          function
%   Sech2(s)     = 1 / cosh(s)^2                  Sech-squared      distribution
%   T(s,nu)      = (1+s^2/nu)^(-nu/2-1/2)  nu>0   Student's t       distribution 
%
% Examples:
% >> s=0.3;   pot = @(s) potLogistic(s,'EP');
% >> s=2.7;   nu=3; pot = @(s) potT(s,'VB',nu);
% >> s=[1;2]; al=4; pot = @(s) [potGauss(s(1),'VB'); potExpPow(s(2),'VB',al)];
%
% Another way of concatenating potentials is provided by pot/potCat.m, which
% allows concatenation by:
% >> pot = @(s,varargin) potCat(s,varargin{:},{@potT,@potGauss},{1:5,6:9});
%
%   See also PENFUNCTIONS.M.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 10