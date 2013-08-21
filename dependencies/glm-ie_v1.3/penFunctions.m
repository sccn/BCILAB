% A penalty function pen(s) is a function IR^q -> IR.
%
%   [p,dp,d2p] = pen(s,psi) where psi contains additional parameters
%
% The return arguments have the following meaning:
%   p   = pen(s),
%   dp  = d   pen(s) / ds,
%   d2p = d^2 pen(s) / ds^2,
% and the following dimensions:
%   p     function value   [1x1] or vector of individual function values [qx1], 
%   dp    gradient vector  [qx1],
%   d2p   Hessian diagonal [qx1] or full matrix [qxq].
%
% Currently, we have implemented in pen/pen<NAME>.m
%   Absolute value:          penAbs(s)               = abs(s),     => potLaplace
%   Smoothed absolute value: penAbsSmooth(s,ep)      = sqrt(s^2+ep),
%   Linear pen. on negative: penNegLin(s)            = max(-s,0),
%   Logarithmic:             penLogSmooth(s,nu)      = log(s^2+nu),=> potT
%   Power:                   penPow(s,d)             = abs(s)^d,   => potExpPow
%   Smoothed power:          penPowSmooth(s,d,ep)    = sqrt(s^2+ep)^d,
%   Quadratic:               penQuad(s)              = s^2/2,      => potGauss
%   Quadr. pen on negative:  penNegQuad(s)           = min(s,0)^2/2,
%   Zero:                    penZero(s)              = 0, and
%
%   Derived from a potential for use with Variational bounding (VB) inference:  
%     penVB(s,pot,tau,z) = tau*b*(r-s) -log(pot(tau*r)), r=sign(s)*sqrt(s^2+z)
%   Norm penalty derived from symmetric potential to be used with VB inference:
%     penVBNorm(s,pot,tau,z,G) = -log( pot(tau*r) ), r=sqrt( G*(s^2+z) )
%
% Examples:
% >>  s=0.3; d=1.5; penPow(s,d)
% >>  s=-4;  pot='potLaplace'; tau=2; z=1.3; penVB(s,pot,tau,z)
%
%   See also POTFUNCTIONS.M.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 09