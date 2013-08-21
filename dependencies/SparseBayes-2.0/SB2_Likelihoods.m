% SB2_LIKELIHOODS  Convenience function to encapsulate likelihood types
% 
% LIKELIHOOD = SB2_LIKELIHOODS(TYPE)
% 
% OUTPUT ARGUMENTS:
% 
%	LIKELIHOOD	A structure containing an enumeration of possible
%				likelihood types, and specification of the type in use.
% 
% INPUT ARGUMENTS:
% 
%	TYPE		A string. One of:
%					'Gaussian'
%					'Bernoulli'
%					'Poisson'
% 
% EXAMPLE:
%
%	LIKELIHOOD	= SB2_Likelihoods('Gaussian');
%
% NOTES:
% 
% This function is "for internal use only" and exists purely for convenience
% of coding and to facilitate potential future expansion.
%

%
% Copyright 2009, Vector Anomaly Ltd
%
% This file is part of the SPARSEBAYES library for Matlab (V2.0).
%
% SPARSEBAYES is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.
%
% SPARSEBAYES is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along
% with SPARSEBAYES in the accompanying file "licence.txt"; if not, write to
% the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
% MA 02110-1301 USA
%
% Contact the author: m a i l [at] m i k e t i p p i n g . c o m
%
function LIKELIHOOD = SB2_Likelihoods(likelihood_)

% For convenience in code, and later expandability,
% we "enum" the available likelihood models
%
% Continuous
LIKELIHOOD.Gaussian		= 1;
% Discrete
LIKELIHOOD.Bernoulli	= 2;
LIKELIHOOD.Poisson		= 3;

% Feel free to add your own ... and work out the maths :-)

% Determine the likelihood to be used
% 
LIKELIHOOD.InUse		= 0;
%
switch upper(likelihood_(1:5))
 case 'GAUSS',
  LIKELIHOOD.InUse	= LIKELIHOOD.Gaussian;
  %
 case 'BERNO',
  LIKELIHOOD.InUse	= LIKELIHOOD.Bernoulli;
  %
 case 'POISS',
  LIKELIHOOD.InUse	= LIKELIHOOD.Poisson;
  %
 otherwise
  error('Unknown likelihood ''%s'' specified\n', likelihood_)
end

