function varargout = sim_simulateSources(varargin)
% Project source amplitudes through a forward model to generate channel data.
%
% This function is a wrapper for sim_fwdProj() and maintained for backwards-
% compatibility. This function may be deprecated in a future release.
% 
% Author: Tim Mullen, SCCN/INC/UCSD, 2013

varargout = sim_fwdProj(varargin{:});
