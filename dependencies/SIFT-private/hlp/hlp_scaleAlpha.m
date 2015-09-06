function hlp_scaleAlpha(handles,factor)
% increment or decrement the FaceAlpha property for a set of objects by a given factor
% Inputs: 
%   handles:    an array of object handles
%   factor :    the scaling factor
%
% Author: Tim Mullen, SCCN/INC, UCSD 2012

arrayfun(@(x) set(x,'FaceAlpha',get(x,'FaceAlpha')*factor),handles);
drawnow;