function data = hlp_rewrite(data,varargin)
% Rewrite (replace) an input data structure into an output data structure.
% Output = hlp_rewrite(Input,Rules...)
%
% In:
%   Input : some data structure to be rewritten
%
%   Rules... : a comma-separated list of rewrite rules (oldval,newval,oldval,newval,oldval,newval, ...)
%
% Out:
%   Output : the input data structure, but rewritten according to the rules (where they matched).
%
%
% Notes:
%   * No two oldval's may be equal. 
%   * If there is no match, data remains unchanged.
%
% Examples:
%   % rewrite short forms of some string into corresponding long forms
%   hlp_rewrite(myinput, 'hp','highpass', 'lp','lowpass', 'bp','bandpass')
% 
%   % rewrite true to 'on' and false to 'off'
%   hlp_rewrite(myinput, true,'on', false,'off')
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-26

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

old = varargin(1:2:end);
new = varargin(2:2:end);

if iscellstr(old)
    % mapping from strings
    match = strcmp(data,old);
    if any(match)
        data = new{match}; end
else
    % mapping from general structures
    match = cellfun(@(x)isequalwithequalnans(x,data),old); %#ok<DISEQN>
    if any(match)
        data = new{match}; end        
end
