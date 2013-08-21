function res = hlp_split(str,delims)
% Split a string according to some delimiter(s).
% Result = hlp_split(String,Delimiters)
%
% In:
%   String : a string (char vector)
%
%   Delimiters : a vector of delimiter characters (includes no special support for escape sequences)
%
% Out:
%   Result : a cell array of (non-empty) non-Delimiter substrings in String
%
% Examples:
%   % split a string at colons and semicolons; returns a cell array of four parts
%   hlp_split('sdfdf:sdfsdf;sfdsf;;:sdfsdf:',':;')
% 
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-05

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

pos = find(diff([0 ~sum(bsxfun(@eq,str(:)',delims(:)),1) 0]));
res = cell(~isempty(pos),length(pos)/2);
for k=1:length(res)
    res{k} = str(pos(k*2-1):pos(k*2)-1); end
