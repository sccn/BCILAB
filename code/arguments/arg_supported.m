function result = arg_supported(func)
% Check if the argument system is supported by a given function.
% Result = arg_supported(Function)
%
% In:
%   Function : A function whose support for the argument system shall be checked.
%
% Out:
%   Result : Whether the function supports the argument system (i.e., uses arg_define).
%
% See also:
%   arg_define
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-02-03

% Copyright (C) Christian Kothe, SCCN, 2014, christian@sccn.ucsd.edu
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

if ischar(func)
    func = str2func(func); end %#ok<NASGU>

try
    [conout,result] = evalc('arg_report(''supported'',func,repmat({''__garbage__''},1,100))'); %#ok<ASGLU>
catch
    result = false;
end
