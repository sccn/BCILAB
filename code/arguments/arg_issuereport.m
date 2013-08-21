function arg_issuereport(payload)
% Internal function to issue a report to a requesting function. 
% This is implemented by throwing an exception. Used by arg_define() and declare_properties()
%
% In:
%   Payload : the payload to report
%
% See also:
%   arg_define, arg_report

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

global tracking;

% first obtain a report ticket
if ~isfield(tracking,'arg_sys')
    tracking.arg_sys = struct(); end
if ~isfield(tracking.arg_sys,'tickets')
    tracking.arg_sys.tickets = java.util.concurrent.LinkedBlockingDeque();
    for k=50000:-1:1
        tracking.arg_sys.tickets.addLast(k); end
end
ticket = tracking.arg_sys.tickets.removeLast();

% ... and store the payload
tracking.arg_sys.reports{ticket} = payload;

% now throw the exception
error('BCILAB:arg:report_args','This (internal) exception is destined to be caught by arg_report(); please do not interfere with it. Ticket=%.0f',ticket);
