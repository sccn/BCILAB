function arg_issuereport(payload)
% Internal function to yield a report to a requesting function. 
% This is implemented by throwing an exception. Used mainly by arg_define() and declare_properties().
%
% In:
%   Payload : the payload to report
%
% See also:
%   arg_report, arg_define

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

% first obtain a report "ticket"
% these tickets are used to uniquely retrieve a given report after it has been issued,
% in particular when multiple reports are in flight in parallel (e.g., from timers)
try
    ticket = tracking.arg_sys.tickets.removeLast();
catch %#ok<CTCH>
    max_inflight_tickets = 10000;
    % initialize data structures if necessary
    if ~isfield(tracking,'arg_sys')
        tracking.arg_sys = struct(); end
    if ~isfield(tracking.arg_sys,'tickets')
        tracking.arg_sys.tickets = java.util.concurrent.LinkedBlockingDeque();
        for t=max_inflight_tickets:-1:1
            tracking.arg_sys.tickets.addLast(t); end
    else
        warning('BCILAB ran out of report tickets; if this happens it means that tickets did not get reclaimed properly elsewhere in the system, perhaps due to an error.');
        for t=max_inflight_tickets:-1:1
            tracking.arg_sys.tickets.addLast(t); end
    end
    ticket = tracking.arg_sys.tickets.removeLast();
end

% store the payload
tracking.arg_sys.reports{ticket} = payload;

% now throw the exception
error('BCILAB:arg:report_args','This (internal) exception is destined to be caught by arg_report(); please do not interfere with it. Ticket=%05u',ticket);
