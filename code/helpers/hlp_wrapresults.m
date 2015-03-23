function a = hlp_wrapresults(f,varargin)
% Wraps all outputs produced by function f for the given arguments into a cell array.
% Results = hlp_wrapresults(Function, Arguments...)
%
% In:
%   Function     : some function handle to execute
%
%   Arguments... : list of arguments to pass to the function
%
% Out:
%   Results : cell array of all function results
%
%
% Notes: 
%   It is not (currently) possible to efficiently determine the number of out-args for a 
%   varargout function; in this case, at most 10 outputs are supported.
%
% Examples:
%   % wrap both outputs of a particualr sort() call into a cell array
%   results = hlp_wrapresults(@sort,[1 4 2 1 3 6 0])
%
% See also:
%   hlp_getresult
%
%					         	 Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

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
dp;

% find out how many arguments f can maximally produce
len = nargout(f);
% for varargout functions, we assume at most 10
if len < 0
    len = 10; end
% ... then invoke the function
while len >= 1
    try
        % get the appropriate number of results from f()
        [a{1:len}] = f(varargin{:});
        return;
    catch e
        % got an exception, check if it is outarg-related
        if ~any(strcmp(e.identifier,{'MATLAB:TooManyOutputs','MATLAB:maxlhs','MATLAB:unassignedOutputs'}))
            % it isn't: propagate the error
            settings = dbstatus;
            if any(strcmp({settings.cond},'error'))
                % if in dbstop if error mode, we're re-running the line to get the debugger to stop at the right place
                fprintf('hlp_wrapresults: caught error while in "dbstop if error" mode; re-running offending code to get break point...\n');
                fprintf('  error traceback was: %s\n',hlp_handleerror(e));
                [a{1:len}] = f(varargin{:});
            else
                rethrow(e);
            end
        end
        % but if it was, we need to retry with fewer out-arguments
        len = len-1;
    end
end

% len = 0: f produces no outputs
f(varargin{:});
a = {};
