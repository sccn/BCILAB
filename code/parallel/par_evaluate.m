function result = par_evaluate(task__)
% Internal: Task processing function of the worker.
% Result = par_evaluate(Task)
%
% Both the task and the result are assumed to be in serialized form, and this function deserializes
% the task, evaluates it, and re-serializes the results.
%
% In:
%   Task : task record; A sequence of bytes which represents a cell array of the form: 
%          {Number, Function, Arugment1, Arguments2, ...} where Number is an arbitrary task 
%          number, Function is the function to evaluate, and ArgumentX are arguments to the function.
%
% Out:
%   Result : result record; A sequence of bytes which represents a cell array of the form:
%            {Number, ResultData} where Number is the corresponding task number and ResultData is 
%            the first output of the function.
%
%            If an exception happened during processing, result represents {Number, LastError} where
%            LastError is an error record as returned by lasterror.
%
% Notes:
%   Depending on the severity of the error (especially when the serialization itself fails), the 
%   error record may contain an inaccurate number or just a string as error
%   record.
%
% See also:
%   par_worker
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-08-29

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

try
    % deserialize the task (names mangled as the task may involve an "eval" call).
    task__ = hlp_deserialize(task__);
    % evaluate Function(Argument1,Argument2, ...) and construct the result record
    result = {task__{1},task__{2}(task__{3:end})};
    % serialize the result record
    result = hlp_serialize(result);
catch e
    try
        % try to display and serialize the error
        fprintf('exception during task processing: %s\n',e.message);
        if exist('hlp_handleerror','file')
            hlp_handleerror(e);
        else
            for k = 1:length(e.stack)
                fprintf('  %s (%i)\n',e.stack(k).name,e.stack(k).line); end
        end
        if iscell(task__)
            result = hlp_serialize({task__{1},e});
        else
            result = hlp_serialize({NaN,e});
        end
    catch
        % fall back to minimal reporting
        try
            fprintf('exception during error serialization: %s\n',e.message);
            result = hlp_serialize({NaN,['Error during error serialization: ' e.message]});
        catch
            fprintf('serialization failed completely.\n');
            result = uint8(200);
        end
    end
end
