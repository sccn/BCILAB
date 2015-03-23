function par_accept_results(tag, payload)
% Accept a result from the scheduler and store it in a result table.
% par_accept_result(Tag, Payload)
%
% This function stores the results produced by the scheduler in a global result table for later
% pickup by par_endschedule. The reason is that the scheduler, being a Java class, cannot hold very
% large (or very many) results since the MATLAB Java memory is limited to a maximum of 64GB.
% Therefore we need to ferry the results out of Java as they are being produced.
%
% In:
%   Tag: tag string that uniquely identifies the result; should begin with 'tag__'.
%
%   Payload: payload string of the same format as returned by the Scheduler

global tracking;
% atomically store payload in global results table
tracking.parallel.results.(tag) = payload;