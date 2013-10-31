function tf = par_haveupdate(current_file,reference_file)
% Return true if a code update is available.
% Result = par_haveupdate(CurrentFile,ReferenceFile);
%
% This is only for expert use (when workers need to update themselves after a code change on the
% master).
%
% In:
%   CurrentFile : name of the file that is currently executing
%
%   ReferenceFile : name of a replacement file that is possibly newer than the CurrentFile
%
% Out:
%   Result : whether a newer ReferenceFile is available
%
% See also:
%   par_worker
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-08-26

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

tf = false;

try
    current_file = env_translatepath(current_file);
    reference_file = env_translatepath(reference_file);
catch
end

if ~exist(current_file,'file')
    fprintf('The currently executing code (%s) is non-existent; cannot check for updates.\n',current_file);
elseif exist(reference_file,'file')
        % reference file present: could potentially update: compare file dates
        ref_info = dir(reference_file);
        cur_info = dir(current_file);
        if isempty(ref_info)
            fprintf('No file info for ReferenceFile (%s) available. Cannot check for updates.\n',current_file);
            return;
        end
        if isempty(ref_info)
            fprintf('No file info for CurrentFile (%s) available. Cannot check for updates.\n',current_file);
            return;
        end
        % compare file dates
        tf = ref_info.datenum > cur_info.datenum;
    end
end
