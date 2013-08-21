function quicklog(logfile,msg,varargin)
% log something to the screen & to a logfile

% Copyright (C) Christian Kothe, SCCN, 2011, christian@sccn.ucsd.edu
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
    msg = sprintf([hlp_hostname '/' datestr(now) '@' hlp_getcaller ': ' msg '\n'],varargin{:});
    fprintf(msg);
    if ~isempty(logfile)
        try
            % try to create / open the log file
            if ~exist(logfile,'file')
                try
                    % re-create
                    fid = fopen(logfile,'w+');
                    if fid == -1
                        error('Error creating file'); end
                catch
                    % failed: try to create directories
                    try
                        io_mkdirs(logfile,{'+w','a'});
                    catch
                        disp(['Could not create logfile directory: ' fileparts(logfile)]);
                    end
                    % try to create file
                    try
                        fid = fopen(logfile,'w+');
                        if fid == -1
                            error('Error creating file'); end
                    catch
                        disp(['Could not create logfile ' logfile]);
                    end
                end
            else
                % append
                fid = fopen(logfile,'a');
                if fid == -1
                    error('Error creating file'); end
            end
        catch
            disp(['Could not open logfile ' logfile]);
        end
        % write message
        try
            fprintf(fid,msg);
        catch
            disp(['Could not write to logfile ' logfile]);
        end
        % close file again
        try
            fclose(fid);
        catch
        end
    end
catch
    disp('Invalid logging parameters.');
end
