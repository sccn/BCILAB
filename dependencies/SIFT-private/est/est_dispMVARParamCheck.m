function [checkcode checkstr] = est_dispMVARParamCheck(ALLEEG,params,display)
% display results of sanity checks on MVAR parameters
%
% Input:
%
%   ALLEEG             EEG data structure.
%   params             struct of MVAR parameters as returned from
%                      est_fitMVAR()
%
%
% Output:
%
%   checkcode:         'ok'      - all checks passed
%                      'error'   - parameters result in critical error
%                      'warning' - possible problem with chosen params
%
% See Also: est_checkMVARParams(), est_fitMVAR()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen, 2011, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if nargin<3
    display = true;
end

checkcode = 'ok';
checkstr = '';

% check that parameters are OK
try
    [infostring warnstring errstring] = est_checkMVARParams(ALLEEG,params);
catch err
    %         errordlg2(err.message,'Error in MVAR configuration!');
    checkcode = err;
    return;
end

if params.verb>0
    % display summary on command line
    for str=1:length(infostring)
        if ~isempty(errstring{str})
            checkstr = [checkstr sprintf('ERROR: \t')];
        elseif ~isempty(warnstring{str})
            checkstr = [checkstr sprintf('WARNING: \t')];
        else
            checkstr = [checkstr sprintf('OK: \t')];
        end
        checkstr = [checkstr sprintf(infostring{str})];
        checkstr = [checkstr sprintf(warnstring{str})];
        checkstr = [checkstr sprintf(errstring{str})];
        checkstr = [checkstr sprintf('\n')];
    end
end

if display
    fprintf(checkstr);
end

if ~all(ismember_bc(errstring,''))
    checkcode = 'error';
elseif ~all(ismember_bc(warnstring,''))
    checkcode = 'warning';
end