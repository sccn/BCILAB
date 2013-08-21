function self = install(self)
%mlunit/install installs mlUnit. The current version of this method only
%renames the assert method of mlUnit, if a built-in assert already exists
%(which is the case for MATLAB 6.6).
%
%  EXAMPLE
%  =======
%         install(mlunit);

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: install.m 270 2007-04-02 21:02:33Z thomi $

which_assert = which('assert', '-all');
assert_built_in = 0;
for i = 1:length(which_assert)
    if (strfind(which_assert{i}, 'built-in'))
        assert_built_in = 1;
        break;
    end;
end;
if (assert_built_in)
    fprintf(1, 'mlUnit has found a built-in assert functions.\n');
    fprintf(1, 'Trying to rename mlUnit assert.m...');
    try
        if (exist('assert', 'file'))
            movefile('assert.m', 'mlunit_assert.m', 'f');
            fprintf(1, 'ok.\n');
        else
            fprintf(1, 'not found.\n');
        end;
    catch
        fprintf(1, 'error.\n');
    end;
end;
