function self = run(self, dock)
%mlunit/run executes the the graphical user interface of mlUnit.
%
%  EXAMPLE
%  =======
%  Run in window mode:
%         run(mlunit);
%  Run in docked mode:
%         run(mlunit, 1);
%
%  See also MLUNIT, GUI_TEST_RUNNER/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: run.m 160 2007-01-03 21:56:21Z thomi $

if (nargin == 1)
    dock = 0;
end;

run(gui_test_runner, '', dock);