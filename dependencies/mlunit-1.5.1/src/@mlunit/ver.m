function ver(self) %#ok
%mlunit/ver prints the version string of mlUnit to the standard output.
%
%  EXAMPLE
%  =======
%         ver(mlunit);
%
%  See also MLUNIT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: ver.m 226 2007-01-21 15:20:53Z thomi $

fprintf(1, 'mlUnit Version 1.5\n');