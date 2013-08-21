function self = mlunit
%mlunit contructor.
%  The constructor creates an object of the class mlunit.
%
%  Class Info / Example
%  ====================
%  The class mlunit is a shortcut to run the graphical user interface of
%  mlUnit.
%         Example: run(mlunit);
%
%  See also GUI_TEST_RUNNER.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: mlunit.m 160 2007-01-03 21:56:21Z thomi $

self = class(struct([]), 'mlunit');