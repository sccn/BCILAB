function result = quickif(condition,ontrue,onfalse)
% Select one of two possible values dependent on a condition
% Result = quickif(Condition,OnTrue,OnFalse)
%
% In:
%   Condition : conditional expression
%
%   OnTrue : value to return if condition evaluates to true
%
%   OnFalse : value to return if condition evaluates to false
%
% Out:
%   Result : either OnTrue or OnFalse
%
% Notes:
%   Note that both branches will always be evaluated, regardless of the value of Condition, since
%   quickif() is a regular function call. This is purely a convenience function that allows to use
%   conditional expressions inside other expressions.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-10-16

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
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

if condition
    result = ontrue;
else
    result = onfalse;
end
