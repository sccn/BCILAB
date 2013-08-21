function s = structpack(names)
%STRUCTPACK        Copies workspace variables into the fields of a structure.
%   S = STRUCTPACK(NAMES) takes a cell array NAMES of strings containing
%   the names of variables in the current workspace.  It returns a
%   structure with field names and contents copied from those variables.
%   For example,
%                    X = 1; Y = 2; S = STRUCTPACK({'X','Y'});
%   will return a structure S such that (S.X = 1) and (S.Y = 2).
%
%   See also STRUCTUNPACK.

for n = 1:length(names)
	nameexist = evalin('caller', ['exist(''' names{n} ''');']);
	if (nameexist ~= 1), 
		error(['The variable ' names{n} ' was not found in the workspace.']); 
	end;
	
	s.(names{n}) = evalin('caller', names{n});
end
