function structure = structindex(structure, index)
%STRUCTINDEX       Indexes each field in a structure of arrays.
%   STRUCTOUT = STRUCTINDEX(STRUCTIN, INDEX) takes a structure STRUCTIN
%   in which every field is a cell or numeric array and returns a
%   structure with the same fields with the data rearranged according to:
%                STRUCTOUT.FIELD1 = STRUCTIN.FIELD1(INDEX).
fields = fieldnames(structure); 
for f = 1:length(fields)
	structure.(fields{f}) = structure.(fields{f})(index);
end
