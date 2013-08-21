function struct12 = structcat(dim, struct1, struct2)
%STRUCTCAT         Concatenates two structures field by field.
%   STRUCT12 = STRUCTCAT(DIM, STRUCT1, STRUCT2) takes two structures with
%   with identical fields and returns a new structure in which each field
%   contains the concatenation of the corresponding fields in STRUCT1 and
%   STRUCT2 along dimension DIM.

% Argument checking.
fields = fieldnames(struct1);      fields2 = fieldnames(struct2);
if ((length(fields) ~= length(fields2)) || (~all(strcmp(fields,fields2))))
	error('Field names do not match.');
end

for f = 1:length(fields)
	struct12.(fields{f}) = cat(dim, struct1.(fields{f}), struct2.(fields{f}));
end