function data=change_row_to_column(data)
% Helper routine to transform 1d arrays into column vectors that are needed
% by other routines in Chronux
%
% Usage: data=change_row_to_column(data)
% 
% Inputs:
% data -- required. If data is a matrix, it is assumed that it is of the
% form samples x channels/trials and it is returned without change. If it
% is a vector, it is transformed to a column vector. If it is a struct
% array of dimension 1, it is again returned as a column vector. If it is a
% struct array with multiple dimensions, it is returned without change
% Note that the routine only looks at the first field of a struct array.
% 
% Ouputs:
% data (in the form samples x channels/trials)
%

dtmp=[];
if isstruct(data);
   C=length(data);
   if C==1;
      fnames=fieldnames(data);
      eval(['dtmp=data.' fnames{1} ';'])
      data=dtmp(:);
   end
else
  [N,C]=size(data);
  if N==1 || C==1;
    data=data(:);
  end;
end;
