function [N,C]=check_consistency(data1,data2,sp)
% Helper routine to check consistency of data dimensions
% Usage: [N,C]=check_consistency(data1,data2,sp)
% Inputs:
% data1 - first dataset
% data2 - second dataset
% sp - optional argument to be input as 1 when one of the two data sets is
% spikes times stored as a 1d array.
% Outputs:
% Dimensions of the datasets - data1 or data2 (note that 
%    routine stops with an error message if dimensions don't match - [N,C]
%    N left empty for structure arrays
N1=[]; N2=[];
if nargin < 3 || isempty(sp); sp=0; end;
if isstruct(data1);
    C1=length(data1);
else
    [N1,C1]=size(data1);
end;
if isstruct(data2);
    C2=length(data2);
else
    [N2,C2]=size(data2);
end;
if C1~=C2; error('inconsistent dimensions'); end;
if sp==0;
   if ~isstruct(data1) && ~isstruct(data2);
      if N1~=N2; error('inconsistent dimensions'); end;
   end;
end;
N=N1; C=C1;