function data=extractdatac(data,Fs,t)
% Extract segements of continuous data between t(1) and t(2)
% Usage: data=extractdatac(data,Fs,t)
%
% Input:
% data: continous data in the form samples x channels or a single vector
% Fs: sampling frequency
% t   : time as a 2d vector [t(1) t(2)]
% Note that sampling frequency and t have to be in consistent units
% Output:
% data: data between t(1) and t(2)

if nargin < 3; 
    error('need all three arguments');
end;
if t(1) < 0 || t(2)<=t(1);
    error('times cannot be negative and t(2) has to greater than t(1)');
end;
data=change_row_to_column(data);
N=size(data,1);
tt=(0:N-1)/Fs;
indx=find(tt>=t(1) & tt<t(2));
data=data(indx,:);