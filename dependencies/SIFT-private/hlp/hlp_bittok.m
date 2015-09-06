function idx = hlp_bittok(x,bit)
% segment logical vector x and return start-stop indices of all sequences 
% of value=bit. 'bit' can be 0 or 1
% 
% Input: 
%   X:      a binary vector
%   bit:    0 or 1 (the bit-sequence to scan for)
% Output:   
%   idx:    matrix where each row contains [start end] indices of the
%           sequence of bits
%
% Example:
%
% x = [0 1 1 0 0 1 0 0 1 1 1 1];
% bit = 1;
% hlp_bittok(x,bit)
% ans =
% 
%      2     3
%      6     6
%      9    12
%
% Author: Tim Mullen (C) 2011, SCCN/INC/UCSD

if ~(bit==0 || bit==1)
    error('bit must be either 0 or 1');
end

if ~all(x==1 | x==0)
    error('x must be a logical (binary) sequence');
end
   
x = x(:)';

x = num2str(x);
x(isspace(x))=[];

[startidx endidx] = regexp(x,sprintf('[%d]*%d',bit),'start','end');

if isempty(startidx)
    idx = [];
else
    idx = [startidx; endidx]';
end
