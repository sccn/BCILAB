function option = setopt(opt_par,default);
%    option = setopt(opt_par,default); 
%
%    SETOPT can modify a default option vector with user specified options.
%
%    Input: 
%       opt_par : Users desired option vector
%       default : Program default option vector
%
%    Output: 
%       option : New option vector
%
%    Example:
%       opt_par = [1 2 3 4];
%       default = [1 1 1 1];
%       option = setopt(opt_par,default)
%       option = 1     2     3     4
%
%Author: Jan Erik Odegard  <odegard@ece.rice.edu>

if (nargin < 2) 
  error('You need two inputs');
end;
len = length(opt_par);
option = zeros(size(default));
option(1:len) = opt_par(1:len);
option = option + (option == 0).*default; 
