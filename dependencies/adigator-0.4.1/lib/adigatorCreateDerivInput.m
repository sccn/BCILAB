function x = adigatorCreateDerivInput(xsize,derivinfo)
% function x = adigatorCreateDerivInput(xsize,derivinfo)
% 
% This function creates an input,x, which has derivatives that is to be 
% used with the function source-to-derivative source transformation
% function, adigator.
%
% Input Information:
% xsize = [mx my]
%       where mx and my are the row and column dimensions, respectively, of
%       the input x.    
%
% derivinfo = string name of variable of differentiation, assuming the
%           input is the variable of differentiation
%
%                           OR
%
% derivinfo = N by 1 structure array containing information pertaining to 
%           the derivatives of the input x, where N is the number of 
%           variables which x contains non-zero derivatives with respect to
%
% Each element of the derivinfo structure array must have the following
% fields:
%
%     .vodname = 'derivstr'
%              where 'derivstr' is a string which uniquely defines the
%              variable of differentiation.
%
%     .vodsize = [mvod nvod]
%              where mvod and nvod are the row and column dimensions of
%              the variable of differentiation.
%
%     .nzlocs = Nnz by 2 integer array of the locations of the possible
%             nonzeros in the unrolled Jacobian of x with respect to the
%             variable of differentiation.
% 
% If derivinfo is a string name, then it will be assumed that the variable
% of differentiation is the input, x. Thus this is the same as defining the
% string name in .vodname together with .vodsize = xsize, .nzlocs =
% [(1:mx*my).' (1:mx*my).'].
%
% Simple Example:
% If we are to differentiate a function, y = myfun(x,a), with respect to the
% input x, where x is a vector of length 3, and a is some auxillory input,
% we would do so as follows: (Assuming input a is previously defined)
%
% % Define the input size
% xsize = [3 1];
%
% % Define the variable of differentiation (in this case it is our input)
% derivinfo.vodname = 'x';    % Call var of diff x
% derivinfo.vodsize = xsize;  % Var of diff and input are equal
% derivinfo.nzlocs  = [1 1;...
%                      2 2;...
%                      3 3];  % Define non-zero locations (along the
%                             % diagonal of the Jacobian of x wrt x)
%
% % Create the input variable for x
% x = adigatorCreateDerivInput(xsize,derivinfo);
%
% % Transform myfun file into myderiv file
% adigator('myfun','myderiv',x,a);
%
% % Can now evaluate the derivative file (Note: Generated file input is
% % different than the original file)
% % Create x input
% xn.f  = rand(5,1);  % Evaluate at random function values
% xn.dx = ones(5,1);  % Possible non-zero elements of dx/dx
% yn = myfun(xn,a);
% % Build Jacobian (Assuming output is a vector)
% Jy_x = ...
% sparse(yn.dx_location(:,1),yn.dx_location(:,2),y.dx,y.dx_size(1),y.dx_size(2));
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% See also adigator adigatorCreateAuxInput adigatorOptions
if isequal(size(xsize),[1 2]) && isequal(xsize,floor(xsize))
  func.size = xsize;
else
  error('first input xsize must be integer array of size 1 by 2')
end
if isinf(xsize(1)) && isinf(xsize(2))
  error('only one dimension of the input may be vectorized');
elseif isinf(xsize(1))
  xvectflag = 2;
elseif isinf(xsize(2))
  xvectflag = 1;
else
  xvectflag = 0;
end

if ischar(derivinfo)
  vodname = strtrim(derivinfo);
  deriv = struct('vodname',[],'vodsize',[],'nzlocs',[]);
  deriv.vodname = vodname;
  deriv.vodsize = xsize;
  if xvectflag
    deriv.nzlocs  = [(1:xsize(xvectflag)).' (1:xsize(xvectflag)).'];
  else
    deriv.nzlocs  = [(1:xsize(1)*xsize(2)).' (1:xsize(1)*xsize(2)).'];
  end
elseif isstruct(derivinfo) && isfield(derivinfo,'vodname') && isfield(derivinfo,'vodsize') && isfield(derivinfo,'nzlocs')
  Nvod = length(derivinfo);
  deriv = struct('vodname',cell(Nvod,1),'vodsize',cell(Nvod,1),'nzlocs',cell(Nvod,1));
  for Vcount = 1:Nvod
    if ischar(derivinfo(Vcount).vodname)
      deriv(Vcount).vodname = strtrim(derivinfo(Vcount).vodname);
    else
      error('error in deriv %1.0f: vodname field must be a string name',Vcount)
    end
    if isequal(size(derivinfo(Vcount).vodsize),[1 2]) && ...
        isequal(derivinfo(Vcount).vodsize,floor(derivinfo(Vcount).vodsize))
      vodsize = derivinfo(Vcount).vodsize;
      deriv(Vcount).vodsize = vodsize;
      if xvectflag && (all(isinf(vodsize)) || ~any(isinf(vodsize)))
        error(['error in deriv %1.0f: if input is vectorized, then one ',...
          'dimension of the vod must be vectorized'],Vcount)
      elseif isinf(vodsize(1))
        dvectflag = 2;
      elseif isinf(vodsize(2))
        dvectflag = 1;
      else
        dvectflag = 0;
      end
    else
      error(['error in deriv %1.0f: vodsize field must be integer array ',...
        'of size 1 by 2'],Vcount)
    end
    if size(derivinfo(Vcount).nzlocs,2) == 2 && ~isempty(derivinfo(Vcount).nzlocs)...
        && isequal(derivinfo(Vcount).nzlocs,floor(derivinfo(Vcount).nzlocs))
      deriv(Vcount).nzlocs = derivinfo(Vcount).nzlocs;
      if (~xvectflag && any(deriv(Vcount).nzlocs(:,1) > prod(xsize) |...
          deriv(Vcount).nzlocs(:,2) > prod(deriv(Vcount).vodsize))) ||...
          (xvectflag && any(deriv(Vcount).nzlocs(:,1) > xsize(xvectflag) |...
          deriv(Vcount).nzlocs(:,2) > deriv(Vcount).vodsize(dvectflag)))
        error(['error in deriv %1.0f: nzlocs outside of defined size of ',...
          'unrolled Jacobian'],Vcount)
      end
      if xvectflag
        testjac = sparse(deriv(Vcount).nzlocs(:,1),deriv(Vcount).nzlocs(:,2),...
          1:size(deriv(Vcount).nzlocs,1),xsize(xvectflag),deriv(Vcount).vodsize(dvectflag));
      else
        testjac = sparse(deriv(Vcount).nzlocs(:,1),deriv(Vcount).nzlocs(:,2),...
          1:size(deriv(Vcount).nzlocs,1),prod(xsize),prod(deriv(Vcount).vodsize));
      end
      if ~isequal(nonzeros(testjac),(1:size(deriv(Vcount).nzlocs,1)).')
        warning(['deriv %1.0f nzlocs defined out of order - deriv inputs  ',...
          'to generated deriv function must be in proper order'],Vcount)
        deriv(Vcount).nzlocs = deriv(Vcount).nzlocs(nonzeros(testjac),:);
      end
    else
      error(['error in deriv %1.0f: vodsize field must be integer array ',...
        'of size Nnz by 2'],Vcount)
    end
  end
else
  error(['second input derivinfo must be character string or structure ',...
    'array with fields: vodname, vodsize, nzlocs'])
end

x = cada(0,func,deriv);