% Resampling matrix
%
% Signal resampling in 2d and 3d.
%
% A = matResample(sz,x,g,method,deriv)
%  sz     size of the image
%  x      new coordinates
%  method is the interpolation method nearest, linear, cubic   [default nearest]
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 February 08

function W = matResample(sz,x,method,deriv)

  ndims = numel(sz);
  if ndims==2
    [i2,i1]=meshgrid(1:sz(2),1:sz(1)); i1=i1(:)'; i2=i2(:)';       % target grid
  else
    [i2,i1,i3] = meshgrid(1:sz(2),1:sz(1),1:sz(3));
    i1=i1(:)'; i2=i2(:)'; i3=i3(:)';
  end

  if nargin<3, method = 'nearest'; end                      % set default values
  if nargin<4, deriv  = 0; end
  if strfind(lower(method),'lin')                         % set window functions
    nn =  0:1; k = @(x) klin(x); dk = @(x) dklin(x);
  elseif strfind(lower(method),'cub')
    nn = -1:2; k = @(x) kcub(x); dk = @(x) dkcub(x);
  else
    nn =  0:1; k = @(x) knea(x); dk = @(x) dknea(x);
  end
  N = prod(sz); M = size(x,2);
  W = sparse(M,N);                                 % interpolation weight matrix
  for n1 = nn,   j1 = floor(x(1,:))+n1;       % linear interpolation from source
    for n2 = nn, j2 = floor(x(2,:))+n2;      % .. neighbours; iterate over those
      if ndims==2
        ok = 0<j1 & j1<=sz(1)  &  0<j2 & j2<=sz(2);     %   indices inside image
        i = i1(ok) + (i2(ok)-1)*sz(1);
        j = j1(ok) + (j2(ok)-1)*sz(1);
        if deriv==1
          w = dk(x(1,ok)-j1(ok));
        else
          w =  k(x(1,ok)-j1(ok));
        end
        if deriv==2
          w = w.*dk(x(2,ok)-j2(ok));
        else
          w = w.* k(x(2,ok)-j2(ok));
        end
        W = W + sparse(i,j,w,M,N);                       % add upp contributions
      else
        for n3=0:1, j3 = floor(x(3,:))+n3;
          ok = 0<j1 & j1<=sz(1) & 0<j2 & j2<=sz(2) ...
                                & 0<j3 & j3<=sz(3);     %   indices inside image
          i = i1(ok) + (i2(ok)-1)*sz(1) + (i3(ok)-1)*prod(sz(1:2));
          j = j1(ok) + (j2(ok)-1)*sz(1) + (j3(ok)-1)*prod(sz(1:2));
          
          if deriv==1
            w = dk(x(1,ok)-j1(ok));
          else
            w =  k(x(1,ok)-j1(ok));
          end
          if deriv==2
            w = w.*dk(x(2,ok)-j2(ok));
          else
            w = w.* k(x(2,ok)-j2(ok));
          end
          if deriv==3
            w = w.*dk(x(3,ok)-j3(ok));
          else
            w = w.* k(x(3,ok)-j3(ok));
          end
          W = W + sparse(i,j,w,M,N);                     % add upp contributions
        end
      end
    end
  end


% Robert G. Keys, Cubic Convolution Interpolation for Digital Image Processing,
% IEEE ASSP, 29:6, December 1981, p. 1153-1160.
function y = kcub(x)
  y = zeros(size(x)); x = abs(x);
  q = (x<=1);          % Coefficients:  1.5, -2.5,  0.0, 1.0
  y(q) =            (( 1.5 * x(q) - 2.5) .* x(q)      ) .* x(q) + 1.0;
  q = (1<x & x<=2);    % Coefficients: -0.5,  2.5, -4.0, 2.0
  y(q) =            ((-0.5 * x(q) + 2.5) .* x(q) - 4.0) .* x(q) + 2.0;

function y = dkcub(x)
  y = sign(x); x = abs(x);
  q = (x<=1);          % Coefficients:  1.5, -2.5,  0.0, 1.0
  y(q) = y(q) .*  ( 4.5 * x(q) - 5.0) .* x(q);
  q = (1<x & x<=2);    % Coefficients: -0.5,  2.5, -4.0, 2.0
  y(q) = y(q) .* ((-1.5 * x(q) + 5.0) .* x(q) - 4.0);
  y(x>2) = 0;

function y =  klin(x)
y = max(1-abs(x),0);
  

function y = dklin(x)
  y = -sign(x);
  y(abs(x)>1) = 0;

function y =  knea(x)
  y = zeros(size(x));
  y(abs(x)<=1/2) = 1;

function y = dknea(x)
  y = zeros(size(x));