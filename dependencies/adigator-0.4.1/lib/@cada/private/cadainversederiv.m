function nzlocs = cadainversederiv(x,ytemp,Vcount,derivstr,DPFLAG)
% This function does the Derivative Calculations for Matrix Inverse, thus,
% mldivide, mrdivide, and inv call this function for matrix derivatives.
% Letting y = inv(x), the derivative of the inverse is calculated using the
% fact that I = y*x, thus dI = dy*x + y*x = 0 => dy = -y*dx*y
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if DPFLAG
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);
end

N = x.func.size(1);

dxind = x.deriv(Vcount).nzlocs;
nv  = ADIGATOR.VAROFDIFF(Vcount).usize;
nzx = size(dxind,1);
% -------------------- X has Derivatives ---------------------------- %
% we need to compute dy = -y*dx*y,
% let dA = y*dx, we can compute dA by
% dA = y*reshape(dx,N,N*nv) then reshape to size
% [N*N,nv]
dx = sparse(dxind(:,1),dxind(:,2),1:nzx,N*N,nv);
dxReshape = reshape(dx,N,N*nv);
dA = reshape(ytemp*dxReshape,N*N,nv);

% We now have dA, let dy = dA*y, we can compute the transpose of dy
% by dy' = y'*dA'
% First get dA transpose
ATranMap = zeros(N,N); ATranMap(:) = 1:N*N;
ATranMap = ATranMap.';
dATranReshape = reshape(dA(ATranMap(:),:),N,N*nv);
% Now can get dy'
dyTran = reshape(ytemp.'*dATranReshape,N*N,nv);
% Now get dy from dy'
yTranMap = ATranMap;
dy = dyTran(yTranMap(:),:);
% we now have -dy really, but has same non-zero locations
[yrows,ycols] = find(dy);
if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end

if ~isempty(yrows)
  nzlocs = [yrows ycols];
  if DPFLAG
    % --------------------- Derivative Printing --------------------- %
    % We need to print out all of the calculations that we just did
    % above
    % 1. Print out Calculations to get to dA
    % 1a. Project x derivative variable into a matrix corresponding to
    %     dxReshape
    TD1 = ['cada',NDstr,'td1']; % dx' in file
    if N*N*nv > 250 && nzx < .6*N*N*nv
      SPflag = 1;
      % Do projection sparsely - dx is larger than 250 and has at least
      % 40% zeros
      dxrInds = zeros(nzx,2);
      [dxrInds(:,1),dxrInds(:,2)] = find(reshape(dxReshape,N,N*nv));
      TDind1 = cadaindprint(dxrInds(:,1));
      TDind2 = cadaindprint(dxrInds(:,2));
      % Print out the Projection
      fprintf(fid,[indent,TD1,' = sparse(',TDind1,',',TDind2,',',...
        x.deriv(Vcount).name,',%1.0d,%1.0d);\n'],N,N*nv);
    else
      SPflag = 0;
      % Project into zeros using a linear index
      dxIndLin = sub2ind([N*N,nv],dxind(:,1),dxind(:,2));
      % Note that dxreshape and dx share the same LINEAR index.
      TDind1 = cadaindprint(dxIndLin(:));
      % Print out the projection
      fprintf(fid,[indent,TD1,' = zeros(%1.0d,%1.0d);\n'],N,N*nv);
      fprintf(fid,[indent,TD1,'(',TDind1,') = ',...
        x.deriv(Vcount).name,';\n']);
    end
    % 1b. We can now print out calculations for dAReshape =
    % y*dxReshape, or equivalently, x\dxReshape, reshape to dA here as
    % well
    if SPflag && nnz(dA) > 0.6*N*N*nv
      fprintf(fid,[indent,TD1,' = reshape(full(',x.func.name,'\\',TD1,...
        '),%1.0d,%1.0d);\n'],N*N,nv);
      SPflag = 0;
    else
      fprintf(fid,[indent,TD1,' = reshape(',x.func.name,'\\',TD1,...
        ',%1.0d,%1.0d);\n'],N*N,nv);
    end
    % TD1 is now dA
    % 2. Print out Calculations to get dy
    % 2a. print out dATranReshape
    TDind1 = cadaindprint(ATranMap(:));
    fprintf(fid,[indent,TD1,' = reshape(',TD1,'(',TDind1,',:),',...
      '%1.0d,%1.0d);\n'],N,N*nv);
    % TD1 is now dATranReshape
    % 2b. Print out calculations to get dyTranReshape =
    % y'*dATranReshape = x'\dATransposeReshape
    fprintf(fid,[indent,TD1,' = ',x.func.name,'.''\\',TD1,';\n']);
    % 2c. Need to get a linear reference which will reference off of
    % dyTranReshape into the non-zeros of dy.
    [yTrows,yTcols] = find(dyTran);
    dyTIndLin = sub2ind([N*N,nv],yTrows,yTcols);
    dyTran = sparse(yTrows,yTcols,dyTIndLin,N*N,nv);
    dy = dyTran(yTranMap(:),:);
    % Note that dyTranReshape and dyTran share the same linear index.
    dyIndLin = nonzeros(dy);
    TDind1 = cadaindprint(dyIndLin(:));
    % Can now print out dy
    if SPflag
      fprintf(fid,[indent,derivstr,' = -full(',TD1,'(',TDind1,'));\n']);
    else
      fprintf(fid,[indent,derivstr,' = -',TD1,'(',TDind1,');\n']);
    end
  end
else
  nzlocs = [];
end