function nzlocs = cadamtimesderivvec(x,y,xtemp,ytemp,Vcount,derivstr,DPFLAG,zvec)
% This function does matrix multiplication derivative calculations when in
% the vectorized mode.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if DPFLAG
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);
end
[xMrow, xNcol] = size(xtemp);
[yMrow, yNcol] = size(ytemp);
% Let y be non-vectorized, Y be vectorized
nv = ADIGATOR.VAROFDIFF(Vcount).usize;
if zvec == 1
  % Z = X*y, where x is N by n (N vectorized dim)
  dxind = x.deriv(Vcount).nzlocs;
  if ~isempty(dxind)
    % cancel any derivs if possible
    xrows = dxind(:,1); xcols = dxind(:,2);
    nzx = length(xrows);
    nx = xMrow*xNcol;
    dx = sparse(xrows,xcols,1:nzx,nx,nv);
    dz = ytemp.'*dx;
    [zrows, zcols] = find(dz);
    if size(zrows,2) > 1; zrows = zrows.'; zcols = zcols.'; end
    nzlocs = [zrows, zcols];
    if DPFLAG
      % To print this out we need to do dZ = dX*A
      % dX: N by nzx
      % dZ: N by nzz
      % A : nzx by nzz
      nzz = length(zrows);
      % Need to build A by referencing off of y
      iAasgn = zeros(nzx,nzz);
      jAasgn = repmat(1:nzz,nzx,1);
      iYref  = zeros(nzx,nzz);
      jYref  = zeros(nzx,nzz);
      for Dcount = 1:nzz
        iasgn = nonzeros(dx(:,zcols(Dcount)));
        iAasgn(iasgn,Dcount) = iasgn;
      end

      jAasgn = jAasgn.*logical(iAasgn);
      iYref(logical(iAasgn)) = xrows(nonzeros(iAasgn));
      jYref(logical(iAasgn)) = zrows(nonzeros(jAasgn));
      
      RefInd  = cadaindprint(sub2ind([yMrow,yNcol],nonzeros(iYref),nonzeros(jYref)));
      TF1 = ['cada',NDstr,'tf1'];
      if nzx*nzz > 250 && nnz(iAasgn) < .6*nzx*nzz
        AsgnInd1 = cadaindprint(nonzeros(iAasgn));
        AsgnInd2 = cadaindprint(nonzeros(jAasgn));
        fprintf(fid,[indent,TF1,' = sparse(',AsgnInd1,',',AsgnInd2,',',...
          y.func.name,'(',RefInd,'),%1.0f,%1.0f);\n'],nzx,nzz);
        fprintf(fid,[indent,derivstr,' = full(',x.deriv(Vcount).name,'*',TF1,');\n']);
      else
        AsgnInd = cadaindprint(sub2ind([nzx,nzz],nonzeros(iAasgn),nonzeros(jAasgn)));
        fprintf(fid,[indent,TF1,' = zeros(%1.0d,%1.0d);\n'],nzx,nzz);
        fprintf(fid,[indent,TF1,'(',AsgnInd,') = ',y.func.name,'(',RefInd,');\n']);
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'*',TF1,';\n']);
      end
    end
  else
    nzlocs = [];
  end
else
  % Z = x*Y, where Y is n by N (N vectorized dim)
  dyind = y.deriv(Vcount).nzlocs;
  if ~isempty(dyind)
    % cancel any derivs if possible
    yrows = dyind(:,1); ycols = dyind(:,2);
    nzy = length(yrows);
    ny = yMrow*yNcol;
    dy = sparse(yrows,ycols,1:nzy,ny,nv);
    dz = xtemp*dy;
    [zrows, zcols] = find(dz);
    if size(zrows,2) > 1; zrows = zrows.'; zcols = zcols.'; end
    nzlocs = [zrows, zcols];
    if DPFLAG
      % To print this out we need to do dZ = dY*A 
      % dY: N by nzy
      % dZ: N by nzz
      % A : nzy by nzz
      nzz = length(zrows);
      % Need to build A by referencing off of y
      iAasgn = zeros(nzy,nzz);
      jAasgn = repmat(1:nzz,nzy,1);
      iXref  = zeros(nzy,nzz);
      jXref  = zeros(nzy,nzz);
      for Dcount = 1:nzz
        iasgn = nonzeros(dy(:,zcols(Dcount)));
        iAasgn(iasgn,Dcount) = iasgn;
      end

      jAasgn = jAasgn.*logical(iAasgn);
      iXref(logical(iAasgn)) = zrows(nonzeros(jAasgn));
      jXref(logical(iAasgn)) = yrows(nonzeros(iAasgn));
      
      RefInd  = cadaindprint(sub2ind([xMrow,xNcol],nonzeros(iXref),nonzeros(jXref)));
      TF1 = ['cada',NDstr,'tf1'];
      if nzy*nzz > 250 && nnz(iAasgn) < .6*nzy*nzz
        AsgnInd1 = cadaindprint(nonzeros(iAasgn));
        AsgnInd2 = cadaindprint(nonzeros(jAasgn));
        fprintf(fid,[indent,TF1,' = sparse(',AsgnInd1,',',AsgnInd2,',',...
          x.func.name,'(',RefInd,'),%1.0f,%1.0f);\n'],nzy,nzz);
        fprintf(fid,[indent,derivstr,' = full(',y.deriv(Vcount).name,'*',TF1,');\n']);
      else
        AsgnInd = cadaindprint(sub2ind([nzy,nzz],nonzeros(iAasgn),nonzeros(jAasgn)));
        fprintf(fid,[indent,TF1,' = zeros(%1.0d,%1.0d);\n'],nzy,nzz);
        fprintf(fid,[indent,TF1,'(',AsgnInd,') = ',x.func.name,'(',RefInd,');\n']);
        fprintf(fid,[indent,derivstr,' = ',y.deriv(Vcount).name,'*',TF1,';\n']);
      end
    end
  else
    nzlocs = [];
  end
end