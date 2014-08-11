function nzlocs = cadamtimesderiv(x,y,xtemp,ytemp,Vcount,derivstr,DPFLAG,caller)
% This function does the Derivative calculations for Matrix Multiplication,
% where it may be called from mtimes, mldivide, or mrdivide. We note that
% the multiplication a*db for matrices is relatively simple given the
% unrolled form in which we store derivatives, thus we use the following
% derivative rules. Allowing that dx' is the derivative of the transpose of
% x, x'
%
% If called from mtimes: 
%   z  = x*y
%   dz = dx*y + x*dy = (y'*dx')' + x*dy
%
% If called from mldivide: 
%   z  = x\y and this will only be called with y having derivatives, thus
%   dz = x\dy
%   see mldivide for more details
%
% If called from mrdivide:
%   z  = x/y and this will only be called with x having derivatives, thus
%   dz = dx/y = (y'\dx')'
%   see mrdivide for more details
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR
if DPFLAG
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);
end

xMrow = x.func.size(1); xNcol = x.func.size(2);
yMrow = y.func.size(1); yNcol = y.func.size(2);
switch caller
  case 'mtimes'
    FMrow = xMrow;  FNcol = yNcol;
  case 'mldivide'
    FMrow = xNcol;  FNcol = yNcol;
  case 'mrdivide'
    FMrow = xMrow;  FNcol = yMrow;
  otherwise
   error('??? invalid caller specified') 
end


dxind = x.deriv(Vcount).nzlocs; dyind = y.deriv(Vcount).nzlocs;
nzx = size(dxind,1); nzy = size(dyind,1);
nv = ADIGATOR.VAROFDIFF(Vcount).usize;
if ~isempty(dxind) && ~isempty(dyind)
  % ---------------------X and Y have Derivatives---------------------- %
  % derivatives of Z from Y can be found by
  % X*reshape(DY,yMrow,yNcol*nv) where DY is the projection of DY using
  % the subs index in dyind.
  % derivatives of Z from X can be found using a similiar approach,
  % though we say Z.' = Y.'*X.' and find the derivatives of Z.'
  if ~strcmp(caller,'mtimes')
    error('?? caller can only be mtimes if both variables have derivatives')
  end
  
  % Get DZ from Y first.
  dy  = sparse(dyind(:,1),dyind(:,2),1:nzy,yMrow*yNcol,nv);
  dzy = reshape(xtemp*reshape(dy,yMrow,yNcol*nv),FMrow*FNcol,nv);
  [dzy_rows, dzy_cols] = find(dzy);
  if size(dzy_rows,2) > 1;dzy_rows = dzy_rows.';dzy_cols = dzy_cols.';end
  
  % Get DZ from X.
  xTranMap    = zeros(xMrow,xNcol);
  xTranMap(:) = 1:xMrow*xNcol;
  xTranMap    = xTranMap.';
  dx          = sparse(dxind(:,1),dxind(:,2),1:nzx,xMrow*xNcol,nv);
  dxTran      = dx(xTranMap(:),:); % DX.'
  dzxTran     = reshape(ytemp.'*reshape(dxTran,xNcol,xMrow*nv),FMrow*FNcol,nv);
  zTranMap    = zeros(FNcol,FMrow); zTranMap(:) = 1:FMrow*FNcol;
  zTranMap    = zTranMap.';
  dzx         = dzxTran(zTranMap(:),:);
  [dzx_rows, dzx_cols] = find(dzx);
  if size(dzx_rows,2) > 1;dzx_rows = dzx_rows.';dzx_cols = dzx_cols.';end
  if ~isempty(dzy_rows) || ~isempty(dzx_rows)
    % Can Now build DZ
    dz = sparse([dzy_rows;dzx_rows],[dzy_cols;dzx_cols],...
      [-ones(size(dzy_rows));2*ones(size(dzx_rows))],FMrow*FNcol,nv);
    % DZ with val < 2 corresponds to Y, DZ with val > 0 corresponds with
    % X
    [zrows,zcols] = find(dz);
    if size(zrows,2) > 1; zrows = zrows.'; zcols = zcols.'; end
    nzlocs = [zrows,zcols];
    
    if DPFLAG
      % ----------------Derivative Printing---------------------------- %
      nz = length(zrows);
      % Print out the Calculations for DZX first. -
      TD1 = ['cada',NDstr,'td1'];
      if ~isempty(dzx_rows)
        % We need to project the derivatives of DX into a matrix of the
        % form of reshape(dxTran,xNcol,xMrow*nv).
        TD2 = ['cada',NDstr,'td2'];
        if xMrow*xNcol*nv > 250 && nzx < .6*xMrow*xNcol*nv
          SPXflag = 1;
          % Do the projection sparsely - dxTran is larger than 250 and
          % has at least 40% zeros.
          xTranInds = zeros(nzx,3);
          [xTranInds(:,2), xTranInds(:,3), xTranInds(:,1)] = ...
            find(reshape(dxTran,xNcol,xMrow*nv));
          xTranInds = sortrows(xTranInds);
          
          % 2nd and 3rd column of xTranInds now correspond with the row
          % and column index to project dx(as is in file) into dxTran
          % reshaped
          TDind1 = cadaindprint(xTranInds(:,2));
          TDind2 = cadaindprint(xTranInds(:,3));
          % Print out the Projection
          fprintf(fid,[indent,TD2,' = sparse(',TDind1,',',TDind2,',',...
            x.deriv(Vcount).name,',%1.0d,%1.0d);\n'],xNcol,xMrow*nv);
        else
          SPXflag = 0;
          % Just project into zeros using a linear index - linear index
          % of dxTran and dxTran-reshaped are the same
          xTranInds = zeros(nzx,2);
          [xTranInds(:,2),~,xTranInds(:,1)] = find(dxTran(:));
          xTranInds = sortrows(xTranInds);
          % 2nd column of xTranInds is now the linear index to project
          % dx(as is in file) into reshaped dxTran
          TDind1 = cadaindprint(xTranInds(:,2));
          fprintf(fid,[indent,TD2,' = zeros(%1.0d,%1.0d);\n'],...
            xNcol,xMrow*nv);
          fprintf(fid,[indent,TD2,'(',TDind1,') = ',...
            x.deriv(Vcount).name,';\n']);
        end
        % Can now print out Y.'*DX.' - where DX.' is in file as TD2.
        fprintf(fid,[indent,TD2,' = ',y.func.name,'.''*',TD2,';\n']);
        % Need to get the mapping from dzxTran to dzx.
        [dzxTrows,dzxTcols] = find(dzxTran);
        dzxTran  = sparse(dzxTrows,dzxTcols,...
          sub2ind([FMrow*FNcol,nv],dzxTrows,dzxTcols),FMrow*FNcol,nv);
        dzx      = dzxTran(zTranMap(:),:);
        dzxTinds = nonzeros(dzx);
        Dind1    = cadaindprint(dzxTinds(:));
        % dzxTinds are now the Linear Reference index off of dzxTran (in
        % the file) to the nonzeros in dzx.
        %
        % need to find which parts of DZ correspond to DZX
        dzxInds = 1:nz;
        dzxInds = dzxInds(nonzeros(dz)>0);
        if length(dzxInds) == nz
          if SPXflag
            fprintf(fid,[indent,TD1,' = full(',TD2,'(',Dind1,'));\n']);
          else
            fprintf(fid,[indent,TD1,' = ',TD2,'(',Dind1,');\n']);
          end
          fprintf(fid,[indent,TD1,' = ',TD1,'(:);\n']);
        else
          Dind2 = cadaindprint(dzxInds);
          fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],nz);
          fprintf(fid,[indent,TD1,'(',Dind2,') = ',TD2,'(',Dind1,');\n']);
        end
      end
      
      if ~isempty(dzy_rows)
        % Print out the Calculations for DZY
        TD2 = ['cada',NDstr,'td2'];
        % Need to project the vector of DY(in the file) into the matrix
        % of the form reshape(DY,yMrow,yNcol*nv)
        if yMrow*yNcol*nv > 250 && nzy < .6*yMrow*yNcol*nv
          SPYflag = 1;
          % Do the projection sparsely - dy is larger than 250 and has at
          % least %40 zeros.
          dyInds = zeros(nzy,2);
          [dyInds(:,1),dyInds(:,2)] = find(reshape(dy,yMrow,yNcol*nv));
          TDind1 = cadaindprint(dyInds(:,1));
          TDind2 = cadaindprint(dyInds(:,2));
          % Print out the Projection
          fprintf(fid,[indent,TD2,' = sparse(',TDind1,',',TDind2,',',...
            y.deriv(Vcount).name,',%1.0d,%1.0d);\n'],yMrow,yNcol*nv);
        else
          SPYflag = 0;
          % Project into zeros using a linear index.
          dyInds = sub2ind([yMrow*yNcol,nv],dyind(:,1),dyind(:,2));
          TDind1 = cadaindprint(dyInds);
          % Print out the Projection
          fprintf(fid,[indent,TD2,' = zeros(%1.0d,%1.0d);\n'],...
            yMrow,yNcol*nv);
          fprintf(fid,[indent,TD2,'(',TDind1,') = ',...
            y.deriv(Vcount).name,';\n']);
        end
        % Can now print out X*DY - where DY is TD2
        fprintf(fid,[indent,TD2,' = ',x.func.name,'*',TD2,';\n']);
        fprintf(fid,[indent,TD2,' = ',TD2,'(:);\n']);
        % TD2 in the file now corresponds with DZY (except reshaped), so
        % we need to know what are the nonzero indices of DZY which we
        % need to reference off of TD2 (even though it is reshaped,
        % linear indexing will still be the same)
        dyInds = sub2ind([FMrow*FNcol,nv],dzy_rows,dzy_cols);
        Dind1  = cadaindprint(dyInds(:));
        % dyInds are now the linear reference off of TD2 which give the
        % nonzeros in DZY - just need to figure out the mapping of DZY
        % into DZ now.
        dzyInds = 1:nz;
        dzyInds = dzyInds(nonzeros(dz) < 2);
        if length(dzyInds) == nz
          if SPYflag
            if ~isempty(dzx_rows)
              fprintf(fid,[indent,TD1,' = ',TD1,' + full(',TD2,'(',Dind1,'));\n']);
            else
              fprintf(fid,[indent,TD1,' = full(',TD2,'(',Dind1,'));\n']);
            end
          else
            if ~isempty(dzx_rows)
              fprintf(fid,[indent,TD1,' = ',TD1,' + ',TD2,'(',Dind1,');\n']);
            else
              fprintf(fid,[indent,TD1,' = ',TD2,'(',Dind1,');\n']);
            end
          end
        else
          Dind2 = cadaindprint(dzyInds(:));
          fprintf(fid,[indent,TD1,'(',Dind2,') = ',...
            TD1,'(',Dind2,') + ',TD2,'(',Dind1,');\n']);
        end
      end
      fprintf(fid,[indent,derivstr,' = ',TD1,';\n']);
    end
  else
    nzlocs = [];
  end
elseif ~isempty(dxind)
  % ---------------------X has Derivatives----------------------------- %
  % we can compute DZ' = Y'*DX' where DX' is reshaped.
  
  % Get DZ from X.
  xTranMap = zeros(xMrow,xNcol); xTranMap(:) = 1:xMrow*xNcol;
  xTranMap = xTranMap.';
  dx       = sparse(dxind(:,1),dxind(:,2),1:nzx,xMrow*xNcol,nv);
  dxTran   = dx(xTranMap(:),:); % DX.'
  switch caller
    case 'mtimes'
      dzxTran  = reshape(ytemp.'*reshape(dxTran,xNcol,xMrow*nv),FMrow*FNcol,nv);
    case 'mrdivide'
      dzxTran  = reshape(ytemp.'\reshape(dxTran,xNcol,xMrow*nv),FMrow*FNcol,nv);
    case 'mldivide'
      error('X should not have derivatives if called from mldivide')
  end
  
  zTranMap = zeros(FNcol,FMrow); zTranMap(:) = 1:FMrow*FNcol;
  zTranMap = zTranMap.';
  dzx      = dzxTran(zTranMap(:),:);
  [dzx_rows, dzx_cols] = find(dzx);
  if size(dzx_rows,2) > 1;dzx_rows = dzx_rows.';dzx_cols = dzx_cols.';end
  
  if ~isempty(dzx_rows)
    nzlocs = [dzx_rows,dzx_cols];
    if DPFLAG
      % ------------------Derivative Printing-------------------------- %
      % We need to project the derivatives of DX into a matrix of the
      % form of reshape(dxTran,xNcol,xMrow*nv).
      TD1 = ['cada',NDstr,'td1'];
      if xMrow*xNcol*nv > 250 && nzx < .6*xMrow*xNcol*nv
        SPXflag = 1;
        % Do the projection sparsely - dxTran is larger than 250 and
        % has at least 40% zeros.
        xTranInds = zeros(nzx,3);
        [xTranInds(:,2), xTranInds(:,3), xTranInds(:,1)] = ...
          find(reshape(dxTran,xNcol,xMrow*nv));
        xTranInds = sortrows(xTranInds);
        % 2nd and 3rd column of xTranInds now correspond with the row
        % and column index to project dx(as is in file) into dxTran
        % reshaped
        TDind1 = cadaindprint(xTranInds(:,2));
        TDind2 = cadaindprint(xTranInds(:,3));
        % Print out the Projection
        fprintf(fid,[indent,TD1,' = sparse(',TDind1,',',TDind2,',',...
          x.deriv(Vcount).name,',%1.0d,%1.0d);\n'],xNcol,xMrow*nv);
      else
        SPXflag = 0;
        % Just project into zeros using a linear index - linear index
        % of dxTran and dxTran-reshaped are the same
        xTranInds = zeros(nzx,2);
        [xTranInds(:,2),~,xTranInds(:,1)] = find(dxTran(:));
        xTranInds = sortrows(xTranInds);
        % 2nd column of xTranInds is now the linear index to project
        % dx(as is in file) into reshaped dxTran
        TDind1 = cadaindprint(xTranInds(:,2));
        fprintf(fid,[indent,TD1,' = zeros(%1.0d,%1.0d);\n'],...
          xNcol,xMrow*nv);
        fprintf(fid,[indent,TD1,'(',TDind1,') = ',...
          x.deriv(Vcount).name,';\n']);
      end
      % Can now print out Y.'*DX.' - where DX.' is in file as TD2.
      if strcmp(caller,'mrdivide')
        fprintf(fid,[indent,TD1,' = ',y.func.name,'.''\\',TD1,';\n']);
      else
        fprintf(fid,[indent,TD1,' = ',y.func.name,'.''*',TD1,';\n']);
      end
      % Need to get the mapping from dzxTran to dzx.
      [dzxTrows,dzxTcols] = find(dzxTran);
      dzxTran  = sparse(dzxTrows,dzxTcols,...
        sub2ind([FMrow*FNcol,nv],dzxTrows,dzxTcols),FMrow*FNcol,nv);
      dzx      = dzxTran(zTranMap(:),:);
      dzxTinds = nonzeros(dzx);
      Dind1    = cadaindprint(dzxTinds(:));
      % dzxTinds are now the Linear Reference index off of dzxTran (in
      % the file) to the nonzeros in dzx.
      fprintf(fid,[indent,TD1,' = ',TD1,'(:);\n']);
      if SPXflag
        fprintf(fid,[indent,derivstr,' = full(',TD1,'(',Dind1,'));\n']);
      else
        fprintf(fid,[indent,derivstr,' = ',TD1,'(',Dind1,');\n']);
      end
    end
  else
    nzlocs = [];
  end
elseif ~isempty(dyind)
  % ---------------------Y has Derivatives----------------------------- %
  % We say DZ = X*DY - where DY is reshaped.
  % If mldivide: DZ = X\DY
  % Get DZ from Y first.
  dy  = sparse(dyind(:,1),dyind(:,2),1:nzy,yMrow*yNcol,nv);
  switch caller
    case 'mtimes'
      dzy = reshape(xtemp*reshape(dy,yMrow,yNcol*nv),FMrow*FNcol,nv);
    case 'mldivide'
      dzy = reshape(xtemp\reshape(dy,yMrow,yNcol*nv),FMrow*FNcol,nv);
    case 'mrdivide'
      error('Y should not have derivatives if called from mrdivide')
  end
  
  [dzy_rows, dzy_cols] = find(dzy);
  if size(dzy_rows,2) > 1;dzy_rows = dzy_rows.';dzy_cols = dzy_cols.';end
  
  if ~isempty(dzy_rows)
    nzlocs = [dzy_rows,dzy_cols];
    if DPFLAG
      % ------------------Derivative Printing-------------------------- %
      % Print out the Calculations for DZY
      TD1 = ['cada',NDstr,'td1'];
      % Need to project the vector of DY(in the file) into the matrix
      % of the form reshape(DY,yMrow,yNcol*nv)
      if yMrow*yNcol*nv > 250 && nzy < .6*yMrow*yNcol*nv
        SPYflag = 1;
        % Do the projection sparsely - dy is larger than 250 and has at
        % least %40 zeros.
        dyInds = zeros(nzy,2);
        [dyInds(:,1),dyInds(:,2)] = find(reshape(dy,yMrow,yNcol*nv));
        TDind1 = cadaindprint(dyInds(:,1));
        TDind2 = cadaindprint(dyInds(:,2));
        % Print out the Projection
        fprintf(fid,[indent,TD1,' = sparse(',TDind1,',',TDind2,',',...
          y.deriv(Vcount).name,',%1.0d,%1.0d);\n'],yMrow,yNcol*nv);
      else
        SPYflag = 0;
        % Project into zeros using a linear index.
        dyInds = sub2ind([yMrow*yNcol,nv],dyind(:,1),dyind(:,2));
        TDind1 = cadaindprint(dyInds(:));
        % Print out the Projection
        fprintf(fid,[indent,TD1,' = zeros(%1.0d,%1.0d);\n'],...
          yMrow,yNcol*nv);
        fprintf(fid,[indent,TD1,'(',TDind1,') = ',...
          y.deriv(Vcount).name,';\n']);
      end
      % Can now print out X*DY - where DY is TD2
      if strcmp(caller,'mldivide')
        fprintf(fid,[indent,TD1,' = ',x.func.name,'\\',TD1,';\n']);
      else
        fprintf(fid,[indent,TD1,' = ',x.func.name,'*',TD1,';\n']);
      end
      % TD2 in the file now corresponds with DZY (except reshaped), so
      % we need to know what are the nonzero indices of DZY which we
      % need to reference off of TD2 (even though it is reshaped,
      % linear indexing will still be the same)
      dyInds = sub2ind([FMrow*FNcol,nv],dzy_rows,dzy_cols);
      Dind1  = cadaindprint(dyInds(:));
      % dyInds are now the linear reference off of TD2 which give the
      % nonzeros in DZY
      fprintf(fid,[indent,TD1,' = ',TD1,'(:);\n']);
      if SPYflag
        fprintf(fid,[indent,derivstr,' = full(',TD1,'(',Dind1,'));\n']);
      else
        fprintf(fid,[indent,derivstr,' = ',TD1,'(',Dind1,');\n']);
      end
    end
  else
    nzlocs = [];
  end
end