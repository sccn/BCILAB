function [yName,yInds] = cadaRepDers(xName,xInds,ySize,Vcount,DPFLAG)
% This function is used to repmat derivatives of scalar variables.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
fid    = ADIGATOR.PRINT.FID;
indent = ADIGATOR.PRINT.INDENT;
NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);

nv = ADIGATOR.VAROFDIFF(Vcount).usize;
dx = sparse(xInds(:,1),xInds(:,2),1:size(xInds,1),1,nv);
dy = repmat(dx,ySize,1);

[yrows,ycols,yind] = find(dy);

if size(yrows,2) > 1
  yrows = yrows.'; ycols = ycols.';
end
yInds = [yrows,ycols];

if DPFLAG
  yName = ['cada',NDstr,'tempd',ADIGATOR.VAROFDIFF(Vcount).name];
  Dind1 = cadaindprint(yind);
  fprintf(fid,[indent,yName,' = ',xName,'(',Dind1,');\n']);
else
  yName = xName;
end