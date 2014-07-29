function NewVar = cadaPrintReMap(x,OverVar,VarCount)
% NewVar = cadaPrintReMap(x,OverVar,varID)
% Function prints out ReMapping of x to OverVar, where OverVar is what
% the variable needs to be, and x is what the variable currently is.
%
% Note: this is also called if want to remove overmapping (i.e. x is the
% overmapped variable, and OverVar is what it needs to be)
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
fid = ADIGATOR.PRINT.FID; indent = ADIGATOR.PRINT.INDENT;
if ~isa(OverVar,'cada') && isempty(OverVar)
  NewVar = x;
  return
end
NewVar.id   = VarCount;
NewVar.func = OverVar.func;


[funcstr,DPflag] = cadafuncname(VarCount);
if ~VarCount
  funcstr = x.func.name;
end
NewVar.func.name = funcstr;

% 
% Dbstuff = dbstack; CallingFile = Dbstuff(2).file;
% fprintf(fid,[indent,'%% Remap ',funcstr,' from ',CallingFile,'\n']);

%---Function Sizes--- 
xMrow = x.func.size(1);       xNcol = x.func.size(2);
oMrow = OverVar.func.size(1); oNcol = OverVar.func.size(2);
% Vectorized Stuff
vecerr = 0;
if (isinf(xMrow) && isinf(oMrow)) || (isinf(xMrow) && ~oMrow) ...
    || (~xMrow && isinf(oMrow))
  if isinf(xNcol) || isinf(oNcol); vecerr=1; end
  vecDim = ['size(',funcstr,',1)'];
  zvec = 1; xMrow = 1; oMrow = 1;
elseif (isinf(xNcol) && isinf(oNcol)) || (isinf(xNcol) && ~oNcol) ...
    || (~xNcol && isinf(oNcol))
  if isinf(xMrow) || isinf(oMrow); vecerr=1; end
  zvec = 2; xNcol = 1; oNcol = 1;
  vecDim = ['size(',funcstr,',2)'];
else
  if isinf(xNcol) || isinf(oNcol) || isinf(xMrow) || isinf(oMrow)
    vecerr = 1;
  end
  zvec = 0;
end
if vecerr
  error(['If a variable may take on different sizes depending upon ',...
    'loop iterations or different branches of a conditional block ',...
    'then it must not switch between vectorized to non-vectorized'])
end

% xref/oref used to remap derivatives if var is changing sizes
if xMrow*xNcol > oMrow*oNcol
  nf = xMrow*xNcol;
  xref = zeros(xMrow,xNcol); xref(:) = 1:xMrow*xNcol;
  oref = xref(1:oMrow,1:oNcol);
else
  nf = oMrow*oNcol;
  oref = zeros(oMrow,oNcol); oref(:) = 1:oMrow*oNcol;
  xref = oref(1:xMrow,1:xNcol);
end

%---Derivatives---
NewVar.deriv = OverVar.deriv;
for Vcount = 1:ADIGATOR.NVAROFDIFF
  if ~isempty(OverVar.deriv(Vcount).nzlocs)
    if VarCount
      derivstr = cadadername(funcstr,Vcount,VarCount);
    else
      derivstr = x.deriv(Vcount).name;
    end
    NewVar.deriv(Vcount).name = derivstr;
    if DPflag && ~ADIGATOR.EMPTYFLAG
      nzover = size(OverVar.deriv(Vcount).nzlocs,1);
      if ~isempty(x.deriv(Vcount).nzlocs)
        orows = oref(OverVar.deriv(Vcount).nzlocs(:,1));
        ocols = OverVar.deriv(Vcount).nzlocs(:,2);
        xrows = xref(x.deriv(Vcount).nzlocs(:,1));
        xcols = x.deriv(Vcount).nzlocs(:,2);
        if ~isequal(xrows,orows) && ~isequal(xcols,ocols)
          % Derivatives are not the same, need to remap
          nzx = size(x.deriv(Vcount).nzlocs,1);
          nv = ADIGATOR.VAROFDIFF(Vcount).usize;
          dover = sparse(orows,ocols,2*ones(nzover,1),nf,nv);
          dx = sparse(xrows,xcols,ones(size(xrows)),nf,nv);
          dc = dover+dx;
          dlocs = nonzeros(dc);
          dclocs = find(dlocs==3);
          nzd = length(dclocs);
          % ----------------Print out re-Mapped Derivatives-------------- %
          if ~(nzd == nzx && nzx == nzover)
            % Derivatives are not the same
            if nzd == nzx
              % the common derivatives are the same as all of the derivatives in
              % x -
              Dind1 = cadaindprint(dclocs);
              TempDerivName = sprintf(['cada%1.0ftempd',ADIGATOR.VAROFDIFF(Vcount).name,'%1.0f'],ADIGATOR.DERNUMBER);
              fprintf(fid,[indent,TempDerivName,' = ',derivstr,';\n']);
              if zvec
                fprintf(fid,[indent,derivstr,' = zeros(',vecDim,',%1.0f);\n'],nzover);
                fprintf(fid,[indent,derivstr,'(:,',Dind1,') = ',TempDerivName,';\n']);
              else
                fprintf(fid,[indent,derivstr,' = zeros(%1.0f,1);\n'],nzover);
                fprintf(fid,[indent,derivstr,'(',Dind1,',1) = ',TempDerivName,';\n']);
              end
            elseif nzd == nzover
              % the common derivatives are the same as all of the derivatives in
              % the overmapped variable -
              Dind1 = cadaindprint(dclocs);
              if zvec
                fprintf(fid,[indent,derivstr,' = ',derivstr,'(:,',Dind1,');\n']);
              else
                fprintf(fid,[indent,derivstr,' = ',derivstr,'(',Dind1,',1);\n']);
              end
            else
              % neither of the above - we need dover(doverlocs) = dx(dxlocs)
              oind = 1:nzover;
              doverlocs = dlocs(dlocs>1);
              doverlocs = oind(doverlocs==3);
              xind = 1:nzx;
              dxlocs = dlocs(dlocs==1 | dlocs==3);
              dxlocs = xind(dxlocs==3);
              Dind1 = cadaindprint(dxlocs);
              TempDerivName = sprintf(['cada%1.0ftempd',ADIGATOR.VAROFDIFF(Vcount).name,'%1.0f'],ADIGATOR.DERNUMBER);
              if zvec
                fprintf(fid,[indent,TempDerivName,' = ',derivstr,'(:,',Dind1,');\n']);
                Dind2 = cadaindprint(doverlocs);
                fprintf(fid,[indent,derivstr,' = zeros(',vecDim,',%1.0f);\n'],nzover);
                fprintf(fid,[indent,derivstr,'(:,',Dind2,') = ',TempDerivName,';\n']);
              else
                fprintf(fid,[indent,TempDerivName,' = ',derivstr,'(',Dind1,',1);\n']);
                Dind2 = cadaindprint(doverlocs);
                fprintf(fid,[indent,derivstr,' = zeros(%1.0f,1);\n'],nzover);
                fprintf(fid,[indent,derivstr,'(',Dind2,',1) = ',TempDerivName,';\n']);
              end
            end
          end
        end
      else
        % Current variable doesnt have derivatives, just make it a zero
        % vector
        if zvec
          fprintf(fid,[indent,derivstr,' = zeros(',vecDim,',%1.0f);\n'],nzover);
        else
          fprintf(fid,[indent,derivstr,' = zeros(%1.0f,1);\n'],nzover);
        end
      end
    end
  elseif ~isempty(x.deriv(Vcount).nzlocs)
    % New var doesnt have derivs but x does - we want to remove deriv
    % field.
    if VarCount
      derivstr = cadadername(funcstr,Vcount,VarCount);
    else
      derivstr = x.deriv(Vcount).name;
    end
    NewVar.deriv(Vcount).name = [];
    NewVar.deriv(Vcount).nzlocs = [];
    derfieldloc = strfind(derivstr,'.');
    if ~isempty(derfieldloc) && DPflag
      derfield = derivstr(derfieldloc(end)+1:end);
      varname  = derivstr(1:derfieldloc(end)-1);
      if strcmp(derfield(1),'d')
        fprintf(fid,[indent,varname,' = rmfield(',varname,',''',derfield,''');\n']);
      end
    end
  end
end

% -----------------Print out re-mapped Functions if needed--------------- %
if ~ADIGATOR.EMPTYFLAG
  if xMrow > oMrow % remove rows
    if xNcol > oNcol % remove cols
      fprintf(fid,[indent,funcstr,' = ',funcstr,'(1:%1.0f,1:%1.0f);\n'],oMrow,oNcol);
    elseif xNcol < oNcol % add cols
      error('This shouldnt happen');
    else % cols same
      fprintf(fid,[indent,funcstr,' = ',funcstr,'(1:%1.0f,:);\n'],oMrow);
    end
  elseif xMrow < oMrow % add rows
    if xNcol > oNcol % remove cols
      error('This shouldnt happen');
    else % add rows/cols same
      if zvec == 2
        fprintf(fid,[indent,funcstr,'(%1.0f,:) = 0;\n'],oMrow);
      else
        fprintf(fid,[indent,funcstr,'(%1.0f,%1.0f) = 0;\n'],oMrow,oNcol);
      end
    end
  else % rows same
    if  xNcol > oNcol % remove cols
      fprintf(fid,[indent,funcstr,' = ',funcstr,'(:,1:%1.0f);\n'],oNcol);
    elseif xNcol < oNcol % add cols
      if zvec == 1
        fprintf(fid,[indent,funcstr,'(:,%1.0f) = 0;\n'],oNcol);
      else
        fprintf(fid,[indent,funcstr,'(%1.0f,%1.0f) = 0;\n'],oMrow,oNcol);
      end
    end
  end
end
NewVar = class(NewVar,'cada');