function z = cadabinaryarraymath(x,y,xzeroflag,yzeroflag,callerstr)
% Binary math function - any binary mathematical functions should call this
% function in order to print the derivative. The derivative of each calling
% function must be defined in the getdzdx/getdzdy sub-functions.
% xzeroflag => if x(i) = 0, then [dz/dy](i,:) = 0
% yzeroflag => if y(i) = 0, then [dz/dx](i,:) = 0;
%
% Copyright 2011-2014 Matthew J. Weinstein, Michael Patterson and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.EMPTYFLAG
  z = cadaEmptyEval(x,y);
  return
end

NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);
xscalarflag = 0;
yscalarflag = 0;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~ Parse Inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
if isa(x,'cada') && isa(y,'cada')
  % Both Inputs are Symbolic
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  yMrow = y.func.size(1); yNcol = y.func.size(2);
  if isinf(xMrow); xvec = 2; elseif isinf(xNcol); xvec =1; else xvec = 0;end
  if isinf(yMrow); yvec = 2; elseif isinf(yNcol); yvec =1; else yvec = 0;end
elseif isa(x,'cada')
  % y is numeric input
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  if isinf(xMrow); xvec = 2; elseif isinf(xNcol); xvec =1; else xvec = 0;end
  yvec = 0;
  [yMrow,yNcol] = size(y);
  ytemp.id = [];
  ytemp.func = struct('name',[],'size',[yMrow yNcol],'zerolocs',[],...
    'value',y);
  if PFLAG
    if yMrow*yNcol == 1
      ytemp.func.name = num2str(y,16);
    else
      ytemp.func.name = cadamatprint(y);
    end
  end
  ytemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  y = ytemp;
  y = class(y,'cada');
else
  % x is numeric input
  yMrow = y.func.size(1); yNcol = y.func.size(2);
  if isinf(yMrow); yvec = 2; elseif isinf(yNcol); yvec =1; else yvec = 0;end
  xvec = 0;
  [xMrow,xNcol] = size(x);
  xtemp.id = [];
  xtemp.func = struct('name',[],'size',[xMrow xNcol],'zerolocs',[],...
    'value',x);
  if PFLAG
    if xMrow*xNcol == 1
      xtemp.func.name = num2str(x,16);
    else
      xtemp.func.name = cadamatprint(x);
    end
  end
  xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  x = xtemp;
  x = class(x,'cada');
end
%% ~~~~~~~~~~~~~~~~~~~~~~~~ Function Sizing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
if (xMrow == yMrow && xNcol == yNcol)
  FMrow = yMrow; FNcol = yNcol;
  if FMrow == 1 && FNcol == 1; xscalarflag = 1; yscalarflag = 1; end
elseif (xMrow == 1 && xNcol == 1)
  FMrow = yMrow; FNcol = yNcol; xscalarflag = 1;
elseif (yMrow == 1 && yNcol == 1)
  FMrow = xMrow; FNcol = xNcol; yscalarflag = 1;
elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
  % In printing run of a FOR loop - Sizes could be off due to IF statement
  if xMrow <= yMrow && xNcol <= yNcol
    y = cadaRemoveRowsCols(y,[xMrow xNcol]);
    yMrow = xMrow; yNcol = xNcol;
  elseif yMrow <= xMrow && yNcol <= xNcol
    x = cadaRemoveRowsCols(x,[yMrow yNcol]);
    xMrow = yMrow; xNcol = yNcol;
  else
    error('Inputs are not of compatible sizes');
  end
  FMrow = xMrow; FNcol = xNcol;
else
  keyboard
  error('Inputs are not of compatible sizes');
end
if isinf(FMrow); zvec = 2; elseif isinf(FNcol); zvec =1; else zvec = 0;end

if (xvec && ~yvec && cadaCheckForDerivs(y)) || ...
    (yvec && ~xvec && cadaCheckForDerivs(x))
  error(['cannot perform ',callerstr,' if scaler input has derivatives and ',...
    'other input is vectorized']);
end




%% ~~~~~~~~~~~~~~~~~~~~ Build Function Properties ~~~~~~~~~~~~~~~~~~~~~~ %%
z.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
z.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
  'value',[]);

% Get Rid of vectorized dims
if isinf(xMrow); xMrow = 1; elseif isinf(xNcol); xNcol = 1; end
if isinf(yMrow); yMrow = 1; elseif isinf(yNcol); yNcol = 1; end
if isinf(FMrow); FMrow = 1; elseif isinf(FNcol); FNcol = 1; end

% ------------------- Logical Reference Check --------------------------- %
if isfield(x.func,'logicref') && isfield(y.func,'logicref') && ...
    isequal(x.func.logicref,y.func.logicref)
  z.func.logicref = x.func.logicref;
elseif isfield(x.func,'logicref') && yscalarflag
  z.func.logicref = x.func.logicref;
elseif isfield(y.func,'logicref') && xscalarflag
  z.func.logicref = y.func.logicref;
elseif isfield(x.func,'logicref') || isfield(y.func,'logicref')
  errorstr = sprintf(['if performing binary array operation on an object which',...
    ' is the result of a symbolic logical reference, then both inputs ',...
    ' to the array command must be the result of a reference using the',...
    ' same logical reference index:\n',...
    ' Example: ',...
    ' ind = x > y;\nz = x(ind).*y(ind);\nis a valid set of operations.\n',...
    ' z = x(x>y).*y(x>y) is not.']);
  error(errorstr);
end

% --------------Function Numeric Values and Sparsity--------------------- %
if ~isempty(x.func.value) && ~isempty(y.func.value)
  % z is numeric
  callerfunc   = str2func(callerstr);
  z.func.value = callerfunc(x.func.value,y.func.value);
else
  spflag = 0;
  if ~isempty(x.func.value)
    xtemp = logical(x.func.value); spflag = 1;
  elseif ~isempty(x.func.zerolocs)
    xtemp = true(xMrow,xNcol); xtemp(x.func.zerolocs) = false; spflag = 1;
  else
    xtemp = true(xMrow,xNcol);
  end
  if ~isempty(y.func.value)
    ytemp = logical(y.func.value); spflag = 1;
  elseif ~isempty(y.func.zerolocs)
    ytemp = true(yMrow,yNcol); ytemp(y.func.zerolocs) = false; spflag = 1;
  else
    ytemp = true(yMrow,yNcol);
  end
  if spflag == 1
    switch callerstr
      case {'plus','minus'}
        ztemp = or(xtemp,ytemp);
      case 'times'
        ztemp = and(xtemp,ytemp);
      case {'rdivide','power','atan2'}
        ztemp = xtemp;
    end
    z.func.zerolocs = find(~ztemp(:));
    if length(z.func.zerolocs) == FMrow*FNcol
      z.func.zerolocs = []; z.func.value = zeros(FMrow,FNcol);
    end
  end
end

%% ~~~~~~~~~~~~~~~~~~~~ Build Derivative Properties ~~~~~~~~~~~~~~~~~~~~ %%
z.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod;
  dxind = x.deriv(Vcount).nzlocs;
  dyind = y.deriv(Vcount).nzlocs;
  % RepMat scalars if needed
  if ~isempty(dxind) && xscalarflag && ~yscalarflag
    [x.deriv(Vcount).name,dxind] =...
        cadaRepDers(x.deriv(Vcount).name,dxind,yMrow*yNcol,Vcount,DPFLAG);
  elseif ~isempty(dyind) && yscalarflag && ~xscalarflag
    [y.deriv(Vcount).name,dyind] =...
        cadaRepDers(y.deriv(Vcount).name,dyind,xMrow*xNcol,Vcount,DPFLAG);
  end
  % ---------------Use Function Sparsity to Cancel Derivatives----------- %
  spdyflag = 0;
  if xzeroflag && xscalarflag && ~isempty(dyind)
    % --x is scalar 0, remove all dy--
    if xtemp == 0; spdyind = []; else spdyind = dyind; end
  elseif xzeroflag && (~isempty(x.func.zerolocs) ||...
      ~isempty(x.func.value)) && ~isempty(dyind)
    % --x is sparse or numeric, remove elements of dy if possible--
    fx = xtemp(dyind(:,1)); nzlocsy = logical(fx(:).*dyind(:,1));
    spdyind = dyind(nzlocsy,:);
    if size(spdyind,1) < size(dyind,1); spdyflag = 1; end
  else
    % --leave dy as is--
    spdyind = dyind;
  end
  spdxflag = 0;
  if yzeroflag && yscalarflag && ~isempty(dxind)
    %--y is scalar 0, remove all dx--
    if ytemp == 0; spdxind = []; else spdxind = dxind; end
  elseif yzeroflag && (~isempty(y.func.zerolocs) ||...
      ~isempty(y.func.value)) && ~isempty(dxind)
    % --y is sparse or numeric, remove elements of dx if possible--
    fy = ytemp(dxind(:,1)); nzlocsx = logical(fy(:).*dxind(:,1));
    spdxind = dxind(nzlocsx,:);
    if size(spdxind,1) < size(dxind,1); spdxflag = 1; end
  else
    % --leave dx as is--
    spdxind = dxind;
  end
  if ~isempty(spdxind) && ~isempty(spdyind)
    %% ---------------------X and Y have Derivatives-------------------- %%
    derivstr = cadadername(funcstr,Vcount);
    z.deriv(Vcount).name = derivstr;
    if zvec
      [z.deriv(Vcount).nzlocs,dzxind,dzyind,xindflag,yindflag] = ...
        cadaunion(spdxind,spdyind,z.func.size(zvec),ADIGATOR.VAROFDIFF(Vcount).usize);
    else
      [z.deriv(Vcount).nzlocs,dzxind,dzyind,xindflag,yindflag] = ...
        cadaunion(spdxind,spdyind,FMrow*FNcol,ADIGATOR.VAROFDIFF(Vcount).usize);
    end
    nz = size(z.deriv(Vcount).nzlocs,1);
    nzx = size(spdxind,1); nzy = size(spdyind,1);
    % ------------------------Derivative Printing------------------------ %
    if DPFLAG == 1
      if spdxflag;dxindin = (1:size(dxind,1)).';dxindin = dxindin(nzlocsx);end
      if spdyflag;dyindin = (1:size(dyind,1)).';dyindin = dyindin(nzlocsy);end
      TD1 = ['cada',NDstr,'td1'];
      
      % ----------------------- Print dZ = [dZ/dX]*dX ------------------- %
      % Get dX
      if spdxflag
        Dind2 = cadaindprint(dxindin);
        if xvec
          DXstr = [x.deriv(Vcount).name,'(:,',Dind2,')'];
        else
          DXstr = [x.deriv(Vcount).name,'(',Dind2,')'];
        end
      else
        DXstr = x.deriv(Vcount).name;
      end
      % Get y and x that [dZ/dX] is function of
      switch callerstr
        case {'plus','minus'}
          Xstr = []; Ystr = [];
        otherwise
          % Get y
          if yscalarflag || ...
              (yvec==2 && nzx == yNcol && isequal(spdxind(:,1),(1:yNcol).'))
            Ystr = y.func.name;
          elseif ~yvec && nzx == yMrow*yNcol && isequal(spdxind(:,1),(1:yMrow*yNcol).')
            Ystr = [y.func.name,'(:)'];
          elseif yvec ==1 && nzx == yMrow && isequal(spdxind(:,1),(1:yMrow).')
            Ystr = [y.func.name,'.'''];
          else
            TF1 = ['cada',NDstr,'tf1'];
            TFind1 = cadaindprint(spdxind(:,1));
            if yvec == 2
              fprintf(fid,[indent,TF1,' = ',y.func.name,'(:,',TFind1,');\n']);
              Ystr = TF1;
            elseif yvec == 1
              fprintf(fid,[indent,TF1,' = ',y.func.name,'(',TFind1,',:).'';\n']);
              Ystr = TF1;
            else
              fprintf(fid,[indent,TF1,' = ',y.func.name,'(',TFind1,');\n']);
              Ystr = [TF1,'(:)'];
            end
          end
          switch callerstr
            case {'times','rdivide'}
              Xstr = [];
            otherwise
              if xscalarflag || ...
                  (xvec==2 && nzx == xNcol && isequal(spdxind(:,1),(1:xNcol).'))
                Xstr = x.func.name;
              elseif ~xvec && nzx == xMrow*xNcol && isequal(spdxind(:,1),(1:xMrow*xNcol).')
                Xstr = [x.func.name,'(:)'];
              elseif xvec ==1 && nzx == xMrow && isequal(spdxind(:,1),(1:xMrow).')
                Xstr = [x.func.name,'.'''];
              else
                TF2 = ['cada',NDstr,'tf2'];
                TFind1 = cadaindprint(spdxind(:,1));
                if xvec == 2
                  fprintf(fid,[indent,TF2,' = ',x.func.name,'(:,',TFind1,');\n']);
                  Xstr = TF2;
                elseif xvec == 1
                  fprintf(fid,[indent,TF2,' = ',x.func.name,'(',TFind1,',:).'';\n']);
                  Xstr = TF2;
                else
                  fprintf(fid,[indent,TF2,' = ',x.func.name,'(',TFind1,');\n']);
                  Xstr = [TF2,'(:)'];
                end
              end
          end
      end
      % get DZ assignment that [dZ/dX]*dX goes into
      if xindflag
        DZXstr = TD1;
      elseif zvec
        fprintf(fid,[indent,TD1,' = zeros(size(',x.deriv(Vcount).name,',1),%1.0f);\n'],nz);
        Dind1 = cadaindprint(dzxind);
        DZXstr = [TD1,'(:,',Dind1,')'];
      else
        fprintf(fid,[indent,TD1,' = zeros(%1.0f,1);\n'],nz);
        Dind1 = cadaindprint(dzxind);
        DZXstr = [TD1,'(',Dind1,')'];
      end
      if strcmp(callerstr,'power')
        % Guard Against x^0
        if yscalarflag
          fprintf(fid,[indent,'cadaconditional1 = ',y.func.name,';\n']);
          fprintf(fid,[indent,'if cadaconditional1\n']);
          fprintf(fid,[indent,'    ',DZXstr,' = ',getdzdx(Xstr,Ystr,DXstr,callerstr),';\n']);
          fprintf(fid,[indent,'end\n']);
        else
          TD2 = ['cada',NDstr,'td2'];
          if zvec
            fprintf(fid,[indent,TD2,' = zeros(size(',x.deriv(Vcount).name,',1),%1.0f);\n'],nzx);
          else
            fprintf(fid,[indent,TD2,' = zeros(%1.0f,1);\n'],nzx);
          end
          TF3 = ['cada',NDstr,'tf3'];
          fprintf(fid,[indent,TF3,' = ',Ystr,'~= 0;\n']);
          if isempty(regexp(Xstr,'\w$','once'))
            TF4 = ['cada',NDstr,'tf4'];
            fprintf([indent,TF4,' = ',Xstr,';\n']);
            Xstr = [TF4,'(',TF3,')'];
          else
            Xstr = [Xstr,'(',TF3,')']; %#ok<AGROW>
          end
          if isempty(regexp(Ystr,'\w$','once'))
            TF5 = ['cada',NDstr,'tf5'];
            fprintf([indent,TF5,' = ',Ystr,';\n']);
            Ystr = [TF5,'(',TF3,')'];
          else
            Ystr = [Ystr,'(',TF3,')']; %#ok<AGROW>
          end
          % Get [dZ/dX]*dX and print assignment to dZ
          fprintf(fid,[indent,TD2,'(',TF3,') = ',getdzdx(Xstr,Ystr,DXstr,callerstr),';\n']);
          fprintf(fid,[indent,DZXstr,' = ',TD2,';\n']);
        end
      else
        % Get [dZ/dX]*dX and print assignment to dZ
        fprintf(fid,[indent,DZXstr,' = ',getdzdx(Xstr,Ystr,DXstr,callerstr),';\n']);
      end
      % -------------------Print DZ = DZ + [dZ/dY]*dY-------------------- %
      % Get dY
      if spdyflag
        Dind2 = cadaindprint(dyindin);
        if yvec
          DYstr = [y.deriv(Vcount).name,'(:,',Dind2,')'];
        else
          DYstr = [y.deriv(Vcount).name,'(',Dind2,')'];
        end
      else
        DYstr = y.deriv(Vcount).name;
      end
      % Get y and x that [dZ/dY] is function of
      switch callerstr
        case {'plus','minus'}
          Xstr = [];  Ystr = [];
        otherwise
          % Get x
          if xscalarflag || ...
              (xvec==2 && nzy == xNcol && isequal(spdyind(:,1),(1:xNcol).'))
            Xstr = x.func.name;
          elseif ~xvec && nzy == xMrow*xNcol && isequal(spdyind(:,1),(1:xMrow*xNcol).')
            Xstr = [x.func.name,'(:)'];
          elseif xvec ==1 && nzy == xMrow && isequal(spdyind(:,1),(1:xMrow).')
            Xstr = [x.func.name,'.'''];
          else
            TF1 = ['cada',NDstr,'tf1'];
            TFind1 = cadaindprint(spdyind(:,1));
            if xvec == 2
              fprintf(fid,[indent,TF1,' = ',x.func.name,'(:,',TFind1,');\n']);
              Xstr = TF1;
            elseif xvec == 1
              fprintf(fid,[indent,TF1,' = ',x.func.name,'(',TFind1,',:).'';\n']);
              Xstr = TF1;
            else
              fprintf(fid,[indent,TF1,' = ',x.func.name,'(',TFind1,');\n']);
              Xstr = [TF1,'(:)'];
            end
          end
          switch callerstr
            case 'times'
              Ystr = [];
            otherwise
              % Get y
              if yscalarflag || ...
                  (yvec==2 && nzy == yNcol && isequal(spdyind(:,1),(1:yNcol).'))
                Ystr = y.func.name;
              elseif ~yvec && nzy == yMrow*yNcol && isequal(spdyind(:,1),(1:yMrow*yNcol).')
                Ystr = [y.func.name,'(:)'];
              elseif yvec ==1 && nzy == yMrow && isequal(spdyind(:,1),(1:yMrow).')
                Ystr = [y.func.name,'.'''];
              else
                TF2 = ['cada',NDstr,'tf2'];
                TFind1 = cadaindprint(spdyind(:,1));
                if yvec == 2
                  fprintf(fid,[indent,TF2,' = ',y.func.name,'(:,',TFind1,');\n']);
                  Ystr = TF2;
                elseif yvec == 1
                  fprintf(fid,[indent,TF2,' = ',y.func.name,'(',TFind1,',:).'';\n']);
                  Ystr = TF2;
                else
                  fprintf(fid,[indent,TF2,' = ',y.func.name,'(',TFind1,');\n']);
                  Ystr = [TF2,'(:)'];
                end
              end
          end
      end
      % get DZ assignment that [dZ/dY]*dY goes into
      if yindflag
        DZYstr = TD1;
      elseif zvec
        Dind1 = cadaindprint(dzyind);
        DZYstr = [TD1,'(:,',Dind1,')'];
      else
        Dind1 = cadaindprint(dzyind);
        DZYstr = [TD1,'(',Dind1,')'];
      end

      % Get [dZ/dY]*dY and print assignment to dZ
      fprintf(fid,[indent,DZYstr,' = ',DZYstr,' + ',getdzdy(Xstr,Ystr,DYstr,callerstr),';\n']);
      
      fprintf(fid,[indent,derivstr,' = ',TD1,';\n']);
    end
  elseif ~isempty(spdxind)
    %% -------------------------X has Derivatives----------------------- %%
    derivstr = cadadername(funcstr,Vcount);
    z.deriv(Vcount).name = derivstr; z.deriv(Vcount).nzlocs = spdxind;
    nzx = size(spdxind,1);
    if DPFLAG == 1
      if spdxflag;dxindin = 1:size(dxind,1);dxindin = dxindin(nzlocsx);end
      
      % ----------------------- Print dZ = [dZ/dX]*dX ------------------- %
      % Get dX
      if spdxflag
        Dind2 = cadaindprint(dxindin);
        if xvec
          DXstr = [x.deriv(Vcount).name,'(:,',Dind2,')'];
        else
          DXstr = [x.deriv(Vcount).name,'(',Dind2,')'];
        end
      else
        DXstr = x.deriv(Vcount).name;
      end
      % Get y and x that [dZ/dX] is function of
      switch callerstr
        case {'plus','minus'}
          Xstr = []; Ystr = [];
        otherwise
          % Get y
          if yscalarflag || ...
              (yvec==2 && nzx == yNcol && isequal(spdxind(:,1),(1:yNcol).'))
            Ystr = y.func.name;
          elseif ~yvec && nzx == yMrow*yNcol && isequal(spdxind(:,1),(1:yMrow*yNcol).')
            Ystr = [y.func.name,'(:)'];
          elseif yvec ==1 && nzx == yMrow && isequal(spdxind(:,1),(1:yMrow).')
            Ystr = [y.func.name,'.'''];
          else
            TF1 = ['cada',NDstr,'tf1'];
            TFind1 = cadaindprint(spdxind(:,1));
            if yvec == 2
              fprintf(fid,[indent,TF1,' = ',y.func.name,'(:,',TFind1,');\n']);
              Ystr = TF1;
            elseif yvec == 1
              fprintf(fid,[indent,TF1,' = ',y.func.name,'(',TFind1,',:).'';\n']);
              Ystr = TF1;
            else
              fprintf(fid,[indent,TF1,' = ',y.func.name,'(',TFind1,');\n']);
              Ystr = [TF1,'(:)'];
            end
          end
          switch callerstr
            case {'times','rdivide'}
              Xstr = [];
            otherwise
              if xscalarflag || ...
                  (xvec==2 && nzx == xNcol && isequal(spdxind(:,1),(1:xNcol).'))
                Xstr = x.func.name;
              elseif ~xvec && nzx == xMrow*xNcol && isequal(spdxind(:,1),(1:xMrow*xNcol).')
                Xstr = [x.func.name,'(:)'];
              elseif xvec ==1 && nzx == xMrow && isequal(spdxind(:,1),(1:xMrow).')
                Xstr = [x.func.name,'.'''];
              else
                TF2 = ['cada',NDstr,'tf2'];
                TFind1 = cadaindprint(spdxind(:,1));
                if xvec == 2
                  fprintf(fid,[indent,TF2,' = ',x.func.name,'(:,',TFind1,');\n']);
                  Xstr = TF2;
                elseif xvec == 1
                  fprintf(fid,[indent,TF2,' = ',x.func.name,'(',TFind1,',:).'';\n']);
                  Xstr = TF2;
                else
                  fprintf(fid,[indent,TF2,' = ',x.func.name,'(',TFind1,');\n']);
                  Xstr = [TF2,'(:)'];
                end
              end
          end
      end
      if strcmp(callerstr,'power') && isempty(y.func.value)
        % Guard Against x^0
        if yscalarflag
          fprintf(fid,[indent,'cadaconditional1 = ',y.func.name,';\n']);
          fprintf(fid,[indent,'if cadaconditional1\n']);
          fprintf(fid,[indent,'    ',derivstr,' = ',getdzdx(Xstr,Ystr,DXstr,callerstr),';\n']);
          fprintf(fid,[indent,'else\n']);
          if zvec
            fprintf(fid,[indent,'    ',derivstr,' = zeros(size(',x.deriv(Vcount).name,',1),%1.0f);\n'],nzx);
          else
            fprintf(fid,[indent,'    ',derivstr,' = zeros(%1.0f,1);\n'],nzx);
          end
          fprintf(fid,[indent,'end\n']);
        else
          if zvec
            fprintf(fid,[indent,derivstr,' = zeros(size(',x.deriv(Vcount).name,',1),%1.0f);\n'],nzx);
          else
            fprintf(fid,[indent,derivstr,' = zeros(%1.0f,1);\n'],nzx);
          end
          TF3 = ['cada',NDstr,'tf3'];
          fprintf(fid,[indent,TF3,' = ',Ystr,'~= 0;\n']);
          if isempty(regexp(Xstr,'\w$','once'))
            TF4 = ['cada',NDstr,'tf4'];
            fprintf([indent,TF4,' = ',Xstr,';\n']);
            Xstr = [TF4,'(',TF3,')'];
          else
            Xstr = [Xstr,'(',TF3,')']; %#ok<AGROW>
          end
          if isempty(regexp(Ystr,'\w$','once'))
            TF5 = ['cada',NDstr,'tf5'];
            fprintf([indent,TF5,' = ',Ystr,';\n']);
            Ystr = [TF5,'(',TF3,')'];
          else
            Ystr = [Ystr,'(',TF3,')']; %#ok<AGROW>
          end
          % Get [dZ/dX]*dX and print assignment to dZ
          fprintf(fid,[indent,derivstr,'(',TF3,') = ',getdzdx(Xstr,Ystr,DXstr,callerstr),';\n']);
        end
      else
        % Get [dZ/dX]*dX and print assignment to dZ
        fprintf(fid,[indent,derivstr,' = ',getdzdx(Xstr,Ystr,DXstr,callerstr),';\n']);
      end
    end
  elseif ~isempty(spdyind)
    %% -------------------------Y has Derivatives----------------------- %%
    derivstr = cadadername(funcstr,Vcount);
    z.deriv(Vcount).name   = derivstr;
    z.deriv(Vcount).nzlocs = spdyind;
    nzy = size(spdyind,1);
    if DPFLAG == 1
      if spdyflag;dyindin = 1:size(dyind,1);dyindin = dyindin(nzlocsy);end
      % -------------------Print DZ = DZ + [dZ/dY]*dY-------------------- %
      % Get dY
      if spdyflag
        Dind2 = cadaindprint(dyindin);
        if yvec
          DYstr = [y.deriv(Vcount).name,'(:,',Dind2,')'];
        else
          DYstr = [y.deriv(Vcount).name,'(',Dind2,')'];
        end
      else
        DYstr = y.deriv(Vcount).name;
      end
      % Get y and x that [dZ/dY] is function of
      switch callerstr
        case {'plus','minus'}
          Xstr = [];  Ystr = [];
        otherwise
          % Get x
          if xscalarflag || ...
              (xvec==2 && nzy == xNcol && isequal(spdyind(:,1),(1:xNcol).'))
            Xstr = x.func.name;
          elseif ~xvec && nzy == xMrow*xNcol && isequal(spdyind(:,1),(1:xMrow*xNcol).')
            Xstr = [x.func.name,'(:)'];
          elseif xvec ==1 && nzy == xMrow && isequal(spdyind(:,1),(1:xMrow).')
            Xstr = [x.func.name,'.'''];
          else
            TF1 = ['cada',NDstr,'tf1'];
            TFind1 = cadaindprint(spdyind(:,1));
            if xvec == 2
              fprintf(fid,[indent,TF1,' = ',x.func.name,'(:,',TFind1,');\n']);
              Xstr = TF1;
            elseif xvec == 1
              fprintf(fid,[indent,TF1,' = ',x.func.name,'(',TFind1,',:).'';\n']);
              Xstr = TF1;
            else
              fprintf(fid,[indent,TF1,' = ',x.func.name,'(',TFind1,');\n']);
              Xstr = [TF1,'(:)'];
            end
          end
          switch callerstr
            case 'times'
              Ystr = [];
            otherwise
              % Get y
              if yscalarflag || ...
                  (yvec==2 && nzy == yNcol && isequal(spdyind(:,1),(1:yNcol).'))
                Ystr = y.func.name;
              elseif ~yvec && nzy == yMrow*yNcol && isequal(spdyind(:,1),(1:yMrow*yNcol).')
                Ystr = [y.func.name,'(:)'];
              elseif yvec ==1 && nzy == yMrow && isequal(spdyind(:,1),(1:yMrow).')
                Ystr = [y.func.name,'.'''];
              else
                TF2 = ['cada',NDstr,'tf2'];
                TFind1 = cadaindprint(spdyind(:,1));
                if yvec == 2
                  fprintf(fid,[indent,TF2,' = ',y.func.name,'(:,',TFind1,');\n']);
                  Ystr = TF2;
                elseif yvec == 1
                  fprintf(fid,[indent,TF2,' = ',y.func.name,'(',TFind1,',:).'';\n']);
                  Ystr = TF2;
                else
                  fprintf(fid,[indent,TF2,' = ',y.func.name,'(',TFind1,');\n']);
                  Ystr = [TF2,'(:)'];
                end
              end
          end
      end
      % Get [dZ/dY]*dY and print assignment to dZ
      fprintf(fid,[indent,derivstr,' = ',getdzdy(Xstr,Ystr,DYstr,callerstr),';\n']);
    end
  end
end
%% --------------------------Function Printing ------------------------- %%
if PFLAG
  switch callerstr
    case 'plus'
      fprintf(fid,[indent,funcstr,' = ',x.func.name,' + ',y.func.name,';\n']);
    case 'minus'
      fprintf(fid,[indent,funcstr,' = ',x.func.name,' - ',y.func.name,';\n']);
    case 'times'
      fprintf(fid,[indent,funcstr,' = ',x.func.name,'.*',y.func.name,';\n']);
    case 'rdivide'
      fprintf(fid,[indent,funcstr,' = ',x.func.name,'./',y.func.name,';\n']);
    case 'power'
      fprintf(fid,[indent,funcstr,' = ',x.func.name,'.^',y.func.name,';\n']);
    otherwise
      fprintf(fid,[indent,funcstr,' = ',callerstr,'(',x.func.name,',',y.func.name,');\n']);
  end
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id z.id],1) = ADIGATOR.VARINFO.COUNT;
z = class(z,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
end
%% ~~~~~~~~~~~~~~~~~~~~~~~ Partial of Z wrt X ~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
function dz = getdzdx(x,y,dx,callerstr)

switch callerstr
  case {'plus','minus'}
    dz = dx;
  case 'times'
    dz = [y,'.*',dx];
  case 'rdivide'
    dz = [dx,'./',y];
  case 'power'
    dz = [y,'.*',x,'.^(',y,'-1).*',dx];
  case 'atan2'
    dz = [y,'./(',y,'.^2+',x,'.^2).*',dx];
  otherwise
    error(['no derivative rule defined for function: ',callerstr])
end
end
%% ~~~~~~~~~~~~~~~~~~~~~~~ Partial of Z wrt Y ~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
function dz = getdzdy(x,y,dy,callerstr)

switch callerstr
  case 'plus'
    dz = dy;
  case 'minus'
    dz = ['-',dy];
  case 'times'
    dz = [x,'.*',dy];
  case 'rdivide'
    dz = ['-',x,'./',y,'.^2.*',dy];
  case 'power'
    dz = ['log(',x,').*',x,'.^',y,'.*',dy];
  case 'atan2'
    dz = ['-',x,'./(',y,'.^2+',x,'.^2).*',dy];
  otherwise
    error(['no derivative rule defined for function: ',callerstr])
end
end