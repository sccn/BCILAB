function y = cadaunarymath(x,zeroflag,callerstr)
% Unary math function - any unary mathematical functions should call this
% function in order to print the derivative. The derivative of each calling
% function must be defined in the getdydx sub-function.
%
% Copyright 2011-2014 Matthew J. Weinstein, Michael Patterson and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.EMPTYFLAG
  y = cadaEmptyEval(x);
  return
end

NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);
xMrow = x.func.size(1); xNcol = x.func.size(2);
if isinf(xMrow); xvec = 2; elseif isinf(xNcol); xvec =1; else xvec = 0;end

% ---------------------- Build Function Properties -----------------------%
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',[xMrow xNcol],'zerolocs',...
  [],'value',[]);
if isfield(x.func,'logicref')
  y.func.logicref = x.func.logicref;
end
if ~isempty(x.func.zerolocs) && zeroflag
  y.func.zerolocs = x.func.zerolocs;
end
if ~isempty(x.func.value)
  callerfun = str2func(callerstr);
  y.func.value = callerfun(x.func.value);
end

% ------------------Build Derivative Properties-----------------------%
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
    if ~isempty(x.deriv(Vcount).nzlocs)
        derivstr = cadadername(funcstr,Vcount);
        y.deriv(Vcount).name    = derivstr;
        y.deriv(Vcount).nzlocs  = x.deriv(Vcount).nzlocs;
        % -----------------------Derivative Printing----------------------%
        if DPFLAG
          if strcmp(callerstr,'uminus')
            fprintf(fid,[indent,derivstr,' = -',x.deriv(Vcount).name,';\n']);
          elseif strcmp(callerstr,'uplus')
            fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
          else
            if xMrow*xNcol == 1 || ...
                (xvec==2 && isequal(x.deriv(Vcount).nzlocs(:,1),(1:xNcol).'))
              Xstr = x.func.name;
            elseif xvec==1 && isequal(x.deriv(Vcount).nzlocs(:,1),(1:xMrow).')
              Xstr = [x.func.name,'.'''];
            elseif ~xvec && isequal(x.deriv(Vcount).nzlocs(:,1),(1:xMrow*xNcol).')
              Xstr = [x.func.name,'(:)'];
            else
              TF1 = ['cada',NDstr,'tf1'];
              TFind1 = cadaindprint(x.deriv(Vcount).nzlocs(:,1));
              if xvec==2
                fprintf(fid,[indent,TF1,' = ',x.func.name,'(:,',TFind1,');\n']);
                Xstr = TF1;
              elseif xvec==1
                fprintf(fid,[indent,TF1,' = ',x.func.name,'(',TFind1,',:).'';\n']);
                Xstr = TF1;
              else
                fprintf(fid,[indent,TF1,' = ',x.func.name,'(',TFind1,');\n']);
                Xstr = [TF1,'(:)'];
              end
            end
            dYdXstr = getdydx(Xstr,callerstr);
            fprintf(fid,[indent,derivstr,' = ',...
              dYdXstr,'.*',x.deriv(Vcount).name,';\n']);
          end
        end
    end
end


% --------------------------Function Printing ----------------------------%
if PFLAG == 1
    fprintf(fid,[indent,funcstr,' = ',callerstr,'(',x.func.name,');\n']);
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
end

function dydx = getdydx(x,callerstr)
switch callerstr
  case 'abs'
    dydx = ['sign(',x,')'];
  case 'exp'
    dydx = ['exp(',x,')'];
  case 'sqrt'
    dydx = ['(1/2)./sqrt(',x,')'];
  case 'log'
    dydx = ['1./',x];
  case 'log10'
    dydx = ['1./log(10)./',x];
  case 'sin'
    dydx = ['cos(',x,')'];
  case 'cos'
    dydx = ['-sin(',x,')'];
  case 'tan'
    dydx = ['sec(',x,').^2'];
  case 'cot'
    dydx = ['-csc(',x,').^2'];
  case 'sec'
    dydx = ['sec(',x,').*tan(',x,')'];
  case 'csc'
    dydx = ['-csc(',x,').*cot(',x,')'];
  case 'asin'
    dydx = ['1./sqrt(1-',x,'.^2)'];
  case 'acos'
    dydx = ['-1./sqrt(1-',x,'.^2)'];
  case 'atan'
    dydx = ['1./(1+',x,'.^2)'];
  case 'acot'
    dydx = ['-1./(1+',x,'.^2)'];
  case 'asec'
    dydx = ['1./',x,'./','sqrt(',x,'.^2-1)'];
  case 'acsc'
    dydx = ['-1./',x,'./','sqrt(',x,'.^2-1)'];
  case 'sind'
    dydx = ['180./pi.*cos(',x,')'];
  case 'cosd'
    dydx = ['-180./pi.*sin(',x,')'];
  case 'tand'
    dydx = ['180./pi.*sec(',x,').^2'];
  case 'cotd'
    dydx = ['-180./pi.*csc(',x,').^2'];
  case 'secd'
    dydx = ['180./pi.*sec(',x,').*tan(',x,')'];
  case 'cscd'
    dydx = ['-180./pi.*csc(',x,').*cot(',x,')'];
  case 'asind'
    dydx = ['180./pi./sqrt(1-',x,'.^2)'];
  case 'acosd'
    dydx = ['-180./pi./sqrt(1-',x,'.^2)'];
  case 'atand'
    dydx = ['180./pi./(1+',x,'.^2)'];
  case 'acotd'
    dydx = ['-180./pi./(1+',x,'.^2)'];
  case 'asecd'
    dydx = ['180./pi./',x,'./','sqrt(',x,'.^2-1)'];
  case 'acscd'
    dydx = ['-180./pi./',x,'./','sqrt(',x,'.^2-1)'];
  case 'sinh'
    dydx = ['cosh(',x,')'];
  case 'cosh'
    dydx = ['sinh(',x,')'];
  case 'tanh'
    dydx = ['sech(',x,').^2'];
  case 'coth'
    dydx = ['-csch(',x,').^2'];
  case 'sech'
    dydx = ['-sech(',x,').*tanh(',x,')'];
  case 'csch'
    dydx = ['-csch(',x,').*coth(',x,')'];
  case 'asinh'
    dydx = ['1./sqrt(',x,'.^2+1)'];
  case 'acosh'
    dydx = ['1./sqrt(',x,'.^2-1)'];
  case 'atanh'
    dydx = ['1./(1-',x,'.^2)'];
  case 'acoth'
    dydx = ['1./(1-',x,'.^2)'];
  case 'asech'
    dydx = ['-1./',x,'./sqrt(1-',x,'.^2)'];
  case 'acsch'
    dydx = ['-1./abs(',x,')./sqrt(1+',x,'.^2)'];
  case 'erf'
    dydx = ['2./sqrt(pi).*exp(-',x,'.^2)'];
  otherwise
    error(['No derivative rule defined for function: ',callerstr]);
end
end