function z = isequalwithequalnans(x,y,varargin)
% CADA overloaded version of function ISEQUALWITHEQUALNANS (arrays)

global ADIGATOR
if ADIGATOR.EMPTYFLAG
  z = cadaEmptyEval(x,y);
  return
end
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;

% ----------------------------Parse Inputs------------------------------- %
if PFLAG
  InputFuncStr = cell(1,nargin);
end

ValEqualFlag = 1;
% ValEqualFlag
%   0: definitely false
%   1: definitely true
%   2: we dont know
if isa(x,'cada')
  RowSize = x.func.size(1);
  ColSize = x.func.size(2);
  if ~isempty(x.func.value)
    Value = x.func.value;
  else
    ValEqualFlag = 2;
  end
  if PFLAG; InputFuncStr{1} = [x.func.name,',']; end
elseif isnumeric(x)
  RowSize = size(x,1);
  ColSize = size(x,2);
  Value = x;
  if PFLAG; InputFuncStr{1} = [cadamatprint(x),',']; end
else
  error('invalid input to ISEQUALWITHEQUALNANS');
end
for Icount = 2:nargin
  if Icount == 2
    a = y;
  else
    a = varargin{Icount-2};
  end
  if isa(a,'cada')
    aMrow = a.func.size(1);
    aNcol = a.func.size(2);
    if isempty(a.func.value)
      ValEqualFlag = 2;
    elseif ValEqualFlag == 1 && ~isequalwithequalnans(a.func.value,Value)
      ValEqualFlag = 0;
    end
    if PFLAG; InputFuncStr{Icount} = [a.func.name,',']; end
  elseif isnumeric(a)
    aMrow = size(a,1);
    aNcol = size(a,2);
    if ValEqualFlag == 1 && ~isequalwithequalnans(a,Value)
      ValEqualFlag = 0;
    end
    if PFLAG; InputFuncStr{Icount} = [cadamatprint(a),',']; end
  else
    error('invalid input to ISEQUALWITHEQUALNANS');
  end
  if (aMrow ~= RowSize || aNcol ~= ColSize)
    ValEqualFlag = 0;
  end
end

% ----------------------Build Function Properties--------------------------
z.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
z.func = struct('name',funcstr,'size',[1 1],'zerolocs',[],...
  'value',[],'logical',[]);

%----------------------Function Numeric Values (if exist)-----------------%
if ValEqualFlag < 2
  z.func.value = ValEqualFlag;
end

% ----------------------------Function Printing ------------------------- %
if PFLAG == 1
  InputFuncStr = cell2mat(InputFuncStr);
  InputFuncStr = ['(',InputFuncStr(1:end-1),')'];
  fprintf(fid,[indent,funcstr,' = isequalwithequalnans',InputFuncStr,';\n']);
end


% ------------------------Build Derivative Properties----------------------
z.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

ADIGATOR.VARINFO.LASTOCC([x.id y.id z.id],1) = ADIGATOR.VARINFO.COUNT;
z = class(z,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
