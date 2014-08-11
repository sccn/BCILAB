function adigatorfilenames = adigatorGenFiles4gpops2(setup,varargin)
% function adigatorfilenames = adigatorGenFiles4gpops2(setup)
%
% This function is for use with the GPOPS2 optimal control software. By
% inputting the setup (the input to the gpops2 function), this function
% will generate the necessary derivative files which the user can then
% supply to gpops2 in order to supply adigator derivatives. These
% derivative files must only be re-generated if the user changes their
% continuous or endpoint functions.
%
% ------------- Required Information in Input Setup Structure ----------- %
% This function does not require that the entire gpops2 setup structure be
% built, but rather only certain fields, the required fields are as
% follows:
%
% setup.endpfunction and setup.contfunction - this tells the function which
%   files to differentiate.
%
% setup.auxdata (if the endpoint or continuous function has auxilliary
%   inputs) - need these in order to generate the derivative files
%
% setup.bounds.phase.state.lower,
% setup.bounds.phase.control.lower,
% setup.bounds.phase.integral.lower (if any integrals), and
% setup.bounds.parameter.lower - This function uses the second dimension of
%   these in order to determine the input sizes for the endpoint and
%   continuous functions, thus arrays must be assigned to these for all
%   phases of the problem. The values assigned do not matter, but only the
%   dimension.
%
% ------------- Optional Information in Input Setup Structure ----------- %
%
% setup.derivatives.derivativelevel - This should be either 'first' or
%   'second', if 'second' this function will generate the Hessian files as
%   well as the Gradient files. If this is not specified, both Gradient and
%   Hessian files will be generated.
%
% ------------------- Generated Derivative File Names ------------------- %
% The generated derivative files will be named according to the user's file
% names with the string 'ADiGatorGrd' and 'ADiGatorHes' appended to the end
% for the first and second derivative files respectively. So, if the user's
% endpoint and continuous function files are named 'myEndp' and 'myCont',
% this function will generate the files 'myEndpADiGatorGrd',
% 'myContADiGatorGrd', 'myEndpADiGatorHes' and 'myContADiGatorHes'. These
% file names will also be given in the structure output where they will be
% assigned to the fields .EndpGrd, .ContGrd, .EndpHes, and .ContHes,
% respectively.

% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0


% ----------------------------------------------------------------------- %
%                       Extract Info from Setup                           %
% ----------------------------------------------------------------------- %
if isfield(setup,'derivatives');
  if isfield(setup.derivatives,'derivativelevel');
    if strcmpi(setup.derivatives.derivativelevel,'first');
      derivativelevel = 1;
    elseif strcmpi(setup.derivatives.derivativelevel,'second');
      derivativelevel = 2;
    else
      derivativelevel = 2;
    end
  else
    derivativelevel = 2;
  end
else
  derivativelevel = 2;
end

nPhases     = size(setup.bounds.phase,2);
numstate    = zeros(nPhases,1);
numcontrol  = zeros(nPhases,1);
numintegral = zeros(nPhases,1);
for p = 1:nPhases;
  numstate(p) = size(setup.bounds.phase(p).state.lower,2);
  if isfield(setup.bounds.phase,'control') && ~isempty(setup.bounds.phase(p).control.lower);
    numcontrol(p) = size(setup.bounds.phase(p).control.lower,2);
  end
  if isfield(setup.bounds.phase,'integral') && ~isempty(setup.bounds.phase(p).integral.lower);
    numintegral(p) = size(setup.bounds.phase(p).integral.lower,2);
  end
end
if isfield(setup.bounds,'parameter');
  numparameter = size(setup.bounds.parameter.lower,2);
else
  numparameter = 0;
end

nPV = 2*sum(numstate) + 2*nPhases + numparameter + sum(numintegral);
nCV = max(numstate) + max(numcontrol) + 1 + numparameter;
np  = numparameter;

if isfield(setup,'auxdata');
  pinput.auxdata = setup.auxdata;
  cinput.auxdata = setup.auxdata;
end

% ----------------------------------------------------------------------- %
%                            First Derivative                             %
% ----------------------------------------------------------------------- %
pinput.phase = struct('initialstate',[],'finalstate',[],'initialtime',[],'finaltime',[],'integral',[]);
cinput.phase = struct('state',[],'control',[],'parameter',[]);
pointCount = 0;
for p = 1:nPhases;
  nx  = numstate(p);
  nu  = numcontrol(p);
  ni  = numintegral(p);
  % ------------------------ Point Variables -------------------------- %
  % Point Variables are all assigned to one variable of differentiation,
  % we will call this 'v'.
  if nx > 0;
    % Initial State
    pinput.phase(p).initialstate = adigatorCreateDerivInput([1 nx],...
      struct('vodname','v','vodsize',[1 nPV],...
      'nzlocs',[(1:nx).' pointCount+(1:nx).']));
    pointCount = pointCount+nx;
    % Final State
    pinput.phase(p).finalstate = adigatorCreateDerivInput([1 nx],...
      struct('vodname','v','vodsize',[1 nPV],...
      'nzlocs',[(1:nx).' pointCount+(1:nx).']));
    pointCount = pointCount+nx;
  end
  % Initial Time
  pinput.phase(p).initialtime = adigatorCreateDerivInput([1 1],...
    struct('vodname','v','vodsize',[1 nPV],...
    'nzlocs',[1 pointCount+1]));
  pointCount = pointCount+1;
  % Final Time
  pinput.phase(p).finaltime = adigatorCreateDerivInput([1 1],...
    struct('vodname','v','vodsize',[1 nPV],...
    'nzlocs',[1 pointCount+1]));
  pointCount = pointCount+1;
  if ni > 0;
    % Integral
    pinput.phase(p).integral = adigatorCreateDerivInput([1 ni],...
      struct('vodname','v','vodsize',[1 nPV],...
      'nzlocs',[(1:ni).' pointCount+(1:ni).']));
    pointCount = pointCount+ni;
  end
  
  % --------------------- Continuous Variables ------------------------ %
  % For continuous variables we have a different variable of
  % differentiation on each phase, we will call these V
  if nx > 0;
    % State
    cinput.phase(p).state = adigatorCreateDerivInput([Inf nx],...
      struct('vodname','V','vodsize',[Inf nCV],...
      'nzlocs',[(1:nx).' (1:nx).']));
  end
  if nu > 0;
    % Control
    cinput.phase(p).control = adigatorCreateDerivInput([Inf nu],...
      struct('vodname','V','vodsize',[Inf nCV],...
      'nzlocs',[(1:nu).' nx+(1:nu).']));
  end
  % Time
  cinput.phase(p).time = adigatorCreateDerivInput([Inf 1],...
    struct('vodname','V','vodsize',[Inf nCV],...
    'nzlocs',[1 nx+nu+1]));
  if np > 0;
    % Continuous Parameter
    cinput.phase(p).parameter = adigatorCreateDerivInput([Inf np],...
      struct('vodname','V','vodsize',[Inf nCV],...
      'nzlocs',[(1:np).' nx+nu+1+(1:np).']));
  end
end

% Point Parameter
if np > 0;
  pinput.parameter = adigatorCreateDerivInput([1 np],...
    struct('vodname','v','vodsize',[1 nPV],...
    'nzlocs',[(1:np).' pointCount+(1:np).']));
end

% ------------------------ Generate Files ----------------------------- %
if nargin == 2
  ADopts = varargin{1};
  ADopts.comments  = 1;
  ADopts.overwrite = 1;
else
  ADopts = adigatorOptions('comments',1,'overwrite',1);
end
if isfield(setup,'displaylevel') && ~isempty(setup.displaylevel);
  ADopts.echo = logical(setup.displaylevel);
end
PointFunc = setup.functions.endpoint;
if ~ischar(PointFunc);
  PointFunc = func2str(PointFunc);
end
ContFunc = setup.functions.continuous;
if ~ischar(ContFunc);
  ContFunc = func2str(ContFunc);
end
PointDeriv1 = [PointFunc,'ADiGatorGrd'];
ContDeriv1  = [ContFunc,'ADiGatorGrd'];
adigatorfilenames.EndpGrd = PointDeriv1;
adigatorfilenames.ContGrd = ContDeriv1;
adigator(PointFunc,{pinput},PointDeriv1,ADopts);
adigator(ContFunc,{cinput},ContDeriv1,ADopts);

% ----------------------------------------------------------------------- %
%                            Second Derivative                            %
% ----------------------------------------------------------------------- %
if derivativelevel == 2;
  for p = 1:nPhases;
    nx  = numstate(p);
    nu  = numcontrol(p);
    ni  = numintegral(p);
    % ------------------------ Point Variables -------------------------- %
    if nx > 0;
      % Initial State
      pinput.phase(p).initialstate = struct(...
        'f',pinput.phase(p).initialstate,'dv',ones(nx,1));
      % Final State
      pinput.phase(p).finalstate = struct(...
        'f',pinput.phase(p).finalstate,'dv',ones(nx,1));
    end
    % Initial Time
    pinput.phase(p).initialtime = struct(...
      'f',pinput.phase(p).initialtime,'dv',1);
    % Final Time
    pinput.phase(p).finaltime = struct(...
      'f',pinput.phase(p).finaltime,'dv',1);
    if ni > 0;
      % Integral
      pinput.phase(p).integral = struct(...
        'f',pinput.phase(p).integral,'dv',ones(ni,1));
    end
    
    % -------------------- Continuous Variables ------------------------- %
    if nx > 0;
      % State
      cinput.phase(p).state = struct('f',cinput.phase(p).state,...
        'dV',adigatorCreateAuxInput([Inf,nx],ones(1,nx)));
    end
    if nu > 0;
      % Control
      cinput.phase(p).control = struct('f',cinput.phase(p).control,...
        'dV',adigatorCreateAuxInput([Inf,nu],ones(1,nu)));
    end
    % Time
    cinput.phase(p).time = struct('f',cinput.phase(p).time,...
      'dV',adigatorCreateAuxInput([Inf,1],1));
    if np > 0;
      % Continuous Parameter
      cinput.phase(p).parameter = struct('f',cinput.phase(p).parameter,...
        'dV',adigatorCreateAuxInput([Inf,np],ones(1,np)));
    end
  end
  if np > 0;
    % Point Parameter
    pinput.parameter = struct('f',pinput.parameter,'dv',ones(np,1));
  end
  
  % ------------------------ Generate Files ----------------------------- %
  if nargin == 2
    ADopts2 = varargin{1};
    ADopts2.comments  = 0;
    ADopts2.overwrite = 1;
  else
    ADopts2 = adigatorOptions('comments',0,'overwrite',1);
  end
  if isfield(setup,'displaylevel') && ~isempty(setup.displaylevel);
    ADopts2.echo = logical(setup.displaylevel);
  end
  PointDeriv2 = [PointFunc,'ADiGatorHes'];
  ContDeriv2  = [ContFunc,'ADiGatorHes'];
  adigatorfilenames.EndpHes = PointDeriv2;
  adigatorfilenames.ContHes = ContDeriv2;
  adigator(PointDeriv1,{pinput},PointDeriv2,ADopts2);
  adigator(ContDeriv1,{cinput},ContDeriv2,ADopts2);
end