function fit=locfit(varargin)

% Smoothing noisy data using Local Regression and Likelihood.
%
% arguments still to add: dc maxit
%
%  Usage: fit = locfit(x,y)   % local regression fit of x and y.
%         fit = locfit(x)     % density estimation of x.
%
%  Smoothing with locfit is a two-step procedure. The locfit()
%  function evaluates the local regression smooth at a set of points
%  (can be specified through an evaluation structure). Then, use
%  the predict() function to interpolate this fit to other points.
%
%  Additional arguments to locfit() are specified as 'name',value pairs, e.g.:
%  locfit( x, 'alpha',[0.7,1.5] , 'family','rate' , 'ev','grid' , 'mg',100 ); 
%
%
%  Data-related inputs:
%
%    x is a vector or matrix of the independent (or predictor) variables.
%      Rows of x represent subjects, columns represent variables.
%      Generally, local regression would be used with 1-4 independent
%      variables. In higher dimensions, the curse-of-dimensionality,
%      as well as the difficulty of visualizing higher dimensional
%      surfaces, may limit usefulness.
%
%    y is the column vector of the dependent (or response) variable.
%      For density families, 'y' is omitted. 
% NOTE: x and y are the first two arguments. All other arguments require
%        the 'name',value notation.
%
%    'weights' Prior weights for observations (reciprocal of variance, or
%           sample size). 
%    'cens' Censoring indicators for hazard rate or censored regression.
%           The coding is '1' (or 'TRUE') for a censored observation, and
%           '0' (or 'FALSE') for uncensored observations. 
%    'base' Baseline parameter estimate. If a baseline is provided,
%           the local regression model is fitted as
%                        Y_i = b_i + m(x_i) + epsilon_i,
%           with Locfit estimating the m(x) term. For regression models,
%           this effectively subtracts b_i from Y_i. The advantage of the
%           'base' formulation is that it extends to likelihood
%           regression models. 
%    'scale' A scale to apply to each variable. This is especially
%           important for multivariate fitting, where variables may be
%           measured in non-comparable units. It is also used to specify
%           the frequency for variables with the 'a' (angular) style.
%     'sty' Character string (length d) of styles for each predictor variable.
%           n denotes `normal'; a denotes angular (or periodic); l and r
%           denotes one-sided left and right; c is conditionally parametric.
% 
%
%  Smoothing Parameters and Bandwidths:
%  The bandwidth (or more accurately, half-width) of the smoothing window
%  controls the amount of smoothing. Locfit allows specification of constant
%  (fixed), nearest neighbor, certain locally adaptive variable bandwidths,
%  and combinations of these. Also related to the smoothing parameter
%  are the local polynmial degree and weight function.
%
%    'nn' 'Nearest neighbor' smoothing parameter. Specifying 'nn',0.5
%         means that the width of each smoothing neighborhood is chosen
%         to cover 50% of the data.
%
%     'h' A constant (or fixed) bandwidth parameter. For example, 'h',2
%         means that the smoothing windows have constant half-width
%         (or radius) 2. Note that h is applied after scaling.
%
%   'pen' penalty parameter for adaptive smoothing. Needs to be used
%         with care.
%
%  'alpha' The old way of specifying smoothing parameters, as used in
%         my book. alpha is equivalent to the vector [nn,h,pen].
%         If multiple componenents are non-zero, the largest corresponding
%         bandwidth is used. The default (if none of alpha,nn,h,pen
%         are provided) is [0.7 0 0].
%
%   'deg' Degree of local polynomial. Default: 2 (local quadratic).
%         Degrees 0 to 3 are supported by almost all parts of the
%         Locfit code. Higher degrees may work in some cases. 
% 
%  'kern' Weight function, default = 'tcub'. Other choices are
%         'rect', 'trwt', 'tria', 'epan', 'bisq' and 'gauss'.
%         Choices may be restricted when derivatives are
%         required; e.g. for confidence bands and some bandwidth
%         selectors. 
% 
%    'kt' Kernel type, 'sph' (default); 'prod'. In multivariate
%         problems, 'prod' uses a simplified product model which
%         speeds up computations. 
% 
%  'acri' Criterion for adaptive bandwidth selection.
% 
%
%  Derivative Estimation.
%  Generally I recommend caution when using derivative estimation
%  (and especially higher order derivative estimation) -- can you
%  really estimate derivatives from noisy data? Any derivative
%  estimate is inherently more dependent on an assumed smoothness
%  (expressed through the bandwidth) than the data. Warnings aside...
%
%  'deriv' Derivative estimation. 'deriv',1 specifies the first derivative
%         (or more correctly, an estimate of the local slope is returned.
%         'deriv',[1 1] specifies the second derivative. For bivariate fits
%         'deriv',2 specifies the first partial derivative wrt x2.
%         'deriv',[1 2] is mixed second-order derivative.
% 
%  Fitting family.
%  'family' is used to specify the local likelihood family.
%         Regression-type families are 'gaussian', 'binomial',
%           'poisson', 'gamma' and 'geom'. If the family is preceded
%           by a q (e.g. 'qgauss', or 'qpois') then quasi-likelihood is
%           used; in particular, a dispersion estimate is computed.
%           Preceding by an 'r' makes an attempt at robust (outlier-resistant)
%           estimation. Combining q and r (e.g. 'family','qrpois') may
%           work, if you're lucky.
%         Density estimation-type families are 'dens', 'rate' and 'hazard'
%           (hazard or failure rate). Note that `dens' scales the output
%           to be a statistical density estimate (i.e. scaled to integrate
%           to 1). 'rate' estimates the rate or intensity function (events
%           per unit time, or events per unit area), which may be called
%           density in some fields.
%         The default family is 'qgauss' if a response (y argument) has been
%         provided, and 'dens' if no response is given.
%    'link' Link function for local likelihood fitting. Depending on the
%           family, choices may be 'ident', 'log', 'logit',
%           'inverse', 'sqrt' and 'arcsin'. 
% 
%  Evaluation structures.
%    By default, locfit chooses a set of points, depending on the data
%    and smoothing parameters, to evaluate at. This is controlled by
%    the evaluation structure.
%      'ev' Specify the evaluation structure. Default is 'tree'.
%           Other choices include 'phull' (triangulation), 'grid' (a grid
%           of points), 'data' (each data point), 'crossval' (data,
%           but use leave-one-out cross validation), 'none' (no evaluation
%           points, effectively producing the global parametric fit).
%           Alternatively, a vector/matrix of evaluation points may be
%           provided. 
%           (kd trees not currently supported in mlocfit)
%     'll' and 'ur' -- row vectors specifying the upper and lower limits
%           for the bounding box used by the evaluation structure.
%           They default to the data range. 
%     'mg' For the 'grid' evaluation structure, 'mg' specifies the
%           number of points on each margin. Default 10. Can be either a
%           single number or vector. 
%    'cut' Refinement parameter for adaptive partitions. Default 0.8;
%           smaller values result in more refined partitions. 
%    'maxk' Controls space assignment for evaluation structures. For the
%           adaptive evaluation structures, it is impossible to be sure
%           in advance how many vertices will be generated. If you get
%           warnings about `Insufficient vertex space', Locfit's default
%           assigment can be increased by increasing 'maxk'. The default
%           is 'maxk','100'. 
%
%    'xlim' For density estimation, Locfit allows the density to be
%           supported on a bounded interval (or rectangle, in more than
%           one dimension). The format should be [ll;ul] (ie, matrix with
%           two rows, d columns) where ll is the lower left corner of
%           the rectangle, and ur is the upper right corner.
%           One-sided bounds, such as [0,infty), are not supported, but can be
%           effectively specified by specifying a very large upper
%           bound. 
% 
%      'module' either 'name' or {'name','/path/to/module',parameters}.
% 
%  Density Estimation
%      'renorm',1  will attempt to renormalize the local likelihood
%           density estimate so that it integrates to 1. The llde
%           (specified by 'family','dens') is scaled to estimate the
%           density, but since the estimation is pointwise, there is
%           no guarantee that the resulting density integrates exactly
%           to 1. Renormalization attempts to achieve this.
%
%  The output of locfit() is a Matlab structure:
%
% fit.data.x (n*d)
% fit.data.y (n*1)
% fit.data.weights (n*1 or 1*1)
% fit.data.censor (n*1 or 1*1)
% fit.data.baseline (n*1 or 1*1)
% fit.data.style (string length d)
% fit.data.scales (1*d)
% fit.data.xlim (2*d)
%
% fit.evaluation_structure.type (string)
% fit.evaluation_structure.module.name (string)
% fit.evaluation_structure.module.directory (string)
% fit.evaluation_structure.module.parameters (string)
% fit.evaluation_structure.lower_left (numeric 1*d)
% fit.evaluation_structure.upper_right (numeric 1*d)
% fit.evaluation_structure.grid (numeric 1*d)
% fit.evaluation_structure.cut (numeric 1*d)
% fit.evaluation_structure.maxk
% fit.evaluation_structure.derivative
%
% fit.smoothing_parameters.alpha = (nn h pen) vector
% fit.smoothing_parameters.adaptive_criterion (string)
% fit.smoothing_parameters.degree (numeric)
% fit.smoothing_parameters.family (string)
% fit.smoothing_parameters.link (string)
% fit.smoothing_parameters.kernel (string)
% fit.smoothing_parameters.kernel_type (string)
% fit.smoothing_parameters.deren 
% fit.smoothing_parameters.deit
% fit.smoothing_parameters.demint
% fit.smoothing_parameters.debug
%
% fit.fit_points.evaluation_points (d*nv matrix)
% fit.fit_points.fitted_values (matrix, nv rows, many columns)
% fit.fit_points.evaluation_vectors.cell
% fit.fit_points.evaluation_vectors.splitvar
% fit.fit_points.evaluation_vectors.lo
% fit.fit_points.evaluation_vectors.hi
% fit.fit_points.fit_limits (d*2 matrix)
% fit.fit_points.family_link (numeric values)
% fit.fit_points.kappa (likelihood, degrees of freedom, etc)
%
% fit.parametric_component
%
%
%  The OLD format:
%
%    fit{1} = data.
%    fit{2} = evaluation structure.
%    fit{3} = smoothing parameter structure.
%    fit{4}{1} = fit points matrix.
%    fit{4}{2} = matrix of fitted values etc.
%           Note that these are not back-transformed, and may have the
%           parametric component removed.
%           (exact content varies according to module).
%    fit{4}{3} = various details of the evaluation points.
%    fit{4}{4} = fit limits.
%    fit{4}{5} = family,link.
%    fit{5} = parametric component values.
%



% Minimal input validation    
if nargin < 1
   error( 'At least one input argument required' );
end

xdata = double(varargin{1});
d = size(xdata,2);
n = size(xdata,1);
if ((nargin>1) && (~ischar(varargin{2})))
  ydata = double(varargin{2});
  if (any(size(ydata) ~= [n 1])); error('y must be n*1 column vector'); end;
  family = 'qgauss';
  na = 3;
else
  ydata = 0;
  family = 'density';
  na = 2;
end;
if mod(nargin-na,2)==0
  error( 'All arguments other than x, y must be name,value pairs' );
end


wdata = ones(n,1);
cdata = 0;
base  = 0;
style = 'n';
scale = 1;
xl = zeros(2,d);

alpha = [0 0 0];
deg = 2;
link = 'default';
acri = 'none';
kern = 'tcub';
kt = 'sph';
deren = 0;
deit  = 'default';
demint= 20;
debug = 0;

ev = 'tree';
ll = zeros(1,d);
ur = zeros(1,d);
mg = 10;
maxk = 100;
deriv=0;
cut = 0.8;
mdl = struct('name','std', 'directory','', 'parameters',0 );

while na < length(varargin)
    inc = 0;
    if (varargin{na}=='y')
        ydata = double(varargin{na+1});
        family = 'qgauss';
        inc = 2;
        if (any(size(ydata) ~= [n 1])); error('y must be n*1 column vector'); end;
    end
    if (strcmp(varargin{na},'weights'))
        wdata = double(varargin{na+1});
        inc = 2;
        if (any(size(wdata) ~= [n 1])); error('weights must be n*1 column vector'); end;
    end
    if (strcmp(varargin{na},'cens'))
        cdata = double(varargin{na+1});
        inc = 2;
        if (any(size(cdata) ~= [n 1])); error('cens must be n*1 column vector'); end;
    end
    if (strcmp(varargin{na},'base')) % numeric vector, n*1 or 1*1.
        base = double(varargin{na+1});
        if (length(base)==1); base = base*ones(n,1); end;
        inc = 2;
    end
    if (strcmp(varargin{na},'style')) % character string of length d.
        style = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'scale')) % row vector, length 1 or d.
        scale = varargin{na+1};
        if (scale==0)
          scale = zeros(1,d);
          for i=1:d
            scale(i) = sqrt(var(xdata(:,i)));
          end;
        end;
        inc = 2;
    end;
    if (strcmp(varargin{na},'xlim')) % 2*d numeric matrix.
        xl = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'alpha')) % row vector of length 1, 2 or 3.
        alpha = [varargin{na+1} 0 0 0];
        alpha = alpha(1:3);
        inc = 2;
    end
    if (strcmp(varargin{na},'nn')) % scalar
        alpha(1) = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'h')) % scalar
        alpha(2) = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'pen')) % scalar
        alpha(3) = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'acri')) % string
        acri = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'deg')) % positive integer.
        deg = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'family')) % character string.
        family = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'link')) % character string.
        link = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'kern')) % character string.
        kern = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'kt')) % character string.
        kt = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'ev')) % char. string, or matrix with d columns.
        ev = varargin{na+1};
        if (isnumeric(ev)); ev = ev'; end;
        inc = 2;
    end;
    if (strcmp(varargin{na},'ll')) % row vector of length d.
        ll = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'ur')) % row vector of length d.
        ur = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'mg')) % row vector of length d.
        mg = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'cut')) % positive scalar.
        cut = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'module')) % string.
        mdl = struct('name',varargin{na+1}, 'directory','', 'parameters',0 );
        inc = 2;
    end;
    if (strcmp(varargin{na},'maxk')) % positive integer.
        maxk = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'deriv')) % numeric row vector, up to deg elements.
        deriv = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'renorm')) % density renormalization.
        deren = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'itype')) % density - integration type.
        deit = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'mint')) % density - # of integration points.
        demint = varargin{na+1};
        inc = 2;
    end;
    if (strcmp(varargin{na},'debug')) % debug level.
        debug = varargin{na+1};
        inc = 2;
    end;
    if (inc==0)
      disp(varargin{na});
      error('Unknown Input Argument.');
    end;
    na=na+inc;
end


fit.data.x = xdata;
fit.data.y = ydata;
fit.data.weights = wdata;
fit.data.censor = cdata;
fit.data.baseline = base;
fit.data.style = style;
fit.data.scales = scale;
fit.data.xlim = xl;

fit.evaluation_structure.type = ev;
fit.evaluation_structure.module = mdl;
fit.evaluation_structure.lower_left = ll;
fit.evaluation_structure.upper_right = ur;
fit.evaluation_structure.grid = mg;
fit.evaluation_structure.cut = cut;
fit.evaluation_structure.maxk = maxk;
fit.evaluation_structure.derivative = deriv;

if (alpha==0); alpha = [0.7 0 0]; end;

fit.smoothing_parameters.alpha = alpha;
fit.smoothing_parameters.adaptive_criterion = acri;
fit.smoothing_parameters.degree = deg;
fit.smoothing_parameters.family = family;
fit.smoothing_parameters.link = link;
fit.smoothing_parameters.kernel = kern;
fit.smoothing_parameters.kernel_type = kt;
fit.smoothing_parameters.deren = deren;
fit.smoothing_parameters.deit = deit;
fit.smoothing_parameters.demint = demint;
fit.smoothing_parameters.debug = debug;

[fpc pcomp] = mexlf(fit.data,fit.evaluation_structure,fit.smoothing_parameters);
fit.fit_points = fpc;
fit.parametric_component = pcomp;

return



