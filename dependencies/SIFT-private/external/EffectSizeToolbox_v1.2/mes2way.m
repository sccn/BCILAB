function [stats,varargout]=mes2way(X,group,esm,varargin)
% ** function [stats,varargout]=mes2way(X,group,esm,varargin)
% computes measures of effect size for samples in a factorial two-way
% design and places the results in output structure stats. A summary table
% including the effect size measures and two-way ANOVA results will be
% displayed in the command window. All input parameters except X, group and
% esm are optional and must be specified as parameter/value pairs in any
% order, e.g. as in
%      mes2way(X,g,'eta2','confLevel',.9);
%
% For information on assumptions of the underlying model and related
% information please see the notes below TABLE 1 further down as well as
% the documentation accompanying this code.
%
%                           INPUT ARGUMENTS
%                           ---------------
% stats=mes2way(X,group,esm)  computes the effect size measure(s) 
%   specified in input variable esm (see TABLE 1) from input array X. X
%   must be a single-column vector. NaNs or Infs will be eliminated in
%   'pairwise' fashion, that is, without concurrent elimination of other
%   values (which would be the case in 'listwise' fashion). group must be a
%   two-column array with the same number of rows as X. The numbers in the
%   columns of group correspond to the levels of the two factors. 
%   PLEASE NOTE: the actual values in group may be arbitrary, but rank
%   order, not the order in which they appear in group, determines the
%   assignment to factor levels. For example, if group is
%     3.1  1
%     3.1  2
%     1.5  1
%     1.5  2
%   the first two data points of X are assigned to the SECOND level of the
%   first factor (because 3.1 > 1.5). This is particularly relevant for all
%   calculations based on contrasts (see below). Input variable esm must be
%   a character array (e.g. 'eta2') or a cell array (e.g.
%   {'eta2','partialeta2'}); see TABLE 1 below.
% stats=mes2way(...,'fName',{'genotype','treat'})  assigns names to the two 
%   factors which will be displayed in the summary table of results; by 
%   default they will be given generic names 'fac1' and 'fac2'
% stats=mes2way(...,'isDep',[1 0])  assumes that the samples in X are 
%   dependent (repeated measures) with respect to the first factor. Other
%   legal input values are [0 1] (dependence with respect to the second
%   factor), [1 1] (completely within-subjects design) and [0 0]
%   (completely between-subjects design, the default). If there is
%   within-subjects dependence along one or both factors, data must be
%   balanced. Furthermore, it is assumed that input data X are sorted
%   according to subjects; that is, within each group data for subjects are
%   listed in the same order.
% stats=mes2way(...,'nBoot',n)  computes confidence intervals of the 
%   statistic in question via bootstrapping n times. n should be on the
%   order of thousands, otherwise bootstrapping will lead to inaccurate
%   results. By default or if the value of nBoot is set to zero or
%   nonsensical values (<0, infinity, nan) confidence intervals will be 
%   computed analytically (where possible).
% stats=mes2way(...,'confLevel',0.9)  computes 90 % confidence  
%   intervals of the statistic in question (95 % ci are the default; any 
%   value may be specified)
% stats=mes2way(...,'cWeight',c)  allows specification of contrast weights
%   for the computation of standardized contrasts (see TABLE 1). The
%   required size of input array c depends on the kind of comparison (see
%   e.g. Kline 2004 for definitions and the documentation accompanying this
%   code for concrete examples):
%   MAIN COMPARISON CONTRAST: c must be a single-column array (comparison
%     of levels of the first factor) or a single-row array (ditto for
%     second factor)
%   SIMPLE COMPARISON CONTRAST: c must match the analysis design, i.e. in a
%     2x3 analysis c must have two rows and three columns, but all elements
%     except the row or the column of interest must be NaN
%  INTERACTION CONTRAST: c must match the analysis design and it should
%     also be doubly centered, that is, row sums and column sums must all
%     be zero
%
% -------------------------------------------------------------------------
% TABLE 1: ANALYSES TO BE SPECIFIED IN INPUT VARIABLE esm
% -------------------------------------------------------------------------
% esm             QUANTITIY COMPUTED ((*) require input of contrast weights)
% 'psi'           unstandardized contrast (*)
% 'g_psi'         standardized contrast (*)  
% 'eta2'          eta squared
% 'partialeta2'   partial eta squared 
% 'omega2',       omega squared
% 'partialomega2' partial omega squared
%
% Notes:
% i. Parameters in TABLE 1 marked by (*) need contrast weights to be
%  computed. The other parameters will be computed for main and interaction
%  effects and contrasts (if specified).
% ii. 'g_psi' is one possible twoway equivalent of Hedges' g, namely,
%  contrast divided by the square root of the pooled within-conditions
%  variance. Its value is identical for dependent and independent data.
% iii. Fixed factors and interactions between the main factors are assumed.
% iv. Data may be unbalanced in a completely between subjects design but 
%  must be balanced in any other design.
% iv. If data is unbalanced it is assumed that the imbalance is due to
%  random loss of data, not that unequal cell size is part of the effect.
%  Accordingly, type III errors are computed and an unweighted means
%  analysis (using harmonic means) is performed
%
%                         OUTPUT ARGUMENTS
%                         ----------------
% The results of the computations are placed into fields of output
% structure stats with names corresponding to the analyses requested. For
% example, stats=effectsz_oneway(X,{'eta2','partialeta2'}) will result in
% fields
%   .eta2 
%   .partialeta2
% and so on. 
% PLEASE NOTE: if contrast weights were specified, all output arguments
% will be single- or two-column arrays holding results in the following
% (row) order:
% 1. main effect of factor 1
% 2. main effect of factor 2
% 3. interaction effect
% 4. contrast-related effect
% If a parameter is not defined for main/interaction effects or for
% contrasts, the corresponding entries will be NaN.
% Where applicable, confidence intervals will be placed in fields
%   .eta2Ci 
%   .partialeta2Ci
% and the computations underlying confidence intervals (either of 'none',
% 'bootstrap', 'approximate analytical', or 'exact analytical') will be
% placed in fields
%   .eta2CiType 
%   .partialeta2CiType
% Additional fields:
%   .n (sample sizes)
%   .confLevel (confidence level)
%   .contrast (containing fields .weight, the contrast weights and .type,
%   the type of comparison)
% stats will also have fields .isDep and .nBoot with the same values as the
% corresponding input arguments, thus providing potentially crucial post-
% analysis information. The second, optional output argument is the summary
% table of results.

% -------------------------------------------------------------------------
% Version 1.2, March 2012
% Code by Harald Hentschke (University of Tübingen) and 
% Maik Stüttgen (University of Bochum)
% For additional information see Hentschke and Stüttgen, 
% Eur J Neurosci 34:1887-1894, 2011
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                       STRUCTURE OF CODE
% The code is composed of two major parts: In PART I input arguments are
% checked and pre-processed, 'constants' set up, and all kinds of other
% preparatory and householding jobs accomplished. PART II contains the
% computations proper. Right at the beginning, summed squares, F and p
% values are computed in local function prepcomp, which also accomplishes
% assembly of bootstrapped data if requested. Effect sizes and confidence
% intervals are then computed and the results placed in structure array
% 'stats' as described above. The results are additionally collected in
% cell array 'table' which is displayed in the command window (and may also
% be retrieved as an output argument)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ------- PART I: CHECKS OF INPUT ARGS & OTHER PREPARATORY WORKS ----------
% -------------------------------------------------------------------------
% ----- default values & varargin -----
% standard values of optional input arguments
fName=[];
isDep=[0 0];
nBoot=0;
cWeight=[];
confLevel=.95;

% check variable number of input args:
% - even number (parameter AND value specified)?
% - convert to proper upper/lowercase spelling as above
% - check whether valid parameter was given
% - overwrite defaults by values, if specified
nvarg=numel(varargin);
v=who;
if nvarg
  if nvarg/2==round(nvarg/2)
    for g=1:2:nvarg
      % check which optional input parameter is given, ignoring case, and
      % impose spelling used internally
      ix=find(strcmpi(varargin{g},v));
      if ~isempty(ix)
        varargin{g}=v{ix};
      else
        error(['invalid optional input parameter (' varargin{g} ') was specified']);
      end
    end
    % finally, the assignment of values to parameters
    pvpmod(varargin);
  else
    error('optional input parameters must be specified as parameter/value pairs, e.g. as in ''par1'',1')
  end
end

% ----- a few 'constants':
alpha=1-confLevel;
% number of factors
nFactor=2;
% minimal number of bootstrapping runs below which a warning will be issued
% (and bootstrapping will be skipped)
minNBootstrap=1000;
% this is an internal switch dictating the way summed squares are computed
% (relevant only to unbalanced designs); currently valid options are
% 'unweighted' and 'Type III'. If set to 'unweighted', an 'unweighted means
% analysis' is performed, that is, grand mean and marginal means are means
% of the corresponding cells NOT weighted by the number of samples in them,
% and 'overall' cell size is the harmonic mean of all cell sizes, which is
% plugged into the traditional formulae for balanced designs. It seems not
% to be used a lot in the literature, so by default the 'Method 1/Type III'
% way according to Overall & Spiegel 1969 is pursued, which is also the
% default in Matlab.
ssType='Type III';

if nargin<3
  error('input arguments X, group and esm are mandatory');
end
% Variable list_analyses below is a list of all possible analyses:
% - column 1 holds the shortcut as strings
% - column 2 holds 1 for quantities REQUIRING a contrast
% - column 3 is a short description of the analyses
% This is the 'master list' of all analyses against which input argument
% esm will be checked (see) below
list_analysis={...
'psi'            1, 'unstandardized contrast';...
'g_psi',         1, 'standardized contrast';...
'eta2',          0, 'eta squared';...
'omega2',        0, 'omega squared';...
'partialeta2',   0, 'partial eta squared';...
'partialomega2', 0, 'partial omega squared';...
};

% in case esm is a char array convert it to a cell array because the code
% below relies on it
if ischar(esm)
  esm={esm};
end

% the tag column, extracted
listTag=cat(1,list_analysis{:,2});
% indices to currently requested analyses
listIx=ismember_bc(list_analysis(:,1),esm);
% unique tags of these analyses
uTag=unique_bc(listTag(listIx));
% catch a few silly errors:
% - typos or non-existent analysis type
if any(~ismember_bc(esm,list_analysis(:,1)))
  error('An illegal type of analysis was specified');
end

% check whether contrast weight(s) were specified, whether they are
% required, and set according flags (exact checks of contrast weights will
% be done further below). Combine all this information in a struct
contrast.weight=cWeight;
if isempty(contrast.weight)
  contrast.type='none';
  contrast.do=false;
  if any(uTag)
    error('part of the requested analyses require contrast weights as input args')
  end
else
  contrast.do=true;
  % kind of contrast unknown at this point (see further below)
  contrast.type='unspec';
end

% --- check input X and group
[nRowX nColX]=size(X);
[nRowG nColG]=size(group);
if nColX~=1 
  error('X must be a single column-array');
end
% - number of rows not matching
if nRowG~=nRowX
  error('input variables group and X must have the same number of rows');
end
if nColG~=nFactor
  error('input variable group must have two columns representing the levels of the two factors');
end
% - undefined elements in group
if any(~isfinite(group))
  error('input variable group contains NaNs or Infs');
end
% - nans and infs in X
badIx=find(~isfinite(X));
if ~isempty(badIx)
  warning('eliminating NaNs and Infs from input variable X');
  X(badIx)=[];
  group(badIx,:)=[];
  nRowX=nRowX-numel(badIx);
  nRowG=nRowG-numel(badIx);
end

% determine how many factors and levels there are, on the occasion
% replacing the values in variable group (i.e. the levels of the factors)
% by integers 1,2,... because function dummyvar used below requires this
for g=1:nFactor
  [factor(g).level,nada,intG]=unique_bc(group(:,g));
  group(:,g)=intG;
  factor(g).nLevel=numel(factor(g).level);
  % factor(g).levelName=cellstr([repmat('level',factor(g).nLevel,1) int2str((1:factor(g).nLevel)')]);
end
% assign labels to factors
if iscell(fName) && numel(fName)==2
  [factor.name]=deal(fName{:});
else
  % if fName is anything but the default empty matrix, warn
  if ~isempty(fName)
    warning('check names for factors - setting generic names')
  end
  [factor.name]=deal('fac1','fac2');
end
% number of individual groups (cells) (nonredundant combinations of levels of all
% factors)
nGroup=prod([factor.nLevel]);
% ****************************** NOTE: ************************************
% the structure/layout of all intermediate and final results hinges on the
% order of groups established below by unique, namely, the rank order of
% the (numeric) group labels!
% *************************************************************************
[uGroup,aIx,bIx]=unique_bc(group,'rows');
% any empty cell?
if nGroup~=size(uGroup,1)
  error('there is at least one empty group (cell), which is not allowed in factorial designs');
end
% if group is unsorted in the sense that samples from any group do not
% form a contiguous block it appears possible that the data are messed
% up, so better generate an error here and force the user to rethink (and
% sort) the data
isContiguousGroups=sum(any(diff(group),2))==nGroup-1;
if ~isContiguousGroups
  error([mfilename ' expects samples from each group to form contiguous blocks. Please sort your data accordingly']);
end
% *** for each group, determine indexes & count samples. Arrange the
% results in 2D cell array groupIx according to the analysis design: first
% factor = rows, second factor = columns ***
% (use output from unique, linear indexing and transposition to embed
% properly)
groupIx=cell([factor.nLevel])';
nSample=zeros([factor.nLevel])';
for gIx=1:nGroup
  tmpIx=find(bIx==gIx);
  groupIx{gIx}=tmpIx;
  nSample(gIx)=numel(tmpIx);
end
% don't forget to transpose
groupIx=groupIx';
nSample=nSample';

% in contrast to the 2D layout of groupIx and nSample most variables (i.e.
% summed squares, etc.) will have a 1D layout, generated by accessing
% elements of groupIx via linear indexing  (the 2nd dim is reserved for
% bootstrapped data). So, in order to compute e.g. marginal means we need
% linear indexes: 1D equivalents of 2D indexes into entire rows and columns
% in groupIx: with the code below, e.g. the first row of s2i2 contains
% indexes into e.g. r.meanGroup (computed in local function prepcomp),
% corresponding to the first row of groupIx, that is, all data of the first
% level of the first factor, and so on. Likewise, the first column will
% contain indexes to the first level of the second factor. (naming: s2i2 as
% a shortcut for 'output of sub2ind, placed in 2D array')
s2i2=reshape(1:nGroup,[factor.nLevel]);

% just to be 100% sure, this is an internal check that elements of groupIx
% are really integers incrementing in steps of 1
for g=1:nGroup
  tmp=unique_bc(diff(groupIx{g}));
  if numel(groupIx{g})>1 && ~(numel(tmp)==1 && tmp==1)
    error('internal: groupIx messed up - tell the programmer');
  end
end

% --- check bootstrapping settings
doBoot=false;
if isfinite(nBoot)
  if nBoot>=minNBootstrap;
    doBoot=true;
  else
    if nBoot~=0
      % warn only if nBoot small but different from zero because zero may
      % be a deliberate input value
      warning('number of bootstrap repetitions is not adequate - not bootstrapping');
    end
    nBoot=0;
  end
end

% --- check other input arguments
if confLevel<=0 || confLevel>=1
  error('input variable ''confLevel'' must be a scalar of value >0 and <1');
end

% deal with contrast weights
if contrast.do
  [cwN1 cwN2]=size(contrast.weight);
  % should there be Infs, convert to NaNs
  contrast.weight(isinf(contrast.weight))=nan;
  if cwN1==factor(1).nLevel && cwN2==1
    % single-column array: main comparison of first level
    contrast.type='main1';
  elseif cwN1==1 && cwN2==factor(2).nLevel
    % single-row array: main comparison of second level
    contrast.type='main2';
  elseif isequal([cwN1 cwN2],[factor.nLevel])
    % layout of contrast.weight mirrors analysis design, but is not a 1 by
    % cwN2 or cwN1 by 1 design
    cwIsFin=isfinite(contrast.weight);
    if all(all(cwIsFin))
      % all values are finite: interaction contrast
      contrast.type='interaction';
    elseif numel(find(any(cwIsFin,1)))==1 && numel(find(any(cwIsFin,2)))==factor(1).nLevel
      % only one column contains finite values: simple comparison of first
      % factor at specific level of second factor
      contrast.type='simple1';
    elseif numel(find(any(cwIsFin,2)))==1 && numel(find(any(cwIsFin,1)))==factor(2).nLevel
      % the inverse
      contrast.type='simple2';
    else
      error('check matrix of contrast weights - it contains NaNs or Infs in improper places');
    end
  else
    error('check layout of contrast weights');
  end
end
% final checks: values
if contrast.do
  if any(strcmp(contrast.type,{'simple1','simple2','main1','main2'}))
    if nansum(nansum(abs(contrast.weight)))-2>eps(2)
      warning('contrast weights for simple or main comparison are not a standard set: the standardized mean difference g_psi computed with this set does not represent the difference between the averages of two subsets of means')
    end
  elseif strcmp(contrast.type,{'interaction'})
    if ~sum(sum(abs(contrast.weight)))
      error('contrast weights are all zero');
    elseif any(abs(sum(contrast.weight,1))>eps(1)) || any(abs(sum(contrast.weight,2))>eps(1))
      warning('array of contrast weights must be doubly centered, otherwise the resulting contrast is NOT independent of the main effects');
    elseif sum(sum(abs(contrast.weight)))-4>eps(4)
      warning('the sum of the absolute values of contrast weights is unequal 4: the standardized mean difference g_psi computed with this set does not represent the difference between a pair of simple comparisons');
    end
  end
end
% check design (isDep) and sample sizes
isDep=logical(isDep(:)');
if ~isequal(size(isDep), [1 2])
  error('check input variable ''isDep'' - it must be a two-element array'),
end
% in repeated measures designs check whether sample sizes match
if any(isDep)
  if numel(unique_bc(nSample))~=1
    error('in a design with within-subjects factors all cell sizes must be equal');
  end
end

% -------------------------------------------------------------------------
% ------- PART II: THE WORKS ----------
% -------------------------------------------------------------------------

% create a few standard fields
stats.isDep=isDep;
stats.nBoot=nBoot;
stats.confLevel=confLevel;
stats.n=nSample;
stats.contrast.type=contrast.type;
stats.contrast.weight=contrast.weight;

% preparatory computations of SS, MS etc.
r=prepcomp(X,group,groupIx,factor,s2i2,nSample,isDep,doBoot,nBoot,contrast,ssType,alpha);

% if partialeta2 or partialomega2 or both are requested AND analytical CI
% shall be computed AND the data are independent, compute CI for partial
% eta2 here because those of partial omega2 can be derived from them (and
% computing them independently of each other would be a waste).
% Note on the side: in twoway designs, exact confidence intervals can be
% computed (via noncentral F) only for the partialled versions of eta2 and
% omega2. This is because the denominator in the formula for e.g. eta2 is
% SStot=SSerr+SS1+SS2+SS12, meaning that eta2 depends on several
% noncentralities. Thus, even if it is mentioned only in passing (or not at
% all) in papers, if analytical ci for what the authors term eta squared or
% omega squared are specified, in all likelihood these are the ci for the
% partialled MES
if any(ismember_bc({'partialeta2','partialomega2'},esm)) && ~doBoot && ~any(isDep)
  % exact analytical confidence intervals for partial eta2 (Smithson 2003, p.43 [5.6])
  % - main effects
  for fi=1:2
    tmp=ncpci(r.factor(fi).F,'F',[r.factor(fi).df,r.dfErr],'confLevel',confLevel);
    pEta2Ci(fi,1:2)=tmp./(tmp+r.factor(fi).df+r.dfErr+1);
  end
  % - interaction
  tmp=ncpci(r.factor12.F,'F',[r.factor12.df,r.dfErr],'confLevel',confLevel);
  pEta2Ci(3,:)=tmp./(tmp+r.factor12.df+r.dfErr+1);
  % - contrast
  if contrast.do
    tmp=ncpci(r.FPsi,'F',[1,r.dfErr],'confLevel',confLevel);
    pEta2Ci(4,:)=tmp./(tmp+1+r.dfErr+1);
  end
  % information on kind of CI
  pEta2CiType=repmat('exact analytical',3+contrast.do,1);
end

% now computations of effect size parameters 
nEs=numel(esm);
% a temporary container to hold es and ci for the table of results to be
% displayed; column order is [es ci_lo ci_up] en block for each es in the
% order in which they are computed; row order is factor1, factor2,
% factor1*factor2, contrast
esStore=repmat(nan,[3+contrast.do nEs*3]);

for tti=1:nEs
  curEs=esm{tti};
  es=repmat(nan,3+contrast.do,1);
  % If analytical confidence intervals for the statistic in question cannot
  % be computed and no bootstrapping is requested, ci must be set to nan.
  % Do this here for all statistics to avoid redundant assignments. ci will
  % be overwritten with meaningful values if i) there is a formula for
  % analytical confidence intervals implemented within the respective case
  % in the switch statement, or ii) confidence intervals based on
  % bootstrapping are computed right after the switch statement
  ci=repmat(nan,3+contrast.do,2);
  % a similar argument applies to variable ciType, which indicates on which
  % method the computation of ci is based (which may differ between main
  % and interaction effects and contrasts, hence we need more than one row
  if doBoot
    ciType=repmat('bootstrap',3+contrast.do,1);
  else
    ciType=repmat('none',3+contrast.do,1);
  end
  switch curEs
    case 'psi'
      % rows 1-3 are nans because they are reserved for main and
      % interaction effects; then the contrast
      es=cat(1,repmat(nan,3,nBoot+1),r.psi);
      if ~doBoot
        ci(4,:)=r.ciPsi;
        ciType=char(repmat('none',3,1),'exact analytical');
      end
      
    case 'g_psi'
      % rows 1-3 are nans because this mes is only defined for contrasts;
      % then the contrast divided by standardizer, the square root of the
      % pooled within-conditions variance, making this one possible
      % equivalent of Hedges' g
      es=cat(1,repmat(nan,3,nBoot+1),r.psi./sqrt(r.msErr));
      if ~doBoot
        % independent samples: exact confidence intervals 
        if ~any(isDep)
          % relation of ncp to g_psi: Kline formula 6.14, p. 177
          ci(4,:)=sqrt(r.denomSsPsi)*ncpci(r.tPsi,'t',r.df4nc_Psi,'confLevel',confLevel);
          ciType=char(repmat('none',3,1),'exact analytical');
        end
      end
      
    case 'eta2'
      % main effects
      es(1,1:nBoot+1)=r.factor(1).ss./r.ssTot;
      es(2,:)=r.factor(2).ss./r.ssTot;
      % interaction
      es(3,:)=r.factor12.ss./r.ssTot;
      % contrast
      if contrast.do
        es(4,:)=r.ssPsi./r.ssTot;
      end
      % in contrast to oneway analyses exact CI are not computable for
      % eta2
      
    case 'partialeta2'
      % main effects
      % (note the second term in the denominator, which is the
      % design-dependent error of the effect)
      es(1,1:nBoot+1)=r.factor(1).ss./(r.factor(1).ss+r.factor(1).ssEffErr);
      es(2,:)=r.factor(2).ss./(r.factor(2).ss+r.factor(2).ssEffErr);
      % interaction
      es(3,:)=r.factor12.ss./(r.factor12.ss+r.factor12.ssEffErr);
      % contrast
      if contrast.do
        % (note second term in denominator)
        es(4,:)=r.ssPsi./(r.ssPsi+r.ssEffErr_Psi);
      end
      % exact analytical ci apply only to completely between subjects
      % design
      if ~doBoot && ~any(isDep)
        ci=pEta2Ci;
        ciType=pEta2CiType;
      end
      
    case 'omega2'
      if ~any(isDep)
        % main effects
        es(1,1:nBoot+1)=(r.factor(1).ss-r.factor(1).df.*r.msErr)./(r.ssTot+r.msErr);
        es(2,:)=(r.factor(2).ss-r.factor(2).df.*r.msErr)./(r.ssTot+r.msErr);
        % interaction
        es(3,:)=(r.factor12.ss-r.factor12.df.*r.msErr)./(r.ssTot+r.msErr);
        if contrast.do
          % contrasts: df(effect)=1 and MS=SS
          es(4,:)=(r.ssPsi-1*r.msErr)./(r.ssTot+r.msErr);
        end
        % ci: see corresponding comments for eta2
      end
      
    case 'partialomega2'
      if ~any(isDep)
        % (simplified formula 6.31, p. 186, Kline 2004)
        tmp=prod([factor.nLevel])*r.hCellSz;
        % main effects
        es(1,1:nBoot+1)=r.factor(1).df*(r.factor(1).F-1)./(r.factor(1).df*(r.factor(1).F-1)+tmp);
        es(2,:)=r.factor(2).df*(r.factor(2).F-1)./(r.factor(2).df*(r.factor(2).F-1)+tmp);
        % interaction
        es(3,:)=r.factor12.df*(r.factor12.F-1)./(r.factor12.df*(r.factor12.F-1)+tmp);
        if contrast.do
          % contrast: df(effect)=1
          es(4,:)=1*(r.FPsi-1)./(1*(r.FPsi-1)+tmp);
        end
        if ~doBoot
          % exact analytical confidence intervals according to Fidler &
          % Thompson 2001, p. 593, based on those for partial eta2
          tmp=(r.ssTot*(1-pEta2Ci))/r.dfErr;
          % main effects, interaction and contrast all in one step as the
          % only free variable is df
          tmp2=cat(1,r.factor.df,r.factor12.df,find(contrast.do))*[1 1];
          ci=(r.ssTot*pEta2Ci-tmp2.*tmp)./(r.ssTot+tmp);
          ciType=pEta2CiType;
        end
      end
      
    otherwise
      error(['internal error: unrecognized test not caught by input checks: ' curEs]);
  end
  
  % *********************************************************************
  % If data were NOT bootstrapped, all computations are done at this
  % point and the results can be placed into appropriate fields of output
  % variable stats. If they were bootstrapped, confidence intervals must
  % be computed and the ES statistics extracted from the first element of
  % variable es.
  % *********************************************************************
  if doBoot
    % determine confidence intervals from array of effect size measures
    % generated from bootstrapped data
    ci=prctile(es(:,2:end)',[alpha/2  1-alpha/2]'*100)';
    % retain first column; this is the es computed from real data
    es=es(:,1);
  end
  % place es and ci in esStore
  esStore(:,(tti-1)*3+1:tti*3)=[es ci];
  
  % finally, use dynamic fields to store currently computed measures in
  % output variable stats
  stats.(curEs)=es;
  stats.([curEs 'Ci'])=ci;
  stats.([curEs 'CiType'])=ciType;
end

% assemble table for display and optional output:
% the concatenated es and ci as cell array
esStore=num2cell(esStore);
% 'fillers' of empty cells
filla=cell(1,nEs*3);
% insert ci titles
ciTi=cell(2,nEs);
[ciTi{1,:}]=deal('ci_lo');
[ciTi{2,:}]=deal('ci_up');
esm=cat(1,esm,ciTi);
esm=esm(:)';
% column order:
% source | SS | df | MS | F | p | effect sizes (in the order requested)
table={...
    'SOURCE',                'SS',              'df',           'MS',              'F',              'p',              esm{1:end}};
caseStr={'between-subjects','within-subjects'};
caseCondit=[false true];
% loop over between/within subjects cases
for caseIx=1:2
  st=caseStr{caseIx};
  if any(isDep==caseCondit(caseIx))
    table=cat(1,table,...
      {caseStr{caseIx},        '---',             '---',          '---',             '---',            '---',          filla{1:end}});
    for g=1:2
      if isDep(g)==caseCondit(caseIx)
        tmp=['- ' factor(g).name];
        table=cat(1,table,...
          {tmp,              r.factor(g).ss(1), r.factor(g).df, r.factor(g).ms(1), r.factor(g).F(1), r.factor(g).p(1), esStore{g,:}});
      end
    end
    if ~any(isDep) && caseIx==1 || (any(isDep) && caseIx==2) 
      tmp=['- ' factor(1).name '*' factor(2).name];
      table=cat(1,table,...
        {tmp,                r.factor12.ss(1),  r.factor12.df, r.factor12.ms(1), r.factor12.F(1), r.factor12.p(1), esStore{3,:}});
    end
  end
end
if contrast.do
  table=cat(1,table,...
    {'contrast',               '---',             '---',          '---',             '---',            '---',          filla{1:end}},...
    {['- ' contrast.type ],  r.ssPsi(1),        1,              r.msPsi(1),        r.FPsi(1),        r.pPsi(1),        esStore{4,:}});
end

table=cat(1,table,...
  {'within-cells',           r.ssErr(1),        r.dfErr,        r.msErr(1),        [],               [],        filla{1:end}});

switch sum(isDep)
  case 1
    % mixed within-subjects
    % - index to between-subjects factor
    bi=find(~isDep);
    % - same for within-subjects 
    wi=find(isDep);
    tmp1=['- subj within ' factor(bi).name];
    tmp2=['- ' factor(wi).name '*subj within ' factor(bi).name];    
    table=cat(1,table,...
      {tmp1,              r.ssSubjWithinNonRM(1), r.dfSubjWithinNonRM, r.msSubjWithinNonRM(1), [], [], filla{1:end}},...
      {tmp2,              r.ssSubjRM(1),          r.dfSubjRM,          r.msSubjRM(1),          [], [], filla{1:end}});
  
  case 2
    % completely within-subjects 
    table=cat(1,table,...
      {'- subjects',         r.ssSubj(1),         r.dfSubj,        r.msSubj(1), [], [], filla{1:end}});
    for g=1:2
      tmp=['- ' factor(g).name '*subj'];
      table=cat(1,table,...                      
        {tmp,              r.factor(g).ssSubj(1), r.factor(g).dfSubj, r.factor(g).msSubj(1), [], [], filla{1:end}});
    end
    tmp=['- ' factor(1).name '*' factor(2).name  '*subj'];
    table=cat(1,table,...
      {tmp,              r.factor12.ssSubj(1), r.factor12.dfSubj, r.factor12.msSubj(1), [], [], filla{1:end}});
end

table  

if nargout>1
  varargout{1}=table;
end
  
% ========================= LOCAL FUNCTIONS ===============================
% ========================= LOCAL FUNCTIONS ===============================

function r=prepcomp(x,group,groupIx,factor,s2i2,nSample,isDep,doBoot,nBoot,contrast,ssType,alpha)
% performs preparatory computations for 2-way analyses on single column
% array x with group assignment coded by input var group

% note on terminology: 'group' is equivalent to 'cell'

% number of factors and groups (cells)
nFactor=numel(factor);
nGroup=prod([factor.nLevel]);

% ------------------ START BY BOOTSTRAPPING (IF REQUESTED) ----------------
if doBoot
  % *********************************************************************
  % Here, individual, groupwise resampled instances of variable X are 
  % generated, expanding X in the second dimension: the first column
  % contains the original data, the others will be filled up with sampled
  % (with replacement) data.
  % *********************************************************************
  % variable bootSamX generated below shall contain indices to be used for
  % randomly sampling (with replacement) the original data (*** requires
  % groupIx to be integers incrementing in steps of 1! ***)
  switch sum(isDep)
    case 0
      % completely within-subjects: resample data within each group
      % independently of data points picked in other groups
      bootSamX=repmat(nan,size(x,1),nBoot);
      for g=1:nGroup
        bootSamX(groupIx{g},1:nBoot)=groupIx{g}(1)-1+ceil(nSample(g)*rand([nSample(g) nBoot]));
      end
    case 1
      % mixed within-subjects: the innermost loop applies a random sampling
      % index (offset-corrected) to all groups within a given column of
      % groupIx; this is done column by column of groupIx. This amounts to
      % a within-subjects design along the first factor; if the
      % within-subjects axis is along the second factor, we have to
      % temporarily transpose groupIx and nSample 
      if isequal(isDep,[false true])
        groupIx=groupIx';
        nSample=nSample';
      end
      bootSamX=repmat(nan,size(x,1),nBoot);
      for cIx=1:size(groupIx,2)
        % generic random sampling index for groups in current column
        partSamIx=ceil(nSample(1,cIx)*rand([nSample(1,cIx) nBoot]));
        for rIx=1:size(groupIx,1)
          bootSamX(groupIx{rIx,cIx},:)=partSamIx+groupIx{rIx,cIx}(1)-1;
        end          
      end
      % retranspose
      if isequal(isDep,[false true])
        groupIx=groupIx';
        nSample=nSample';
      end
      clear partSamIx;
    case 2
      % completely within-subjects: generate random sampling index for
      % first group and replicate for all other groups (which are identical
      % in size)
      bootSamX=repmat(ceil(nSample(1)*rand([nSample(1) nBoot])),[nGroup,1]);
      % add offset due to arrangement of data in a single-column array
      for g=1:nGroup
        bootSamX(groupIx{g},:)=bootSamX(groupIx{g},:)+groupIx{g}(1)-1;
      end
  end
  % from second column onwards fill with sampled data
  x(:,2:nBoot+1)=x(bootSamX);
  % delete resampling index 
  clear bootSam* 
end

%  --------------------- COMPUTATIONS PROPER ------------------------------
% within-groups SS, the individual groups' means and SS_error
r.ssErr=0;
for g=nGroup:-1:1
  % means of groups
  r.meanGroup(g,:)=sum(x(groupIx{g},:))/nSample(g);
  % within-groups SS
  r.ssErr=r.ssErr+sum((x(groupIx{g},:)-repmat(r.meanGroup(g,:),nSample(g),1)).^2);
end
% within-groups df
r.dfErr=sum(sum(nSample-1));
% within-groups MS
r.msErr=r.ssErr/r.dfErr;

% unweighted arithmetic mean of all cell means
r.meanGrand=mean(r.meanGroup);
% cell size: harmonic mean 
r.hCellSz=nGroup/sum(1./nSample(:));
% cell size: arithmetic mean 
r.aCellSz=mean(nSample(:));

% for both factors, compute df, marginal means and SS_a, or SS_effect, the
% between-groups SS (or effect sum of squares)
switch ssType

  case 'unweighted'
    % ** unweighted means method using simple formulae for balanced data
    % and harmonic mean of cell size to accomodate unbalanced data
    for fi=1:nFactor
      % permute if fi==2
      s2i2=permute(s2i2,mod((1:nFactor)+fi,nFactor)+1);
      r.factor(fi).df=factor(fi).nLevel-1;
      r.factor(fi).ss=0;
      for lix=factor(fi).nLevel:-1:1
        % unweighted marginal means
        r.factor(fi).margMean(lix,:)=mean(r.meanGroup(s2i2(lix,:),:));
        r.factor(fi).ss=r.factor(fi).ss+...
          factor(mod(fi,nFactor)+1).nLevel*r.hCellSz*(r.factor(fi).margMean(lix,:)-r.meanGrand).^2;
      end
      r.factor(fi).ms=r.factor(fi).ss/r.factor(fi).df;
      % re-permute if fi==2
      s2i2=permute(s2i2,mod((1:nFactor)+fi,nFactor)+1);
    end
    % interaction terms
    r.factor12.df=r.factor(1).df*r.factor(2).df;
    r.factor12.ss=0;
    for lix1=factor(1).nLevel:-1:1
      for lix2=factor(2).nLevel:-1:1
        r.factor12.ss=r.factor12.ss+...
          r.hCellSz*(r.meanGroup(s2i2(lix1,lix2),:)-...
          r.factor(1).margMean(lix1,:)-...
          r.factor(2).margMean(lix2,:)+...
          r.meanGrand).^2;
      end
    end
    r.factor12.ms=r.factor12.ss/r.factor12.df;
  
  case 'Type III'
    % ** Method 1/Type III SS
    % i. create design matrix and bring it in proper shape:
    dm=dummyvar(group);
    % temporary helper var used for indexing
    tmp=[0 factor.nLevel];
    for fi=1:nFactor
      % within each factor, subtract last level's column from other columns
      % in dm
      dm(:,(1:tmp(fi+1))+tmp(fi))=dm(:,(1:tmp(fi+1))+tmp(fi))-...
        repmat(dm(:,sum(tmp(fi:fi+1))),1,tmp(fi+1));
      % df of main factor
      r.factor(fi).df=factor(fi).nLevel-1;
      % marginal mean: permute groupIx and/or s2i2 if fi==2
      groupIx=permute(groupIx,mod((1:nFactor)+fi,nFactor)+1);
      s2i2=permute(s2i2,mod((1:nFactor)+fi,nFactor)+1);
      for lix=factor(fi).nLevel:-1:1
        % unweighted marginal mean (=mean of cell means across rows or
        % columns regardless of cell size)
        r.factor(fi).margMean(lix,:)=mean(r.meanGroup(s2i2(lix,:),:));
        % % marginal mean from raw scores - not recommended
        % r.factor(fi).margMean(lix,:)=mean(x(cat(1,groupIx{lix,:}),:));
      end
      % re-permute if fi==2
      groupIx=permute(groupIx,mod((1:nFactor)+fi,nFactor)+1);
      s2i2=permute(s2i2,mod((1:nFactor)+fi,nFactor)+1);
    end
    % within each factor, delete columns corresponding to last level in dm
    dm(:,cumsum(tmp(2:end)))=[];
    % now add interaction columns:
    ct=prod([r.factor.df]);
    for lix1=r.factor(1).df:-1:1
      for lix2=r.factor(2).df:-1:1
        dm(:,ct+sum([r.factor.df]))=dm(:,lix1).*dm(:,lix2+r.factor(1).df);
        ct=ct-1;
      end
    end
    % df of interaction
    r.factor12.df=r.factor(1).df*r.factor(2).df;
    
    % ii. compute regression-related SS
    % the vector of ones which must be added to obtain a constant term from
    % the regressions
    vo=ones(size(x,1),1);
    % regressions:
    % - a helper index, pointing to the column in the design matrix dm
    % corresponding to the last entry for 'main' effects
    tmpIx=sum([r.factor.df]);
    % - a,b,ab
    ss_abab=bbregress(x,[vo dm]);
    % - a,b
    ss_ab=bbregress(x,[vo dm(:,1:tmpIx)]);
    % - a,ab
    ss_aab=bbregress(x,[vo dm(:,[1:r.factor(1).df  tmpIx+1:end])]);
    % - b,ab
    ss_bab=bbregress(x,[vo dm(:,[r.factor(1).df+1:end])]);
    
    % iii. finally, SS and MS associated with main and interaction effects
    r.factor(1).ss=ss_abab-ss_bab;
    r.factor(2).ss=ss_abab-ss_aab;
    r.factor12.ss=ss_abab-ss_ab;
    r.factor(1).ms=r.factor(1).ss/r.factor(1).df;
    r.factor(2).ms=r.factor(2).ss/r.factor(2).df;
    r.factor12.ms=r.factor12.ss/r.factor12.df;
    clear dm ct vo tmp* ss*
end

% SS_t, the sum of all (balanced data)
r.ssTot=r.factor(1).ss+r.factor(2).ss+r.factor12.ss+r.ssErr;
% corresponding df, ms
r.dfTot=r.factor(1).df+r.factor(2).df+r.factor12.df+r.dfErr;
r.msTot=r.ssTot/r.dfTot;

% compute contrasts and associated SS
if contrast.do
  switch contrast.type
    case {'simple1','simple2'}
      % linear index to all non-nan elements in contrast.weight...
      mIx=find(isfinite(contrast.weight(:)));
      % ...to be used for calculation of contrast
      r.psi=sum(repmat(contrast.weight(mIx),1,nBoot+1).*r.meanGroup(mIx,:));
      % denominator for SS, which will later also be used for computation
      % of confidence interval
      r.denomSsPsi=sum(contrast.weight(mIx).^2./nSample(mIx));
      % SS
      r.ssPsi=r.psi.^2./r.denomSsPsi;
      
    case {'main1','main2'}
      % index to factor under consideration (note usage of last letter of
      % contrast.type)
      ix=sscanf(contrast.type(end),'%i');
      % index to other factor
      oix=mod(ix,2)+1;
      % instead of computing contrast from marginal means replicate
      % contrast weights according to design and compute cellwise, then
      % average
      cw=repmat(contrast.weight,circshift([1 factor(oix).nLevel],[0 ix-1]));
      r.psi=sum(repmat(cw(:),1,nBoot+1).*r.meanGroup)/factor(oix).nLevel;
      % denominator for SS
      r.denomSsPsi=sum(contrast.weight.^2./sum(nSample,oix));
      % SS
      r.ssPsi=r.psi.^2./r.denomSsPsi;
      
    case 'interaction'
      r.psi=sum(repmat(contrast.weight(:),1,nBoot+1).*r.meanGroup);
      % denominator for SS
      r.denomSsPsi=sum(contrast.weight(:).^2./nSample(:));
      % SS
      r.ssPsi=r.psi.^2./r.denomSsPsi;
  end
  % as df is always 1 MS=SS, but create field nonetheless
  r.msPsi=r.ssPsi;
end

% compute design-dependent SS, F, p for main and interaction effects as
% well as contrast
switch sum(isDep)
  case 0
    % *****************************
    % completely between-subjects:
    % *****************************
    % effect error SS, F, p for main and interaction effects
    for fi=1:2
      % the effect error SS - generally it depends on the design; it is
      % needed for the computation of partialeta2
      r.factor(fi).ssEffErr=r.ssErr;
      r.factor(fi).F=r.factor(fi).ms./r.msErr;
      r.factor(fi).p=1-fcdf(r.factor(fi).F,r.factor(fi).df,r.dfErr);
    end
    r.factor12.ssEffErr=r.ssErr;
    r.factor12.F=r.factor12.ms./r.msErr;
    r.factor12.p=1-fcdf(r.factor12.F,r.factor12.df,r.dfErr);
      
    % effect error SS, df, ci, F, t, p for contrast
    if contrast.do
      r.ssEffErr_Psi=r.ssErr;
      % the df needed for the computation of the t noncentrality parameter
      % for the contrast (NOT the df of the contrast, which is by
      % definition always 1)
      r.df4nc_Psi=r.dfErr;
      % confidence intervals (original data only)
      r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.msErr(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.dfErr);
      % F (original & bootstrapped, because bootstrapped Fs are needed for
      % the computation of CI of partial eta2 and omega2)
      r.FPsi=r.msPsi./r.msErr;
      % compute t from F and sign of contrast (original data only)
      r.tPsi=sign(r.psi(1)).*sqrt(r.FPsi(1));
      r.pPsi=1-fcdf(r.FPsi(1),1,r.dfErr(1));
    end
    
  case 1
    % *****************************
    % mixed within-subjects:
    % *****************************
    % identify the repeated-measures (RM) factor
    rmFacIx=find(isDep);
    nonRmFacIx=find(~isDep);
    % the code further below assumes the second factor to be
    % repeated-measures (i.e. its levels span the rows of groupIx). If the
    % first factor is RM, we have to temporarily transpose all variables
    % whose 2D layout reflects the factors and levels.
    if rmFacIx==1
      % first factor is repeated-measures: 
      % - marginal mean
      margMean=r.factor(2).margMean;
      % - transpose groupIx, nSample and s2i2
      groupIx=groupIx';
      nSample=nSample';
      s2i2=s2i2';
    else
      % second factor is repeated-measures:
      % - marginal mean
      margMean=r.factor(1).margMean;
    end
    % (note on the side: the quantities are computed in a 'raw' fashion
    % here, not by subtraction, as this may be necessary in a future
    % version of the code allowing for unbalanced data across the
    % non-repeated measures factor)
    nSubj=nSample(:,1);
    % df
    r.dfSubjWithinNonRM=sum(nSample(:,1)-1);
    r.dfSubjRM=r.dfSubjWithinNonRM*(factor(rmFacIx).nLevel-1);
    % initialize
    r.ssSubjWithinNonRM=0;
    r.ssSubjRM=0;
    % loop over levels of non-RM factor
    for rIx=factor(nonRmFacIx).nLevel:-1:1
      % loop over subjects
      for sIx=nSubj(rIx):-1:1
        % index to entries of current subject in current level of non-RM factor
        tmpIx=intersect_bc(sIx:nSubj(rIx):size(x,1),cat(1,groupIx{rIx,:}));
        % mean of those
        tmpMn=mean(x(tmpIx,:));
        % SS of subjects within levels of the non-RM factor (error between)
        r.ssSubjWithinNonRM=r.ssSubjWithinNonRM+factor(rmFacIx).nLevel*...
          (tmpMn-margMean(rIx,:)).^2;
        % loop over RM factor
        for cIx=factor(rmFacIx).nLevel:-1:1
          % interaction term: RM factor x subjects within levels of the
          % non-RM factor (error within)
          r.ssSubjRM=r.ssSubjRM+...
            (x(tmpIx(cIx),:)-tmpMn-r.meanGroup(s2i2(rIx,cIx),:)+margMean(rIx,:)).^2;
        end
      end
    end
    % MS
    r.msSubjWithinNonRM=r.ssSubjWithinNonRM/r.dfSubjWithinNonRM;
    r.msSubjRM=r.ssSubjRM/r.dfSubjRM;
    if rmFacIx==1
      % re-transpose variables transposed above
      groupIx=groupIx';
      nSample=nSample';
      s2i2=s2i2';
    end
    
    % effect error SS, F, t, p
    % - between-subjects factor
    r.factor(nonRmFacIx).ssEffErr=r.ssSubjWithinNonRM;
    r.factor(nonRmFacIx).F=r.factor(nonRmFacIx).ms./r.msSubjWithinNonRM;
    r.factor(nonRmFacIx).p=1-fcdf(r.factor(nonRmFacIx).F,r.factor(nonRmFacIx).df,r.dfSubjWithinNonRM);
    % - within-subjects factor
    r.factor(rmFacIx).ssEffErr=r.ssSubjRM;
    r.factor(rmFacIx).F=r.factor(rmFacIx).ms./r.msSubjRM;
    r.factor(rmFacIx).p=1-fcdf(r.factor(rmFacIx).F,r.factor(rmFacIx).df,r.dfSubjRM);
    % - interaction
    r.factor12.ssEffErr=r.ssSubjRM;    
    r.factor12.F=r.factor12.ms./r.msSubjRM;
    r.factor12.p=1-fcdf(r.factor12.F,r.factor12.df,r.dfSubjRM);
    % - § contrast: 
    if contrast.do
      % if it is an interaction contrast or a single factor-contrast across
      % the within-subjects factor...
      if strcmp(contrast.type,'interaction') || rmFacIx==sscanf(contrast.type(end),'%i')
        r.ssEffErr_Psi=r.ssSubjRM;
        r.df4nc_Psi=r.dfSubjRM;
        % confidence intervals (original data only)
        r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.msSubjRM(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.dfSubjRM);
        % F (original & bootstrapped)
        r.FPsi=r.msPsi./r.msSubjRM;
        % p (original data only)
        r.pPsi=1-fcdf(r.FPsi(1),1,r.dfSubjRM);
      else
        r.ssEffErr_Psi=r.ssSubjWithinNonRM;
        r.df4nc_Psi=r.dfSubjWithinNonRM;
        % confidence intervals (original data only)
        r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.msSubjWithinNonRM(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.dfSubjWithinNonRM);
        % F (original & bootstrapped)
        r.FPsi=r.msPsi./r.msSubjWithinNonRM;
        % p (original data only)
        r.pPsi=1-fcdf(r.FPsi(1),1,r.dfSubjWithinNonRM);
      end
      % compute t from F and sign of contrast
      r.tPsi=sign(r.psi).*sqrt(r.FPsi);
    end

  case 2
    % *****************************
    % completely within-subjects:
    % *****************************
    nSubj=nSample(1);
    r.dfSubj=nSubj-1;
    r.ssSubj=zeros(1,nBoot+1);
    for sIx=nSubj:-1:1
      % means of subjects
      r.meanSubj(sIx,:)=mean(x(sIx:nSubj:end,:));
      % SS_subj, between-subjects SS
      r.ssSubj=r.ssSubj+nGroup*(r.meanSubj(sIx,:)-r.meanGrand).^2;
    end
    r.msSubj=r.ssSubj/r.dfSubj;
    % loop over factors
    for fi=1:nFactor
      % permute if fi==2
      groupIx=permute(groupIx,mod((1:nFactor)+fi,nFactor)+1);
      % df
      r.factor(fi).dfSubj=r.dfSubj*r.factor(fi).df;
      % initialize
      tmpSs=0;
      % loop over levels
      for lix=factor(fi).nLevel:-1:1
        % loop over subjects
        for sIx=nSubj:-1:1
          % index to entries of current subject in current level of current factor
          tmpIx=intersect_bc(sIx:nSubj:size(x,1),cat(1,groupIx{lix,:}));
          % mean of those
          tmpMn=mean(x(tmpIx,:));
          % SS 
          tmpSs=tmpSs+factor(mod(fi,nFactor)+1).nLevel*(tmpMn-r.meanGrand).^2;
        end
      end
      % re-permute groupIx if fi==2
      groupIx=permute(groupIx,mod((1:nFactor)+fi,nFactor)+1);
      % obtain factor x subject term by subtraction
      r.factor(fi).ssSubj=tmpSs-r.factor(fi).ss-r.ssSubj;
      % ms
      r.factor(fi).msSubj=r.factor(fi).ssSubj/r.factor(fi).dfSubj;
    end
    % (interaction between factors) x subjects term:
    % - df
    r.factor12.dfSubj=r.dfSubj*r.factor12.df;
    % - obtain SS by subtraction
    r.factor12.ssSubj=r.ssErr-sum(cat(1,r.factor(:).ssSubj))-r.ssSubj;
    % - ms
    r.factor12.msSubj=r.factor12.ssSubj/r.factor12.dfSubj;
    % effect error SS, F, p
    for fi=1:2
      r.factor(fi).ssEffErr=r.factor(fi).ssSubj;
      r.factor(fi).F=r.factor(fi).ms./r.factor(fi).msSubj;
      r.factor(fi).p=1-fcdf(r.factor(fi).F,r.factor(fi).df,r.factor(fi).dfSubj);
    end
    r.factor12.ssEffErr=r.factor12.ssSubj;
    r.factor12.F=r.factor12.ms./r.factor12.msSubj;
    r.factor12.p=1-fcdf(r.factor12.F,r.factor12.df,r.factor12.dfSubj);
    
    % effect error SS, df, ci, F, t, p for contrast
    if contrast.do
      if strcmp(contrast.type,'interaction')
        r.ssEffErr_Psi=r.factor12.ssSubj;
        r.df4nc_Psi=r.factor12.dfSubj;
        % confidence intervals (original data only)
        r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.factor12.msSubj(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.factor12.dfSubj);
        % F (original & bootstrapped)
        r.FPsi=r.msPsi./r.factor12.msSubj;
        % p (original data only)
        r.pPsi=1-fcdf(r.FPsi(1),1,r.factor12.dfSubj);
      else
        fi=sscanf(contrast.type(end),'%i');
        r.ssEffErr_Psi=r.factor(fi).ssSubj;
        r.df4nc_Psi=r.factor(fi).dfSubj;
        % confidence intervals (original data only)
        r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.factor(fi).msSubj(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.factor(fi).dfSubj);
        % F (original & bootstrapped)
        r.FPsi=r.msPsi./r.factor(fi).msSubj;
        % p (original data only)
        r.pPsi=1-fcdf(r.FPsi(1),1,r.factor(fi).dfSubj);
      end
      % compute t from F and sign of contrast (original data only)
      r.tPsi=sign(r.psi(1)).*sqrt(r.FPsi(1));
    end
end

function ssReg=bbregress(dep,idep)
% 'bare bones' regression which performs only computations essential in the
% context of ANOVA computations cast as a linear regression, and returns
% the regression summed squares
% - coefficients
c=idep\dep;
% - expected values
e=idep*c;
% - regression SS
ssReg=sum((e-repmat(mean(dep),size(dep,1),1)).^2);