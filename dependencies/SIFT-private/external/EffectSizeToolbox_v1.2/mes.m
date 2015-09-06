function stats=mes(X,Y,esm,varargin)
% ** function stats=mes(X,Y,esm,varargin)
% computes measures of effect size or related quantities between two
% samples (2-sample analyses) or one sample and a null value (1-sample
% analyses). The output consists of effect size measure(s), t statistics,
% and 95% confidence intervals where applicable. All input parameters
% except X, Y and esm are optional and must be specified as parameter/value
% pairs in any order, e.g. as in
%      mes(X,Y,'glassdelta','isDep',1,'nBoot',3000);
% 
%                           INPUT ARGUMENTS
%                           ---------------
% stats=mes(X,Y,esm)  computes the effect size measure(s) specified in 
%   input variable esm (see table at bottom) between samples X and Y or
%   between sample X and null value Y. X and Y may have multiple columns;
%   in this case the numbers of columns in X and Y must be equal as
%   comparisons will be made between matching columns. Furthermore, it is
%   assumed that each row corresponds to one subject/case (see treatment of
%   missing values below). No correction for multiple comparisons is
%   applied.
%   In 2-sample analyses it is assumed that the samples are independent
%   (not paired); if they are dependent (paired) use optional input
%   argument isDep (see below). Both the effect size measures available
%   and their computation are contingent on whether the samples are
%   independent or dependent. 
%   Input variable esm must be a character array (e.g. 'hedgesg') or a cell
%   array (e.g. {'hedgesg','glassdelta'}). 
% stats=mes(...,'isDep',1)  assumes that the samples in X and Y are
%   dependent (paired); accordingly, X and Y must have an equal number of 
%   rows. If this parameter is omitted the assumption is one of independent 
%   samples. 
% stats=mes(...,'missVal','listwise')  omits NaNs and Infs in X and Y, 
%   summarily termed here as missing values, in a within-variable listwise
%   fashion (which is the default): if there is a missing value anywhere in
%   X, the entire row will be omitted from analysis. If X has only one
%   column, only the single data point will be omitted. If the data are
%   unpaired the same applies to Y, independent of foul values in X. In
%   case of paired data deletion of rows is done in truly listwise fashion:
%   matching rows in X and Y will be omitted. If 'missVal' is set to
%   'pairwise' only the individual missing data points will be omitted.
%   Note that in case of paired data matching data points in X or Y will be
%   omitted.
% stats=mes(...,'nBoot',n)  computes confidence intervals of the statistic 
%   in question via bootstrapping n times. n should be on the order of
%   thousands, otherwise bootstrapping will lead to inaccurate results. By
%   default or if the value of nBoot is set to zero or nonsensical values
%   (<0, infinity, nan) confidence intervals will be computed analytically.
%   In cases where bootstrapping is not requested and computation of
%   analytical confidence intervals is not implemented the corresponding
%   fields of output struct stats will contain NaNs.
% stats=mes(...,'exactCi',true)  computes exact analytical confidence 
%   intervals for effect size measures for which both exact and approximate
%   CIs can be computed. Exact confidence intervals are based on iterative
%   determination of noncentrality parameters of noncentral Chi square, t
%   or F distributions. As this can be a very time-consuming process, and
%   as good bias corrections for small sample sizes are implemented for a
%   number of effect size measures, by default approximate analytical CIs
%   will be computed. See the documentation for details. Note that setting
%   this option is without effect if bootstrapping is requested (see input
%   variable 'nBoot' above).
% stats=mes(...,'confLevel',0.90)  computes 90 % confidence intervals (ci)
%   of the statistic in question (95 % ci are the default; any value may be
%   specified). If ci cannot be computed analytically and bootstrappig was
%   not requested, the corresponding fields of output struct stats (see
%   below) will contain NaNs.
% stats=mes(...,'ROCtBoot',1)  computes bootstrap confidence intervals for 
%   the area under the receiver-operating curve according to the 'bootstrap
%   t' method, which is more conservative than the 'bootstrap percentile'
%   method, the default (see the documentation for details). This option
%   will be ignored if bootstrapping is not requested
% stats=mes(...,'trCutoff',1.5)  sets 1.5 as the cutoff value expressed in 
%   standard deviations of the combined data set beyond the grand mean in
%   the computation of tail ratios. For positive values the right tail
%   ratio will be computed, for negative values the left tail ratio (see
%   documentation for further information). Default is 1.
% stats=mes(...,'trMeth','analytic')  determines that the tail ratios shall be computed
%   'analytically', assuming normal distributions, in the computation of
%   tail ratios. Default is 'count', which means that the ratios will be
%   determined by counting the actual data points beyond the cutoff
%   (relatively insensitive to deviations from normality)
% stats=mes(...,'doPlot',1)  will produce very simple plots of the results 
%   (one figure per effect size measure requested)
%
% -------------------------------------------------------------------------
% TABLE 1: ANALYSES TO BE SPECIFIED IN INPUT VARIABLE esm
% -------------------------------------------------------------------------
% esm           QUANTITIY COMPUTED
%        -- measures between one sample and a null value: --
% 'g1'          standardized difference (sample mean - comparison value)
% 'U3_1'        fraction of values below comparison value
%        -- measures between two samples: --
% 'md'          mean difference
% 'hedgesg'     Hedges' g (standardized mean difference)
% 'glassdelta'  Glass's delta (standardized mean difference)
% 'mdbysd'      mean difference divided by std of difference score
% 'requiv'      point-biserial correlation coefficient 
% 'cles'        common language effect size 
% 'U1'          Cohen's U1
% 'U3'          Cohen's U3
% 'auroc'       receiver-operating characteristic: area under curve
% 'tailratio'   tail ratios
% 'rbcorr'      rank-biserial correlation coefficient 
%
% 
%                         OUTPUT ARGUMENTS
%                         ----------------
%   The results of the computations are placed into fields of output
%   argument stats with names corresponding to the analyses requested. For 
%   example, stats=mes(X,Y,{'requiv','hedgesg'}) will result in 
%     stats.requiv 
%     stats.hedgesg
%   Where applicable, confidence intervals will be placed in fields 
%     .requivCi
%     .hedgesgCi
%   and the computations underlying confidence intervals (either of 'none',
%   'bootstrap', 'bootstrap t' (auroc only), 'approximate analytical', or
%   'exact analytical') will be placed in fields
%     .requivCiType
%     .hedgesgCiType
%   Additional fields:
%     .n (sample sizes)
%     .confLevel (confidence level)
%     .tstat (t test statistics) 
%   stats will also have fields .isDep and .nBoot with the same values as
%   the corresponding input arguments, thus providing potentially crucial
%   post-analysis information. 

% -------------------------------------------------------------------------
% Version 1.2, March 2012
% Code by Harald Hentschke (University of Tübingen) and 
% Maik Stüttgen (University of Bochum)
% For additional information see Hentschke and Stüttgen, 
% Eur J Neurosci 34:1887-1894, 2011
% -------------------------------------------------------------------------

% ----- default values & varargin -----
% standard values of optional input arguments (see description above)
isDep=false;
missVal='listwise';
nBoot=0;
ROCtBoot=false;
exactCi=false;
confLevel=.95;
trCutoff=1;
trMeth='count';
doPlot=false;

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
      ix=find(strcmp(lower(varargin{g}),lower(v)));
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
% minimal number of bootstrapping runs below which a warning will be issued
% (and bootstrapping will be skipped)
minNBootstrap=1000;

% -------------------------------------------------------------------------
% ------- PART I: CHECKS OF INPUT ARGS & OTHER PREPARATORY WORKS ----------
% -------------------------------------------------------------------------
% Variable list_analyses below is a list of all possible analyses:
% - column 1 holds the shortcut as strings
% - column 2 contains a numerical tag coding for the kinds of samples 
%   (1=one-sample, 2=two-sample)
% - column 3 is a short description of the analyses
% This is the 'master list' of all analyses against which input argument
% esm will be checked (see) below
list_analysis={...
  'g1',          1, 'standardized difference (sample mean - comparison value)';...
  'U3_1',        1, 'fraction of values below comparison value';...
  'md',          2, 'mean difference (unstandardized)';...
  'hedgesg',     2, 'Hedges'' g (standardized mean difference)';...
  'glassdelta',  2, 'Glass''s delta (standardized mean difference)';...
  'mdbysd',      2, 'mean difference divided by std of difference score';...
  'requiv',      2, 'pointbiserial correlation coefficient';...
  'cles',        2, 'common language effect size';...
  'auroc',       2, 'receiver-operating characteristic';...
  'U1',          2, 'Cohen''s U1';...
  'U3',          2, 'Cohen''s U3';...
  'tailratio',   2, 'tail ratios';...
  'rbcorr',      2, 'rank-biserial correlation coefficient'...
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
% unique tags of these analyses: must be a scalar of value 1 (1-sample
% tests) or 2 (2-sample tests)
uTag=unique_bc(listTag(listIx));
% catch a few silly errors:
% - typos or non-existent analysis type
if ~any(listIx)
  error('An illegal type of analysis was specified');
end
% - mixed one- and two-samples tests
if numel(uTag)>1
  error('A mix of one-sample and two-sample tests was requested')
end
% - programmer's error 
if isempty(intersect_bc(uTag,[1 2]))
  error('internal: uTag different from 1 or 2');
end
% illegal values for missVal
if isempty(find(strcmpi(missVal,{'listwise','pairwise'})))
  error('illegal value for input parameter ''missVal'' (choose ''listwise'' or ''pairwise)''');
end

% --- check input X and Y
[nRowX nColX]=size(X);
[nRowY nColY]=size(Y);
% reshape X and Y in a few selected scenarios
% - if X is a single row array
if nRowX==1 && nColX>1
  warning('input variable X is a single row array - reshaping');
  X=X(:);
  [nRowX nColX]=size(X);
end
% - if X is a single column array and Y a single row array, reshape Y as 
% well
if nColX==1 && nRowY==1 && nColY>1
  warning('input variable Y is a single row-array - reshaping')
  Y=Y(:);
  [nRowY nColY]=size(Y);
end
% (if Y is a single row array while X is not the situation is ambiguous, so
% we'd better let the code produce an error further below)

% scalar expansion: if X has several columns and Y is a scalar, expand it
if nColX>1 && nRowY==1 && nColY==1
  Y=repmat(Y,1,nColX);
  nColY=nColX;
end
% now perform strict checks
if nColX~=nColY
  error('input variables X and Y must have the same number of columns');
end
if isDep
  if nRowX~=nRowY
    error('for paired data input variables X and Y must have the same number of rows');
  end
end
if uTag==1 && nRowY>1
  error('1-sample analyses require input variable y to be a scalar (or a row array with as many columns as x)');
end
  
% deal with foul values:
% - first, convert infs to nans
X(~isfinite(X))=nan;
Y(~isfinite(Y))=nan;
[nanRowX,nanColX]=find(isnan(X));
[nanRowY,nanColY]=find(isnan(Y));
% - depending on how missing values shall be treated, set selected
% values or entire rows to NaN. Differentiate between 1-sample and 2-sample
% analyses
if ~isempty(nanRowX) || ~isempty(nanRowY)
  if uTag==1 
    % 1-sample case (easy)
    if ~isempty(nanRowY)
      error('in 1-sample analyses Y is not allowed to contain NaN or Inf');
    end
    if strcmp(lower(missVal),'listwise')
      X(nanRowX,:)=nan;
    end
  else
    % 2-sample case
    switch lower(missVal)
      case 'listwise'
        if isDep
          nanRowX=union_bc(nanRowX,nanRowY);
          nanRowY=nanRowX;
        end
        X(nanRowX,:)=nan;
        Y(nanRowY,:)=nan;
      case 'pairwise'
        if isDep
          badIx=union_bc(sub2ind([nRowX nColX],nanRowX,nanColX),sub2ind([nRowY nColY],nanRowY,nanColY));
          X(badIx)=nan;
          Y(badIx)=nan;
        end
        % nothing to do in case of unpaired data and pairwise elemination
    end
  end
end
% temp variables not needed anymore
clear badIx nanRowX nanColX nanRowY nanColY;

% --- test-specific checks
% - mostly, we must check for the 'isDep' option: 
% a) for some analyses, data in x and y are usually not considered 'paired'
% (e.g. ROC), so by issuing an error or warning the user is encouraged to
% rethink his/her analysis
% b) the exact procedure/results of bootstrapping and t-test depend on
% whether the data are paired or not, hence a correct specification may be
% essential
if any(ismember_bc(esm,'auroc'))
  if isDep
    error('''auroc'' (receiver-operating characteristic) is not implemented for dependent samples (see parameter ''isDep'')');
  end
end
if any(ismember_bc(esm,'glassdelta'))
  if isDep
    warning('''glassdelta'' does not make sense for dependent samples (see parameter ''isDep'')');
  end
end
if any(ismember_bc(esm,'mdbysd'))
  if ~isDep
    error('''mdbysd'' (mean difference divided by std of difference score) is not defined for independent samples (see parameter ''isDep'')');
  end
end

if strcmpi(trMeth,'analytic')
  doAnalyticalTails=true;
elseif strcmp(lower(trMeth),'count')
  doAnalyticalTails=false;
else
  error('illegal value specified for input parameter ''trMeth''');
end

% --- check bootstrapping settings
doBoot=false;
if isfinite(nBoot)
  if nBoot>=minNBootstrap;
    doBoot=true;
  else
    if nBoot~=0
      % warn only if nBoot small but different from zero because zero may
      % by a deliberate input value
      warning('number of bootstrap repetitions is not adequate - not bootstrapping');
    end
    nBoot=0;
  end
end
if ROCtBoot && ~doBoot
  warning('option ''ROCtBoot'' is only effective if bootstrapping is requested - check input parameter ''nBoot''');
end

% --- check other input arguments
if confLevel<=0 || confLevel>=1
  error('input variable ''confLevel'' must be a scalar of value >0 and <1');
end

% -------------------------------------------------------------------------
% ------- PART II: THE WORKS ----------
% -------------------------------------------------------------------------
% The outermost loop cycles through all columns of X (and matching columns
% of Y), performs bootstrapping (if requested) and then runs the tests

% preallocate results arrays:
% - ttest associated parameters:
allocationNanny=repmat(nan,[1 nColX]);
% value of test statistic 
stats.t.tstat=allocationNanny;
% p values (which will be kept only for original data)
stats.t.p=allocationNanny;
% degrees of freedom
stats.t.df=allocationNanny;
% create a few standard fields
stats.isDep=isDep;
stats.nBoot=nBoot;
stats.confLevel=confLevel;
stats.n=[nRowX; nRowY];
% reorder field names: tstats moves from first place to last
stats=orderfields(stats,[2:numel(fieldnames(stats)) 1]);

% ** effect size and confidence intervals cannot be preallocated here
% because they depend on the exact tests required (will be done at first
% run of loop)

% waitbar only if X (and Y) has more than one column
if nColX>1
  wbH=waitbar(0,'Computing...','name','Please wait');
end

% loop over columns of X (and Y)
for g=1:nColX
  % pick non-nan values in columns
  x=X(~isnan(X(:,g)),g);
  y=Y(~isnan(Y(:,g)),g);
  n1=size(x,1);
  n2=size(y,1);

  % -----------------------------------------------------------------------
  % ------------------ BOOTSTRAPPING (IF REQUESTED) -----------------------
  % -----------------------------------------------------------------------
  if doBoot
    % *********************************************************************
    % Here, x, representing individual columns of variable X, is expanded
    % in the second dimension: the first column contains the original data,
    % the others will be filled up with sampled (with replacement) data. y
    % will be expanded in the same (or similar) manner in case of 2-sample
    % tests; in case of 1-sample tests its (scalar) value will be
    % replicated to nBoot+1 columns 
    % *********************************************************************
    % indices to be used for randomly sampling (with replacement) the
    % original data
    bootSamX=ceil(n1*rand([n1 nBoot]));
    % if data are paired use same indices for y; if they are unpaired
    % compute additional one
    if isequal(uTag,2) && ~isDep
      bootSamY=ceil(n2*rand([n2 nBoot]));
    end
    % from second column onwards fill with sampled data
    x(:,2:nBoot+1)=x(bootSamX);
    if isequal(uTag,2)
      % do the same with y in case of 2-sample tests
      if isDep
        y(:,2:nBoot+1)=y(bootSamX);
      else
        y(:,2:nBoot+1)=y(bootSamY);
      end
    elseif isequal(uTag,1)
      % expand y in case of 1-sample tests 
      y(:,2:nBoot+1)=y;
    end
    % delete index and temporary variables
    clear bootSam* tmp*
  end
    
  % -----------------------------------------------------------------------
  % ----------------- ES COMPUTATIONS -------------------------------------
  % -----------------------------------------------------------------------
  % In this section, first some some preparatory computations will be
  % performed, notably terms needed for t statistics and several mes. Then,
  % the FOR loop below will sequentially work on all types of analysis as 
  % listed in variable esm.

  % preparatory computations: 
  % - means, variances, number of samples
  m1=mean(x,1);
  s1=var(x,0,1); % var and std are by default divided by n-1
  m2=mean(y,1);
  s2=var(y,0,1);
  % exand n1 and n2
  n1=repmat(n1,[1 nBoot+1]);
  n2=repmat(n2,[1 nBoot+1]);

  % - degrees of freedom, t statistics, sd of diff score & pooled var
  if isDep || isequal(uTag,1)
    % df in dependent case: n=n1-1=n2-1
    df=n1-1;
    % standard deviation of difference score
    if uTag==1
      % take into account that in 1-sample analyses if x is a matrix y is a
      % row array of identical values
      stdD=std(x-y(1));
    else
      stdD=std(x-y);      
    end
    % standard error
    seD=stdD./sqrt(n1);
    % t statistic
    tst=(m1-m2)./seD;
  else
    % independent case: n=[n1+n2-2]
    df=n1+n2-2;
    % pooled (within-groups) variance
    sP=((n1-1).*s1 + (n2-1).*s2)./(n1+n2-2);
    % standard error 
    seP=sqrt(sP.*((n1+n2)./(n1.*n2)));
    % t statistic
    tst=(m1-m2)./seP;
  end
  % p value
  p=2*(1-tcdf(abs(tst),df(1)));
  
  % in output struct stats.t keep only first element of t statistic 
  % computed from original data
  stats.t.tstat(g)=tst(1);
  % ditto for p and sd
  stats.t.p(g)=p(1);
  stats.t.df(g)=df(1);

  % - inverse cumulative t distribution:
  % n1, n2 and therefore df are identical for the original data and the
  % sampled (bootstrapped) data, so compute inverse cumulative t
  % distribution only from the former. However, as the computations below
  % rely on ciFac having dimensions [1 by size(x,2)] expand it here
  ciFac=repmat(-tinv(alpha/2,df(1)),[1 size(x,2)]);
  
  % ***********************************************************************
  %                  loop over all ES computations
  % ***********************************************************************
  % PROGRAMMER'S NOTES:
  % -> At this point in the code, we are within the outermost FOR loop 
  %    (mentioned above) which cycles through the columns of X, g being the 
  %    column index (see code line saying 'for g=1:nColX') 
  % -> Variable x represents data from an individual column of X. In case
  %    of 1-sample tests, y is the corresponding null value; in case of
  %    2-sample tests, it contains data from the corresponding column of Y.
  % -> If bootstrapping was requested, x has [nBoot+1] columns: the
  %    original data are in column 1 and the bootstrapped data are in 
  %    columns 2 and above. The same holds for y in case of 2-sample tests. 
  % -> In case of bootstrapping and 1-sample tests, x has [nBoot+1] columns
  %    as described above; y is a [1 by nBoot+1] array full of identical
  %    values (namely, the null value of the current column of interest in
  %    x): it has been expanded above to facilitate computations here.
  % -> To sum up, the columns in x and y always match and correspond to
  %    each other; thus any new code to be inserted below will work if it
  %    takes the potentially 2D layout of x and y, with the standard Matlab
  %    columnar organization of data, into account 
  % -> sample sizes, means and variances of x and y have been precomputed
  %    above (variables n, m, s, respectively) and need not (in fact should 
  %    not) be computed here
  % -> same holds for t statistics (variable tst, see case 'requiv')
  % -> if a formula for analytical confidence intervals is known, the code
  %    must be embedded (for each statistic of course) in 
  %       if ~doBoot
  %       end
  %    (see e.g. case 'hedgesg'). If they cannot be computed analytically,
  %    NaNs are inserted in the results variable ci. If bootstrapping was
  %    requested, the NaNs will be replaced by meaningful values
  % -> nRowX and nColX are the number or rows and columns, respectively, of 
  %    the original input variables. As the loop hops through the columns
  %    of input variable X, x is by definition either a single-column array
  %    or an array with nBoot columns - not with nColX columns! Moreover,
  %    as any NaNs will have been eliminated here from x, n1 is the number
  %    of rows in x, not nRowX. To make a long story short, n1 can be used
  %    for e.g. preallocation or the like, and make sure in your new code
  %    that you do not confuse nBoot and nColX. size(x,2) is the safest
  %    bet. Of course, the same applies to Y and y.
  
  nEs=numel(esm);
  for tti=1:nEs
    curEs=esm{tti};
    % If analytical confidence intervals for the statistic in question
    % cannot be computed and no bootstrapping is requested, ci must be set
    % to nan. Do this here for all statistics to avoid lots of redundant
    % code. ci will be overwritten with meaningful values if i) there is a
    % formula for analytical confidence intervals implemented within the
    % respective case in the switch statement, or ii) confidence intervals
    % based on bootstrapping are computed right after the switch statement
    ci=[nan; nan];
    % a similar argument applies to variable ciType, which indicates on
    % which method the computation of ci is based
    if g==1
      if doBoot
        ciType='bootstrap';
      else
        ciType='none';
      end
    end
    
    switch curEs
      case 'g1'
        es=(m1-y)./sqrt(s1);

      case 'U3_1'
        % number of values in x below (scalar value) y
        es=sum(x<repmat(y,[n1(1) 1]))./n1;
        % count values on threshold half
        es=es+.5*sum(x==repmat(y,[n1(1) 1]))./n1;
        
      case 'md'
        % mean difference (unstandardized), the most basic mes
        % imaginable...
        es=m1-m2;
        if ~doBoot
          if isDep
            % se=standard error of difference scores
            ci=cat(1,es-ciFac.*seD,es+ciFac.*seD);
          else
            % se=standard error of mean difference
            ci=cat(1,es-ciFac.*seP,es+ciFac.*seP);
          end
          if g==1
            % confidence intervals of mean differences computed from
            % central t distributions are not approximate, hence we term 
            % them exact here, although exactness does not, as in the cases
            % of other mes below, imply computation via noncentral
            % distributions
            ciType='exact analytical';
          end
        end
        
      case 'mdbysd'
        % mean difference divided by std of difference score (defined only
        % for dependent data)
        es=(m1-m2)./std(x-y);
        if ~doBoot
          % exact ci (Smithson 2003, p.36)
          ci=ncpci(tst,'t',n1-1,'confLevel',confLevel)'/sqrt(n1);
          if g==1
            ciType='exact analytical';
          end
        end
              
      case 'hedgesg'
        % Hedges' g
        if isDep
          % n=n1=n2
          es=tst.*sqrt(2*stdD.^2./(n1.*(s1+s2)));
        else
          es=(m1-m2)./sqrt(sP);
        end
        % correct for bias due to small n (both dependent and independent
        % data, Kline 2004 (p. 102, 106); Nakagawa & Cuthill 2007)
        biasFac=(1-(3./(4*n1+4*n2-9)));
        es=es.*biasFac;
        if ~doBoot
          if isDep
            % approximate ci for paired data (Nakagawa & Cuthill 2007) -
            % note that n=n1=n2 and that correlation coeff r is needed; no
            % bias correction here
            se=sqrt((2-2*corr(x,y))./n1 + es.^2./(2*n1-1));
            ci=cat(1,es-ciFac.*se,es+ciFac.*se);
            if g==1
              ciType='approximate analytical';
            end
          else
            if exactCi
              % exact ci (Smithson 2003, p. 37), including bias correction
              ci=biasFac*ncpci(tst,'t',n1+n2-2,'confLevel',confLevel)'*sqrt((n1+n2)/(n1*n2));
              if g==1
                ciType='exact analytical';
              end
            else
              % approximate ci (Nakagawa & Cuthill 2007, eq. 17 in table
              % 3), including bias correction
              se=sqrt((n1+n2)./(n1.*n2) + (es.^2./(2*n1+2*n2-4)));
              ci=biasFac*cat(1,es-ciFac.*se,es+ciFac.*se);
              if g==1
                ciType='approximate analytical';
              end
            end
          end
        end
        
      case 'glassdelta'
        es=(m1-m2)./sqrt(s1);
        % analytical confidence intervals: only approximate
        if ~doBoot
          se=sqrt(es^2/(2*n2-2)+(n1+n2)/(n1*n2));
          ci=cat(1,es-ciFac.*se,es+ciFac.*se);
          if g==1
            ciType='aproximate analytical';
          end
        end
        
      case 'requiv'
        % note: an alternative computation would work as follows: assign
        % discrete numbers to the two different groups (say, 0 to x and 1
        % to y), concatenate the group indices and the real data and then
        % plug the resulting two variables into corr (Pearson's
        % correlation). This yields identical results, but is
        % computationally more demanding and also more complicated when x
        % and y are 2D (bootstrapped). 
        es=tst./sqrt(tst.^2 + n1+n2-2);
        if ~doBoot
          % if data are independent & exact analytical ci requested...
          if ~isDep && exactCi
            % exact analytical confidence intervals for partial eta2, of
            % which requiv is a special case (Smithson 2003, p.43 [5.6])
            ci=ncpci(tst^2,'F',[1 stats.t.df],'confLevel',confLevel)';
            % don't forget to sqrt and to take care of sign
            ci=sqrt(ci./(ci+1+stats.t.df+1));
            if es<0
              ci=-1*fliplr(ci);
            end
            if g==1
              ciType='exact analytical';
            end
          else
            % all other cases: approximate CI via Z-transform
            tmp=.5*log((1+es)/(1-es));
            tmp=tmp+[-1; 1]*ciFac./sqrt(n1+n2-3);
            % transform back
            ci=(exp(2*tmp)-1)./(exp(2*tmp)+1);
            if g==1
              ciType='approximate analytical';
            end
          end
        end

      case 'cles'
        % 'For continuous data, it is the probability that a score sampled
        % at random from one distribution will be greater than a score
        % sampled from some other distribution'
        es=normcdf((m1-m2)./sqrt(s1+s2));
        
      case 'U1'
        % line arrays of minimal and maximal values in x
        maxX=max(x);
        minX=min(x);        
        % ditto for y
        maxY=max(y);
        minY=min(y);        
        es=(sum(x>repmat(maxY,n1(1),1))+sum(y>repmat(maxX,n2(1),1))...
          +sum(x<repmat(minY,n1(1),1))+sum(y<repmat(minX,n2(1),1)))./(n1+n2);
        
      case 'U3'
        % we need the medians of both groups
        med1=median(x);
        med2=median(y);
        es=sum(x<repmat(med2,n1(1),1))./n1;
        % count identical values half
        es=es+.5*sum(x==repmat(med2,n1(1),1))./n1;
        % in case both medians are equal, U3 must by definition be 0.5, but
        % the code lines above may have resulted in divergent values if
        % the lower group (x) contains heavily aliased data (e.g. histogram
        % data with one dominant bin). The following two lines correct for
        % this
        tmpIx=med1==med2;
        es(tmpIx)=.5;
        
      case 'auroc'
        % compute area under curve using approach in Bamber 1975 (J Math
        % Psych 12:387-415), eq. 3, also presented in Hanley & McNeil
        % Radiology 143:29-36, 1982 (p. 31)
        es=zeros(1,nBoot+1);
        % loop over elements (rows) in x
        for xIx=1:n1
          tmp=repmat(x(xIx,:),n2(1),1);
          es=es+sum(tmp>y)+0.5*sum(tmp==y);
        end
        es=es./(n1.*n2);
        if ~doBoot
          % analytical confidence intervals: start by computing standard
          % error of AUROC according to Hanley & McNeil (Radiology
          % 143:29-36, 1982), who refer to the classic work by Bamber 1975
          % (J Math Psych 12:387-415)
          se=se_auroc(es,n1,n2);
          % note that normality is assumed in this step
          tmp=norminv(1-alpha/2).*se;
          ci=cat(1,es-tmp,es+tmp);
          if g==1
            ciType='approximate analytical';
          end
        else
          if ROCtBoot
            % this portion of code performs preparatory steps for the
            % 'bootstrap t' method of estimating CI (Obuchowski & Lieber
            % Acad Radiol 5:561-571, 1998), using the approach of Hanley &
            % McNeil (Radiology 143:29-36, 1982) for the estimation of each
            % bootstrapped data set's standard error
            % - compute standard error of original and bootstrapped sets
            se=se_auroc(es,n1,n2);
            % - overwrite the bootstrapped auroc values with the
            % studentized pivot statistic
            es(2:end)=(es(2:end)-es(1))./se(2:end);
            % - replace Infs resulting from bootstrapped cases with zero
            % se by NaNs because Matlab function prctile will eliminate
            % these
            es(isinf(es))=nan;
            % - finally, multiply by se of original data and add auroc of
            % original data. This way, the CI can be directly extracted
            % from es(2:end), as is done with all other effect size
            % measures, too
            es(2:end)=es(2:end)*se(1)+es(1);
            if g==1
              ciType='bootstrap t';
            end
          else
            % nothing to do here: CI will be extracted below, which
            % corresponds to the standard 'percentile bootstrap' method
          end
        end
        clear tmp;
        
      case 'tailratio'
        % total mean
        m_t=(n1.*m1 + n2.*m2)./(n1+n2);
        % total std (compute in two steps)
        std_t=...
          n1.*(m1-m_t).^2 + n2.*(m2-m_t).^2 + (n1-1).*s1 + (n2-1).*s2;
        std_t=sqrt(std_t./(n1+n2-1));
        % (that one would work as well: std(cat(1,x,y)))
        % the cutoff (threshold)
        cutoff=m_t+trCutoff*std_t;
        % now compute the proportions of samples...
        if doAnalyticalTails
          % ...the 'analytical' way (as in Kline p.126) - this has the
          % advantage that in cases where two samples are completely
          % separated we don't get Infinity as result. However, this really
          % requires normality
          if trCutoff>=0
            es=(1-normcdf((cutoff-m1)./sqrt(s1)))./(1-normcdf((cutoff-m2)./sqrt(s2)));
          else
            es=(normcdf((m1-cutoff)./sqrt(s1)))./(normcdf((m2-cutoff)./sqrt(s2)));
          end
        else
          % ...by counting
          if trCutoff>=0
            % right tail ratio
            es=(sum(x>ones(n1(1),1)*cutoff)./n1)./(sum(y>ones(n2(1),1)*cutoff)./n2);
          else
            % left tail ratio
            es=(sum(x<ones(n1(1),1)*cutoff)./n1)./(sum(y<ones(n2(1),1)*cutoff)./n2);
          end
        end
        % note that there may be division by zero, particularly in
        % bootstrapped data, which prevents proper computation of
        % confidence intervals, so replace the infs by nans
        es([false isinf(es(2:end))])=nan;
        
      case 'rbcorr'
        % rank-biserial correlation coefficient
        % § compared to all other analyses this one takes an awfully long
        % time because of function corr and the functions it calls
        % (tiedrank and so on) - possibly this could be sped up
        es=corr(cat(1,x,y),cat(1,zeros(n1(1),1),ones(n2(1),1)),'type','Spearman');
                
      otherwise
        error(['internal: unrecognized test not caught by input checks: ' curEs]);
    end
    
    % *********************************************************************
    % If data were NOT bootstrapped, all computations are done at this
    % point and the results can be placed into appropriate fields of output
    % variable stats. If they were bootstrapped, confidence intervals must
    % be computed and the ES statistics extracted from the first element of
    % variable es.
    % *********************************************************************

    % start by preallocating major results fields of output variable stats
    if g==1
      stats.(curEs)=allocationNanny;
      stats.([curEs 'Ci'])=[allocationNanny; allocationNanny];
    end
    
    if doBoot
      % determine confidence intervals from array of effect size measures
      % generated from bootstrapped data
      ci=prctile(es(2:end),[alpha/2  1-alpha/2]'*100);
      % retain first element; this is the es computed from real data
      es=es(1);
    end
    
    % finally, use dynamic fields to store currently computed measures in
    % output variable stats
    stats.(curEs)(g)=es;
    stats.([curEs 'Ci'])(:,g)=ci;
    if g==1
      stats.([curEs 'CiType'])=ciType;
    end
  end
  
  % update waitbar
  if nColX>1
    waitbar(g/nColX,wbH);
  end
end

% kill waitbar window
if nColX>1
  delete(wbH);
end

% call simple plot routine if requested
if doPlot
  simpleEsPlot(stats,esm);
end


% ========================= LOCAL FUNCTIONS ===============================
% ========================= LOCAL FUNCTIONS ===============================

function se=se_auroc(es,n1,n2)
% ** function se_auroc(es,n1,n2) computes standard error of AUROC according
% to Hanley & McNeil Radiology 143:29-36, 1982
% let's break down the computations somewhat & stick to their terminology
q1=es./(2-es);
q2=2*es.^2./(1+es);
es2=es.^2;
se=sqrt((es-es2+(n1-1).*(q1-es2)+(n2-1).*(q2-es2))./(n1.*n2));


function simpleEsPlot(stats,esm)
% a simple routine for plotting results
nEs=numel(esm);
for tti=1:nEs
  curEs=esm{tti};
  % assign values to generic local variables es and ci
  curEsCi=[curEs 'Ci'];
  es=stats.(curEs);
  ci=stats.(curEsCi);
  % one figure per statistic
  figure(tti)
  clf
  hold on
  if numel(es)==1
    plot(1,es,'ko');
    plot([1 1],ci,'b+');
  else
    plot(es,'ko-');
    plot(ci','b-');
  end
end
  
  
% ======================= REPOSITORY =======================================

% The code below is a 'traditional' and inefficient way to compute the area
% under the curve of the ROC

%         % in the following lines, we generate an array of criterion values,
%         % equally spaced from the minimum to the maximum in [x and y
%         % combined], with at most rocNBin values. The array thus
%         % generated is used for both original and bootstrapped (if any)
%         % data.
%         xyso=sort([x(:,1); y(:,1)]);
%         xysod=unique_bc(sort(diff(xyso)));
%         % zero is not an acceptable bin width
%         if ~xysod(1)
%           xysod(1)=[];
%         end
%         % number of bins: use the lesser of [rocNBin, max-min divided
%         % by smallest acceptable difference between values]
%         nBin=min(rocNBin,round((xyso(end)-xyso(1))/xysod(1)));
%         % with the bins defined like this and function histc, used below,
%         % the last bin contains only values falling exactly on the last
%         % value. Shift the border of the last bin by a tiny amount to the
%         % right to ensure that i) the last bin will contain a count of zero
%         % (and can be discarded), and ii) the last but one bin will contain
%         % the maximal value in the sample(s) without extending far beyond
%         % this value
%         critVal=linspace(xyso(1),xyso(end)+eps(xyso(end)),nBin+1);
%         % now compute the cum histograms - separately for x and y
%         % because their sizes (number of rows) may differ
%         xch=cumsum(histc(x,critVal))./repmat(n1,[nBin+1 1]);
%         ych=cumsum(histc(y,critVal))./repmat(n2,[nBin+1 1]);
%         % the first bin always contains a nonzero value. We need to have a
%         % first bin with zero count (needed for proper calculation of the
%         % area under the curve), so append these
%         xch=cat(1,zeros(1,size(xch,2)),xch);
%         ych=cat(1,zeros(1,size(ych,2)),ych);
%         % close the curves by replacing the useless values in the last bin
%         % by coordinate (defined in x-y space) [1 0]
%         xch(end,:)=ones(1,size(xch,2));
%         ych(end,:)=zeros(1,size(ych,2));
%         % compute area under curve using polyarea 
%         es=polyarea(xch,ych);

% here are the 95% CI based on smax: Bamber, J Math Psychol 1975
% % -> appear unrealistically narrow
% smax=es.*(1-es)./(min(n1,n2)-1);
% tmp=norminv(1-alpha/2).*smax;

