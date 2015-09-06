function [stats,varargout]=mes1way(X,esm,varargin)
% ** function [stats,varargout]=mes1way(X,esm,varargin)
% computes measures of effect size between two or more samples and places
% the results in output structure stats. A summary table including the
% effect size measures and one-way ANOVA results will be displayed in the
% command window. All input parameters except X and esm are optional and
% must be specified as parameter/value pairs in any order, e.g. as in
%      mes1way(X,'eta2','group',g);
%
% For information on assumptions of the underlying model and related
% information please see the notes below TABLE 1 further down as well as
% the documentation accompanying this code.
%
%                           INPUT ARGUMENTS
%                           ---------------
% stats=mes1way(X,esm)  computes the effect size measure(s) specified in 
%   input variable esm (see table at bottom) from input array X. X must
%   have two or more columns, each representing a different group, or may
%   be a single column-array (see below). Input variable esm must be a
%   character array (e.g. 'eta2') or a cell array (e.g.
%   {'eta2','partialeta2'}).
% stats=mes1way(...,'group',g)  allows for an alternative structure of 
%   input data X. X must be a single-column vector and g a single-column
%   vector of arbitrary numbers coding for the different groups.
% stats=mes1way(...,'isDep',1)  assumes that the samples in X are dependent
%   (repeated measures); accordingly, the number of samples in each group
%   of X must be equal. If this parameter is omitted the assumption is one
%   of independent samples. Both the effect size measures available and
%   their computation are contingent on whether the samples are independent
%   or dependent. NOTE: if data is dependent, it is assumed that input data
%   X are sorted according to subjects. That is, if X is a multi-column-
%   array all data from each subject reside in one row. If X is a single-
%   column array, the equivalent assumption is made: that data for subjects
%   are listed in the order in which they occur within each group.
% stats=mes1way(...,'missVal','listwise')  omits NaNs and Infs in X,
%   summarily termed here as missing values, in listwise fashion: if X is a
%   multi-column array and there is a missing value anywhere in X, the
%   entire row of X will be omitted from analysis. If X is in single-column
%   format, positions corresponding to the row(s) in question will be
%   omitted. If set to 'pairwise', only the individual missing data will be
%   omitted. Note: 'listwise' is mandatory (and therefore the default) for
%   dependent data. Likewise, 'pairwise' is mandatory and default for
%   independent unbalanced data. It is also the default for independent,
%   balanced data.
% stats=mes1way(...,'nBoot',n)  computes confidence intervals of 
%   the statistic in question via bootstrapping n times. n should be on the
%   order of thousands, otherwise bootstrapping will lead to inaccurate
%   results. By default or if the value of nBoot is set to zero or
%   nonsensical values (<0, infinity, nan) confidence intervals will be 
%   computed analytically (where possible).
% stats=mes1way(...,'confLevel',0.90)  computes 90 % confidence intervals 
%   of the statistic in question (95 % ci are the default; any value may be 
%   specified)
% stats=mes1way(...,'cWeight',c)  allows specification of contrast weights
%   for the computation of effect size measures like standardized contrasts
%   and eta squared for focussed comparisons (see TABLE 1). Input array c
%   must contain as many columns as there are groups; you may specifiy as
%   many rows (=contrasts) as you wish.
% stats=mes1way(...,'tDenom','msw')  computes, in case of dependent data, F
%   and p values of contrasts, as well as confidence intervals of psi and
%   g_psi, from MS_[between x subject]. The default is 'sd', which means
%   that the values are based on the contrasts' difference score. Please
%   see the documentation for detailed explanation.
% -------------------------------------------------------------------------
% TABLE 1: ANALYSES TO BE SPECIFIED IN INPUT VARIABLE esm
% -------------------------------------------------------------------------
% esm             QUANTITIY COMPUTED ((*) require input of contrast weights)
% 'psi'           unstandardized contrast (*)
% 'g_psi'         standardized contrast (*)  
% 'psibysd'       standardized contrast for dependent data (*)
% 'eta2'          eta squared
% 'partialeta2'   partial eta squared 
% 'omega2',       omega squared
% 'partialomega2' partial omega squared
%
% Notes:
%   i. Parameters marked by (*) need contrast weights to be
%   computed; by definition they do not exist for an omnibus effect. The
%   other parameters will be computed for both the omnibus effect and the
%   specified contrasts.
%  ii. 'g_psi' is the oneway equivalent of Hedges' g, namely, contrast
%   divided by the square root of the pooled within-conditions variance.
%   Its value is identical for dependent and independent data, but for
%   dependent data confidence intervals will be smaller in case of
%   correlations between the groups.
% iii. 'psibysd' is a standardized contrast for dependent data only, 
%   namely, contrast divided by the standard deviation of the difference
%   score. It is the oneway equivalent of mdbysd in mes.m.
% iv. Fixed factors are assumed.
% v. Independent data may be unbalanced.
% vi. In the computation of omega2 and partial omega2 for dependent data,
%   an additive model (no subject x treatment interaction) is assumed
%
%                         OUTPUT ARGUMENTS
%                         ----------------
% The results of the computations are placed into fields of output
% structure stats with names corresponding to the analyses requested. For
% example, stats=mes1way(X,{'eta2','partialeta2'}) will result in fields
%   .eta2 
%   .partialeta2
% PLEASE NOTE: if contrast weights were specified, all output arguments
% will be single- or two-column arrays holding results in the following
% (row) order:
% 1. omnibus effect
% 2. first contrast
% 3. second contrast
% and so on. If a parameter is not defined for the omnibus effect or for
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
%   .cWeight (the contrast weights)
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

% ----- default values & varargin -----
% standard values of optional input arguments
group=[];
isDep=false;
missVal='pairwise';
nBoot=0;
cWeight=[];
confLevel=.95;
tDenom='sd'; 

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
% minimal number of bootstrapping runs below which a warning will be issued
% (and bootstrapping will be skipped)
minNBootstrap=1000;

% -------------------------------------------------------------------------
% ------- PART I: CHECKS OF INPUT ARGS & OTHER PREPARATORY WORKS ----------
% -------------------------------------------------------------------------
% Variable list_analyses below is a list of all possible analyses:
% - column 1 holds the shortcut as strings
% - column 2 holds 1 for quantities REQUIRING a contrast
% - column 3 is a short description of the analyses
% This is the 'master list' of all analyses against which input argument
% esm will be checked (see) below
list_analysis={...
'psi',           1, 'unstandardized contrast';...
'g_psi',         1, 'standardized contrast';...
'psibysd',       1, 'standardized contrast for dependent data';...
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
tmp=find(strcmp('psibysd',esm));
if ~isempty(tmp) && ~isDep
  warning('parameter ''psibysd'' can only be computed for dependent data - eliminating it from list of parameters to be computed')
  esm(tmp)=[];
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
% be done further below)
if isempty(cWeight)
  doContrast=false;
  if any(uTag)
    error('part of the requested analyses require contrast weights as input args')
  end
else
  doContrast=true;
end

% illegal values for missVal
if isempty(find(strcmpi(missVal,{'listwise','pairwise'})))
  error('illegal value for input parameter ''missVal'' (choose ''listwise'' or ''pairwise)''');
end
if isDep && ~strcmpi(missVal,'listwise')
  warning('dependent data - setting input parameter ''missVal'' to ''listwise''');
  missVal='listwise';
end

% --- check and possibly transform input X: if it is a multiple column
% array, transform it into a single column array and set up variable group.
% This may appear counterintuitive as the computations are much easier and
% elegantly done on a matrix, in a vectorized version. However, this holds
% true only for balanced data. In real life more often than not we have to
% deal with unbalanced data, a vectorized computation of which would be
% cumbersome and, in conjunction with bootstrapping, inefficient.
% Therefore, the code is designed to deal with the most general case, and
% requires reshaping multiple column arrays into a single column array.
[nRowX nColX]=size(X);
if nColX>1
  % *** X is a multiple column-array:
  if ~isempty(group)
    % - group specified although not necessary
    warning('input variable ''group'' will be ignored with X having several columns');
    % instead, it will be generated internally
  end
  nGroup=nColX;
  group=repmat(1:nGroup,nRowX,1);
  goodIx=isfinite(X);
  nSample=sum(goodIx);
  % transformation of X and group into single column arrays and generation
  % of groupIx may need to be done after the switch construct, unless the
  % 'pairwise' case below is run
  doTransform=true;
  % missing values (nans and infs) in X
  if any(any(~goodIx))
    switch lower(missVal)
      case 'listwise'
        % omit bad row(s)
        badRowIx=any(~goodIx,2);
        X(badRowIx,:)=[];
        group(badRowIx,:)=[];
        nRowX=size(X,1);
        nSample=repmat(nRowX,1,nColX);
        % (value of doTransform unchanged)
      case 'pairwise'
        % omit bad values, at the same time transforming X and group into
        % single col-format
        X=X(goodIx);
        group=group(goodIx);
        % generate groupIx, the indexes coding for group assignment
        tmp=cumsum([0 nSample]);
        for gIx=nGroup:-1:1
          groupIx{gIx}=tmp(gIx)+1:tmp(gIx+1);
        end
        doTransform=false;
        % nSample unchanged
    end
  end   
  if doTransform
    % transform X and group into single col-format
    X=X(:);
    group=group(:);
    for gIx=nGroup:-1:1
      groupIx{gIx}=(gIx-1)*nSample+1:gIx*nSample;
    end
  end
else
  % *** X is a single column-array:
  if isempty(group)
    % - group not specified 
    error('input variable ''group'' must be specified when X is a single column-array');
  else
    [tmpNRow tmpNCol]=size(group);
    % - group not a single col-array
    if tmpNCol>1
      error('input variable group must be a single-column array');
    end
    % - numel not matching
    if tmpNRow~=nRowX
      error('input variables group and X must have the same number of elements');
    end
    % - undefined elements in group
    if any(~isfinite(group))
      error('input variable group contains NaNs or Infs');
    end
  end
  % determine how many groups there are
  [uGroup,aIx]=unique_bc(group);
  nGroup=numel(uGroup);
  % revert sorting implicitly done by unique because the user probably
  % expects the original order to be maintained
  uGroup=group(sort(aIx));
  % if group is unsorted in the sense that samples from any group do not
  % form a contiguous block it appears possible that the data are messed
  % up, so better generate an error here and force the user to rethink (and
  % sort) the data
  isContiguousGroups=numel(find(diff(group)))==nGroup-1;
  if ~isContiguousGroups
    error([mfilename ' expects samples from each group to form contiguous blocks. Please sort your data accordingly']);
  end
  if isDep && mod(nRowX,nGroup)
    error('input variable X supposedly holds dependent data, but sample sizes are unequal');
  end
  if strcmpi(missVal,'listwise') && mod(nRowX,nGroup)
    warning('listwise elimination of missing data is undefined with unbalanced data - employing pairwise elimination (if any)');
    missVal='pairwise';
  end
  % logical index to missing values (nans and infs) in X
  badIx= ~isfinite(X);
  % if listwise elimination is requested and we have dependent or balanced data ...
  if strcmpi(missVal,'listwise') && ~mod(nRowX,nGroup)
    if ~isempty(any(badIx))
      badIx=reshape(badIx,nRowX/nGroup,nGroup);
      badIx(any(badIx,2),:)=true;
      badIx=badIx(:);
      X(badIx)=[];
      group(badIx)=[];
      nRowX=size(X,1);
    end
    % scalar expansion of variable nSample
    nSample=repmat(nRowX/nGroup,1,nGroup);
    for gIx=nGroup:-1:1
      groupIx{gIx}=(gIx-1)*nSample+1:gIx*nSample;
    end
    % variable group is not needed anymore
    group=[];
  else
    % pairwise elimination requested or independent/unbalanced data 
    if ~isempty(any(badIx))
      X(badIx)=[];
      group(badIx)=[];
    end
    for gIx=nGroup:-1:1
      groupIx{gIx}=find(group==uGroup(gIx));
      nSample(gIx)=numel(groupIx{gIx});
    end
  end
end
% At this point, variables nGroup and nSample hold the number of groups (a
% scalar) and the number of samples in each group (a row array),
% respectively, and X is a single column-array of data devoid of NaNs and
% Infs; group association is coded by variable group

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

[nCo nCoCol]=size(cWeight);
if doContrast
  if nCoCol~=nGroup
    error('number of contrast weights (per row) does not match number of groups in input data');
  end
  if any(any(~isfinite(cWeight)))
    error('contrast weights contain foul data (NaN or Inf)');
  end
  % check whether we're dealing with standard sets
  if any(abs(sum(abs(cWeight),2)-2)>eps(2))
    warning('at least one set of contrast weights is not a standard set: standardized mean differences (g_psi, psibysd) computed with this set do not represent the difference between the averages of two subsets of means')
  end
  % is there more than one contrast weight of zero? This may reflect the
  % user's intention to have one or more conditions 'set aside', in which
  % case the computations will not be valid
  if any(sum(cWeight==0,2)>1)
    warning('at least one set of contrast weights contains more than one zeroes - if you wish to exclude the corresponding groups from analysis you should eliminate them prior to computation')
  end
end
if ~any(ismember_bc(tDenom,{'sd','msw'}))
  error('bad choice for input parameter ''tDenom'' (choose ''sd'' or ''msw'')');
end
% -------------------------------------------------------------------------
% ------- PART II: THE WORKS ----------
% -------------------------------------------------------------------------

% create a few standard fields
stats.isDep=isDep;
stats.nBoot=nBoot;
stats.confLevel=confLevel;
stats.n=nSample;
stats.cWeight=cWeight;

% preparatory computations of SS, MS etc.
r=prepcomp(X,groupIx,nSample,isDep,doBoot,nBoot,cWeight,tDenom);

% if eta2, partialeta2, omega2 or partialomega2 or any combination thereof
% is requested AND analytical CI shall be computed AND the data are
% independent, compute CI for partial eta2 here because those of most other
% ESM can be derived from them (and computing all independently of each
% other would be a waste).
if any(ismember_bc({'eta2','partialeta2','omega2','partialomega2'},esm)) && ~doBoot && ~isDep
  % exact analytical confidence intervals for eta2 and/or partial eta2
  % (Smithson 2003, p.43 [5.6]), omnibus effects of which are identical in
  % oneway independent designs
  pEta2Ci=ncpci(r.F,'F',[r.dfGroup,r.dfErr],'confLevel',confLevel);
  pEta2Ci=pEta2Ci./(pEta2Ci+r.dfGroup+r.dfErr+1);
  for g=1:nCo
    % exact ci for contrasts: these apply only to PARTIAL eta2 because 
    % partial eta2 for a contrast = SSpsi/(SSpsi+SSerr) whereas eta2 for a
    % contrast = SSpsi/SStot = SSpsi/(SSa+SSerr), involving three terms
    pEta2Ci(g+1,:)=ncpci(r.FPsi(g),'F',[1,r.dfErr],'confLevel',confLevel);
    pEta2Ci(g+1,:)=pEta2Ci(g+1,:)./(pEta2Ci(g+1,:)+1+r.dfErr+1);
  end
  % information on kind of CI
  pEta2CiType=repmat('exact analytical',nCo+1,1);
end

% now computations of effect size parameters:
nEs=numel(esm);
% a temporary container to hold [es ci_lo ci_up] in the order in
% which they are computed (for the table of results to be displayed)
esStore=repmat(nan,[nCo+1 nEs*3]);

for tti=1:nEs
  curEs=esm{tti};
  es=repmat(nan,nCo+1,1);
  % If analytical confidence intervals for the statistic in question cannot
  % be computed and no bootstrapping is requested, ci must be set to nan.
  % Do this here for all statistics to avoid redundant assignments. ci will
  % be overwritten with meaningful values if i) there is a formula for
  % analytical confidence intervals implemented within the respective case
  % in the switch statement, or ii) confidence intervals based on
  % bootstrapping are computed right after the switch statement
  ci=repmat(nan,nCo+1,2);
  % a similar argument applies to variable ciType, which indicates on which
  % method the computation of ci is based (which may differ between the
  % omnibus effect and contrasts, hence we need more than one row
  if doBoot
    ciType=repmat('bootstrap',nCo+1,1);
  else
    ciType=repmat('none',nCo+1,1);
  end
  switch curEs
    case 'psi'
      % nan as first row because a contrast is not an omnibus effect; then
      % the contrast(s)
      es=cat(1,repmat(nan,1,nBoot+1),r.psi);
      if ~doBoot
        % note: confidence intervals of contrasts can be computed in an
        % exact fashion from central t distributions and are therefore
        % tagged 'exact analytical' here, although exactness does not, as
        % in the case of other mes below, imply computation via noncentral
        % distributions
        if isDep
          % dependent samples: 
          if strcmpi(tDenom,'sd')
            % exact CI of contrast
            ci(2:nCo+1,:)=repmat(r.psi,1,2) + r.sD/sqrt(nSample(1)) * tinv([alpha/2 1-alpha/2],r.dfSubj);
          elseif strcmpi(tDenom,'msw')
            % exact CI of contrast
            for g=1:nCo
              % SE of contrast
              fac=sqrt(r.msGroupSubj*sum((cWeight(g,:)).^2./nSample));
              ci(g+1,:)=repmat(r.psi(g),1,2) + fac * tinv([alpha/2 1-alpha/2],r.dfGroupSubj);
            end
          end
          ciType=char('none', repmat('exact analytical',nCo,1));
        else
          % independent samples:
          for g=1:nCo
            % SE of contrast
            fac=sqrt(r.msErr*sum((cWeight(g,:)).^2./nSample));
            ci(g+1,:)=repmat(r.psi(g),1,2) + fac * tinv([alpha/2 1-alpha/2],r.dfErr);
          end
          ciType=char('none', repmat('exact analytical',nCo,1));
        end
      end
    
    case 'g_psi'
      % nan as first row because this mes is not defined for an omnibus
      % effect; then the contrasts divided by standardizer proposed by
      % Kline 2004, the square root of the pooled within-conditions
      % variance, making this the equivalent of Hedges' g
      es=cat(1,repmat(nan,1,nBoot+1),r.psi./repmat(sqrt(r.msErr),nCo,1));
      if ~doBoot
        if isDep
          % dependent samples: 
          if strcmpi(tDenom,'sd')
            % exact CI of contrast...
            tmp=repmat(r.psi,1,2) + r.sD/sqrt(nSample(1)) * tinv([alpha/2 1-alpha/2],r.dfSubj);
            % ...leading to approximate CI of g_psi
            ci(2:nCo+1,:)=tmp/sqrt(r.msErr);
          elseif strcmpi(tDenom,'msw')
            % exact CI of contrast...
            for g=1:nCo
              % SE of contrast
              fac=sqrt(r.msGroupSubj*sum((cWeight(g,:)).^2./nSample));
              % ...leading to approximate CI of g_psi
              tmp=repmat(r.psi(g),1,2) + fac * tinv([alpha/2 1-alpha/2],r.dfGroupSubj);
              ci(g+1,:)=tmp/sqrt(r.msErr);
            end
          end
          ciType=char('none', repmat('approximate analytical',nCo,1));
        else
          % independent samples: exact confidence intervals 
          for g=1:nCo
            fac=sqrt(sum((cWeight(g,:)).^2./nSample));
            ci(g+1,:)=fac*ncpci(r.tPsi(g),'t',r.dfErr,'confLevel',confLevel);
          end
          ciType=char('none', repmat('exact analytical',nCo,1));
        end
      end
      
    case 'psibysd'
      % psibysd is defined only for dependent data; as the criterion of
      % dependence has been dealt with in the code above we do not heed it
      % here
      % - nan as first row, then the contrasts divided by standard
      % deviation of difference scores
      es=cat(1,repmat(nan,1,nBoot+1),r.psi./r.sD);
      if ~doBoot
        % as for g_psi above: exact CI of contrast, from which approximate
        % CI for es can be derived
        tmp=repmat(r.psi,1,2) + r.sD/sqrt(nSample(1)) * tinv([alpha/2 1-alpha/2],r.dfSubj);
        ci(2:nCo+1,:)=tmp./repmat(r.sD,1,2);
        ciType=char('none', repmat('approximate analytical',nCo,1));
      end
      
    case 'eta2'
      % eta2 does not differentiate between dependent and independent
      % samples
      es=r.ssGroup./r.ssTot;
      if doContrast
        es=cat(1,es,r.ssPsi./repmat(r.ssTot,nCo,1));
      end
      if ~doBoot
        if ~isDep
          % as explained above, ci for the omnibus effect are computable,
          % those for the contrasts are not
          ci(1,:)=pEta2Ci(1,:);
          ciType=char('exact analytical', repmat('none',nCo,1));
        end
      end
      
    case 'partialeta2'
      if isDep
        % dependent case: SSerror=r.ssGroupSubj
        es=r.ssGroup./(r.ssGroup+r.ssGroupSubj);
        if doContrast
          es=cat(1,es,r.ssPsi./(r.ssPsi+repmat(r.ssGroupSubj,nCo,1)));
        end
      else
        % independent samples: main effect same as eta squared because
        % r.ssTot=r.ssGroup+r.ssErr
        es=r.ssGroup./(r.ssTot);
        % contrast: denominator differs from dependent case
        if doContrast
          es=cat(1,es,r.ssPsi./(r.ssPsi+repmat(r.ssErr,nCo,1)));
        end
        if ~doBoot
          % all ci computed above apply to partial eta2
          ci=pEta2Ci;
          ciType=pEta2CiType;
        end
      end      
                
    case 'omega2'
      if isDep
        % formula according to Kline p. 188 & Table 6.8, setting
        % MSeffect=r.msGroup and MSerror=r.msGroupSubj
        es=(r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj))./...
          (r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj) + ...
          1/nGroup*(r.msSubj-r.msGroupSubj) + r.msGroupSubj);
        if doContrast
          tmp=(1/sum(nSample)*(r.ssPsi-repmat(r.msGroupSubj,nCo,1)))./...
            repmat(r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj) + ...
            1/nGroup*(r.msSubj-r.msGroupSubj) + r.msGroupSubj,nCo,1);
          es=cat(1,es,tmp);      
        end
      else
        es=(r.ssGroup-(nGroup-1).*r.msErr)./(r.ssTot+r.msErr);
        if doContrast
          % contrasts: derived from formula 6.30, p. 186, Kline 2004, and
          % setting df(group)=1 and MS=SS
          es=cat(1,es,(r.ssPsi-repmat(r.msErr,nCo,1))./repmat(r.ssTot+r.msErr,nCo,1));
        end
        if ~doBoot
          % exact analytical confidence intervals according to Fidler &
          % Thompson 2001, p. 593, based on those for partial eta2
          tmp=(r.ssTot*(1-pEta2Ci))/r.dfErr;
          ci(1,:)=(r.ssTot*pEta2Ci(1,:)-r.dfGroup*tmp(1,:))./(r.ssTot+tmp(1,:));
          % as in the case of eta2 and contrasts, we cannot compute ci for
          % contrasts for omega2 for the same reasons stated above
          ciType=char('exact analytical', repmat('none',nCo,1));
        end
      end
      
    case 'partialomega2'
      if isDep
        % formula according to Kline p. 188 & Table 6.8, setting
        % MSeffect=r.msGroup and MSerror=r.msGroupSubj (note that the only
        % difference to the corresponding case of omega2 is one term less
        % in the denominator, namely the subjects' variance component
        es=(r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj))./...
          (r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj) + ...
          r.msGroupSubj);
        if doContrast
          tmp=(1/sum(nSample)*(r.ssPsi-repmat(r.msGroupSubj,nCo,1)))./...
            repmat(r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj) + ...
            r.msGroupSubj,nCo,1);
          es=cat(1,es,tmp);      
        end
      else
        % independent samples: same as omega squared
        es=(r.ssGroup-(nGroup-1).*r.msErr)./(r.ssTot+r.msErr);
        if doContrast
          % contrasts: assembled from Kline p. 186 & Table 6.7
          tmp=(1/sum(nSample)*(r.ssPsi-repmat(r.msErr,nCo,1)))./...
            (1/sum(nSample)*(r.ssPsi-repmat(r.msErr,nCo,1)) + repmat(r.msErr,nCo,1));
          es=cat(1,es,tmp);
        end
        if ~doBoot
          % exact analytical confidence intervals according to Fidler &
          % Thompson 2001, p. 593, based on those for partial eta2
          tmp=(r.ssTot*(1-pEta2Ci))/r.dfErr;
          ci=(r.ssTot*pEta2Ci(1,:)-r.dfGroup*tmp(1,:))./(r.ssTot+tmp(1,:));
          if doContrast
            % same formula applied to CI of contrasts, taking into account df(group)=1 
            ci(2:nCo+1,:)=(r.ssTot*pEta2Ci(2:nCo+1,:)-1*tmp(2:nCo+1,:))./(r.ssTot+tmp(2:nCo+1,:));
          end          
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
% chars to be added to headers of the F and p columns
if isDep && nCo>0 && strcmpi(tDenom,'sd')
  ac=' (**)';
else
  ac='';
end

esm=cat(1,esm,ciTi);
esm=esm(:)';
% column order:
% source | SS | df | MS | F | p | effect sizes (in the order requested)
table={...
  'source',                'SS',            'df',          'MS',            ['F',ac]   ['p',ac] esm{1:end};...
  'between-groups',        r.ssGroup(1),    r.dfGroup,     r.msGroup(1),    r.F(1),    r.p(1),  esStore{1,:};...
  };
for g=1:nCo
  table=cat(1,table,...
    {['- psi' int2str(g)], r.ssPsi(g),      1,             r.msPsi(g),      r.FPsi(g), r.pPsi(g), esStore{1+g,:}});
end
table=cat(1,table,...
  {'within-groups',        r.ssErr(1),      r.dfErr,       r.msErr(1),      [],      [],        filla{1:end}});
if isDep
  table=cat(1,table,...
    {'- subj',             r.ssSubj(1),     r.dfSubj,      r.msSubj(1),     [],      [],        filla{1:end};...
    '- group x subj',      r.ssGroupSubj(1),r.dfGroupSubj, r.msGroupSubj(1),[],      [],        filla{1:end}});
end
table=cat(1,table,...
  {'total',                r.ssTot(1),       r.dfTot,       r.msTot(1),     [],      [],        filla{1:end}});

table  
if isDep && nCo>0 && strcmpi(tDenom,'sd')
  disp(char({...
    '(**) NOTE: F and p values and the confidence intervals of psi and';...
    '     g_psi are derived from the standard deviation of the contrasts'' ';...
    '     difference scores. For an alternative computation set optional ';...
    '     input parameter ''tDenom'' to ''msw''.';...
    '     See the documentation for detailed information.'}));
end

if nargout>1
  varargout{1}=table;
end

% ========================= LOCAL FUNCTIONS ===============================
% ========================= LOCAL FUNCTIONS ===============================

function r=prepcomp(x,groupIx,nSample,isDep,doBoot,nBoot,cWeight,tDenom)
% performs preparatory computations for oneway analyses on single column
% array x with group assignment coded by input var group

nGroup=numel(groupIx);
nCo=size(cWeight,1);
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
  if isDep
    % dependent data: generate random sampling index for first group and
    % replicate for all other groups (which are identical in size)
    bootSamX=repmat(ceil(nSample(1)*rand([nSample(1) nBoot])),[nGroup,1]);
    % add offset due to arrangement of data in a single-column array
    for g=1:nGroup
      bootSamX(groupIx{g},:)=bootSamX(groupIx{g},:)+groupIx{g}(1)-1;
    end
  else
    % independent data: resample data within each group independent of data
    % points picked in other groups
    bootSamX=repmat(nan,sum(nSample),nBoot);
    for g=1:nGroup
      bootSamX(groupIx{g},1:nBoot)=groupIx{g}(1)-1+ceil(nSample(g)*rand([nSample(g) nBoot]));
    end
  end
  % from second column onwards fill with sampled data
  x(:,2:nBoot+1)=x(bootSamX);
  % delete resampling index 
  clear bootSam* 
end

% between-groups df
r.dfGroup=nGroup-1;
% within-groups df
r.dfErr=sum(nSample-1);
% within-subjects df
r.dfSubj=nSample(1)-1;
% group x subjects df 
r.dfGroupSubj=r.dfGroup*r.dfSubj;

% grand mean
r.meanGrand=mean(x);

r.ssGroup=0;
r.ssErr=0;
for g=nGroup:-1:1
  % means of groups
  r.meanGroup(g,:)=sum(x(groupIx{g},:))/nSample(g);
  % SS_a, or SS_effect, between-groups SS
  r.ssGroup=r.ssGroup + nSample(g)*(r.meanGroup(g,:)-r.meanGrand).^2;
  % SS_error, within-groups SS
  r.ssErr=r.ssErr + sum((x(groupIx{g},:)-repmat(r.meanGroup(g,:),nSample(g),1)).^2);
end
% between-groups MS
r.msGroup=r.ssGroup/r.dfGroup;
% within-groups MS
r.msErr=r.ssErr/r.dfErr;

% SS_t, the sum of both
r.ssTot=r.ssGroup+r.ssErr;
% corresponding df
r.dfTot=r.dfGroup+r.dfErr;
r.msTot=r.ssTot/r.dfTot;

% parameters for dependent data
if isDep
  nSubj=nSample(1);
  r.ssSubj=zeros(1,nBoot+1);
  for g=nSubj:-1:1
    % means of subjects
    r.meanSubj(g,:)=mean(x(g:nSubj:end,:));
    % SS_subj, between-subjects SS
    r.ssSubj=r.ssSubj+nGroup*(r.meanSubj(g,:)-r.meanGrand).^2;
  end
  % SS group x subject
  r.ssGroupSubj=r.ssErr-r.ssSubj;
  % MS
  r.msSubj=r.ssSubj/r.dfSubj;
  r.msGroupSubj=r.ssGroupSubj/r.dfGroupSubj;
end
    
% compute contrasts and their SS
if nCo>0
  for g=nCo:-1:1
    % contrasts
    r.psi(g,:)=sum(repmat(cWeight(g,:)',1,nBoot+1).*r.meanGroup);
    % factor needed for computation of SS_psi and t for contrasts
    r.cFac(g,:)=sum(cWeight(g,:).^2./nSample);
    % SS_psi: see Kline 2004 (formula 6.8, p.168)
    r.ssPsi(g,:)=r.psi(g,:).^2./r.cFac(g,:);
    if isDep
      % std of the difference score
      tmp=repmat(cWeight(g,:),nSample(1),1);
      tmp=repmat(tmp(:),1,nBoot+1);
      tmp=reshape(x.*tmp,[nSample(1) nGroup nBoot+1]);
      r.sD(g,:)=permute(std(sum(tmp,2)),[1 3 2]);
    end
  end
  % as df is always 1 MS=SS, but create field nonetheless
  r.msPsi=r.ssPsi;
end

% finally, t, F, p
if isDep
  r.F=r.msGroup./r.msGroupSubj;
  r.p=1-fcdf(r.F,r.dfGroup,r.dfGroupSubj);
  if nCo>0
    if strcmpi(tDenom,'sd')
      % compute p values of contrasts from t statistic (Kline 2004, p.
      % 168): note that this t statistic is based on the std of the
      % contrast's difference score (in the denominator, see line below).
      % This means that i) t takes into account only the variability
      % between the groups/subjects compared in the contrast; specifically,
      % groups with a contrast weight of zero do not play a role, and ii) a
      % potential lack of sphericity among groups in the full data set is
      % not a problem.
      r.tPsi=r.psi./(r.sD/sqrt(nSample(1)));
      r.pPsi=2*(1-tcdf(abs(r.tPsi),r.dfSubj));
      r.FPsi=r.tPsi.^2;
    elseif strcmpi(tDenom,'msw')
      % F and p values are computed the 'classical' way, which takes
      % variability from all groups into account
      r.FPsi=r.msPsi./repmat(r.msGroupSubj,nCo,1);
      r.pPsi=1-fcdf(r.FPsi,1,r.dfGroupSubj);
      r.tPsi=sign(r.psi).*sqrt(r.FPsi);
    end
  end
else
  r.F=r.msGroup./r.msErr;
  r.p=1-fcdf(r.F,r.dfGroup,r.dfErr);
  if nCo>0
    % compute p values of contrasts from independent data on the basis of t
    % statistic (Kline 2004, formula 6.7, p. 168)
    r.tPsi=r.psi./sqrt(repmat(r.msErr,nCo,1).*repmat(r.cFac,1,nBoot+1));
    r.pPsi=2*(1-tcdf(abs(r.tPsi),r.dfErr));    
    % F is by definition the square of t
    r.FPsi=r.tPsi.^2;
    % § note: for the example in Kline 2004, table 6.4, the F values are
    % all correct but the contrast-associated p values do not match
  end
end
