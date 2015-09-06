function stats=mestab(table,varargin)
% ** function stats=mestab(table,varargin) computes basal parameters and
% effect size measures from N by M tables of categorical outcomes (mostly 
% 2 by 2 tables) and according confidence intervals where applicable and
% possible. The exact analyses/parameters need not be specified; as the
% computations involved are trivial all parameters will be computed and
% placed in output variable stats. All input parameters except table are
% optional and must be specified as parameter/value pairs in any order,
% e.g. as in
%      mestab(table,'confLevel',.90)
%
%                           INPUT ARGUMENTS
%                           ---------------
% stats=mestab(table)  computes all implemented parameters and effect size 
%   measure(s) listed below in TABLE 1. table must be minimally a 2 by 2 
%   table of frequencies of occurrence with the following row and column 
%   order:
%                          CHARACTERISTIC 1
%                        | - PRESENT         - NOT PRESENT
%       ------------------------------------------------
%       CHARACTERISTIC 2 |      
%       - PRESENT        |
%       - NOT PRESENT    |
%
%                        or (as a specific example)
% 
%                        | TRUE POSITIVE     TRUE NEGATIVE
%       ------------------------------------------------
%       TEST POSITIVE    |
%       TEST NEGATIVE    |
% 
%   To give another example, the rows are the dichotomous values of the
%   independent variable (e.g. control and treatment group) and the columns
%   those of the dependent variable (e.g. outcome):
%
%                        | OUTCOME NEGATIVE  OUTCOME POSITIVE
%       ------------------------------------------------
%       CONTROL          |
%       TREATMENT        |

%   If table is larger than 2 by 2 the only parameter to be computed is
%   Cramer's V.
% stats=mestab(...,'confLevel',0.90)  computes 90 % confidence intervals of 
%   the statistic in question (95 % ci are the default)
%
%                         OUTPUT ARGUMENTS
%                         ----------------
%   The results of the computations are placed into fields of output
%   argument stats with names corresponding to the analyses listed in
%   TABLE 1, e.g. 
%     .riskDifference
%     .riskRatio
%     ...
%   Where applicable, confidence intervals will be placed in fields 
%     .riskDifferenceCi
%     .riskRatioCi
%     ...
% Additional output fields are:
%     .confLevel (confidence level)
% -------------------------------------------------------------------------
% TABLE 1: PARAMETERS COMPUTED
% -------------------------------------------------------------------------
% FIELD NAME           QUANTITIY COMPUTED, SYNONYMS
%        -- measures for 2 by 2 tables --
% riskDifference       risk difference, proportion difference
% riskRatio            risk ratio, rate ratio
% oddsRatio            odds ratio
% phi                  degree of association
% baseRate             base rate
% sensitivity          sensitivity
% specificity          specificity
% posPredictValue      positive predictive value
% negPredictValue      negative predictive value
% successTreat         binomial effect size display: success rate (treated)
% successCtrl          binomial effect size display: success rate (control)
%        -- measures for M by N tables --
% cramerV              Cramer's V
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Version 1.2, March 2012
% Code by Harald Hentschke (University of Tübingen) and 
% Maik Stüttgen (University of Bochum)
% For additional information see Hentschke and Stüttgen, 
% Eur J Neurosci 34:1887-1894, 2011
% -------------------------------------------------------------------------

% ----- default values & varargin -----
% standard values of optional input arguments
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

% -------------------------------------------------------------------------
% ------- PART I: CHECKS OF INPUT ARGS & OTHER PREPARATORY WORKS ----------
% -------------------------------------------------------------------------
alpha=1-confLevel;
% --- check input table
[nRow nCol]=size(table);
if nRow<2 || nCol<2
  error('table size must at least be 2 by 2');
end
% what kinda table?
if nRow==2 && nCol==2
  is2by2=true;
else
  is2by2=false;
end
% reject foul values
if any(any(~isfinite(table)))
  error('input variable table contains foul data (nan or inf))');
end
% --- check other input arguments
if confLevel<=0 || confLevel>=1
  error('input variable ''confLevel'' must be a scalar of value >0 and <1');
end

% -------------------------------------------------------------------------
% ------- PART II: THE WORKS ----------
% -------------------------------------------------------------------------
if is2by2
  % first, for convenience, assign entries in table to alphabetic letters
  % according to Kline 2004 (p. 146) (all lowercase, though).
  % note linear indexing of a 2D-var
  a=table(1);
  b=table(3);
  c=table(2);
  d=table(4);
  
  % compute basal parameters 
  stats.baseRate=(a+c)/(a+b+c+d);
  stats.sensitivity=a/(a+c);
  stats.specificity=d/(b+d);
  stats.posPredictValue=a/(a+b);
  stats.negPredictValue=d/(c+d);
  
  % now compute values based on proportions. As all computations involved
  % are trivial and variable sizes miniscule let's not worry about their
  % redundancy and proliferation, respectively, and use the same
  % nomenclature as Kline 2004
  p_c=a/(a+b);
  p_t=c/(c+d);
  % will be needed for confidence intervals
  z2tailFactor=norminv(alpha/2)*[1 -1];
  
  % risk difference  
  stats.riskDifference=p_c-p_t;
  % confidence intervals
  stats.riskDifferenceCi=stats.riskDifference + ...
    sqrt(p_c*(1-p_c)/(a+b) + p_t*(1-p_t)/(c+d)) * z2tailFactor;
  % risk ratio
  stats.riskRatio=p_c/p_t;
  % confidence intervals of risk ratio: for better readibility compute in
  % two steps
  tmp=log(stats.riskRatio) + ...
    sqrt((1-p_c)/((a+b)*p_c) + (1-p_t)/((c+d)*p_t)) * z2tailFactor;
  stats.riskRatioCi=exp(tmp);
  % odds ratio
  stats.oddsRatio=a*d/(b*c);
  % ci: same story 
  tmp=log(stats.oddsRatio) + ...
    sqrt( 1/((a+b)*p_c*(1-p_c)) + 1/((c+d)*p_t*(1-p_t))) * z2tailFactor;
  stats.oddsRatioCi=exp(tmp);
  % the measure of association termed phi [equivalent:
  % stats.phi=(a*d-b*c)/sqrt((a+b)*(c+d)*(a+c)*(b+d));]
    
  % ** use code for Cramer's V for confidence intervals
  [stats.phi stats.phiCi chi2]=cramerv(table,nRow,nCol,confLevel,is2by2);
  % finally, binomial effect size display (to quote Randolph & Edmondson
  % 2005, 'What would the correlationally equivalent effect of the
  % treatment be if 50% of the participants had the occurrence and 50% did
  % not and 50% received treatment and 50% did not?')
  stats.successTreat=.5+stats.phi/2;
  stats.successCtrl=.5-stats.phi/2;
else
  % the only thing to be computed here is Cramer's V (and chi2)
  [stats.cramerV stats.cramerVCi chi2]=cramerv(table,nRow,nCol,confLevel,is2by2);
end

% for the sake of being complete, add chi square to stats
stats.chi2=chi2;


% ========================= LOCAL FUNCTIONS ===============================
% ========================= LOCAL FUNCTIONS ===============================

function [es,esCi,chi2]=cramerv(table,nRow,nCol,confLevel,is2by2)
% this function computes Cramer's V, including exact analytical CI
% ** NOTE: in the case of 2 by 2 tables Cramer's V is identical to phi
% except possibly for the sign), which will be taken care of in the last
% lines
colSum=sum(table);
rowSum=sum(table,2);
n=sum(sum(table));
k=min(nRow,nCol);
df=(nRow-1)*(nCol-1);
% expected frequency of occurrence in each cell: product of row and
% column totals divided by total N
ef=(rowSum*colSum)/n;
% chi square stats
chi2=(table-ef).^2./ef;
chi2=sum(chi2(:));
% Cramer's V
es=sqrt(chi2/(n*(k-1)));
% CI (Smithson 2003, p. 40)
ncp=ncpci(chi2,'X2',df,'confLevel',confLevel);
esCi=sqrt((ncp+df)/(n*(k-1)));
% in case we are dealing with 2 by 2 tables heed sign
if is2by2
  if det(table)<0
    es=es*-1;
    esCi=fliplr(esCi)*-1;
  end
end
