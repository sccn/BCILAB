function ci=ncpci(x,fType,df,varargin)
%  ** function ci=ncpci(x,fType,df,varargin)
% iteratively approaches two-sided confidence intervals for the
% noncentrality parameter of a noncentral Chi square (abbreviated X2), F or
% t distribution with degrees of freedom df, given an abscissa value (X2, F
% or t value). This is achieved by varying the X2, F or t noncentrality
% parameter of the corresponding probability distribution function (pdf)
% until the given abscissa value is, within a certain precision, at the
% percentile values required for the confidence interval (2.5th and 97.5th
% percentile for lower and upper 95 % confidence intervals, respectively).
% All input parameters listed below except x, fType and df are
% optional and must be specified as parameter/value pairs, e.g. as in
%          ncpci(x,'t',df,'confLevel',.9)
%
%               >>> INPUT VARIABLES >>>
% NAME          TYPE/DEFAULT      DESCRIPTION
% x             double scalar     X2, F or t value
% fType         char              'X2','F' or 't'
% df            scalar or array   degrees of freedom 
%                                 (F pdf: [numerator denominator])
% confLevel     double, 0.95      confidence level
% prec          double scalar,    precision: iteration will run until the
%                1e-6             estimated percentile is <=prec away from
%                                 the requested percentile
% doAnimate     logical           if true, the iteration process will be
%                                 graphically displayed in a figure window
%                     
%               <<< OUTPUT VARIABLES <<<
% NAME          TYPE/DEFAULT           DESCRIPTION
% ci            2 element array        confidence intervals

% defaults
prec=1e-6;
confLevel=.95;
doAnimate=false;
% replace defaults by input, if any
pvpmod(varargin);

% convert df to cell for automatic expansion of parameters
df=num2cell(df);
% convert confidence level to alpha
alpha=1-confLevel;
% target p values 
pTarget=[1-alpha/2  alpha/2];

% --- error checks, assignments of function handles, etc.
% are we dealing with pdf defined only for positive abscissa values?
isPosPdf=ismember_bc(fType,{'X2','F'});
% if so...
if isPosPdf && x<0
  error('input arg ''x'' is negative but must be positive for X2 and F distributions')
end
% start index for outermost loop below, determining whether lower CI shall
% be computed or not
loopStartIx=1;
switch fType
  case 'X2'
    curPdf=@ncx2pdf;    
    curCdf=@ncx2cdf;
    curInv=@chi2inv;
    % abscissa limits for plots (if doAnimate==true): first row for lower
    % CI, second row for upper CI
    abscissLim=[0 2*x;0 5*x];
    % check: if cdf of x with noncentrality parameter 0 is less than
    % 1-alpha/2 don't even start on the lower CI because the iteration will
    % not converge (that is, there is no lower CI for given values of x and
    % df)
    if chi2cdf(x,df{:})<1-alpha/2
      % lower CI cannot be constructed as it is too close to zero - set to
      % NaN
      ci=nan;
      loopStartIx=2;
    end
    
  case 'F'
    curPdf=@ncfpdf;    
    curCdf=@ncfcdf;
    curInv=@finv;
    abscissLim=[0 2*x;0 5*x];
    % similar check as above
    if fcdf(x,df{:})<1-alpha/2
      % lower CI cannot be constructed as it is too close to zero - set to
      % NaN
      ci=nan;
      loopStartIx=2;
    end
    
  case 't'
    curPdf=@nctpdf;
    curCdf=@nctcdf;
    curInv=@tinv;
    abscissLim=x+[-4 2;-2 4]*sqrt(abs(x));
    
  otherwise
    error('illegal distribution function specified');
end

if prec>.001
  warning('results will be inaccurate - set input parameter ''prec'' to a lower value');
end

if doAnimate
  fh=figure;
  ph0=plot(x,0,'k^');
  hold on
  set(ph0,'markerfacecolor','k','markersize',6);
  ph=[];
  ti={'lower CI','upper CI'};
end

% loop twice: first lower ci (but see above), then upper ci
for iIx=loopStartIx:2
  % determine initial values: there are probably better ways of estimating
  % the limits of ncp for X2 and F pdfs than the guesses below (which work
  % best if the X2/F/t value is small)
  switch fType
    case 'X2'
      if iIx==1
        % lower CI
        ncp=x+curInv(pTarget(iIx),df{:});
      else
        % upper CI
        ncp=5*x;
      end
      
    case 'F'
      if iIx==1
        ncp=x+curInv(pTarget(iIx),df{:});
      else
        ncp=10*x;
      end
      
    case 't'
      % as a rough first approximation, assume that lower/upper limit of
      % ncp is close to corresponding percentiles of central pdfs
      if iIx==1
        ncp=x+curInv(pTarget(iIx),df{:});
      else
        ncp=x-curInv(pTarget(iIx),df{:});
      end
  end
  
  % interval of first estimates: guessed ncp enlarged by x/2 on either side
  ncp=ncp+abs(x)*[-.5 .5];
  % p values of current estimates
  p=curCdf(x,df{:},ncp);
  % deviations of p of current noncentral x pdfs from target p value
  deltaP=p-pTarget(iIx);
  nIter=1;
  if doAnimate
    ph=plotPdf(x,ncp,ph,curPdf,df,iIx,nIter,abscissLim,ti);
  end
  % while desired precision is not reached...
  while ~any(abs(deltaP)<=prec)
    if all(deltaP>0)
      % shift interval to the right by one interval length
      ncp=[ncp(2) ncp(2)+abs(diff(ncp))];
    elseif all(deltaP<0)
      % shift left by one interval length
      ncp=[ncp(1)-abs(diff(ncp)) ncp(1)];
    else
      % halve interval around mean
      ncp=mean(ncp)+.25*abs(diff(ncp))*[-1 1];
    end
    % X2 and F distributions need an extra check: the lower ncp must be >=0
    if isPosPdf
      if ncp(1)<0
        ncp(1)=0;
      end
      % if both values of ncp are zero here the uppper CI is zero, too, so
      % stop here
      if ~any(ncp)
        break
      end
    end
    % p values of current estimates
    p=curCdf(x,df{:},ncp);
    % deviations of p of current nc x pdfs from target
    deltaP=p-pTarget(iIx);
    nIter=nIter+1;
    if doAnimate
      ph=plotPdf(x,ncp,ph,curPdf,df,iIx,nIter,abscissLim,ti);
    end
  end
  % pick border which is closer to the target value
  [nada,ix]=min(abs(deltaP));
  ci(iIx)=ncp(ix);
end
% close figure
if doAnimate
  pause(1)
  close(fh)
end

% ======================== LOCAL FUNCTION =================================
function ph=plotPdf(x,ncp,ph,pdfH,df,iIx,nIter,abscissLim,ti)
% ** function ph=plotPdf(x,ncp,ph,pdfH,df,iIx,nIter,abscissLim,ti)
% If doAnimate==true, plotPdf plots x (first input arg to ncpci) and
% noncentral pdfs with the noncentrality parameter estimates of each
% iteration step 
abscissVal=linspace(abscissLim(iIx,1),abscissLim(iIx,2),200);
if isempty(ph)
  ph(1)=plot(abscissVal,pdfH(abscissVal,df{:},ncp(1)),'-');
  ph(2)=plot(abscissVal,pdfH(abscissVal,df{:},ncp(2)),'-');
  set(ph(1),'color',[.9 .3 .3]);
  set(ph(2),'color',[.3 .3 .9]);
else
  set(ph(1),'xdata',abscissVal,'ydata',pdfH(abscissVal,df{:},ncp(1)));
  set(ph(2),'xdata',abscissVal,'ydata',pdfH(abscissVal,df{:},ncp(2)));  
end
title([ti{iIx} ', iteration # ' int2str(nIter)])
% supposedly, in animation mode we would like to be able to follow the
% iterative process with our eyes, so slow things down
drawnow
pause(.1)
