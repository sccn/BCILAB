function z = lt(x,y)
% CADA overloaded version of function LT - calls cadabinarylogical
z = cadabinarylogical(x,y,'lt');

% global ADIGATOR
% if ADIGATOR.EMPTYFLAG
%   z = cadaEmptyEval(x,y);
%   return
% end
% NUMvod  = ADIGATOR.NVAROFDIFF;
% fid     = ADIGATOR.PRINT.FID;
% PFLAG   = ADIGATOR.PRINT.FLAG;
% indent  = ADIGATOR.PRINT.INDENT;
% 
% % ----------------------------Parse Inputs------------------------------- %
% if isa(x,'cada') && isa(y,'cada')
%   % Both Inputs are Symbolic
%   xMrow = x.func.size(1); xNcol = x.func.size(2);
%   yMrow = y.func.size(1); yNcol = y.func.size(2);
% elseif isa(x,'cada')
%   % y is numeric input
%   xMrow = x.func.size(1); xNcol = x.func.size(2);
%   [yMrow,yNcol] = size(y);
%   ytemp.id = [];
%   ytemp.func = struct('name',[],'size',[yMrow yNcol],'zerolocs',[],...
%     'value',y);
%   if PFLAG
%     if yMrow*yNcol == 1
%       ytemp.func.name = num2str(y,16);
%     else
%       ytemp.func.name = cadamatprint(y);
%     end
%   end
%   ytemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
%   y = ytemp;
% else
%   % x is numeric input
%   yMrow = y.func.size(1); yNcol = y.func.size(2);
%   [xMrow,xNcol] = size(x);
%   xtemp.id = [];
%   xtemp.func = struct('name',[],'size',[xMrow xNcol],'zerolocs',[],...
%     'value',x);
%   if PFLAG
%     if xMrow*xNcol == 1
%       xtemp.func.name = num2str(x,16);
%     else
%       xtemp.func.name = cadamatprint(x);
%     end
%   end
%   xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
%   x = xtemp;
% end
% 
% % ----------------------------Function Sizing---------------------------- %
% if (xMrow == yMrow && xNcol == yNcol)
%   FMrow = yMrow; FNcol = yNcol;
% elseif (xMrow == 1 && xNcol == 1)
%   FMrow = yMrow; FNcol = yNcol;
% elseif (yMrow == 1 && yNcol == 1)
%   FMrow = xMrow; FNcol = xNcol;
% elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
%   % In printing run of a FOR loop - Sizes could be off due to IF statement
%   if xMrow <= yMrow && xNcol <= yNcol
%     y = cadaRemoveRowsCols(y,[xMrow xNcol]);
%     yMrow = xMrow; yNcol = xNcol;
%   elseif yMrow <= xMrow && yNcol <= xNcol
%     x = cadaRemoveRowsCols(x,[yMrow yNcol]);
%     xMrow = yMrow; xNcol = yNcol;
%   else
%     error('Inputs are not of compatible sizes');
%   end
%   FMrow = xMrow; FNcol = xNcol;
% else
%   error('Inputs are not of compatible sizes');
% end
% 
% % ----------------------Build Function Properties--------------------------
% z.id = ADIGATOR.VARINFO.COUNT;
% [funcstr,~] = cadafuncname();
% z.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
%   'value',[]);
% 
% %----------------------Function Numeric Values (if exist)-----------------%
% if ~isempty(x.func.value) && ~isempty(y.func.value)
%   z.func.value = x.func.value < y.func.value;
% end
% 
% % ----------------------------Function Printing ------------------------- %
% if PFLAG == 1
%   fprintf(fid,[indent,funcstr,' = ',x.func.name,' < ',y.func.name,';\n']);
% end
% 
% 
% % ------------------------Build Derivative Properties----------------------
% z.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
% 
% ADIGATOR.VARINFO.LASTOCC([x.id y.id z.id],1) = ADIGATOR.VARINFO.COUNT;
% z = class(z,'cada');
% ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
