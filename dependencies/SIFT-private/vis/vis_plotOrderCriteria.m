function handles = vis_plotOrderCriteria(varargin)
% Visualize the results of model order selection as computed by
% est_selModelOrder()
% 
% Input: 
%
%       IC:             Cell array of structures containing results of
%                       model order selection for one or more datasets 
%                       (see pop_est_selModelOrder())
%
% Optional: 
%       conditions:     Cell array of strings containing names of
%                       conditions for figure labeling
%                       Default: {}
%       icselector:     Cell array of string denoting which information
%                       criteria to plot (must be fields of IC{.}).
%                       Default: all selectors in IC structure
%       minimizer:      'min' - find model order that minimizes info crit.
%                       'elbow' - find model order corresponding to elbow
%                                 of info crit. plot
%                       Default: 'min'
%       prclim:         upper limit for percentile
%                       Default: 90
%                       
% Output:
%
%       Figures will be rendered showing model selection results for each
%       dataset. The handles to the figures are returned.
%
% See Also: est_selModelOrder(), pop_est_selModelOrder()
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% 
% Author: Tim Mullen, 2010, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

g = arg_define([0 1],varargin, ...
    arg_norep({'IC','InformationCriterion'},mandatory,[],'Informaton criteria object. This can also be a cell array of objects.'), ...
    arg({'conditions','TitleString'},[],[],'Figure title strings. If g.IC is a cell array, then this must be a cell array of strings.'), ...
    arg({'icselector','InformationCriteriaToPlot'},[],[],'Information criteria to plot. This is a cell array of strings of the information criteria to plot'), ...
    arg({'minimizer','OptimalModelSelectionMethod','optimalModelSelection'},{'min'},{'elbow','min'},'Method for determining optimal model order. If "min", then the optimal model order is the one that minimizes the information criterion. If "elbow" then the optimal order is the elbow of the function of informaton criterion versus model order.','type','logical'), ...
    arg({'prclim','PercentileLimits'},90,[1 100],'Upper percentile limit for order selection. If PercentileLimits = L, This places a marker at the model order, p, for which L% of all sampled windows indicate an optimal model order of p or lower.') ...
    );

% commit IC variable to workspace
[data g] = hlp_splitstruct(g,{'IC'});
arg_toworkspace(data);
clear data;

% do some input checking
if ~iscell(IC)
    IC = {IC};
end
if ~iscell(g.conditions)
    g.conditions = {g.conditions};
end
if length(g.conditions)<length(IC)
    error('You must specify title strings for all conditions');
end
if isempty(g.icselector)
    g.icselector = IC{1}.selector;
end
if ~iscell(g.minimizer)
    g.minimizer = {g.minimizer};
end

% recursively plot a separate figure for each minimizer
if length(g.minimizer)>1
    for k=1:length(g.minimizer)
        try, g = rmfield(g,'report_args'); catch, end;
        handles{k} = vis_plotOrderCriteria(IC,g,'minimizer',g.minimizer{k});
    end
    return;
end

for cond = 1:length(IC)
    
    if isempty(g.conditions{cond})
        g.conditions{cond} = sprintf('Condition %d',cond);
    end
    
    numinfocriteria = length(g.icselector);
    numWins = size(IC{cond}.(g.icselector{1}).ic,2);


    % plot the results
    % ----------------------
    handles(cond) = figure('name',sprintf('%s - Model Order Selection Results (%s)',g.conditions{cond},fastif(strcmpi(g.minimizer{1},'min'),'min ic','elbow ic')));
    for i=1:numinfocriteria
        allinfocrit(i,:) = mean(IC{cond}.(lower(g.icselector{i})).ic,2);
    end

    if numWins>1
        numrows = floor(sqrt(numinfocriteria))+1;
        numcols = ceil(numinfocriteria/(numrows-1));
    else
        numrows = 1;
        numcols = 1;
    end

    set(gca,'Position',[0 0 1 1]);

    % plot information criteria
    subplot(numrows,numcols,1:numcols); 
    plot(IC{cond}.pmin:IC{cond}.pmax,fastif(diff(size(allinfocrit))==0,allinfocrit',allinfocrit),'linewidth',2);
    set(gca,'xtick',IC{cond}.pmin:IC{cond}.pmax);
    axis auto
    xlabel('model order    ','fontsize',12);
    ylabel('Information criteria (bits)    ','fontsize',12);
    set(gca,'xgrid','on');
    title({'Mean Info. Criteria across sampled windows   ',['Optimal order determined by ' g.minimizer{1} ' of mean curve   ']},'fontsize',12);
    axcopy(gca);

    hold on

    defcolororder = get(0,'defaultaxescolororder');
    
    % make markers at minima
    for i=1:numinfocriteria
        sel = g.icselector{i};
        
        switch g.minimizer{1}
            case 'min'
                [minic popt] = min(allinfocrit(i,:)); %ceil(mean(IC{cond}.(lower(sel)).popt));
            case 'elbow'
                [minic popt] = hlp_findElbow(allinfocrit(i,:));
        end
        
%         minic = allinfocrit(i,popt);
        plot(popt+IC{cond}.pmin-1,minic,'rs','markersize',10,'MarkerEdgeColor','k','markerfacecolor',defcolororder(i,:));
        lmin = vline(popt+IC{cond}.pmin-1);
        set(lmin,'color',defcolororder(i,:),'linestyle','--','linewidth',2);
        legendstr{i} = sprintf('%s (%d)',g.icselector{i},popt+IC{cond}.pmin-1);
    end
    
%     for i=1:numinfocriteria
%         sel = g.icselector{i};
%         popt = ceil(mean(IC{cond}.(lower(sel)).popt));
%         minic = allinfocrit(i,popt-IC{cond}.pmin+1);
%         plot(popt,minic,'rs','markersize',10,'MarkerEdgeColor','k','markerfacecolor',defcolororder(i,:));
%         lmin = vline(popt);
%         set(lmin,'color',defcolororder(i,:),'linestyle','--','linewidth',2);
%         legendstr{i} = sprintf('%s (%d)',g.icselector{i},popt);
%     end

    xlim([IC{cond}.pmin IC{cond}.pmax]);
    yl=ylim;
    ylim([yl(1)-0.001*(yl(2)-yl(1)) yl(2)]);
    
    legend(legendstr,'linewidth',2);
    
    % popt = popt+IC{cond}.pmin-1;
    % plot(popt',minic);

    % popt = popt(:);
    % ylim = get(gca,'Ylim');
    % ylim = ylim(ones(1,numinfocriteria),:);
    % line([popt popt],ylim);
    
    hold off
    % 


    if numWins>1
        % plot histograms of optimal model order across windows
        for i=1:numinfocriteria
            ax=subplot(numrows,numcols,i+numcols);
            xScale = IC{cond}.pmin:IC{cond}.pmax;
            
            switch g.minimizer{1}
                case 'min'
                    bar(xScale,histc(IC{cond}.(lower(g.icselector{i})).popt,xScale),'k');         %  plot histogram
                    popt = round(mean(IC{cond}.(lower(g.icselector{i})).popt));                    %  mean
                    poptstd = round(std(IC{cond}.(lower(g.icselector{i})).popt));                  %  stdev
                    poptprctile = round(prctile(IC{cond}.(lower(g.icselector{i})).popt,g.prclim));   %  upper 95th prctile
                case 'elbow'
                    bar(xScale,histc(IC{cond}.(lower(g.icselector{i})).pelbow,xScale),'k');
                    popt = round(mean(IC{cond}.(lower(g.icselector{i})).pelbow));
                    poptstd = round(std(IC{cond}.(lower(g.icselector{i})).pelbow));
                    poptprctile = round(prctile(IC{cond}.(lower(g.icselector{i})).pelbow,g.prclim));
            end
            axes(ax);
            
            % shade region for stdev
            [patchHandles textHandles] = hlp_vrect([popt-poptstd popt+poptstd], ...
                'patchProperties', ...
                {'FaceColor',defcolororder(i,:), ...
                 'FaceAlpha',0.2,...
                 'EdgeColor',[0.2 0.2 0.2],...
                 'EdgeAlpha',0.2 ...
                 });
             
            % mark mean popt
            lmin = vline(popt,'-',[num2str(popt) '+-' num2str(poptstd)],[-0.01 0.95],gca,defcolororder(i,:));
            set(lmin,'color',defcolororder(i,:),'linestyle','-','linewidth',2);
            
            % mark upper prctile
            lprctile = vline(poptprctile,':',sprintf('%g%%',g.prclim),[-0.01 0.05],gca,defcolororder(i,:));
            set(lprctile,'color',defcolororder(i,:),'linewidth',2);
            
            title([g.icselector{i} '  '],'fontsize',12,'fontweight','bold','color',defcolororder(i,:));
            axcopy(gca);
            xlabel('opt. model order   ','fontsize',12);
            ylabel('histogram count   ','fontsize',12);
            xlim([IC{cond}.pmin-1 IC{cond}.pmax+1]);
            yl=ylim;
            ylim([yl(1)-0.001*(yl(2)-yl(1)) yl(2)]);
        end
    %     axis on
        %set(gca,'yaxislocation','right');
    end


    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
    
end