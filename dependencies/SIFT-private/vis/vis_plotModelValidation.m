function handles = vis_plotModelValidation(whitestats,PCstats,stabilitystats,varargin)
% Visualize the results of model validation as computed by
% pop_est_validateMVAR()
%
% Inputs:
%
%       whitestats:     Cell array of structs each containing results of
%                       tests for whiteness of residuals for a single
%                       dataset/condition
%                       (as output by est_checkMVARWhiteness())
%                       If not available, pass in an empty matrix to skip
%                       plotting this result.
%       PCstats:        Cell array of structs each containing results of
%                       tests for consistency of the model for a single
%                       dataset/condition
%                       (as output by est_checkMVARConsistency())
%                       If not available, pass in an empty matrix to skip
%                       plotting this result.
%       stabilitystats: Cell array of structs each containing results of
%                       tests for stability of the model for a single
%                       dataset/condition
%                       (as output by est_checkMVARStability())
%                       If not available, pass in an empty matrix to skip
%                       plotting this result.
%
%
% Optional:
%
%       whitenessCriteria:   Cell array of whiteness criteria
%       checkWhiteness:      (boolean) plot whiteness results
%       checkConsistency:    (boolean) plot consistency results
%       checkStability:      (boolean) plot stability results
%       windowTimes:        Times (sec) of windows for x-axis labeling.
%                           This MUST be the same length as number of
%                           time windows (e.g.) whitestats{1}.winStartIdx
%       conditions:         Cell array of condition labels
%                           e.g. conditions = {ALLEEG.condition}
%
% Outputs:
%
%       Figures will be rendered showing model validation results for each
%       dataset. The handles to the figures are returned.
%
% See Also: pop_est_validateMVAR(),est_checkMVARWhiteness(),
%           est_checkMVARStability(), est_checkMVARConsistency()
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


% if ~isempty(varargin) && isstruct(varargin{1})
%     varargin = [hlp_struct2varargin(varargin{1}) varargin];
% end

whitenessCriteria = {};

% find the names of all whiteness tests stored in whitestats
if ~isempty(whitestats{1})
    fn = fieldnames(whitestats{1});
    for i=1:length(fn)
        if isfield(whitestats{1}.(fn{i}),'fullname')
            whitenessCriteria = [whitenessCriteria fn{i}];
        end
    end
end

num_conds = max([length(whitestats),length(PCstats),length(stabilitystats)]);


g = finputcheck(varargin, ...
    {'whitenessCriteria'   'cell'  whitenessCriteria   whitenessCriteria; ...
    'checkWhiteness',     'boolean'   []          true; ...
    'checkConsistency'    'boolean'   []          true; ...
    'checkStability'      'boolean'   []          true; ...
    'windowTimes'         'real'      []          [];   ...
    'conditions'          'cell'      {}          {}; ...
    },'mode','ignore','quiet');

if isempty(whitestats{1}),         g.checkWhiteness    = false;    end
if isempty(PCstats{1}),            g.checkConsistency  = false;    end
if isempty(stabilitystats{1}),     g.checkStability    = false;    end
if isempty(g.conditions),          g.conditions        = cell(1,num_conds); end

numrows = sum([g.checkWhiteness g.checkConsistency g.checkStability]);
numcols = 1;

for cond = 1:num_conds
    
    if isempty(g.conditions{cond})
        g.conditions{cond} = sprintf('Condition %d',cond);
    end
    
    % plot results
    handles(cond) = figure('Name',sprintf('%s - Model Validation Results',g.conditions{cond}));
    curplot=1;
    
    if g.checkWhiteness
        % Plot results of residual whiteness checks
        if ~iscell(whitestats), whitestats = {whitestats}; end
        
        ax=subplot(numrows,numcols,curplot);
        for i = 1:length(g.whitenessCriteria)
            wcstr = lower(hlp_variableize(g.whitenessCriteria{i}));
            wc = whitestats{cond}.(wcstr);
            pvals(i,:) = wc.pval;
            if size(pvals,2)>1
                lgnd{i} = sprintf('%s (%d/%d white)',wc.fullname, sum(wc.w),length(wc.w));
            else
                lgnd{i} = sprintf('%s (%swhite)',wc.fullname, fastif(wc.w,'','not '));
            end
        end
        
        if size(pvals,2)>1
            % more than one window -- make lineplot
            if ~isempty(g.windowTimes)
                abscissa = g.windowTimes;
            else
                abscissa = 1:length(whitestats{cond}.winStartIdx);
            end
            plot(ax,abscissa,pvals','Marker','.');
            xlabel(ax,fastif(isempty(g.windowTimes),'Window number','Time (sec)'));
            
            legend(ax,lgnd);
            set(ax, ...
                'Xlim',[abscissa(1)-abs(diff(abscissa(1:2))) abscissa(end)+abs(diff(abscissa(1:2)))], ...
                'Ylim',[max(0,min(pvals(:))-0.5), min(1,max(pvals(:))+0.5)]);
            axcpystr = 'legend(''''show'''');';
        else
            % single window -- create bar plots
            h = bar(ax,pvals);
            ch = get(h,'Children');
            
            set(ax,'xticklabel',lgnd);
            colors = [[1 0 0];[0 0 1];[0 1 0];[0 0 0];[1 0 1];[0 1 1]];
            colors = colors(1:length(pvals),:);
            set(ch,'FaceVertexCData',colors);
            %             set(gca,'xlim',[0 length(length(g.whitenessCriteria))+1],'Ylim',[max(0,min(pvals(:))-0.5), min(1,max(pvals(:))+0.5)]);
            axcpystr = '';
        end
        
        hl=hline(whitestats{cond}.alpha,'b','1-p_{port}',[1.01 -0.01],ax);
        set(hl,'linestyle','--','linewidth',2);
        ylabel({'Whiteness Significance'});
        
        curplot=curplot+1;
        
        if ismember_bc('acf',lower(g.whitenessCriteria))
            hl=hline(1-whitestats{cond}.alpha,'r','1-P_{acf}',[1.01 -0.01],ax);
            set(hl,'linestyle','--','linewidth',2);
        end
        
        axcopy(ax,axcpystr);
    end
    
    
    
    if g.checkConsistency
        % Plot results of consistency checks
        
        if ~iscell(PCstats), PCstats = {PCstats}; end
        
        ax=subplot(numrows,numcols,curplot);
        if length(PCstats{cond}.PC)>1
            % more than one window -- make lineplot
            meanpc = mean(PCstats{cond}.PC);
            stdpc  = std(PCstats{cond}.PC);
            if ~isempty(g.windowTimes)
                abscissa = g.windowTimes;
            else
                abscissa = 1:1:length(PCstats{cond}.winStartIdx);
            end
            plot(ax,abscissa,PCstats{cond}.PC,'Marker','.');
            xlabel(ax,fastif(isempty(g.windowTimes),'Window number','Time (sec)'));
            set(ax,...
            'Xlim',[abscissa(1)-abs(diff(abscissa(1:2))) abscissa(end)+abs(diff(abscissa(1:2)))], ...
            'Ylim',[min(PCstats{cond}.PC)-50 min(max(PCstats{cond}.PC)+50,100)]);
            % make legend
            text(0.98,0.9,sprintf('Mean PC: %0.2f+/-%0.3g%%',meanpc,stdpc), ...
                    'units','normalized','horizontalalignment','right', ...
                    'edgecolor','k','backgroundcolor','w','parent',ax);
            % make a small histogram on right side of plot
            %             axpos = get(ax,'Position');
            %             axhist = axesRelative(ax, 'Position',[1.01 0 0.1 1], 'Units','Normalized');   %axes('Position',[axpos(1)+axpos(3)+0.01 axpos(2) 0.05 axpos(4)]);
            %             hist(axhist,PCstats{cond}.PC,10);
            %             vline(mean(PCstats{cond}.PC),':r');
            %             set(axhist,'View',[90 90]);
            %             set(axhist,'xdir','rev');
        else
            % single window -- make barplot
            bar(ax,PCstats{cond}.PC);
            % make legend
            text(0.98,0.9,sprintf('(%0.2f%% Consistent)',PCstats{cond}.PC), ...
                'units','normalized','horizontalalignment','right', ...
                'edgecolor','k','backgroundcolor','w','parent',ax);
%             legend(sprintf('(%0.2f%% Consistent)',PCstats{cond}.PC));
            set(ax,'Ylim',[min(0,PCstats{cond}.PC)-1 max(PCstats{cond}.PC,100)+1]);
            set(ax,'xtick',[]);
        end
        
        ylabel('Percent Consistency');
        axcopy(ax);
        curplot = curplot+1;
    end
    
    if g.checkStability
        % Plot results of stability checks
        
        if ~iscell(stabilitystats), stabilitystats = {stabilitystats}; end
        
        % plot stability results
        ax=subplot(numrows,numcols,curplot);
        if length(stabilitystats{cond}.stability)>1
            % more than one window -- make lineplot
            %boxplot(real(lambda)');
            maxlambda = max(real(stabilitystats{cond}.lambda),[],2);
            if ~isempty(g.windowTimes)
                abscissa = g.windowTimes;
            else
                abscissa = 1:length(stabilitystats{cond}.winStartIdx);
            end
            plot(ax,abscissa,maxlambda,'Marker','.');
            xlabel(ax,fastif(isempty(g.windowTimes),'Window number','Time (sec)'));
            
            set(ax,'Xlim',[abscissa(1)-abs(diff(abscissa(1:2))) abscissa(end)+abs(diff(abscissa(1:2)))], ...
                   'Ylim',[1.2*min(maxlambda(:)) max(0.01,1.3*max(maxlambda(:)))]);
        else
            % single window -- make barplot
            maxlambda = max(real(stabilitystats{cond}.lambda(:)));
            bar(ax,maxlambda);
            set(ax,'xtick',[]);
        end
        
        
        %     set(gca,'ylim',[max(0,0.7*min(abs(lambda(:)))) max(1.3,1.3*max(abs(lambda(:))))]);
        %     axis auto
        hl=hline(0,[],[],[],ax);
        set(hl,'linestyle','--','linewidth',2);
        ylabel(ax,{'Stability Index','(should be < 0)'});
        numstable = sum(stabilitystats{cond}.stability);
        % make legend
        text(0.98,0.9,sprintf('(%d/%d stable)',numstable,length(stabilitystats{cond}.stability)), ...
                'units','normalized','horizontalalignment','right', ...
                'edgecolor','k','backgroundcolor','w','parent',ax);
%         legend(ax,sprintf('(%d/%d stable)',numstable,length(stabilitystats{cond}.stability)));
        
        axcopy(ax);
    end
end


try
    icadefs; set(handles, 'color', BACKCOLOR); 
catch
end;
