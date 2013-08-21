function showclust(spikes, useassigns, show)
%    SHOWCLUST  temporary utility to show clusters
%       SHOWCLUST(SPIKES, [USEASSIGNS], [SHOW]);

if (nargin < 2)
	useassigns = spikes.hierarchy.assigns;
end
if (nargin < 3)
    show = unique(useassigns);
end
show = reshape(show, 1, []);

smallWindow = 0.01;

tmin = (size(spikes.waveforms,2) - spikes.threshT + 1)./spikes.Fs;   
if (isfield(spikes,'options') && isfield(spikes.options, 'refractory_period'))
    tref = spikes.options.refractory_period;
else
    tref = max(0.002, tmin*1.5);
end
tminl = tmin - 1/spikes.Fs;

clf;
ylims = [min(spikes.waveforms(:)) max(spikes.waveforms(:))];
tlims = [0 max(spikes.spiketimes)];
cmap = colormap;
for clust = 1:length(show)
    members = find(useassigns == show(clust));
    memberwaves = spikes.waveforms(members,:);
    membertimes = sort(spikes.spiketimes(members));
    subplot(length(show),3, 3 * (clust-1) + 1);
    [n,x,y] = histxt(memberwaves);
    imagesc(x,y,n); axis xy;
    if (clust < length(show))
        set(gca,'XTickLabel',{});
    end
    set(gca, 'YLim', ylims, 'YTickLabel', {}, 'Color', cmap(1,:));
	
    if (show(clust) ~= 0),  clustname = ['Cluster# ' num2str(show(clust))];
    else                    clustname = 'Outliers';
    end
	hy = ylabel({clustname, ['N = ' num2str(size(members,1))]});
	
    subplot(length(show),3,3 * (clust-1) + 2);
    [a, scores] = isiQuality(membertimes, membertimes, tmin, smallWindow, tref, spikes.Fs);
    isis = sort(diff(membertimes));   isis = isis(isis <= smallWindow);
	isis = round(isis*spikes.Fs)/spikes.Fs;
	smalltimes = linspace(0,smallWindow,smallWindow*spikes.Fs+1);
	if (~isempty(isis)), n = histc(isis,smalltimes);  else  n = zeros(length(smalltimes));  end;
    plot(smalltimes,n);    ylim = get(gca, 'YLim');
	patch([0 tref  tref  0]', [0 0 ylim(2) ylim(2)]', -[0.1 0.1 0.1 0.1], [0.8 0.8 0.8], 'EdgeColor', 'none');
	patch([0 tminl tminl 0]', [0 0 ylim(2) ylim(2)]', -[0.1 0.1 0.1 0.1], [0.6 0.6 0.6], 'EdgeColor', 'none');
	set(gca,'Xlim',[0 0.01]);
    if (clust < length(show)),  set(gca,'XTickLabel',{});
    else                        xlabel('ISI (sec)');
    end
    
    subplot(length(show),3,3 * (clust-1) + 3);
    [n,x] = hist(membertimes,1:10:max(tlims));   bar(x,n,1.0);  shading flat;
    set(gca,'Xlim',tlims);
    if (clust < length(show)),  set(gca,'XTickLabel',{});
    else                        xlabel('t (sec)');
    end
    if (clust == 1), title('Spike Times');  end;
	
	% shift text to make more readable
	set(hy, 'Units', 'char');
	hy_pos = get(hy, 'Position');    hy_pos = hy_pos + [-6*rem(clust,2),0,0];   set(hy, 'Position', hy_pos);
end
