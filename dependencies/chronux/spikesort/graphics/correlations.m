function correlations(spikes, useassigns, show)
%    CORRELATIONS  another temporary utility to show clusters
%       CORRELATIONS(SPIKES, [USEASSIGNS], [SHOW]);

if (nargin < 2),
    if (isfield(spikes.hierarchy, 'assigns')),  useassigns = spikes.hierarchy.assigns;
    elseif (isfield(spikes.overcluster, 'assigns')), useassigns = spikes.overcluster.assigns;
    else useassigns = ones(size(spikes.waveforms,1),1);
    end
end
if (nargin < 3)
    list = unique(useassigns);
    show = list(1:min(5,length(list)));
end
    
K = length(show);
maxlag = 0.025;  % (msec);

for r = 1:K
	for c = r:K  % upper right triangle only (b/c symmetric)
		subplot(K,K,c+(r-1)*K);
		selectrow = find(useassigns == show(r));
		selectcol = find(useassigns == show(c));
		if ((length(selectrow) > 1) && (length(selectcol) > 1))
			[cross,lags] = pxcorr(spikes.spiketimes(selectrow), spikes.spiketimes(selectcol), 1000, maxlag);
			if (r == c),  cross(lags == 0) = 0;  end;  % blank out autocorr peak
			bar(lags,cross,1.0);  shading flat;
			set(gca, 'XLim', [-maxlag, maxlag]);
		else
			cla;  % show blank if <= 1 pts
		end
		if (r == c),  ylabel(sprintf('Cluster #%d', show(r)));  end;
		if (r == 1),  title(sprintf('Cluster #%d', show(c)));  end;
	end
end
