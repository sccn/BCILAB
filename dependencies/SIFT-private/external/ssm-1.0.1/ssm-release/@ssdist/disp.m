function disp(d, concise)

%@SSDIST/DISP Display of state space non-Gaussian distribution.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, concise = false; end

dmask       = zeros(size(d.ssmat));
for i = 1 : size(d.diagmask, 2)
    dmask(d.diagmask(:, i), d.diagmask(:, i))  = i;
end
ngdvmask    = dmask(d.ssmat.dmmask);
dvec        = d.ssmat.dvec(ngdvmask == 0, :);

if concise
    fprintf(1, '\tStationary part:\n\n');
    printmat(d.ssmat.mat, [], d.ssmat.dmmask, dmask);
    if ~isempty(dvec)
        fprintf('\tDynamic part:\n\n');
        printmat(dvec);
    end
    if size(dvec, 2) > 1 && any(d.ssmat.dvec(ngdvmask ~= 0) ~= 0)
        fprintf('\tGaussian approximation part:\n\n');
        printmat(d.ssmat.dvec(ngdvmask ~= 0, :));
    end
else
    %% Stationary and constant part %%
    line    = sprintf('Stationary and constant part (%d * %d)', size(d.ssmat.mat, 1), size(d.ssmat.mat, 2));
    fprintf(1, '\t%s\n\t%s\n\n', line, repmat('-', 1, length(line)));
    printmat(d.ssmat.mat, d.ssmat.mmask, d.ssmat.dmmask, dmask);

    %% Dynamic part %%
    if ~isempty(dvec)
        line    = sprintf('Dynamic part (%d * %d)', size(dvec, 1), size(dvec, 2));
        fprintf(1, '\t%s\n\t%s\n\n', line, repmat('-', 1, length(line)));
        printmat(dvec);
    end

    %% Non-Gaussian part %%
    line    = sprintf('Non-Gaussian distributions (%d)', length(d.type));
    fprintf(1, '\t%s\n\t%s\n', line, repmat('-', 1, length(line)));
    for i = 1 : length(d.type)
        fprintf(1, '\t[%d] ', i);
        if d.type(i), fprintf(1, 'Type:               additive noise\n');
        else fprintf(1, 'Type:               exponential family\n'); end
        if isa(d.matf{i}, 'function_handle'), fprintf(1, '\t\tMatrix function:    %s\n', func2str(d.matf{i})); else fprintf(1, '\t\tMatrix function:    Not set\n'); end
        if isa(d.logpf{i}, 'function_handle'), fprintf(1, '\t\tLog prob function:  %s\n', func2str(d.logpf{i})); else fprintf(1, '\t\tLog prob function:  Not set\n'); end
        if size(dvec, 2) > 1 && any(d.ssmat.dvec(ngdvmask == i) ~= 0)
            fprintf(1, '\t\tGaussian approximation:\n');
            printmat(d.ssmat.dvec(ngdvmask == i), [], [], [], [], 2);
        end
        fprintf(1, '\n');
    end
end




