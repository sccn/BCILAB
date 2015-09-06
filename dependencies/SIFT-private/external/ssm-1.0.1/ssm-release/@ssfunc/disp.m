function disp(f, concise)

%@SSFUNC/DISP Display of State space nonlinear function.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, concise = false; end

fmask       = zeros(size(f.ssmat));
for i = 1 : size(f.horzmask, 2)
    fmask(f.vertmask(:, i), f.horzmask(:, i))  = i;
end
fdvmask     = fmask(f.ssmat.dmmask);
dvec        = f.ssmat.dvec(fdvmask == 0, :);

if concise
    fprintf(1, '\tStationary part:\n\n');
    printmat(f.ssmat.mat, [], f.ssmat.dmmask, [], fmask);
    if ~isempty(dvec)
        fprintf('\tDynamic part:\n\n');
        printmat(dvec);
    end
    if size(dvec, 2) > 1 && any(f.ssmat.dvec(fdvmask ~= 0) ~= 0)
        fprintf('\tLinear approximation part:\n\n');
        printmat(f.ssmat.dvec(fdvmask ~= 0));
    end
else
    %% Stationary and constant part %%
    line    = sprintf('Stationary and constant part (%d * %d)', size(f.ssmat.mat, 1), size(f.ssmat.mat, 2));
    fprintf(1, '\t%s\n\t%s\n\n', line, repmat('-', 1, length(line)));
    printmat(f.ssmat.mat, f.ssmat.mmask, f.ssmat.dmmask, [], fmask);

    %% Dynamic part %%
    if ~isempty(dvec)
        line    = sprintf('Dynamic part (%d * %d)', size(dvec, 1), size(dvec, 2));
        fprintf(1, '\t%s\n\t%s\n\n', line, repmat('-', 1, length(line)));
        printmat(dvec);
    end

    %% Nonlinear function part %%
    line    = sprintf('Nonlinear functions (%d)', length(f.f));
    fprintf(1, '\t%s\n\t%s\n', line, repmat('-', 1, length(line)));
    for i = 1 : length(f.f)
        fprintf(1, '\t[%d] ', i);
        if isa(f.f{i}, 'function_handle'), fprintf(1, 'Function:   %s\n', func2str(f.f{i})); else fprintf(1, 'Function:   Not set\n'); end
        if isa(f.df{i}, 'function_handle'), fprintf(1, '\t\tDerivative: %s\n', func2str(f.df{i})); else fprintf(1, '\t\tDerivative: Not set\n'); end
        if size(dvec, 2) > 1 && any(f.ssmat.dvec(fdvmask == i) ~= 0)
            fprintf(1, '\t\tLinear approximation:\n');
            printmat(f.ssmat.dvec(fdvmask == i), [], [], [], [], 2);
        end
        fprintf(1, '\n');
    end
end


