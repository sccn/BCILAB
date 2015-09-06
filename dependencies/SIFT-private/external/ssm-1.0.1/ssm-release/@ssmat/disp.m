function disp(A, concise)

%@SSMAT/DISP Display of state space matrix.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, concise = false; end

if concise
    if isempty(A.dmmask)
        printmat(A.mat);
    else
        fprintf(1, '\tStationary part:\n\n');
        printmat(A.mat, [], A.dmmask);
        fprintf(1, '\tDynamic part:\n\n');
        printmat(A.dvec);
    end
else
    %% Stationary part %%
    line    = sprintf('Stationary part (%d * %d)', size(A.mat, 1), size(A.mat, 2));
    fprintf(1, '\t%s\n\t%s\n\n', line, repmat('-', 1, length(line)));
    printmat(A.mat, A.mmask, A.dmmask);

    %% Dynamic part %%
    if ~isempty(A.dmmask)
        line    = sprintf('Dynamic part (%d * %d)', size(A.dvec, 1), size(A.dvec, 2));
        fprintf(1, '\t%s\n\t%s\n\n', line, repmat('-', 1, length(line)));
        printmat(A.dvec);
    end
end

