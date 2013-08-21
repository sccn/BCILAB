function [MInew ord] = pop_miorder(EEG)

if nargin < 1
    help pop_miorder;
    return;
end;

% check whether ica data is availible
% ---------------------------------------------------
EEG.icaact = eeg_getica(EEG);
if ~isempty( EEG.icaact )
    [dummy MI] = minfo(EEG.icaact);
    MI = (MI + MI') / 2;
    MI = MI - diag(diag(MI));
    [MInew ord] = arrminf2(MI)
    groups = find_blocks(MInew)
    no_clusters = length(groups)/3;
    cluster_indices = zeros(no_clusters,2);
    fprintf('Dependent ICs after block diagonalization:')
    fprintf('\n')
    if no_clusters == 0
        fprintf('None')
    else
        
        for i = 1:no_clusters
            cluster_indices(i,:) = groups(3*i-2:3*i-1);
            str = '';
            for j = 1:cluster_indices(i,2)
                strnew = num2str(ord(cluster_indices(i,1)+j-1));
                str = strcat(str, strnew ,', ');
                
            end
            disp(str)
            fprintf('\n')
        end
    end
else
    error('You must run ICA first');
end;
