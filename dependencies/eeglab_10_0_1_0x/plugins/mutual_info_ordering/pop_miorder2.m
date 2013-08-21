function [MInew ord] = pop_miorder2(EEG)

if nargin < 1
    help pop_miorder;
    return;
end;

% check whether ica data is availible
% ---------------------------------------------------
EEG.icaact = eeg_getica(EEG);
if ~isempty( EEG.icaact )
    MI = minfo(EEG.icaact);
    %MI = (MI + MI') / 2;
    MI = MI - diag(diag(MI));
    [MInew ord] = arrminf2(MI);
    maxi = max(max(MInew));
    thmax_ratio = 0.6;
    
    groups = find_blocks2(MInew,maxi*thmax_ratio,1)
    no_clusters = length(groups)/3;
    %cluster_indices = zeros(no_clusters,2);
    fprintf('Dependent ICs after block diagonalization:')
    fprintf('\n')
    if no_clusters == 0
        fprintf('None')
    else
        
        for i = 1:no_clusters
            temp_group =  groups(3*i-2:3*i-1);
            
            disp(num2str(ord(temp_group(1):temp_group(2))));
            
            fprintf('\n')
        end
    end
else
    error('You must run ICA first');
end;
