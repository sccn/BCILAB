function needed = isneeded(curmethod,allmethods)
% check whether curmethod will need to be calculated

if iscell(curmethod) && length(curmethod)>1
    for i=1:length(curmethod)
        needed(i) = isneeded(curmethod{i},allmethods);
    end
    return;
end

% depedency graph -- each of the measures on the right depend on the
% measure on the left (dependency arrow goes from left to right
% e.g., 'Coh',{'iCoh'} means iCoh depends on Coh (Coh-->iCoh)
% If a new measure is added, make sure to update this depedency table
% (add the method descriptor to the appropriate RHS and make a LHS entry
% for the method)
dependencies = {...
                'dDTF',     {};
                'ffDTF2',   {'dDTF'};
                'ffDTF',    {};
                'nDTF',     {};
                'GGC',      {};
                'iCoh',     {};
                'Coh',      {'iCoh'};
                'S',        {'Coh','GGC'};
                'pCoh',     {'dDTF'};
                'mCoh',     {};
                'RPDC',     {};
                'GPDC',     {};
                'nPDC',     {};
                'DTF',      {'nDTF','ffDTF','ffDTF2'};
                'G',        {'mCoh','pCoh'};
                'PDC',      {'DTF','G'};
                'Rinv',     {'RPDC'}};
                
if isempty(curmethod)
    needed =0;
    return;
elseif ismember_bc(curmethod,allmethods)
    needed = 1;
    return;
else
    % recursively determine whether any of the other methods the user 
    % wants depend on curmethod
    idx=strmatch(curmethod,dependencies(:,1),'exact');
    if isempty(idx) || isempty(dependencies(idx,2)), needed=0; return; end
    children=dependencies(idx,2);
    if any(isneeded(children{1},allmethods))
        needed = 1;
        return;
    end
end

% if we get here, this particular curmethod is not needed
needed = 0;