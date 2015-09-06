function EEG = hlp_averageConnectivity(ALLEEG)
% average connectivity across datasets stored in ALLEEG(:).CAT.Conn
% all Conn objects must have identical fields and
% all connectivity matrices must be of identical dimension
% all fields/content other than Conn are copied from ALLEEG(1) to EEG
%
% Author: Tim Mullen, SCCN/INC/UCSD 2012

EEG = ALLEEG(1);

num_datasets = length(ALLEEG);

connmethods = hlp_getConnMethodNames(EEG.CAT.Conn);

for i=1:length(connmethods)
    % for each connectivity method ...
    m = connmethods{i};
    
    % ... compute average connectivity across datasets
    for k = 2:num_datasets
        EEG.CAT.Conn.(m) = EEG.CAT.Conn.(m) + ALLEEG(k).CAT.Conn.(m);
    end
    EEG.CAT.Conn.(m) = EEG.CAT.Conn.(m)/num_datasets;
    
end

