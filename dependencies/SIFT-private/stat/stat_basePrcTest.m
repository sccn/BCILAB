function Stats = stat_basePrcTest(Conn,Baseline,connmethods,tail,alpha)


% get baseline period from original dataset
baseConn = hlp_filterConns(Conn,'trange',Baseline,'method',{'time','getrange'});

if isempty(connmethods)
    % get the fieldnames of all connectivity matrices
    connmethods = hlp_getConnMethodNames(Conn);  
end

% for each method...
for m=1:length(connmethods)
    
    % ... get conn matrix
    C   = baseConn.(connmethods{m});
    sz  = size(C);
    
    % ... reshape to push time samples into a new last dimension 
    % (this is the dimension across which prctile will be taken) 
    C = reshape(C,[sz(1) sz(2) sz(3) 1 sz(4)]);
    
    % ... replicate matrix for the number of total time points
    C = repmat(C,[1 1 1 length(Conn.erWinCenterTimes) 1]);
    
    % ... C is now a distribution matrix and we can get p-values for deviation
    % from baseline distribution
    Stats.(connmethods{m}).pval = stat_surrogate_pvals(C,Conn.(connmethods{m}),tail);
    
    if exist('alpha','var')
        % also compute percentile thresholds
        Stats.(connmethods{m}).thresh = prctile(C,alpha,ndims(C));
    end
    
end