% prevent warnings caused by eeglab calling "clear functions"
warning off MATLAB:lang:cannotClearExecutingFunction
eeglab; 

if ~isdeployed
    close(gcf); end
