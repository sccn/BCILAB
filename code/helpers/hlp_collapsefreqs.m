function bp = hlp_collapsefreqs(spectrum,method,fidx)
% collapse spectrum over frequencies
% Tim Mullen, SCCN/INC, UCSD 2013

switch method
    case 'sum'
        bp = sum(spectrum(:,fidx),2);
    case 'integrate'
        bp = trapz(fidx(2)-fidx(1),spectrum(:,fidx),2);
    case 'mean'
        bp = mean(spectrum(:,fidx),2);
    case 'max'
        bp = max(spectrum(:,fidx),[],2);
    otherwise
        error('Unsupported collapse method: %s',method);
end