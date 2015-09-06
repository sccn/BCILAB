function Conn = stat_threshConn(Conn,Stat,ThreshMethod,ThreshValue,remNonStatMethods)
% threshold connectivity values according to Stats structure
% Author: Tim Mullen, 2011

if nargin<3
    error('You must supply at least 3 arguments');
end

if nargin<4 && strcmpi(ThreshMethod,'pval')
    error('You must supply a p-value threshold');
end

if nargin<5
    remNonStatMethods = false;
end

statconnmethods = hlp_getConnMethodNames(Stat);
connmethods     = hlp_getConnMethodNames(Conn);

methodsToSkip   = ~ismember_bc(connmethods,statconnmethods);

if remNonStatMethods
    % remove all methods that we don't have stats for
    Conn = rmfield(Conn,connmethods(methodsToSkip));
end

for m=1:length(connmethods)
    if methodsToSkip(m)
        % skip thresholding for methods
        % that don't have associated Stats
        continue;
    end
    
    switch ThreshMethod
        case 'pval'
            % zero out any connections with p-value greater than
            % threshold
            Conn.(connmethods{m})(Stat.(connmethods{m}).(ThreshMethod) > ThreshValue) = 0;
        case 'thresh'
            % zero out any connections with amplitude less than
            % stat threshold
            
            Conn.(connmethods{m})( Conn.(connmethods{m}) < (Stat.(connmethods{m}).(ThreshMethod)) ) = 0;
            
        otherwise
            error('ThreshMethod must be ''pval'' or ''thresh''');
    end
    
end
