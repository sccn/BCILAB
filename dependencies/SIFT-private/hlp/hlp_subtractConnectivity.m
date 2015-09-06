function Conn = hlp_subtractConnectivity(Conn1,Conn2)
    % return a Connectivity object containing the difference in
    % connectivity between two datasets
    
    Conn = Conn1;
    
    connmethods = hlp_getConnMethodNames(Conn1);
    
    for c=1:length(connmethods)
        Conn.(connmethods{c}) = Conn1.(connmethods{c}) - Conn2.(connmethods{c});
    end
    
    otherfields = setdiff_bc(fieldnames(Conn1),connmethods);
    for f=1:length(otherfields)
        Conn.(otherfields{f}) = Conn1.(otherfields{f});
    end
    
    