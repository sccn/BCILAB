function s = firstcaps(s)
% Capitalize the first letter of a string
% Author: Tim Mullen, SCCN/INC/UCSD, 2011
    
    % find the first non-whitespace char...
    fc   = regexp(s,'\S','once');
    %...and make upper-case
    s(fc) = upper(s(fc));
end
