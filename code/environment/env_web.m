function env_web(url)
% Replacement for the 'web' command that also works in deployed mode.
% env_web(URL)
%
% In:
%   URL : the web page or file to display
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2012-04-15

try
    if ~isdeployed || hlp_matlab_version < 714        
        web(url);
    else
        % web doesn't work in deployed mode in more recent MATLABs
        if ispc
            system(['start ' url]);
        elseif ismac
            system(['open ' url]);
        elseif isunix
            system(['xdg-open ' url]);
        else
            error('Unsupported OS.');
        end
    end
catch e
    disp(['Cannot open browser on your system; reason: ' e.message]);
end
