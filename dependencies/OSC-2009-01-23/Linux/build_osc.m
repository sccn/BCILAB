if isunix && ~exist('osc_recv')
	% build OSC if necessary...
    try
        osc_make; 
    catch
        output = evalc('try,osc_make;catch,end');
        if ~isempty(strfind(r,'cannot find -llo'))
            disp('Could not compile OSC (OpenSoundControl) external interface library for your machine.');
            disp('For this, you need install (perhaps using your package manager) the liblo library,');
            disp('or build and install it from dependencies/OSC-2009-01-23/Linux/liblo-0.26');
            error('Cannot compile OSC');
        end
    end
end
