function cl = hlp_bm_updateLights(BMopts,hlights)
% update lighting for a brainmovie
% returns handles to camlight objects

if nargin<2
    hlights = []; end

if isempty(hlights) || ~ishandle(hlights(1))
    lighting(BMopts.vars.hbrainax,BMopts.facelighting);
        
    % delete existing lights
    hlights = findobj(BMopts.vars.hbrainax,'type','light');
    delete(hlights)
    hlights = [];
end

if ~isempty(BMopts.theme)
    cl=feval(BMopts.theme.lightingfcn,BMopts.vars.hbrainax,BMopts.theme.name,hlights);
    set(cl,'tag','brain_camlight','parent',BMopts.vars.hbrainax);
else
    cl(1) = camlight(BMopts.vars.hbrainax,'left');
    cl(2) = camlight(BMopts.vars.hbrainax,'right');
    set(cl,'tag','brain_camlight','parent',BMopts.vars.hbrainax);
end