function hlp_refreshEEGLABMenus(EEG,mnu_tag,state)
% refresh checkmark state of SIFT/EEGLAB menus
% this function automatically determines which operations have been
% completed (using hlp_checkeegset()) and sets the 'checked' property of
% the corresponding menu items.
%
% Author: Tim Mullen, 2013, SCCN/INC, UCSD

if nargin<2
    checks = hlp_checkeegset('supported_checks');
    mnu_tag  = [];
else
    narginchk(3,3);
    checks    = {mnu_tag};
    mnu_state = 1;
end
% 1) for each check, if check passes then tick the appropriate menu item
% (tag should be mnu_<check_name>

if isempty(mnu_tag)
    % determine which menus should be checked
    mnu_state = true(1,length(checks));
    for k=1:length(EEG)
        mnu_state = mnu_state ...
                    & cellfun(@(chk) isempty(hlp_checkeegset(EEG(k),chk)),checks);
    end
end

siftrootmnu  = findobj(findobj(0,'tag','EEGLAB'),'type','uimenu','tag','siftroot');
siftmenus    = get(siftrootmnu(1),'children');
for k=1:length(mnu_state)
    mnutag  = checks{k};
    mnu     = unique_bc(findobj(siftmenus,'tag',mnutag));
    set(mnu,'checked',fastif(mnu_state(k),'on','off'));
end