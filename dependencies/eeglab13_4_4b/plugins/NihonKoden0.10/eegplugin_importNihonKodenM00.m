% eegplugin_importNihonKodenM00() - EEGLAB plugin for importing Nihon Koden
%                                   .m00 text-format data.
%
% Usage:
%   >> eegplugin_nihonkodenM00import(fig, trystrs, catchstrs);
%
% Author: Makoto Miyakoshi SCCN, INC, UCSD
%
% See also: pop_importNihonKodenM00() importNihonKodenM00()

% History:
% 07/17/2013 ver 1.0 by Makoto. Created for my collaboration with Eishi Asano.

% Copyright (C) 2013 Makoto Miyakoshi, SCCN, INC, UCSD.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function eegplugin_nihonkodenM00import(fig, trystrs, catchstrs)
    
    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'import data');
    
    % menu callbacks
    % --------------
    comcnt = [ trystrs.no_check '[EEG LASTCOM] = pop_importNihonKodenM00('''');' catchstrs.new_and_hist ];
    
    % create menus
    % ------------
    uimenu( menu, 'label', 'From Nihon Koden .m00 text files', 'callback', comcnt);