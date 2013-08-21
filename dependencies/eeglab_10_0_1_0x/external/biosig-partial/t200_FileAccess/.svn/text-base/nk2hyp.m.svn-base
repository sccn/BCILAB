function [hyp,g1,g2] = nk2hyp(fn)
% NK2HYP extracts the hypnogram from NihonKohden data 
%
% [hyp,g] = nk2hyp(fn)
%
% fn 	filename
% hyp	hypnogram in 'W1234R' encoding
% g	hypnogram in 0=W,1,2,3,4,5=R encoding


%	$Id$
%	Copyright (C) 2010 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


[p,f,e] = fileparts(fn); 

fn  = fullfile(p,[f,'.dat']); 
fid = fopen(fn,'r');
dat = fread(fid,inf,'uint8'); 
fclose(fid); 

pos = length(dat);
ix  = [];
nn  = [];
while all(dat(pos+[-1,0])==255)
	ix  = [pos,ix];
	pos = pos-62;
	nn  = [[1,256]*dat(pos+[1:2]),nn];
end; 
if (nn(1)~=1) || any(diff(nn)~=1)
	% sanity check
	error('nk2hyp failed');
end; 
g2 = dat(ix-3)'; 
g1 = dat(ix-2)'; 
clear dat;

g2(g2=='W')=0;
g2(g2=='R')=5;
g1(g1=='W')=0;
g1(g1=='R')=5;

H = 'W12345R';
hyp = H(g2+1); 
