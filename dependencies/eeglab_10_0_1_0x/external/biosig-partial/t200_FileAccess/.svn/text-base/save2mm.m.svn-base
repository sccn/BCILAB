function [HDR] = save2mm(fn,MM,montage);
% SAVE2MM  saves Matrix into MatrixMarket format 
%
%       HDR = save2mm(filename,M,comment);  
%
% filename	destination file 
% M		Matrix
%
% see also: SLOAD, getMontage, regress_eog, save2gdf
%

%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.


% 	$Id$
%	Copyright (C) 2010 by Alois Schloegl <a.schloegl@ieee.org>		
%       This file is part of the biosig project http://biosig.sf.net/



        [I,J,V] = find(MM); 
        fid = fopen(fn,'w+'); 
        fprintf(fid,'%%%%MatrixMarket matrix coordinate real general\n');
        fprintf(fid,'%% generated on %04i-%02i-%02i %02i:%02i:%02.0f\n',clock);

        if ischar(montage) m = montage; else m = '? (user specified)'; end;  
        fprintf(fid,'%% Spatial Filter for %s \n',m);
        fprintf(fid,'%i %i %i\n',size(MM),length(V));

        for k = 1:length(V),
                fprintf(fid,'%2i %2i %f\n',I(k),J(k),V(k));
        end;
        fclose(fid);        
        
        HDR.Calib = MM; 
        HDR.FileName = fn; 
        HDR.TYPE = 'MatrixMarket'; 