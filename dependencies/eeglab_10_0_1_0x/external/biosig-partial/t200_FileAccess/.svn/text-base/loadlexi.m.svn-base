function [s,HDR]=loadlexi(filename); 
% LOADLEXI loads LEXICORE EEG data
%
%  [data,HDR]=loadlexi(filename); 
%
%

%	$Id$
%	Copyright (C) 2009 Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% BioSig is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 3 of the License, or (at your option) any later version.


warning('LOADLEXI is experimental and not well tested. ')
%% The solution is based on a single data file, with the comment:
%% "It was used to create a QEEG. It might have been collected 
%% on an older Lexicore - I don't know.   That was 4 years ago [in 2005]." 


HDR.TYPE = 'LEXICORE';
fid = fopen(filename,'r','ieee-le');
HDR.H1   = fread(fid,[1,128],'uint8')';
HDR.data = fread(fid,[24,inf],'int16')';
fclose(fid); 

[HDR.NRec]=size(HDR.data,1);
HDR.NS = 20; 
HDR.SPR = 1; 
HDR.QEEG.status = HDR.data(:,21:24); 
s = HDR.data(:,1:20); 

%% unkwown parameters
HDR.SampleRate = NaN; 
HDR.FLAG.UCAL = 1;	% data is not scaled
HDR.PhysDimCode = zeros(1,HDR.NS);
HDR.Cal = ones(1,HDR.NS);
HDR.Off = zeros(1,HDR.NS); 
HDR.Calib = sparse([HDR.Off;diag(HDR.Cal)]); 
HDR.Label = repmat({' '},HDR.NS,1);
HDR.EVENT.TYP = [];
HDR.EVENT.POS = [];
