function [INI,S,E] = detectmuscle(S, iter, Mode)
% Muscle detection with inverse filtering
% Artifacts are indicated with NaN's. 
%
% [INI,S,E] = detectmuscle(S [, iter [,1]])% [INI,S,E] = detectmuscle(S, Fs, Mode)  with Mode>1% [INI,S,E] = detectmuscle(S, arg2, Mode)%
% iter		number of iterations [default:1]
% INI.MU	mean of S
% INI.InvFilter	coefficients of inverse filter 
% S		outlier replaced by NaN
% E		isnan(E) indicates muscle artifact
% Mode	1: [default] inverse filtering
%	2: based on gradient, range and amplitude [BrainVision method]
% 	3: slope > 11 uV/sample [1]% 	4: beta2 > 0.9 uV^2/Hz [1] %
%
% References: 
% [1]	Van de Velde, M., Van Erp, G., Cluitmans, P., 1998. 
%	Detection of muscle artefact in the normal human awake EEG. 
% 	Electroencephalography and Clinical Neurophysiology 107 (2), 149-158.
%

%	$Id: detectmuscle.m 2202 2009-10-27 12:06:45Z schloegl $%	Copyright (C) 2003,2008,2009 by Alois Schloegl <a.schloegl@ieee.org>	%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
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


if nargin<2,
        iter=[];
end;
INI = [];

% Muscle Detection
if Mode==1, 
	%% TODO: Validation 
	if isempty(iter) iter = 1; end; 
	%% inverse filter   
	TH = 5; 
	[se,INI.MU] = sem(S);
	if TH*se<abs(INI.MU),
		[S] = center(S);
	end;

	while iter>0,
        	[AR,RC,PE]= lattice(S',10);
		INI.InvFilter = ar2poly(AR);
        	E = zeros(size(S));
        	for k = 1:size(S,2),
                	E(:,k) = filter(INI.InvFilter(k,:),1,S(:,k));
	        end;
        	INI.V  = std(E);
	        INI.TH = INI.V * TH; 
        	iter   = iter-1;
        
	        for k = 1:size(S,2),
        	        S(E(:,k)>(INI.TH(k)) | E(:,k)<(-INI.TH(k)),k) = NaN;
	        end;
	end; 
	% the following part demonstrates a possible correction 	for k = 1:size(S,2),
                E(:,k) = filter(INI.InvFilter(k,:),1,S(:,k));
	end;
	% isnan(E), returns Artifactselse 

elseif Mode==2,
	Fs = iter; 
	% Criterion for bad gradient:
	% Maximal allowed voltage step / sampling point: 100.00 µV
	% Mark as bad before event: 100.00 ms
	% Mark as bad after event: 100.00 ms
	ix1 = [abs(diff(S,[],1))>100;zeros(1,size(S,2))];
	
	% Criterion for bad max-min:
	% Maximal allowed absolute difference: 150.00 µV
	% Interval Length: 100.00 ms
	% Mark as bad before event: 100.00 ms
	% Mark as bad after event: 100.00 ms
	
	ix2 = abs(filter([-1;zeros(Fs/10-1,1);1],1,S))>150;
	ix2 = [ix2(Fs/20+1:end,:);zeros(Fs/20+1,size(S,2))]; 

	% Criterion for bad amplitude:
	% Minimal allowed amplitude: -75.00 µV
	% Maximal allowed amplitude: 75.00 µV
	% Mark as bad before event: 100.00 ms
	% Mark as bad after event: 100.00 ms
	ix3 = abs(S)>75;
	
	ix = filtfilt(ones(Fs/10,1),1,double(ix1|ix2|ix3))>0;
	S(ix) = NaN;	
	INI=[]; 


elseif Mode==3,
	Fs = iter; 
	E  = filter(ones(1,Fs),Fs,abs(diff(S))/Fs*250);
	Threshold = 11; % uV/sample [1]
	S(E>Threshold) = NaN;


elseif Mode==4,
	%% TODO: free butter 
	Fs = iter; 
	if Fs>= 250, 
		[A,B]=butter(5,25*2/Fs,'high');
	else 	
		[A,B]=butter(5,[25,125]*2/Fs);
	end; 
	E = filtfilt(ones(1,Fs),Fs,filtfilt(A,B,S).^2);
	Threshold = 0.9; % uV^2/Hz ??? [1]
 	S(E>Threshold) = NaN; 

end; 

if nargout<3, return; end;

