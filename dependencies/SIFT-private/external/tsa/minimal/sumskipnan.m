function [o,count,SSQ,S4M] = sumskipnan(x,DIM)
% SUMSKIPNAN adds all non-NaN values. 
%
% All NaN's are skipped; NaN's are considered as missing values. 
% SUMSKIPNAN of NaN's only  gives O; and the number of valid elements is return. 
% SUMSKIPNAN is also the elementary function for calculating 
% various statistics (e.g. MEAN, STD, VAR, RMS, MEANSQ, SKEWNESS, 
% KURTOSIS, MOMENT, STATISTIC etc.) from data with missing values.  
% SUMSKIPNAN implements the DIMENSION-argument for data with missing values.
% Also the second output argument return the number of valid elements (not NaNs) 
% 
% Y = sumskipnan(x [,DIM])
% [Y,N,SSQ] = sumskipnan(x [,DIM])
% 
% DIM	dimension
% Y	resulting sum
% N	number of valid (not missing) elements
% SSQ	sum of squares
%
% the function FLAG_NANS_OCCURED() returns whether any value in x
%  is a not-a-number (NaN)
%
% features:
% - can deal with NaN's (missing values)
% - implements dimension argument. 
% - compatible with Matlab and Octave
%
% see also: FLAG_NANS_OCCURED, SUM, NANSUM, MEAN, STD, VAR, RMS, MEANSQ, 
%      SSQ, MOMENT, SKEWNESS, KURTOSIS, SEM


%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; If not, see <http://www.gnu.org/licenses/>.

%	$Id: sumskipnan.m 5599 2009-03-10 12:52:20Z schloegl $
%    	Copyright (C) 2000-2005,2009 by Alois Schloegl <a.schloegl@ieee.org>	
%       This function is part of the NaN-toolbox
%       http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/NaN/


global FLAG_NANS_OCCURED;

if nargin<2,
        DIM = [];
end;

% an efficient implementation in C of the following lines 
% could significantly increase performance 
% only one loop and only one check for isnan is needed
% An MEX-Implementation is available in sumskipnan.cpp
%
% Outline of the algorithm: 
% for { k=1,o=0,count=0; k++; k<N} 
% 	if ~isnan(i(k)) 
% 	{ 	o     += x(k);
%               count += 1;
%		tmp    = x(k)*x(k)
%		o2    += tmp;
%		o3    += tmp.*tmp;
%       }; 

if isempty(DIM),
        DIM = min(find(size(x) > 1));
        if isempty(DIM), DIM = 1; end;
end
if (DIM<1) DIM = 1; end; %% Hack, because min([])=0 for FreeMat v3.5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-float data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isa(x,'float') %%%% || ~flag_implicit_skip_nan(), %%% skip always NaN's
	x = double(x); 
	o = sum(x,DIM);
	if nargin>1
		sz = size(x);	
		count = repmat(sz(DIM),sz([1:DIM-1,DIM+1:end]));
		if nargin>2
			x = x.*x; 
			SSQ = sum(x,DIM);
			if nargin>3
				S4M = sum(x.*x,DIM);
			end; 
		end; 
	end; 	
	return; 
end; 	

%% mex and oct files expect double
x = double(x); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use Matlab-MEX function when available  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%if 0,	
try
	
	%% using sumskipnan_mex.mex

	%% !!! hack: FLAG_NANS_OCCURED is an output argument, reserve memory !!!
	if isempty(FLAG_NANS_OCCURED),
		FLAG_NANS_OCCURED = logical(0);  % default value 
	end;
	if (nargout<3),
		[o,count] = sumskipnan_mex(x,DIM,FLAG_NANS_OCCURED);
		return; 
	elseif (nargout==3),
		[o,count,SSQ] = sumskipnan_mex(x,DIM,FLAG_NANS_OCCURED);
		return; 
	elseif (nargout==4) && isreal(i),
		[o,count,SSQ,S4M] = sumskipnan_mex(x,DIM,FLAG_NANS_OCCURED);
		return; 
	end; 	
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use OCT function when available  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% no performance analysis showing a distinct advantage has been performed.  
if 0 %%if exist('OCTAVE_VERSION','builtin'),

	%% for Octave, there might be an sumskipnan_oct.oct
	try
	if isreal(x)
		%% using sumskipnan_oct.mex
		if (nargout<3),
			[o,count] = sumskipnan_oct(x,DIM);
			FLAG_NANS_OCCURED = any(count(:)<size(x,DIM));
			return; 
		elseif (nargout==3),
			[o,count,SSQ] = sumskipnan_oct(x,DIM);
			FLAG_NANS_OCCURED = any(count(:)<size(x,DIM));
			return; 
		elseif (nargout==4) && isreal(i),
			[o,count,SSQ,S4M] = sumskipnan_oct(x,DIM);
			FLAG_NANS_OCCURED = any(count(:)<size(x,DIM));
			return; 
		end; 	
	else	%% if ~isreal(x)
		if (nargout<3),
			[or,count] = sumskipnan_oct(real(x),DIM);
			[oi,counti] = sumskipnan_oct(imag(x),DIM);
			if any(count(:)-counti(:))
				error('number of NaN differ between real and imag part\n');
			end;
			o = or + i*oi;
			FLAG_NANS_OCCURED = any(count(:)<size(x,DIM));
			return; 
		elseif (nargout<4),
			[or,count,SSQ] = sumskipnan_oct(real(x),DIM);
			[oi,counti,SSQi] = sumskipnan_oct(imag(x),DIM);
			if any(count(:)-counti(:))
				error('number of NaN differ between real and imag part\n');
			end;
			o = or + oi;
			SSQ = SSQ+SSQi; 
			FLAG_NANS_OCCURED = any(count(:)<size(x,DIM));
			return; 
		end; 	
	end;	%% if 
	end; 	%% try
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fall-back in case no OCT/MEX support is avialable   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% performance tweak: some tests have shown that x*ones(:,1) is faster than sum(x,2)
FLAG	 = (length(size(x))<3); 
if FLAG, FLAG = DIM; end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% count non-NaN's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>1,
        if FLAG~=2,
                count = sum(x==x,DIM); 
        else
                count = real(x==x)*ones(size(x,2),1);
        end;
	FLAG_NANS_OCCURED = any(count(:)<size(x,DIM));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replace NaN's with zero 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(x~=x) = 0;	

if FLAG~=2, 
        o = sum(x,DIM);
else 
	o = x*ones(size(x,2),1);
end;

if nargout>2,
        x = real(x).^2 + imag(x).^2;
        if FLAG~=2,
                SSQ = sum(x,DIM);
        else
                SSQ = x*ones(size(x,2),1);
        end;
        if nargout>3,
                S4M = sumskipnan(x.^2,DIM);
        end;
end;

