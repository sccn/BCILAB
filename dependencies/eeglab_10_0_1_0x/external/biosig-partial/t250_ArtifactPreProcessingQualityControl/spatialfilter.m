function [W,Label] = spatialfilter(arg1,Mode,s2)
% Spatial filter provides different spatial filters. 
%
% W = spatialfilter(NS,...)	
% W = spatialfilter(s,...)	
% W = spatialfilter(HDR,...)% W = spatialfilter([.], Mode)% W = spatialfilter(s,'CSP',s2)	
% 
% NS	number of channels
% s	data matrix (one column per channel)	% Mode	'Mono': monopolar
%	'CAR':	common average reference
%	'Laplace': Hjorth's Laplacian
%	'bipolar': all posible combinations of bipolar channels 
%	'PCA': Principle Component Analysis 
%	'ICA': Independent Component Analysis
%	'CSP': Common Spatial patterns
%	'ALL': Combination of all 

%	$Id: spatialfilter.m 2202 2009-10-27 12:06:45Z schloegl $%	Copyright (C) 2008 by Alois Schloegl <a.schloegl@ieee.org>	%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
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

s = []; 
if any(size(arg1)>1)
	s = arg1; 
	NS = size(arg1,2); 
elseif isscalar(arg1) 
	NS = arg1; 
elseif isstruct(arg1) && isfield(arg1,'NS')
	NS = arg1.NS; 
end; 

if strcmpi(Mode,'bipolar')
	W = sparse(NS,NS*(NS-1)/2); 
	n = 0; 
	for k1 = 1:NS,
	for k2 = 1+k1:NS,
		n = n+1;
		W([k1,k2],n) = [1,-1];
		Label{n} = sprintf('Bip #%i-#%i',k1,k2);
	end;
	end;

elseif strcmpi(Mode,'Mono')
	W = speye(NS); 
	Label = cellstr(int2str([1:NS]'));
	for k=1:NS,Label{k} = ['Mono #',int2str(k)]; end;
	  
elseif strcmpi(Mode,'CAR')
	W = [eye(NS) - 1/NS, ones(NS,1)/NS]; 
	Label = cellstr(int2str([1:NS]'));
	for k=1:NS,Label{k} = ['CAR #',int2str(k)]; end;  
	Label{NS+1} = 'CAR REF'; 
	
elseif strcmpi(Mode,'ICA')
	W = []; 
	Label = {};
	warning('ICA not supported yet');
	%for k=1:size(W,2),Label{k} = ['ICA #',int2str(k)]; end;  

elseif strcmpi(Mode,'PCA')
	if isempty(s)
		error('PCA requires data matrix s')		
	else
		[W,D] = eig(s'*s);
		for k=1:size(W,2),Label{k} = ['PCA #',int2str(k)]; end;  
	end; 	 

elseif strcmpi(Mode,'CSP')
	if isempty(s) || nargin<3
		error('CSP requires data matrix s and s2')		
	else
		[W,D] = eig(covm(s,'E'),covm(s2,'E'));
		for k=1:size(W,2),Label{k} = ['CSP #',int2str(k)]; end;  
	end; 	 

elseif strcmpi(Mode,'ALL')
	[W1,L1] = spatialfilter(NS,'Mono');
	[W2,L2] = spatialfilter(NS,'bipolar');
	[W3,L3] = spatialfilter(NS,'CAR');
	W = [W1,W2,W3];
	Label = [L1(:); L2(:); L3(:)];

	if ~isempty(s)
		[W1,L1] = spatialfilter(s,'PCA');
		[W2,L2] = spatialfilter(s,'ICA');
		W     = [W,W1,W2];
		Label = [Label; L1(:); L2(:)];
	end; 

elseif strcmpi(Mode,'Laplace')
	error('Laplace not supported yet');
	
else
	error(sprintf('%s not supported yet',Mode));
end; 	
