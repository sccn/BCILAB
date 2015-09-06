function [ARF,PE,argsout] = mvar_vieiramorf(varargin)
% Algorithm: Vieira-Morf
%
% Description:
%
% Unconstrained maximum entropy VAR modeling via
% the Vieira-Morf lattice filter algorithm [2].
%
% This algorithm is a good choice when you have 
% substantially more data than parameters to fit.
% 
% Author Credits:
%
% This implementation is based on the mvar()
% function from Alois Schloegl's BIOSIG and TSA 
% toolboxes [3][4]
% 
% References and Code:
%
% [1] A. Schloegl, Comparison of Multivariate Autoregressive Estimators. Signal processing, Elsevier B.V.
% [2] S.L. Marple. Digital Spectral Analysis with Applications. Prentice Hall, 1987
% [3] http://biosig.sourceforge.net/
% [4] http://hci.tugraz.at/schloegl/matlab/tsa/
% 
% ------------------------------------------------------------------------
% INPUT:
%  Data	 Multivariate data series 
%  p     Model order
%  Mode	 determines estimation algorithm 
%
% OUTPUT:
%  AR    multivariate autoregressive model parameter
%  RC    reflection coefficients (= -PARCOR coefficients)
%  PE    remaining error variance
%
% All input and output parameters are organized in columns, one column 
% corresponds to the parameters of one channel.
%
%  Partial Correlation Estimation: Vieira-Morf [2] using unbiased covariance estimates.
%         In [1] this mode was used and (incorrectly) denominated as Nutall-Strand. 
%
%
% REFERENCES:
%  [1] A. Schlogl, Comparison of Multivariate Autoregressive Estimators.
%       Signal processing, Elsevier B.V. (in press). 
%       available at http://dx.doi.org/doi:10.1016/j.sigpro.2005.11.007
%  [2] S.L. Marple "Digital Spectral Analysis with Applications" Prentice Hall, 1987.
%
%
% A multivariate inverse filter can be realized with 
%   [AR,RC,PE] = mvar(Data,P);
%   e = mvfilter([eye(size(AR,1)),-AR],eye(size(AR,1)),Data);
%  
% see also: MVFILTER, MVFREQZ, COVM, SUMSKIPNAN, ARFIT2

%	$Id: mvar.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1996-2006 by Alois Schloegl <a.schloegl@ieee.org>	
%       This is part of the TSA-toolbox. See also 
%       http://hci.tugraz.at/schloegl/matlab/tsa/
%       http://octave.sourceforge.net/
%       http://biosig.sourceforge.net/
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

%   Revised 24-10-10 by Tim Mullen      - removed all modes except unbiased
%                                         vieira-morf
%                                       - optimized some computations

g = arg_define(varargin, ...
            arg_norep({'data','Data'},mandatory,[],'Data matrix. Format must be [channels x time x trials]'), ...
            arg_nogui({'morder','ModelOrder','Pmax','p'},10,[],'Maximum model order') ...
            );
        
arg_toworkspace(g);
% make data [time x chans x trials]
data = permute(data,[2 1 3]);
% make 2D [(pnts+padding)*trials x chans)]
% with ModelOrder+2 NaNs between trials
data = nanpad(data,g.morder);
  
[N,M] = size(data);

if isempty(morder)
    morder=max([N,M])-1;
end

[C] = zeros(M,'double');
[ARF ARB RCF RCB] = deal(zeros(M,M*morder,'double'));
PE = zeros(M,M*morder+M,'double');
[C(:,1:M),n] = covm(data,'M');
PE(:,1:M)    = C(:,1:M)./n;

%%%%% Partial Correlation Estimation: Vieira-Morf Method [2] with unbiased covariance estimation
%===== In [1] this algorithm is denoted "Nutall-Strand with unbiased covariance" =====%

F = data;
B = data;

PEF = PE(:,1:M);
PEB = PE(:,1:M);
for K = 1:morder,
        [D,n]	= covm(F(K+1:N,:),B(1:N-K,:),'M');
        D       = D./n;


        ARF(:,K*M+(1-M:0)) = D / PEB;	
        ARB(:,K*M+(1-M:0)) = D'/ PEF;	

        tmp        = F(K+1:N,:) - B(1:N-K,:)*ARF(:,K*M+(1-M:0)).';
        B(1:N-K,:) = B(1:N-K,:) - F(K+1:N,:)*ARB(:,K*M+(1-M:0)).';
        F(K+1:N,:) = tmp;

        for L = 1:K-1,
                tmp                    = ARF(:,L*M+(1-M:0)) - ARF(:,K*M+(1-M:0))*ARB(:,(K-L)*M+(1-M:0));
                ARB(:,(K-L)*M+(1-M:0)) = ARB(:,(K-L)*M+(1-M:0)) - ARB(:,K*M+(1-M:0))*ARF(:,L*M+(1-M:0));
                ARF(:,L*M+(1-M:0))     = tmp;
        end

        RCF(:,K*M+(1-M:0)) = ARF(:,K*M+(1-M:0));
        RCB(:,K*M+(1-M:0)) = ARB(:,K*M+(1-M:0));

        [PEF,n] = covm(F(K+1:N,:),F(K+1:N,:),'M');
        PEF     = PEF./n;

        [PEB,n] = covm(B(1:N-K,:),B(1:N-K,:),'M');
        PEB     = PEB./n;

        PE(:,K*M+(1:M)) = PEF;        
end

if any(ARF(:)==inf),
    % Test for matrix division bug. 
    % This bug was observed in LNX86-ML5.3, 6.1 and 6.5, but not in SGI-ML6.5, LNX86-ML6.5, Octave 2.1.35-40; Other platforms unknown.
    p = 3;
    FLAG_MATRIX_DIVISION_ERROR = ~all(all(isnan(repmat(0,p)/repmat(0,p)))) | ~all(all(isnan(repmat(nan,p)/repmat(nan,p))));

    if FLAG_MATRIX_DIVISION_ERROR, 
        %fprintf(2,'### Warning MVAR: Bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.\n');
        warning('MVAR: bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.');

        %%%%% Workaround 
        ARF(ARF==inf)=NaN;
        RCF(RCF==inf)=NaN;
    end
end

if nargout>2
    argsout.RC = RCF;
end