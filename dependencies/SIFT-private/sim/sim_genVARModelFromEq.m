
function A = sim_genVARModelFromEq(expr,morder)
% 
% Generate VAR model coefficient matrix from text-based system of equations
% The model can then be realized using the tvarsim() function modified from 
% ARfit [2].
% If the expr specification requests time-varying coefficients, then this
% function will return inline function placeholders for coefficients, which
% can be evaluated to actual coefficient matrices by sim_genTVARcoeffs().
%
% Inputs:
%
%   expr:       A cell vector containing each equation as a string (one
%               equation per cell element). See Examples for format.
%   morder:     The model order
% 
% Outputs:
%
%   A:          VAR[p] model coefficients in format A=[A1,A2, ... Ap] where
%               p = morder. Ai is the M x M coefficient matrix for lag i.
%
% Example1: generate a static VAR[3] model from a (text-based) system of equations
%
% expr = { ...
%     'x1(t) = 0.9*x1(t-1)  + 0.3*x2(t-2)  + e1(t)' ...
%     'x2(t) = 1.3*x2(t-1)  + -0.8*x2(t-2) + e2(t)' ...
%     'x3(t) = 0.3*x1(t-2)  + 0.6*x2(t-1)  + e3(t)' ...
%     'x4(t) = -0.7*x4(t-3) + -0.7*x1(t-3) + 0.3*x5(t-3) + e4(t)' ...
%     'x5(t) = 1*x5(t-1)    + -0.4*x5(t-2) + 0.3*x4(t-2) + e5(t)' ...
%     };
% A = sim_genVARModelFromEq(expr,3)
% A =
%     0.9000         0         0         0         0         0    0.3000         0         0         0         0         0         0         0         0
%          0    1.3000         0         0         0         0   -0.8000         0         0         0         0         0         0         0         0
%          0    0.6000         0         0         0    0.3000         0         0         0         0         0         0         0         0         0
%          0         0         0         0         0         0         0         0         0         0   -0.7000         0         0   -0.7000    0.3000
%          0         0         0         0    1.0000         0         0         0    0.3000   -0.4000         0         0         0         0         0
%
% % Now generate some data from the model (this requires arsim.m or tvarsim.m)
%
% M   = size(A,1);        % number of variables
% Nl  = 1000;             % length of each trial
% Nr  = 100;              % number of trials
% C = eye(M);             % specify the covariance matrix
% data = zeros(M,Nl,Nr);  
% for tr=1:Nr             % simulate data from VAR model
%    data(:,:,tr) = arsim(zeros(1,M),A,C,Nl)';
% end
% eegplot(data,'srate',1); % visualize the simulated data
%
%
% Example 2: generate a time-varying VAR[2] model
% 
%
% See Also: tvarsim(), est_fitMVAR(), sim_genTVARcoeffs()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1 
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% [2] Schneider T, Neumaier A (2001) Algorithm 808: ARfit---a matlab package
%   for the estimation of parameters and eigenmodes of multivariate 
%   autoregressive models. ACM Transactions on Mathematical Software 27:58-65
%   http://www.gps.caltech.edu/~tapio/arfit/
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
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


Nvars = length(expr);

if nargin>1
    % create cell array of zeros
    A = repmat({0},Nvars,Nvars*morder);
end

fnstr ='{[\S]+}';

for i = 1:Nvars
    eq = expr{i};
    
    % Identify all inline functions
    % Specifically, we extract any string f(x) encapsulated by 
    % curly-brackets '{f(x)}'
    [matchstr] = regexp(eq,'{[\S]+}','match');
    
    % replace each function with an identifier '$i' designating that the
    % ith equation will go here
    for fi=1:length(matchstr)
        eq = strtrim(strrep(eq,matchstr{fi},sprintf('$%d',fi)));
    end
    
    % remove all whitespace
    eq(isspace(eq))=[];
%     
%     for fi=1:length(matchstr)
%         eq = [eq(1:matchstart(fi)-1) sprintf('$%d',fi) eq(matchend(fi)+1:end)];
%     end
    
    terms = regexp(eq,'([=+*])','split');
    
    row = str2double(regexp(terms{1},'(\d+)','match'));
    for k = 2:2:length(terms)-1
        vars = regexp(terms{k+1},'(\d+)','match');
        vars = str2num(char(vars));
        col = vars(1)+(vars(2)-1)*Nvars;
        
        if strcmpi(terms{k}(1),'$')
            % term is a function expression
            A{row,col} = matchstr{str2double(terms{k}(2))}(2:end-1);
        elseif ~isnan(str2double(terms{k}))
            % term is a number
            A{row,col} = str2double(terms{k});
        else
            % term is an expression to be directly evaluated
            A{row,col} = eval(terms{k});
        end
    end
end

% if A is all numeric, convert to standard array
if all(cellfun(@(x)isnumeric(x),A))
    A = cell2mat(A);
end
