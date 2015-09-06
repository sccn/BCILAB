
function [A stable] = sim_genTVARcoeffs(varargin)
%
% Generate coefficients for a time-varying VAR model from a prototype
% specification produced by sim_genVARModelFromEq().
%
% Inputs                Information
% ----------------------------------------------------------------------
% Aproto:               Prototype structure returned by sim_genVARModelFromEq().
%
% ModelOrder:           The VAR model order                             
%                       Input Data Type: real number (double)                            
%                                                                       
% NumSamplesToSimulate: Number of samples to simulate (for each trial)  
%                       Input Data Type: string                         
%                                                        
% ----------------------------------------------------------------------
% Optional              Information                                     
% ----------------------------------------------------------------------
% NumSamplesToDiscard:  Number of samples to discard in VAR sim         
%                       Input Range  : Unrestricted                     
%                       Default value: 1000                             
%                       Input Data Type: real number (double)            
%                                                                       
% Verbose:              Verbose                                         
%                       Input Range  : Unrestricted                     
%                       Default value: 1                                
%                       Input Data Type: boolean   
%                    
% ----------------------------------------------------------------------
% Outputs               Information
% ----------------------------------------------------------------------
%
% A:                    Cell array of dimension Nl+ndisc+p containing VAR
%                       coefficient matrices for each time point
%
%
% See Also: tvarsim(), sim_genVARModelFromEq()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% [2] Schneider T, Neumaier A (2001) Algorithm 808: ARfit---a matlab package
%   for the estimation of parameters and eigenmodes of multivariate 
%   autoregressive models. ACM Transactions on Mathematical Software 27:58-65
%   http://www.gps.caltech.edu/~tapio/arfit/
%
% Author: Tim Mullen 2011, SCCN/INC, UCSD. 
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


arg_define([0 3],varargin, ...
          arg_norep('Aproto',[],[],'Prototype VAR specification'), ...
          arg('ModelOrder',[],[1 Inf],'The VAR model order','type','denserealdouble'), ...
          arg({'Nl','NumSamplesToSimulate'},500,[1 Inf],'Number of samples to simulate (for each trial)'), ...
          arg({'ndisc','NumSamplesToDiscard'},1000,[],'Number of ''burn-in'' samples to discard'), ...
          arg({'checkStability','CheckStability'},true,[],'Check whether process is stable'), ...
          arg({'verb','VerbosityLevel','Verbose'},int32(0),{int32(0) int32(1) int32(2)},'Verbose','type','int32') ...
          );

if isempty(Aproto)
    error('SIFT:sim_GenTVARcoeffs:badParam','Aproto must be specified');
end
if isempty(ModelOrder)
    error('SIFT:sim_GenTVARcoeffs:badParam','ModelOrder must be specified');
end
if isempty(Nl)
    error('SIFT:sim_GenTVARcoeffs:badParam','NumSamplesToSimulate must be specified');
end

%% convert prototype time-varying VAR specification into coefficients

A = repmat({Aproto},1,Nl+ndisc+ModelOrder);

if nargout>1
    stable = ones(1,Nl+ndisc+ModelOrder);
end
    

if iscell(Aproto)
    % Initiate progress bar
    % ---------------------------------------------------------------------
    if verb==1
        fprintf('Generating coefficients: ');
    elseif verb==2
        waitbarTitle = sprintf('Generating coefficients...');
        multiWaitbar(waitbarTitle, ...
                     'Color', [1.0 0.4 0.0], ...
                     'CanCancel','off');
    end
    % ---------------------------------------------------------------------
    
    % find all the function handles
    numidx = find(cellfun(@(x)isnumeric(x),Aproto));
    funidx = setdiff_bc(1:numel(Aproto),numidx);

    % funidx = find(cellfun(@(x)strcmpi(class(x),'function_handle'),Aproto));
    
    % construct inline objects
    for fun = 1:length(funidx)
        fx{fun} = inline(Aproto{funidx(fun)});
    end

    % total number of samples to simulate
    nt = Nl+ndisc+ModelOrder;
    % evaluate each function handle
    for t=1:nt;

        % update verbose progress
        % -----------------------------------------------------------------
        if verb==2
            cancel = multiWaitbar(waitbarTitle,t/nt);
            if cancel && hlp_confirmWaitbarCancel(waitbarTitle)
                multiWaitbar(waitbarTitle,'Close');
                return;
            else
                multiWaitbar(waitbarTitle,'ResetCancel',true);
            end
        elseif verb==1
            if ~mod(t,1000)
                fprintf('%d/%d - ',t,nt);
            end
        end
        % -----------------------------------------------------------------
        
        for fun = 1:length(funidx)
            A{t}{funidx(fun)} = double(fx{fun}(t));
        end

        A{t} = cell2mat(A{t});
        
        if checkStability
            % check VAR stability

            A1 = A{t};
            M = size(A1,1);
            % ModelOrder = size(A1,2)/M;

            A_hat = [[A1(:,1:M*(ModelOrder-1)); eye(M*(ModelOrder-1))] [A1(:,M*(ModelOrder-1)+1:end); zeros(M*(ModelOrder-1),M)]];

            lambda = eig(A_hat);

            if any(abs(lambda)>1)
                if verb, fprintf('System is unstable at sample %d! MaxLambda=%0.5g\n',t,max(abs(lambda))); end
                if nargout>1, stable(t) = 0; end
            end
        end
    end
    
    
    % clean up
    % ---------------------------------------------------------------------
    if verb==1
        fprintf('\n');
    elseif verb==2
        multiWaitbar(waitbarTitle,'Close');
    end
    % ---------------------------------------------------------------------
    
end




