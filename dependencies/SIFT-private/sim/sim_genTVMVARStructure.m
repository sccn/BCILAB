function A = sim_genTVMVARStructure(Aproto,Nl,ndisc,p)
% generate time-varying MVAR structure from a prototype structure
%
% Inputs:
%
%       Aproto:   Prototype structure returned by sim_genVARModelFromEq().
%       Nl:       Number of sample points
%       ndisc:    Number of 'startup' points we will discard when generating data
%       p:        MVAR model order
%
% Outputs:
%
%       A:        Cell array of dimension Nl+ndisc+p containing VAR
%                 coefficient matrices for each time point
%
% See Also: sim_genVARModelFromEq()
%
% References: 
%
%
% Author: Tim Mullen, May 2011, SCCN/INC, UCSD. 
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

    A = repmat({Aproto},1,Nl+ndisc+p);

    if iscell(Aproto)
        % find all the function handles
        numidx = find(cellfun(@(x)isnumeric(x),Aproto));
        funidx = setdiff_bc(1:numel(Aproto),numidx);

        % funidx = find(cellfun(@(x)strcmpi(class(x),'function_handle'),Aproto));


        % construct inline objects
        for fun = 1:length(funidx)
            fx{fun} = inline(Aproto{funidx(fun)});
        end

        % evaluate each function handle
        for t=1:Nl+ndisc+p;

            if ~mod(t,1000)
                fprintf('%d/%d - ',t,Nl+ndisc+p);
            end

            for fun = 1:length(funidx)
                A{t}{funidx(fun)} = double(fx{fun}(t));
            end

            A{t} = cell2mat(A{t});
        end
    end
end