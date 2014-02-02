function spec = override_flags(spec,varargin)
% helper function for the arg_sub* specifiers, to override flags of sub-arguments.
% Spec = override_flags(Spec,Overrides...)
%
% In:
%   Spec : an array of arg_specifier structs
%
%   Overrides... : a sequence of the form {'argname',{'flagname',value,'flagname',value,...},'argname',{...}}
%   
% Out:
%   Spec : the updated specification
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-10-18

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if ~isempty(varargin)    
    % sanity check
    if ~iscellstr(varargin(1:2:end)) || ~all(cellfun('isclass',varargin(2:2:end),'cell'))
        error('The varargin argument must be of the form: {''argname'',{''flagname'',value,''flagname'',value,...},''argname'',{...}}'); end
    % build a remapping table from names to argument indices
    all_names = {spec.names};
    flat_names = [all_names{:}];
    name2index = cumsum(cell2mat(cellfun(@(x)[1 zeros(1,length(x)-1)],all_names,'UniformOutput',false)));
    % for each varargin name...
    for k=1:2:length(varargin)
        % for each affected argument...
        for j=name2index(strcmp(varargin{k},flat_names))
            % more sanity checks
            if ~iscellstr(varargin{k+1}(1:2:end))
                error('The varargin argument must be of the form: {''argname'',{''flagname'',value,''flagname'',value,...},''argname'',{...}}'); end
            % apply each override
            for l=1:2:length(varargin{k+1})
                spec(j).(varargin{k+1}{l}) = varargin{k+1}{l+1}; end
        end
    end
end
