function ALLEEG = hlp_sift_convAlphaToBeta(ALLEEG)
% Upgrade SIFT dataset(s) from alpha-release to beta-release
% 
% Input: 
%       EEGLAB dataset or array of datasets previously processed by SIFT
%       Datasets must contain SIFT datastructures (in 'ALLEEG.CAT')
% Output:
%       Upgraded EEGLAB dataset or array of datasets
%
% 
% Author: Tim Mullen 2013, SCCN/INC, UCSD.
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

if any(arrayfun(@(EEG) ~isfield(EEG,'CAT'),ALLEEG))
    error('All datasets in ALLEEG must contain a SIFT datastructure (''CAT'' field)');
end

for k=1:length(ALLEEG)
    % copy fields from Alpha structure into an empty Beta structure
    args = hlp_struct2varargin(ALLEEG(k).CAT);
    CAT  = hlp_sift_emptyset(args{:});
    % do the same for subfields
    args = hlp_struct2varargin(CAT.Conn);
    CAT.Conn = hlp_sift_emptyconn(args{:});
    args = hlp_struct2varargin(CAT.MODEL);
    CAT.MODEL = hlp_sift_emptyconn(args{:});
    % fill in additional fields
    CAT.pnts  = size(CAT.srcdata,2);
    CAT.trials= size(CAT.srcdata,3);
    if isfield(CAT.configs,'prepData') ....
       && isfield(CAT.configs.prepData,'newtlims') ...
       && ~isempty(CAT.configs.prepData.newtlims)
        nt = CAT.configs.prepData.newtlims;
        CAT.times = linspace(nt(1),nt(end),CAT.pnts);
    else
        CAT.times = linspace(ALLEEG(k).xmin,ALLEEG(k).xmax,CAT.pnts);
    end
    % wipe original configs
    CAT.alpha_configs = CAT.configs;
    CAT.configs       = getfield(hlp_sift_emptyset(),'configs');
    ALLEEG(k).CAT     = CAT;
end
