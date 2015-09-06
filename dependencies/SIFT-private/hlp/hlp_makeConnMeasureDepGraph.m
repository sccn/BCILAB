function dependencies = hlp_makeConnMeasureDepGraph()
% Construct dependency graph for MVAR measures.
% Each of the measures on the right depend on the
% measure on the left (dependency arrow goes from left to right
% e.g., 'Coh',{'iCoh'} means iCoh depends on Coh (Coh-->iCoh)
% On the RHS, only include measures which are *directly* dependent on
% the LHS measure. e.g., if A->B->C, then you should have three entries
%  'A', {'B'};
%  'B', {'C'};
%  'C', {}
% If a new measure is added, you MUST update this depedency table
% (add the method descriptor to the LHS, add any children of this method
% to the RHS and update the RHS of any other entries which this method
% depends on)
% 
% See Also: est_mvtransfer(), est_mvarConnectivity()
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD
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

dependencies = {...
    'dDTF',     {};
    'dDTF08',   {};
    'ffDTF',    {'dDTF'};
    'nDTF',     {};
    'GGC',      {};
    'GGC2',     {};
    'iCoh',     {};
    'Coh',      {'iCoh'};
    'S',        {'Coh','GGC','GGC2','mCoh'};
    'pCoh',     {'dDTF','dDTF08'};
    'pCoh2',    {};
    'mCoh',     {};
    'RPDC',     {};
    'GPDC',     {};
    'nPDC',     {};
    'PDCF',     {};
    'dtf_denom' {'ffDTF','nDTF'};
    'pdc_denom' {};
    'DTF',      {'nDTF','ffDTF','dDTF08','dtf_denom','S','GGC','GGC2'};
    'Sinv',     {'mCoh','pCoh','PDCF'};
    'PDC',      {'DTF','G','nPDC','GPDC','RPDC','PDCF','pCoh2','pdc_denom','Sinv'};
    'Rinv',     {'RPDC','Vpdc'};
    'Vpdc',     {}};