function declare_properties(varargin)
% Declare properties of a function.
% declare_properties('Property',value,'Property',value, ...)
%
% declare_properties must be a top-level statement in the function (like persistent) that must come
% before arg_define. The function must use varargin in its argument list.
%
% This construct can be used to declare any key-value properties for a function which can later be
% retrieved using arg_report('properties',@myfunction). BCILAB looks for certain properties in
% plugin functions (e.g., filters, dataset_editing, machine_learning), to support improved
% integration of these functions into the GUI and the evaluation system.
%
% In:
%   Options...: list of name-value pairs.
%
%               Ordering hints applicable for filters (respected by flt_pipeline) are the following:
%                'depends': expresses that some other filter(s) *must* have been called before running 
%                           the filter whose properties are being declared (this is usually due to some 
%                           required meta-data that is supplied by that other filter)
%                'cannot_follow': expresses that the filter being declared cannot follow after some 
%                                 other filter(s) (this is usually due to destructive
%                                 editing performed by those other filters)
%                'cannot_precede': expresses that the filter being declared cannot precede some 
%                                  other filter(s)
%                'follows': expresses that the filter being declared *prefers* to follow some 
%                           other filter(s) (this is usually for numeric or efficiency reasons)
%                'precedes': expresses that the filter being declared *prefers* to precede some 
%                           other filter(s) (this is usually for numeric or efficiency reasons)
%
%               Optional properties respected by the evaluation system (for improved usability):
%                'independent_channels': specify that this filter does not transfer information
%                                        across channels (e.g. channel selection, FIR filter)
%                                        (allows the online system to auto-infer which data channels 
%                                        are actually required by a given BCI model)
%                'independent_trials': specify that this filter does transfer information
%                                      across second-length or larger time scales (on the order of 
%                                      the duration of a trial or larger); this determines whether
%                                      the filter will be called repeatedly for each partition of the 
%                                      data in a cross-validation and other offline analyses
%
%               Further optional properties include:
%                'name': specify the GUI/human-readable name of this function (auto-determined if
%                        left unspecified)
%                'experimental': whether the function implementation is experimental or otherwise 
%                                in "prototype" stage
%                'deprecated': whether the function is deprecated
%
% Examples:
%   function myfunction(...,varargin)
%
%   % declare a property called 'name', and assign the string 'MyFunction' to it, and another
%   % property called 'price', with some value attached
%   declare_properties('name','MyFunction','price',1999)
%
%   % in a signal processing function, declare the name of the function as it should appear in the 
%   % GUI panel (default is the name of the file without the flt_ / set_ prefix)
%   declare_properties('name','MyFilter')
%
%   % in a signal processing function, declare that the filter depends on two other (previously applied)
%   % filters (here: flt_ica and set_fit_dipoles)
%   declare_properties('depends',{'flt_ica','set_flt_dipoles'});
%   
%   % in a signal processing function, declare a hint that the filter typically follows some other 
%   % filter in the pipeline (if that other filter was enabled; here: flt_resample)
%   declare_properties('follows','flt_resample');
%
%   % in a signal processing function, declare a hint that the filter typically precedes some other 
%   % filters in the pipeline (if that other filter was enabled; here: flt_iir and flt_fir)
%   declare_properties('precedes',{'flt_iir','flt_fir'});
%
%   % in a signal processing function, declare the constraint that the filter cannot follow some other 
%   % filter
%   declare_properties('cannot_follow','set_makepos');
%
%   % in a signal processing function, declare the constraint that the filter cannot follow some other 
%   % filter
%   declare_properties('cannot_precede','set_makepos');
%
%   % in a signal processing function, declare that the contents of each output channel produced 
%   % by the filter depend only on the contents of the corresponding input channel (i.e. there is not
%   % cross-mixing of channels (default is false)
%   declare_properties('independent_channels',true);
%
%   % in a signal processing function, declare that the contents of each output trial depend not
%   % only on data of the respective input trial, but perhaps on other trials, as well (i.e., there is 
%   % cross-mixing over time at trial granularity); the effect of this is that the respective filter
%   % will be executed once per cross-validation fold
%   declare_properties('independent_trials',false);
%   
% See also:
%   arg_report
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-02

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
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

caller_args = evalin('caller','varargin');
if length(caller_args)>1 && isequal(caller_args(end-1:end),{'__arg_report__','properties'})
    % properties are being requested via arg_report; return them as a struct
    caller = hlp_getcaller();
    
    % collect properties as name-value pairs
    if ~iscellstr(varargin(1:2:end))
        error('The given arguments must be name-value pairs.'); end
    nvps = reshape([{'depends',{},'cannot_follow',{},'cannot_precede',{},'follows',{},'precedes',{},'independent_channels',[], ...
        'independent_trials',[],'name',caller,'deprecated',false,'experimental',false} varargin],2,[]);
    % find the indices of the last assignment for each name and convert them to a struct
    [s,inds] = sort(nvps(1,:)); inds(strcmp(s(1:end-1),s(2:end))) = [];
    properties = cell2struct(nvps(2,inds),nvps(1,inds),2);
    
    if isempty(properties.independent_trials) && (strncmp(caller,'flt_',4) || strncmp(caller,'set_',4))
        properties.independent_trials = false;
        % warn about this omission: in practice, this means that this filter and all that come after it
        % in a pipeline will be recomputed for every partition of the data during (nested) cross-validation
        % (perhaps 10x or 100x as often as necessary)
        disp_once('Note: The function "%s" does not declare the property independent_trials; assuming that it is false. This may be a conservative assumption that can have massive performance cost during offline processing. Consider declaring the property explicitly in the declare_properties() clause.',caller);
    end
    
    if isempty(properties.independent_channels) && (strncmp(caller,'flt_',4) || strncmp(caller,'set_',4))
        properties.independent_channels = false;
        % warn about this omission: in practice, this often means that BCILAB has to assume that a given
        % model requires all channels that were present in the source data set -- even if it selects
        % only a subset somewhere in its filter pipeline (e.g., excluding peripheral measures). If such
        % channels are not available in an online stream, the model will fail to work, perhaps
        % unneccessarily.
        disp_once('Note: The function "%s" does not declare the property independent_channels; assuming that it is false. This may be a conservative assumption that can make it unnecessariy difficult to set up online processing. Consider declaring the property explicitly in the declare_properties() clause.',caller);
    end
    
    % report them
    arg_issuereport(properties);
end
