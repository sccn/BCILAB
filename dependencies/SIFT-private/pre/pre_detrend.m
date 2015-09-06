function [data g] = pre_detrend(varargin)
%
% Detrend each channel and trial using a least-squares linear fit.
% Alternately, one can center the data by removing the mean from each
% trial. See [1] for more details.
%
% Inputs:
%
%   data:           raw EEG data
%   SamplingRate:   data sampling rate
%
% Optional:     <'Name',Value> pairs
%
%     VerbosityLevel:   Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                       Possible values: 0,1,2
%                       Default value  : 1
%                       Input Data Type: real number (double)
%
%     DetrendingMethod: Detrending options
%                       Linear: removes the least-squares fit of a straight line from each trial.
%                       Constant: removes the mean from each trial (centering)
%                       Possible values: 'linear','constant'
%                       Default value  : 'linear'
%                       Input Data Type: boolean
% Outputs:
%
%   data:       processed EEG data
%   g:          argument specification structure
%
% See Also: pop_pre_prepData(), pre_prepData(), detrend()
%
% Refences:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
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

SegLenRange = [eps 10];
StepSizeRange = [eps 10];

g = arg_define(varargin, ...
    arg_norep('data',mandatory,[],'Raw data. This should be [num_chans x num_time x num_trials]'), ...
    arg_norep({'srate','SamplingRate'},mandatory,[],'Sampling Rate'), ...
    arg({'verb','VerbosityLevel'},1,{int32(0) int32(1) int32(2)}, ...
    'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
    arg({'method','DetrendingMethod'},{'linear'},{'linear','constant'}, ...
    'Detrending options. Linear: removes the least-squares fit of a straight line from each trial. Constant: removes the mean from each trial (centering)','type','logical'), ...
    arg_subtoggle({'piecewise','Piecewise'},[], ...
    { ...
    arg({'seglength','SegmentLength'},0.33,SegLenRange,'Length of each detrending segment (sec).'), ...
    arg({'stepsize','StepSize'},0.0825,StepSizeRange,'Step size between segment centers (sec). It is recommended to use at least 0.5*SegmentLength (50% overlap).'), ...
    },'Use piecewise detrending. Divide the data into (overlapping) segments and detrend each segment separately. Segment endpoints are merged and stitched together using a cubic spline function to minimize discontinuities at segment intersection points. This is useful for as an alternative to high-pass filtering for removing infraslow oscillations (e.g. SC drift) from the '), ...
    arg({'plot','Plot'},false,[],'Plot results for inspection.') ...
    );

% commit data to workspace
data = g.data;
g=rmfield(g,'data');

[nbchan pnts trials] = size(data);

% set sliding-window parameters
if ~g.piecewise.arg_selection
    % global detrending
    g.piecewise.seglength = pnts/g.srate;
    g.piecewise.stepsize = g.piecewise.seglength;
end
windowing_params = [g.piecewise.seglength g.piecewise.stepsize];

% apply each detrending method (centering,linear detrending, ...)
for i=1:length(g.method)
    
    % initialize verbose output
    m = g.method{i};
    if g.verb && g.piecewise.arg_selection,
        pre = 'Piecewise ';
        post = sprintf(' (using %0.4f sec segment length)',g.piecewise.seglength);
    else
        pre = '';
        post = '';
    end
    if g.verb && strcmpi(m,'mean')
        fprintf('%sCentering data%s...\n',pre,post);    end
    if g.verb && strcmpi(m,'linear')
        fprintf('%sDetrending data%s...\n',pre,post);   end
    
    if g.verb==2
        multiWaitbar('Detrending','Reset','Color',hlp_getNextUniqueColor);
    end
end

% variable intialization
fitlines = zeros(size(data));

% apply the detrending/centering
for ch=1:nbchan
    
    if g.verb==2
        multiWaitbar('Detrending',ch/nbchan);
    end
    
    if g.plot
        % return detrended data as well as fitted curves
        [data(ch,:,:) fitlines(ch,:,:)] = locdetrend_siftmod(squeeze(data(ch,:,:)),g.srate,windowing_params,m);
    else
        % return only detrended data (faster)
        data(ch,:,:) = locdetrend_siftmod(squeeze(data(ch,:,:)),g.srate,windowing_params,m);
    end
    
end

if g.verb==2
    multiWaitbar('Detrending','Close');
end

% plot results, if requested
if g.plot
    eegplot(data+fitlines,'srate',g.srate,'data2',fitlines,'title','Original Data');
    h = gcf;
    ax = findobj(gcf,'tag','eegaxis');
    plts = get(ax,'children');
    legend([plts(end) plts(1)],'original','best local-linear fit');
    eegplot(data,'srate',g.srate,'title','Detrended Data','children',h);
    ax = findobj(gcf,'tag','eegaxis');
    plts = get(ax,'children');
    legend(plts(end),'detrended data');
end



