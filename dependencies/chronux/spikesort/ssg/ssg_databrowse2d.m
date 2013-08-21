function ssg_databrowse2d(spikes, assigns, show)
% SSG_DATABROWSE2D  Feature projection GUI in 2D (work in progress).
%     SSG_DATABROWSE(SPIKES) creates a 2-D databrowse figure.  Call with a
%     spike sorting object SPIKES and an optional assignments vector whose length
%     corresponds to the number of waveforms in the SPIKES structure.  If no
%     assignments vector is specified, the function chooses the first of 
%     {final, local, or none} clustering depending on which have been computed.
%     Once the figure appears, you can:
%         Click on the axis labels to change features.
%         Click on a data point to bring the associated cluster to the front.
%         Double-click outside of the axes to make a density snapshot of the current view.

%   Last Modified By: sbm on Fri Jul 29 16:37:08 2005

if (nargin == 1)
    assigns = [];  show = [];
elseif (nargin == 2)
    if (length(assigns) ~= size(spikes.waveforms, 1))
        error('SSG:assignments_length_mismatch', 'The assignments vector length must match the number of waveforms in SPIKES.');
    end
    show = [];    
elseif (nargin == 3)
    if (~all(ismember(show, assigns)))
        error('SSG:show_request_invalid', 'The requested clusters are not in the assignments list');
    end
else
    error('SSG:invalid_number_args', 'The SSG_DATABROWSE2D function only accepts 1-3 inputs.');
end

if (size(spikes.waveforms, 2) < 2)
    error('SSG:data_dimensions_too_small', 'There are not enough data points per waveform.');
end

ssgtest(spikes, assigns, show, 'xy');