function ssg_databrowse3d(spikes, assigns, show)
% SSG_DATABROWSE3D  Feature projection GUI in 3D (work in progress).
%     SSG_DATABROWSE(SPIKES) creates a 3-D databrowse figure.  Call with a
%     spike sorting object SPIKES and an optional assignments vector whose length
%     corresponds to the number of waveforms in the SPIKES structure.  If no
%     assignments vector is specified, the function chooses the first of 
%     {final, local, or none} clustering depending on which have been computed.
%     Once the figure appears, you can:
%         Click on the axis labels to change features.
%         Click anywhere in the axis to rotate the camera.
%         Double-click outside of the axes to make a density snapshot of the current view.

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
    error('SSG:invalid_number_args', 'The SSG_DATABROWSE3D function only accepts 1-3 inputs.');
end

if (size(spikes.waveforms, 2) < 3)
    error('SSG:data_dimensions_too_small', 'There are not enough data points per waveform.');
end

ssgtest(spikes, assigns, show, 'xyz');