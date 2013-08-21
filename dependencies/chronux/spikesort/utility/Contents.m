% SBM Matlab utilities  
% Last updated: 10-07-2006
%
%
% DATATOOLS        [Matrix manipulation, signal processing and statistics]
%   alignwaveforms    - Finds best alignment of waveforms to a reference
%   buttersimple      - Fixed order Butterworth digital filter design
%   cmtm              - Coherence function estimate via a multitaper method 
%   cmplx2rgb         - Convert complex data to RGB data via HSB space
%   convnsep          - Separable N-dimensional convolution
%   embed             - Lag embedding for scalar time series
%   filterm           - Memory-efficient 1-D filter for non-double data
%   gausskernel       - N-dimensional discretized Gaussian kernels 
%   histnd            - N-D histogram
%   histxt            - Column-by-column histogram
%   histxy            - 2-D density histogram
%   histxyz           - 3-D density histogram
%   knn               - K-Nearest Neighbors (naive algorithm)
%   leadingedges      - Finds 0->nonzero transitions
%   minmax            - Simultaneously find min/max elements of N-D arrays
%   mindist           - Indices of min distances in a pairwise dist mtx
%   padmatrix         - Add rows/columns to a matrix
%   pairdist          - Compute a Euclidean distance matrix
%   pcasvd            - Principal Components Analysis
%   pxcorr            - Cross-correlation estimate for point process data
%   rasterize         - Convert a point process to a binary time series
%   removeoffset      - Remove DC offset
%   resetintegrator   - Running counter with resets
%   rescale           - Rescale the range of a data set
%   sindata           - Generate noisy sinusoidal data
%   smooth3f          - Smooth 3D data (fast algorithm)
%   triggerevents     - Extracts events from a time series
%
% DATATYPES        [Manipulation of Matlab arrays, structs, cell arrays and objects]
%   isvectord         - Returns the orientation of a 1-D vector
%   structcat         - Concatenates two structures field-by-field
%   structindex       - Indexes each field in a structure of arrays
%   structpack        - Copies workspace variables into structure field
%   structunpack      - Copy structure fields to workspace variables
%
% GRAPHICSTYLES    [Stylistic elements for visual consistency]
%   Cdgy              - Dark grey color vector
%   Clgy              - Light gray color vector
%   Cmar              - Maroorn color vector
%   Cmbl              - Muted blue color vector
%   Cmgr              - Muted green color vector
%   Cmrd              - Muted red color vector
%   fancy             - Test OpenGL graphics
%   jetm              - Muted JET colormap
%   maplecolor        - Truecolor mimicking default Maple color scheme
%   nice              - Generic visual style settings
%   saturate          - Add saturation markers to a colormap
%   tintmap           - Tinted monochrome color map
%   top               - Set camera to make a surface look like an image
%   xray              - Inverted grayscale colormap
%
% MATLABTOOLS      [General utilities for Matlab programming]
%   bool2onoff        - Converts boolean values to "on"/"off" string
%   busyfigure        - Disable/restore a UI figure
%   callerfile        - Returns the name of the calling .m file
%   copyaxes          - Duplicate axes in a new figure
%   copyfig           - Duplicate all children of a figure
%   datetime2ind      - Convert date/time strings to serial time indices
%   lasterrid         - Returns ID string associated with previous error
%   matlabversion     - Cached access to current Matlab version
%   onoff2bool        - Converts "on"/"off" string to boolean value
%   printf            - Display formatted text
%   progressBar       - Improved progress indicator
%   ps2pdf            - Converts PostScript to Adobe PDF.
%   report_addpage    - Create multipage PostScript from Matlab figures
%   strmatchre        - Select strings matching regular expressions
%
% PLOTTYPES        [Graphical display of data]
%   anglehist         - Magnitude-weighted histogram of polar data
%   errorarea         - Line with confidence regions
%   linec             - Line with varying color
%   mplot             - Efficient plotting of large numbers of lines
%
% UITOOLS          [User Interface menus and tools]
%   linelabel         - Line plot with interactive labeling
%   UIcolorshift      - Interactive colorbar
%   UIlogmenu         - Toggles logarithmic scaling
%   UImovieplayer     - Show a 3-D matrix as an XY-T movie
%   UIsubzoom         - Expands a subplot to cover entire figure
