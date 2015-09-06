% brainmovie3d() - generate a sequence of images showing event-related coherence,
%               event-related spectral perturbations, and inter-trial coherence
%               of localized EEG waveforms. Uses outputs of timef() and cross().
% Usage:
%   >> brainmovie3d(ersps,itcs,crossfs_amp,crossfs_phase,TIMES,freqs,SELECTED,...
%                         'keyword1',value1,...); % creates files image0001.eps, etc.
%
% Inputs:
% ersps         - Cell array (components,conditions) of ERSP arrays
% (freqs,TIMES)
%                 ERSP = event-related spectral perturbation; returned by
%                 timef()
% itcs          - Cell array (components,conditions) of ITC arrays
% (freqs,TIMES)
%                 ITC = inter-trial coherence; returned by timef()
% crossfs_amp   - Cell array (components,components,conditions) of crossf()
%
%                 amplitude output arrays of size (freqs,TIMES).
% crossfs_phase - Cell array (components,components,conditions) of crossf() phase
%                 output arrays of size (freqs,TIMES). (Only the upper diagonal part
%                 of the matrix is taken into account).
% TIMES         - Array of TIMES returned by timef() or crossf()
% freqs         - Indices into the array of freqs returned by timef() or
% crossf()
%                 (v.henv.g., [1:2] means plot the mean of the first two frequencies).
%                 These indexes determine for which freqs plotting will be
%                 performed.
% SELECTED      - Component indices to plot (default all)
%
% Optional 'keyword' parameters:
% 'latency'   - plot only a subset of latencies. The time point closest to the
%               latency given are plotted. Default = empty, all latencies.
% 'frames'    - vector of frame indices to compute. [1:2] only computes the
%               first two frames.
% 'envelope'  - (2,points,conditions) envelopes of the average data (ERP)
%               in each condition
%               (envelope =  min and max traces of each ERP across all
%               channels and TIMES)
% 'rt'        - cell array of vector containing reaction TIMES of the subject in
%               each conditions. This will plot a small bar which height will vary
%               based on the probability of response (default {} -> ignored)
% 'flashes'   - vector of time indices at which the background flashes.  Specify the color
%               of the flash with a cell array of [1,2] cell arrays.
%               Ex. { { 200 'y' } { 1500 '5' }} will generate two flashes,
%               yellow at 200 ms and red at 1500 ms
%
% Movie ITC, Power and Crossf options:
% 'power'     - ['on'|'off'] vary the size of the component disks according to spectral power
%                                                           {default: on}
% 'itc'       - ['on'|'off'] vary component disk colors according to inter-trial coherence
%							    {default: on}
% 'crossf'    - ['on'|'off'] plot | do not plot coherence   {default: on}
% 'crossfcoh' - ['on'|'off'] vary the width of the connecting arc
%                               according to cross-coherence magnitude {def: on}
% 'crossfphasecolor'   -['on'|'off'] vary the arc color according to coherence {default: on}
% 'crossfphasespeed'   - ['on'|'off'] vary the arc speed according to
%                                      cross-coherence phase {def: off}
% 'crossfphaseunit'    - ['degree'|'radian']. Coherence phase angle unit
% {Default is degree}.
% 'edgeColormap'       - colormap array for arcs {default: hsv(64) with green as 0}
% 'nodeColormap'       - colormap array for disks {default: hot(64)}
% 'nodeSizeDataRange'  - [min max] data range for disk size {default is min/max of data}
% 'nodeSizeLimits'     - [min max] limits for node size {default 0.1 to 1}
% 'nodeColorDataRange' - [min max] data range for disk color {default is min/max of data}
% 'nodeColorLimits'    - [min max] limits for node color {default full colormap 0 to 1}
% 'edgeSizeDataRange'  - [min max] data range for edge size {default is min/max of data}
% 'edgeSizeLimits'     - [min max] limits for edge size {default 0.1 to 1}
% 'edgeColorDataRange' - [min max] data range for edge color {default is
% min/max of data}
% 'edgeColorLimits'    - [min max] limits for edge color {default full colormap 0 to 1}
% 'polarity'  - ['pos'|'posneg'] polarity for ITC and crossf. 'pos' = only positive values
%               'posneg' = positive and negative values.
%
% Movie coordinates and axis options:
% 'magnify'   - integer magnification factor for graphics. Default is 1.
% 'diskscale'   - numeric value that scales the size of disks {default: [1.0]}
% 'xlimaxes'    - x-axis limits axis for the component locations {default: [-1 1]}
% 'ylimaxes'    - y-axis limits axis for the component locations {default: [-1 to 1]}
% 'coordinates' - 3-column array of [x y z] coordinates of the SELECTED components
%                 {default: spaced evenly around the head circle boundary}
% 'square'    - ['on'|'off'] re-square all coordinates (so X and Y width is the same)
%               default is 'on';
% 'project3d' - ['on'|'off'] project disks on each 3-D axis. Default is 'off'.
% 'circfactor'  - (ncomps,ncomps) array of arc curvatures (0=straight; 1=half-round,
%                 positive or negative values give the sense of rotation)
%                 {def: 0s}
% 'envylabel'   - ordinate label for envelope. {Default 'Potential \muV'}
% 'envvert'     - cell array of time indices at which to draw vertical lines.
%                 Can also be a cell array of cell to specify line aspect. For instance
%                 { { 0 'color' 'b' 'linewidth' 2 } {1000 'color' 'r' }} would draw two
%                 lines, one blue thick line at latency 0 and one thin red line at latency 1000.
% 'rthistloc' - location and size of rt histograms in individual axes.
%               [abscissa ordinate width maxheight].
% 'title'       - (string) main movie title
% 'condtitle'   - (string array) condition titles (one condition title per row)
% 'condtitleformat' - list of title properties. Ex: { 'fontize', 12,
% 'fontweight', 'bold' }
% 'plotorder'   - [integer vector] component plot order from 1 to the number of SELECTED
%                 components.
% 'backcolor' - [float array] background color. Default is [1 1 1] (white).
%
% Picture and movie output options:
% 'moviename'  - ['string'] Movie file name. Default is "output.avi".
% 'movieopts'  - [cell] Movie options for avifile function. See "help
% avifile".
% 'framesout'  - ['eps'|'ppm'|'fig'|'tiff'|'none'] Default format for saving frames on disk.
%                Default is 'tiff'.
% 'framefolder' - [string] frames output folder. Default uses current directory.
%               the directory is created if it does not exist.
% 'visible'    - ['on'|'off'] show the images on the screen or keep them hidden {default 'on'}
% 'size'      - [widthcond height] output image size {default [400,400]}
%               widthcond is the width of a single condition plot (in pixels)
% 'view'      - 3D static starting view.  See help view. Default is [1 0 0].
% 'path3d'    - ['on'|'off'|[thetafact phifact]] 'on' activate automatic rotation in 3-D. Use
%               [exttheta extphi] to specify theta and phi multiplicative factor (default is
%               [1 0.75]. Use parameter 'view' to specify starting view point. Default is
% 'stereo'    - [Real] Create a stereo movie. The figure should contain a [left right]
%               display of two identical 3-D plots. The left plot view will follow the
%               given 'path' (see above). The right plot axis will be 3-D rotated by an
%               additional horizontal disparity angle specified by the 'stereo' argument:
%               6 (degrees) suggested. Default is [] = mono display.
%               'off'.
%Outputs to disk:
% imageX      - brainmovie3d() saves an output.avi movie (see 'moviename' option above)
%               and a sequence of image files to disk (image0001.eps, as define in the
%               'framesout' option).
%Example:
%
% % Given ICA activations in array icaact (size ncomps,nframes,ntrials), animate (here)
% % activity at/between two components at 176 points per epoch (from -100 ms to 600 ms
% % re stimulus onset) assuming a 250-Hz sampling rate and 100 output
% frames
%
% >> [ersps{1,1},itcs{1,1},powbase,TIMES,freqs] = ...                          % timef for
%                timef(icaact(1,:),176,[-100 600],'Component
%                1',250,1,32,100); %     1st comp
% >> [ersps{2,1},itcs{2,1},powbase,TIMES,freqs] = ...                          % timef for
%                timef(icaact(2,:),176,[-100 600],'Component 2',250,1,32,100); %     2nd comp
% >> [crossfs_amp{1,2},mcoh,TIMES,freqs,cohboot,crossfs_phase{1,2}] = ...      % crossf for
%      crossf_(icaact(1,:),icaact(2,:),176,[-100 600],'Crossf 1 and
%      2',250,1,32,100); % both
%
% >> brainmovie3d( ersps, itcs, crossfs_amp, crossfs_phase, TIMES, [1:2] );
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 30 Mai 2003
%         Tim Mullen, SCCN/INC, UCSD 2010
%
% Note: Better resolution movies can be generated by .eps -> .ppm -> .avi,
%       (or, under a planned upgrade to brainmovie3d, from Matlab6 to .avi directly).
% >> !/usr/local/bin/convert images*.eps movie.mpg % ImageMagic 'convert' may be
%                                                  % used to generate the movie.

% arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2003
% tim@ucsd.edu,  Tim Mullen, SCCN/INC/UCSD 2010

% This program is free software; you can redistribute it and/or
% modify it.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% $Log: brainmovie3d.m,v $
% Revision 1.24  2008/11/25 23:19:28  arno
% test for cirfactor
%
% Revision 1.23  2008/08/20 22:42:15  arno
% fixed last changes
%
% Revision 1.22  2008/08/20 22:28:37  arno
% Reworded documentation, fixed various problems
%
% Revision 1.21  2008/05/09 23:30:05  arno
% new coordinate format for dipoles
%
% Revision 1.20  2007/03/09 18:58:33  arno
% axes size and coherence thickness
%
% Revision 1.19  2007/03/06 18:56:01  arno
% *** empty log message ***
%
% Revision 1.18  2007/03/02 22:42:56  arno
% trying to fix the function
%
% Revision 1.17  2005/07/26 21:24:47  arno
% add option 'project3d' to function header
%
% Revision 1.16  2004/06/30 01:08:31  arno
% test if width>1280
%
% Revision 1.15  2004/06/29 23:31:19  arno
% recompute size
%
% Revision 1.14  2004/04/13 18:08:43  arno
% dipole projections
%
% Revision 1.13  2003/10/21 21:55:19  arno
% + for dB
%
% Revision 1.12  2003/10/21 19:24:30  arno
% new look for caption ...
%
% Revision 1.11  2003/10/09 01:01:34  arno
% fixing coordinate problem for multiple conditions
%
% Revision 1.10  2003/10/08 23:59:03  arno
% removing s to framesfolder
%
% Revision 1.9  2003/10/08 23:55:12  arno
% fix framefolder
%
% Revision 1.8  2003/10/08 23:21:15  arno
% 3d problem
%
% Revision 1.7  2003/10/08 22:54:16  arno
% updating for 3 conditions
%
% Revision 1.6  2003/09/18 22:47:06  arno
% adding frames output
%
% Revision 1.5  2003/08/06 14:42:19  arno
% *** empty log message ***
%
% Revision 1.4  2003/07/03 16:50:20  arno
% header typo
%
% Revision 1.3  2003/07/03 16:49:44  arno
% test for nan itc
%
% Revision 1.2  2003/07/03 16:07:22  arno
% disable axis redrawing
%
% Revision 1.1  2003/07/03 00:23:01  arno
% Initial revision
%

function [alltimepoints mov g] = brainmovie3d_causal(ALLNODESIZE,ALLNODECOLOR,ALLEDGESIZE,ALLEDGECOLOR,TIMES,FREQS,SELECTED,varargin)

if nargin < 6
    help brainmovie3d_causal;
    return;
end

% create structure for option if necessary
%-----------------------------------------
if ~isempty( varargin ),
    for index=1:length(varargin)
        if iscell(varargin{index})
            varargin{index} = { varargin{index}};
        end
    end
    g=struct(varargin{:});
else
    g= [];
end

if nargin < 7
    SELECTED = 1:size(ALLNODESIZE, 1);
end

nbconditions = size(ALLNODESIZE,2);
nbcomponents = size(ALLNODESIZE,1);


try, g.dipplotopt;          catch, g.dipplotopt = {}; end
try, g.mri;                 catch, g.mri = ''; end
try, g.figurehandle;        catch, g.figurehandle = []; end
try, g.opengl;              catch, g.opengl='on'; end
try, g.head;                catch, g.head=''; end
try, g.visible;             catch, g.visible='on'; end
try, g.square;              catch, g.square='on'; end
try, g.moviename;           catch, g.moviename='output.avi'; end
try, g.movieopts;           catch, g.movieopts={}; end
try, g.rt;                  catch, g.rt={}; end
try, g.modulateNodeColor;    catch, g.modulateNodeColor='on'; end
try, g.modulateNodeSize;     catch, g.modulateNodeSize='on'; end
try, g.modulateEdgeColor;    catch, g.modulateEdgeColor='on'; end
try, g.modulateEdgeSize;     catch, g.modulateEdgeSize='on'; end
try, g.latency;             catch, g.latency=[]; end
try, g.magnify;             catch, g.magnify=1; end
try, g.crossfcoh;           catch, g.crossfcoh='on'; end
try, g.size;                catch, g.size=[400 400]; end
try, g.crossfphasespeed;    catch, g.crossfphasespeed='off'; end
try, g.crossfphaseunit;     catch, g.crossfphaseunit='degree'; end

try, g.nodeSizeDataRange;   catch, g.nodeSizeDataRange  = []; end
try, g.nodeSizeLimits;      catch, g.nodeSizeLimits     = [0.1 1]; end
try, g.nodeColorDataRange;  catch, g.nodeColorDataRange = []; end
try, g.nodeColorLimits;     catch, g.nodeColorLimits    = [0 1]; end

try, g.edgeSizeDataRange;   catch, g.edgeSizeDataRange  = []; end
try, g.edgeSizeLimits;      catch, g.edgeSizeLimits     = [0.1 1]; end
try, g.edgeColorDataRange;  catch, g.edgeColorDataRange = []; end
try, g.edgeColorLimits;     catch, g.edgeColorLimits    = [0 1]; end

try, g.diskscale;           catch, g.diskscale = 1; end
try, g.framefolder;         catch, g.framefolder = ''; end
try, g.caption;             catch, g.caption = true; end
try, g.frames;              catch, g.frames = []; end
try, g.envvert;             catch, g.envvert = {}; end
try, g.events;              catch, g.events = {}; end
try, g.polarity;            catch, g.polarity = 'pos'; end
try, g.framesout;           catch, g.framesout = 'tiff'; end
try, g.condtitle;           catch, g.condtitle = []; end
try, g.condtitleformat;     catch, g.condtitleformat = {'fontsize', 14, 'fontweight', 'bold', 'color' 'w'}; end
try, g.title;               catch, g.title = ''; end                            %% This determines the text for the upper-right corner. Can be multi-line (cell array)
try, g.envylabel;           catch, g.envylabel = '\muV'; end
try, g.plotorder;           catch, g.plotorder = SELECTED; end
try, g.coordformat;         catch, g.coordformat = 'spherical'; end
try, g.stereo;              catch, g.stereo = []; end
try, g.backcolor;           catch, g.backcolor = [0 0 0]; end
try, g.rotationpath3d;      catch, g.rotationpath3d = struct('AngleFactor',1,'PhaseFactor',0.75,'FramesPerCycle',max(1,length(g.latency))); end
try, g.project3d;           catch, g.project3d = 'off'; end
try, g.view;                catch, g.view = [43.6650 30.4420]; end
try, g.cameraMenu;          catch, g.cameraMenu = true; end                         %% TM: whether or not to show the camera menu (lighting etc).
try, g.plotNodeLabels;      catch; g.plotNodeLabels = true; end                     %% TM: Whether or not to plot node labels
try, g.drawmode;            catch, g.drawmode = 'normal'; end                       %% DrawMode for brainmovie ('fast' renders brain more quickly (but poorer) than 'normal'). See Axes Properties for additional details
try, g.footerPanelPlotMode; catch, g.footerPanelPlotMode = {'all','envelope'}; end %% Plot mode for footer panel display (plot all traces and/or envelope)
try, g.envColor;            catch, g.envColor = [1 0 0]; end                       %% TM: evelope color
try, g.makeCompass;         catch, g.makeCompass = false; end                       %% TM: label cardinal directions (posterior,anterior, left, right)
try, g.compassLabels;       catch, g.compassLabels = {'Dorsal ','Anterior ','Right Lat '};   end % labels for Z, Y, X coordinates
try, g.windowLength;        catch; g.windowLength = []; end                        %% length of sliding window (for footer panel display)
try, g.footerPanelData;     catch, g.footerPanelData = []; end
try, g.footerPanelTitle;    catch, g.footerPanelTitle = ''; end                    %% TM: Title for footer panel display
try, g.flashEvents;         catch, g.flashEvents = true; end                    	%% TM: whether or not to flash screen at event times
try, g.LONImesh;            catch, g.LONImesh = []; end                             %% TM
try, g.LONITransparency;    catch, g.LONITransparency = 1; end                      %% TM
try, g.causality;           catch, g.causality = 0; end                             %% TM: added
try, g.nodelabels;          catch, g.nodelabels = {}; end                           %% TM: added
try, g.collapsefun;         catch, g.collapsefun = 'mean'; end                      %% TM: added  (can be 'integrate','mean','max','absmax','peak')
try, g.footerPanelTimes;    catch, g.footerPanelTimes = TIMES; end                  %% TM: added (for envelope plot)
try, g.showLatency;         catch, g.showLatency = 0; end                           %% TM: added
try, g.verb;                catch, g.verb = 0;  end                                 %% TM: added (display progress bar for making movie)
try, g.csf;                 catch, g.csf = []; end                                  %% TM
try, g.cortexTransparency;  catch, g.cortexTransparency = 1; end                    %% TM: transparency of superimposed cortex (1 = don't plot)
try, g.speedy;              catch, g.speedy = true; end                         	%% TM: for fast rendering -- some features disabled
try, g.mode;                catch, g.mode = 'init_and_render'; end                  %% TM: initialization mode. if 'init', then initialize only and current state of BM will be returned in g.vars, pass g in as subsequent input. If 'render', then skip init and use g.vars for rendering. if 'init_and_render' then do both (overwrites g.vars)
try, g.vars;                catch, g.vars = []; end                                 %% TM: a structure to hold initialization variables.
try, g.facelighting;        catch, g.facelighting = 'phong'; end % 'gouraud'        %% TM: facelighting -- phong better but slower than gouraud
try, g.updateLegendText;    catch, g.updateLegendText = true; end                   %% TM: always update the legend text to reflect edge/node size/color min/max
try, g.displayLegendLimitText; catch, g.displayLegendLimitText = true; end                  %% TM: whether or not to show text in the legend indicating the node/edge datalimits
try, g.collapsedFreqs;      catch, g.collapsedFreqs = FREQS; end                    %% TM: optionally provide the frequencies which we have collapsed over (for text display)
try, g.plotimgs;            catch, g.plotimgs = true; end                          %% plot background MRI plates on axes
try, g.theme;               catch, g.theme = hlp_getBrainMovieTheme('theme','classic');  end %% TM: color theme structure
try, g.EdgeColorMappedToDirectionality; catch, g.EdgeColorMappedToDirectionality = false; end

% some parameters for captions
try, g.modulateNodeColor;    catch, g.modulateNodeColor = ''; end
try, g.EdgeColorMapping;    catch, g.EdgeColorMapping = ''; end
try, g.NodeSizeMapping;     catch, g.NodeSizeMapping = ''; end
try, g.EdgeSizeMapping;     catch, g.EdgeSizeMapping = ''; end
try, g.ConnMethod;          catch, g.ConnMethod = 'Connectivity'; end
try, g.nodeColorPolarity;   catch, g.nodeColorPolarity = 'pos'; end
try, g.edgeColorPolarity;   catch, g.edgeColorPolarity = 'pos'; end
try, g.centerDataRange;     catch, g.centerDataRange = true; end
try, g.colorshadow;         catch, g.colorshadow = 1; end
try, g.nodeColormap;        catch,
    colormtmp = hot(64);
    colormtmp(end,3) = (colormtmp(end,3)+colormtmp(end-1,3))/2; % white does not come out when the
    g.nodeColormap = colormtmp;                                    % the figure is printed to ppm
    g.nodeColormap(:,1) =  colormtmp(:,2);
    g.nodeColormap(:,2) =  colormtmp(:,3);
    g.nodeColormap(:,3) =  colormtmp(:,1);
    g.nodeColormap = [ g.nodeColormap; colormtmp(end:-1:1,:)];
    g.nodeColormap = jet(64);
end
try, g.edgeColormap; catch,
    g.edgeColormap = jet(64);
    %g.edgeColormap = hsv(64);
    %g.edgeColormap = [ g.edgeColormap(55:end,:);
    %g.edgeColormap(1:54,:)]; g.edgeColormap = g.edgeColormap(linspace(64, 1, 64),:); % reorganize the colormap
    %g.edgeColormap = hsv(64);
    %g.edgeColormap = [g.edgeColormap(16:end,:); g.edgeColormap(1:5,:)];
end
try, g.xlimaxes; 		catch, g.xlimaxes = [-1 1]; end
try, g.ylimaxes; 		catch, g.ylimaxes = [-1 1]; end
try, g.rthistloc; 	    catch, g.rthistloc(1) = (g.xlimaxes(2)-g.xlimaxes(1))*0.74 + g.xlimaxes(1); % abscicia
    g.rthistloc(3) = (g.xlimaxes(2)-g.xlimaxes(1))*0.1; % width
    g.rthistloc(2) = (g.ylimaxes(2)-g.ylimaxes(1))*0.34 + g.ylimaxes(1); % ordinate
    g.rthistloc(4) = (g.ylimaxes(2)-g.ylimaxes(1))*0.1; % max height
end
try, g.coordinates; catch
    % coordinates around a circle
    g.coordinates = zeros( nbcomponents, 3 );
    count = 0;
    for index = SELECTED
        if length(SELECTED) > 1
            g.coordinates( index,:) = [ cos(count/length(SELECTED)*2*pi) sin(count/length(SELECTED)*2*pi) 0 ] * 0.7;
        else	g.coordinates(index,:) = [ 0.01 0.01 0 ];
        end
        count = count + 1;
    end
end
try, g.circfactor; catch, g.circfactor = ones( nbcomponents, nbcomponents )*0.01; end
if isempty(g.circfactor), g.circfactor = ones( nbcomponents, nbcomponents )*0.01; end
% if ischar(g.rotationpath3d)
%     switch g.rotationpath3d
%         case 'on', g.rotationpath3d = [ 1 0.75];
%         case 'off', g.rotationpath3d = [];
%     end
% else
%     if length(g.rotationpath3d) ~= 2, error('path3d length have to be a string or a 2 element vector'); end
% end

% add defaults ranges
%--------------------
if isempty(g.nodeSizeDataRange)
    g.nodeSizeDataRange = [Inf -Inf];
    for i=1:length(ALLNODESIZE)
        g.nodeSizeDataRange(1) = min(g.nodeSizeDataRange(1), min(ALLNODESIZE{i}(:)));
        g.nodeSizeDataRange(2) = max(g.nodeSizeDataRange(2), max(ALLNODESIZE{i}(:)));
    end
    
    % make 0 in the center of colormap
    if g.centerDataRange
        g.nodeSizeDataRange(1) = -max(g.nodeSizeDataRange);
        g.nodeSizeDataRange(2) =  max(g.nodeSizeDataRange);
    end
    if g.verb, fprintf('Node size data range automatically set to %1.2f to %1.2f\n', g.nodeSizeDataRange(1), g.nodeSizeDataRange(2)); end
end
if isempty(g.nodeColorDataRange)
    g.nodeColorDataRange = [Inf -Inf];
    for i=1:length(ALLNODECOLOR)
        g.nodeColorDataRange(1) = min(g.nodeColorDataRange(1), min(ALLNODECOLOR{i}(:)));
        g.nodeColorDataRange(2) = max(g.nodeColorDataRange(2), max(ALLNODECOLOR{i}(:)));
    end
    
    % make 0 in the center of colormap
    if g.centerDataRange
        g.nodeColorDataRange(1) = -max(g.nodeColorDataRange);
        g.nodeColorDataRange(2) =  max(g.nodeColorDataRange);
    end
    
    if g.verb, fprintf('Node color data range automatically set to %1.2f to %1.2f\n', g.nodeColorDataRange(1), g.nodeColorDataRange(2)); end
end
if isempty(g.edgeSizeDataRange)
    g.edgeSizeDataRange = [Inf -Inf];
    for i=1:length(ALLEDGESIZE(:))
        if ~isempty(ALLEDGESIZE{i})
            g.edgeSizeDataRange(1) = min(g.edgeSizeDataRange(1), min(ALLEDGESIZE{i}(:)));
            g.edgeSizeDataRange(2) = max(g.edgeSizeDataRange(2), max(ALLEDGESIZE{i}(:)));
        end
    end
    
    % make 0 in the center of colormap
    if g.centerDataRange
        g.edgeSizeDataRange(1) = -max(g.edgeSizeDataRange);
        g.edgeSizeDataRange(2) =  max(g.edgeSizeDataRange);
    end
    if g.verb, fprintf('Edge size data range automatically set to %1.2f to %1.2f\n', g.edgeSizeDataRange(1), g.edgeSizeDataRange(2)); end
end
if isempty(g.edgeColorDataRange)
    g.edgeColorDataRange = [Inf -Inf];
    for i=1:length(ALLEDGECOLOR(:))
        if ~isempty(ALLEDGECOLOR{i})
            g.edgeColorDataRange(1) = min(g.edgeColorDataRange(1), min(ALLEDGECOLOR{i}(:)));
            g.edgeColorDataRange(2) = max(g.edgeColorDataRange(2), max(ALLEDGECOLOR{i}(:)));
        end
    end
    
    % make 0 in the center of colormap
    if g.centerDataRange
        g.edgeColorDataRange(1) = -max(g.edgeColorDataRange);
        g.edgeColorDataRange(2) =  max(g.edgeColorDataRange);
    end
    
    if g.verb, fprintf('Edge color data range automatically set to %1.2f to %1.2f\n', g.edgeColorDataRange(1), g.edgeColorDataRange(2)); end
end

% check size of inputs
% --------------------
try
    if ~all(size(ALLNODESIZE) == size(ALLNODECOLOR))
        disp('Error: ERSP and ITC cells array must be the same size'); return;
    end
    if ~isempty(ALLEDGESIZE)
        if ~all(size(ALLEDGESIZE) == size(ALLEDGECOLOR))
            disp('Error: Crossf amplitude and Crossf angle cells array must be the same size'); return;
        end
        if ~(size(ALLEDGESIZE,2) == size(ALLNODESIZE,1))
            disp('Error: number of components different in ERSP and Crossf arrays'); return;
        end
        if ~(size(ALLEDGESIZE,3) == size(ALLNODESIZE,2))
            disp('Error: number of conditions different in ERSP and Crossf arrays'); return;
        end
        if ~(size(ALLEDGESIZE{1,2,1},1) == size(ALLNODESIZE{1,1},1))
            disp('Error: number of frequencies (rows) different in ERSP and Crossf arrays'); return;
        end
        if ~(size(ALLEDGECOLOR{1,2,1},2) == size(ALLNODECOLOR{1,1},2))
            disp('Error: number of time points (columns) different in ERSP and Crossf arrays'); return;
        end
        if ~(size(ALLEDGESIZE{1,2,1},2) == length(TIMES))
            disp('Error: number of time points (columns) different in TIMES and Crossf arrays'); return;
        end
    end
    try tmp = ALLNODESIZE{1,1}; tmp(FREQS,:); catch, disp('Error: unable to access the defined frequencies in ERSPs (out of bounds) '); return; end
    try ALLNODESIZE{SELECTED,1}; catch, disp('Error: unable to access the defined components in ERSPs (out of bounds)'); return; end
catch
    disp('Error accessing one of the variable. Remember: Except for SELECTED, freqs, TIMES and circfactor, all vars are cell arrays. Check also: dimensions and content.'); return;
end

% check structure content
% -----------------------
if ~isempty(g.rt)
    if length(g.rt) ~= nbconditions
        disp('Error: Rt must be either an array of the size of the number of conditions (might be 0 for some conditions)'); return;
    end
end
switch lower(g.visible)
    case {'on', 'off'} ;
    otherwise disp('Error: Visibility must be either ''on'' or ''off'''); return;
end
switch lower(g.square)
    case {'on', 'off'} ;
    otherwise disp('Error: Square must be either ''on'' or ''off'''); return;
end
switch lower(g.modulateNodeSize)
    case {'on', 'off'} ;
    otherwise disp('Error: Power must be either ''on'' or ''off'''); return;
end
switch lower(g.modulateNodeColor)
    case {'on', 'off'} ;
    otherwise disp('Error: Itc must be either ''on'' or ''off'''); return;
end
switch lower(g.modulateEdgeSize)
    case {'on', 'off'} ;
    otherwise disp('Error: Crossf must be either ''on'' or ''off'''); return;
end
switch lower(g.crossfcoh)
    case {'on', 'off'} ;
    otherwise disp('Error: Crossfcoh must be either ''on'' or ''off'''); return;
end
switch lower(g.modulateEdgeColor)
    case {'on', 'off'} ;
    otherwise disp('Error: Crossfphasecolor must be either ''on'' or ''off'''); return;
end
switch lower(g.crossfphasespeed)
    case {'on', 'off'} ;
    otherwise disp('Error: Crossfphasespeed must be either ''on'' or ''off'''); return;
end
switch lower(g.crossfphaseunit)
    case {'degree', 'radian'} ;
    otherwise disp('Error: Crossfphaseunit must be either ''degree'' or ''radian'''); return;
end
% switch g.caption
%     case {'on', 'off'} ;
%     otherwise disp('Error: Caption must be either ''on'' or ''off'''); return;
% end
switch lower(g.polarity)
    case {'pos', 'posneg'} ;
    otherwise disp('Error: Polarity must be either ''pos'' or ''posneg'''); return;
end
if ~isempty(g.envvert),
    if ~iscell(g.envvert) && ~( isstruct(g.envvert{1}) || isnumeric(g.envvert{1}) )
        disp('Error: Invalid type for Envvert.'); return;
    end
end
if ~isempty(g.latency) && ~isnumeric(g.latency)
    disp('Error: Latency must be a vector'); return;
end
if length(g.nodeSizeDataRange) ~= 2
    disp('Error: Scalepower must be a 2-element array'); return;
end
if length(g.edgeSizeDataRange) ~= 2
    disp('Error: Scalecoher must be a 2-element array'); return;
end
if (length(g.diskscale) ~= 1 || g.diskscale < 0)
    disp('Error: Diskscale must be a scalar value >= 0.'); return;
end
if size(g.nodeColormap,2) ~= 3
    disp('Error: Colmapcoh must be a colormap (3 columns)'); return;
end
if size(g.edgeColormap,2) ~= 3
    disp('Error: Colmapcrossf must be a colormap (3 columns)'); return;
end
if size(g.circfactor,1) ~= size(g.circfactor,2)
    disp('Error: Circfactor must be a square matrix'); return;
end
if ~iscell(g.coordinates) && ~isempty(g.circfactor)
    if size(g.circfactor,1) ~= size(g.coordinates,1)
        disp('Error: Circfactor must have the same number of rows as the number of rows of coordinates'); return;
    end
    if nbcomponents ~= size(g.coordinates,1)
        disp('Error: The array of SELECTED components must have length nrows of the array coordinates'); return;
    end
end
if ~ischar(g.envylabel)
    disp('Error: envelope label must be a string'); return;
end
if ~isempty(g.footerPanelData) && isempty(g.footerPanelTimes)
    if (size( g.footerPanelData,1 ) > 2) || (size( g.footerPanelData,2 ) ~= length(g.footerPanelTimes)) || (size( g.footerPanelData,3 ) ~= nbconditions)
        fprintf('Error: Envelope array does not have the right size (%s), i.v.henv. 2 x %d (number of time points) x %d (number of conditions)\n', int2str(size( g.footerPanelData)), length(TIMES), nbconditions); return;
    end
end
if ~isempty(g.condtitle)
    if iscell(g.condtitle), g.condtitle = strvcat(g.condtitle{:}); end
    if ~isempty(g.condtitle) && size( g.condtitle,1 ) ~= nbconditions
        fprintf('Error: The number of rows in the title array(%d) must match the number of conditions (%d)\n', size(g.condtitle,1), nbconditions); return;
    end
end
if length(g.plotorder) ~= length(SELECTED)
    error([ 'Error: ''plotorder'' must be the same size as the number of SELECTED components:' int2str(length(SELECTED)) ]);
end
if max(g.plotorder) > max(SELECTED)
    error([ 'Error: ''plotorder'' must be below the number of SELECTED components:' int2str(max(SELECTED)) ]);
end

SLASH = fastif(isunix,'/','\');

if all(cellfun(@isempty,g.nodelabels))
    g.plotNodeLabels = false;
end
    
if ~isempty(g.framefolder) && ~isdir(g.framefolder)
    error('%s does not exist',g.framefolder);
    %     [tmp1 tmp2] = mkdir(SLASH, g.framefolder(2:end) );
end

% append a slash to framefolder, if necessary
if ~isempty(g.framefolder) && g.framefolder(end) ~= SLASH
    g.framefolder = [g.framefolder SLASH];
end

if g.verb==2
    g.vars.hwaitbar=waitbar(0,'Initializing Brainmovie...');
end

mov = [];
alltimepoints = [];

% other variables
% ---------------
g.factproj = [-88 88 -71]; %[-71 88 -71];
g.projcolor = [0.35 0.35 0.35];
g.rthistcolor  = [1 1 1];
g.resmult = 1;

% latency and color of flash events
g.vars.flashTimes = [];
g.vars.flashColor = {};

        
% prune the list of nodelabels to selected
if g.plotNodeLabels
    g.nodelabels = g.nodelabels(SELECTED);
end

if ismember_bc(lower(g.mode),{'init','init_and_render'})
    
    % create movie
    % ------------
    if ~isempty(g.moviename)
        %     fprintf('A movie is being saved under %s (movie parameters shown below):',g.moviename);
        mov = avifile(g.moviename, g.movieopts{:});
    end
    
    
    currentphase   = zeros( length(SELECTED), length(SELECTED), nbconditions);
    tmp = ALLNODESIZE{1,1};
    nwin = size(tmp,2);
    
    %     % optional resquare of all coordinates
    %     % -----------------------------------
    %     g.magnify = g.magnify/4;
    
    % compute RT distribution
    % -----------------------
    if ~isempty(g.rt)
        RTdist = zeros(nbconditions,nwin);
        for index = 1:nbconditions
            if ~isempty(g.rt{index})
                timestep = (TIMES(2)-TIMES(1))/2;
                for indeximage = 1:nwin
                    RTdist(index, indeximage) = length( intersect_bc( find( g.rt{index} > TIMES(indeximage)-timestep ) , ...
                        find(  g.rt{index} <= TIMES(indeximage)+timestep ) ) );
                end
                RTdist(index,:) = RTdist(index,:)/max(RTdist(index,:));
            end
        end
        RTdist = RTdist/max(RTdist(:));
    end
    
    if ~isempty(g.figurehandle) && ishandle(g.figurehandle)
        % clear the current figure
        clf(g.figurehandle);
    else
        % create a new figure
        g.figurehandle = figure( 'position', [100, 100, ceil(nbconditions*g.size(1)/4)*4, ceil(g.size(2)/4)*4], ...
            'PaperPositionMode', 'auto', 'papertype', 'A1', 'visible',g.visible,'tag','BrainMovieFigure'); %'paperorientation', 'landscape' );
        set(g.figurehandle,'DefaultTextInterpreter','none');
        if g.cameraMenu
            cameramenu('noreset');
        end
    end
    
    % create an axis
    htmp=axes('parent',g.figurehandle);
    axis(htmp,'off');
    
    if 	strcmpi(g.framesout, 'ppm')
        r = 0.8465;
        pos = get(g.figurehandle,'position');
        if floor(pos(3)/r)> 1280
            fact = 1280/(pos(3)/r);
            set(g.figurehandle, 'position', [ 0 0 1280  floor(pos(4)/r*fact) ]);
        else
            set(g.figurehandle, 'position', [ 0 0 floor(pos(3)/r), floor(pos(4)/r) ]);
        end
    end
    
    pos = get(htmp,'position');
    %    left   bottom
    q = [pos(1) pos(2) 0 0];
    %    width  height width  height
    s = [pos(3) pos(4) pos(3) pos(4)];
    
    % compute SELECTED latency point
    % ------------------------------
    if ~isempty(g.latency)
        g.vars.alltimepoints = [];
        for index = 1:length(g.latency)
            [tmp tmptimepoint] = min(abs(g.latency(index)-TIMES));
            g.vars.alltimepoints = [ g.vars.alltimepoints tmptimepoint];
        end
    else
        if isempty(g.frames)
            g.vars.alltimepoints = 1:nwin;
        else
            g.vars.alltimepoints = g.frames;
        end
    end
    
    % initialize the rotation frame index for path3d
    g.vars.path3d_frameindex = g.vars.alltimepoints(1);
    
    % set up axis for figure background
    % ---------------------------------
    g.vars.hFigureBg = axes('position' , [0 0 1 1], 'color',g.backcolor,'xtick', [], 'ytick', [], 'box', 'off','tag','figureBackground','parent',g.figurehandle);
    xlim(g.vars.hFigureBg,[0 1]);
    ylim(g.vars.hFigureBg,[0 1]);
    
    % setup event flash objects
    % -------------------------
    if g.flashEvents
        for index = 1:length(g.events)
            g.vars.flashTimes(index) = g.events{index}{1};  % convert to ms
            g.vars.flashColor{index} = g.events{index}{2};
        end
        g.vars.allFlashIndices = getindex(TIMES,g.vars.flashTimes);
        
        % create a rectangular border out of two patches
        % outer patch (flash colored)
        g.vars.hFlashPatch = patch([0.03 0.8 0.8 0.03], ...  %   [ 0.1 0.84 0.84 0.1 ]
            [0.01 0.01 0.99 0.99], ...
            g.backcolor,'parent',g.vars.hFigureBg,'visible','off');
        % inner patch (background color)
        patch([ 0.05 0.78 0.78 0.05 ], [0.03 0.03 0.9 0.9], ...  %   [ 0.13 0.8 0.8 0.13 ]
            g.backcolor,'parent',g.vars.hFigureBg, ...
            'facecolor' , g.backcolor, 'edgecolor', 'none');
        
        %         g.vars.hFlashPatch = patch([ 0.13 0.84 0.84 0.13 ], [0.8 0.8 0.93 0.93], [0.5 0.5 0.5],'parent',g.vars.hFigureBg);
        
        posf = 0; % used as a counter to preserve color
    end
    
    % draw axes and display images
    % ----------------------------
    ordinate = 0.2;
    max_ordinate = 1-1.4*ordinate;   % makes space at top for figure title
    maxcoordx    = 1.1-1/nbconditions/4;
    coords = g.coordinates;
    g.coordinates = cell(nbconditions);
    
    for tmpcond=1:nbconditions
        
        % plot 3d head
        % ------------
        g.vars.hBrain(tmpcond) = axes('position', [0+maxcoordx/nbconditions*(tmpcond-1), ordinate, ...
            maxcoordx/nbconditions*0.9, max_ordinate].*s+q ,...
            'parent',g.figurehandle,'DrawMode',g.drawmode);
        gr = [ 0.3 0.3 0.3 ];
        g.dipplotopt = [{ 'coordformat' g.coordformat 'gui', 'off', 'cornermri', 'on', ...
            'color', { gr gr gr gr gr gr gr gr gr } } g.dipplotopt];
        
        if ~isempty(g.mri)
            g.dipplotopt = [g.dipplotopt 'mri' g.mri];
        end
        
        % determine whether or not to render MRI image plates
%         g.dipplotopt = [g.dipplotopt 'plotimgs' fastif(g.plotimgs,'on','off')];
            
        if iscell(coords)
            for index = 1:length(coords)
                if size(coords{index},1) == 2 && all(coords{index}(2,:) == 0), coords{index}(2,:) = []; end
                dipstruct(index).posxyz = coords{index};
                dipstruct(index).momxyz = [1 1 1];
                if size(dipstruct(index).posxyz,1) == 2, dipstruct(index).momxyz(2,:) = [1 1 1]; end
                dipstruct(index).component = index;
                dipstruct(index).rv = 0.1;
            end
        else
            for index = 1:size(coords, 1);
                dipstruct(index).posxyz = coords(index,:);
                dipstruct(index).momxyz = [0 0 0];
                dipstruct(index).component = index;
                dipstruct(index).rv = 0.1;
            end
        end
        
        % unfortunately, we have to temporarily change the current axis so
        % dipplot will render into the correct axis
        curax = get(gcf,'currentaxes');
        curfig = gcf;
        set(0,'CurrentFigure',g.figurehandle);
        set(g.figurehandle,'CurrentAxes',g.vars.hBrain(tmpcond));
        %         axes(g.vars.hBrain(tmpcond));
%         CameraPosition = campos(g.vars.hBrain(tmpcond));
        dipplot_sift( dipstruct, 'view', g.view, g.dipplotopt{:});
%         campos(g.vars.hBrain(tmpcond),CameraPosition);
        axis(g.vars.hBrain(tmpcond),'off');
        if ~isempty(curax)
            set(0,'CurrentFigure',curfig);
            set(curfig,'CurrentAxes',curax);
        end
        %         if ~isempty(curax)
        %             axes(curax);  % revert focus back to original axis
        %         end
        
        set(g.vars.hBrain(tmpcond), 'cameraviewanglemode', 'manual','CameraPositionMode','manual'); % disable change size
        set(g.vars.hBrain(tmpcond),'tag',['brain' num2str(tmpcond)]);
        axis(g.vars.hBrain(tmpcond),'vis3d') % same as above (for security)
        set(g.vars.hBrain(tmpcond),'XLimMode','manual','YLimMode','manual','ZLimMode','manual');
%         set(g.vars.hBrain(tmpcond),'Xlim',[-200 200],'Ylim',[-200 200],'Zlim',[-200 200]);
        
        
        view(g.vars.hBrain(tmpcond),g.view);
        camva(g.vars.hBrain(tmpcond),(1/g.magnify)*3.44); % fix the camera zoom
        
%         camproj(g.vars.hBrain(tmpcond),'perspective');
        
        if ~g.plotimgs
            set(findobj('parent', g.vars.hBrain(tmpcond), 'tag', 'img'),'visible','off');
            %             camzoom(g.vars.hBrain(tmpcond),1/(1.2));  % to offset camzoom in dipplot
        end
        
        g.vars.coordinates = g.coordinates;
        
        % get the coordinates of the nodes in MNI space
        for index = 1:length(dipstruct)
            htmp = findobj(g.vars.hBrain(tmpcond), 'tag', [ 'dipole' int2str(index) ]);
            for dipindex = 1:length(htmp)
                tmpstruct = get(htmp(dipindex), 'userdata');
                if isstruct(tmpstruct) % look for dipole location % THIS DOES NOT WORK
                    if isfield(tmpstruct, 'pos3d') && ~all(tmpstruct.pos3d == 0)
                        if length(g.vars.coordinates{tmpcond}) >= index && ~isempty(g.vars.coordinates{tmpcond}{index}) && ~isequal(g.vars.coordinates{tmpcond}{index}, tmpstruct.pos3d)
                            g.vars.coordinates{tmpcond}{index}(2,:) = tmpstruct.pos3d;
                        else
                            g.vars.coordinates{tmpcond}{index} = tmpstruct.pos3d;
                        end
                    elseif isfield(tmpstruct, 'eleccoord') && ~all(tmpstruct.eleccoord == 0)
                        if length(g.vars.coordinates{tmpcond}) >= index && ~isempty(g.vars.coordinates{tmpcond}{index}) && ~isequal(g.vars.coordinates{tmpcond}{index}, tmpstruct.eleccoord)
                            g.vars.coordinates{tmpcond}{index}(2,:) = tmpstruct.eleccoord;
                        else
                            g.vars.coordinates{tmpcond}{index} = tmpstruct.eleccoord;
                        end
                    else
                        tmpstruct
                        error('Field not found in tmpstruct');
                    end
                end
            end
            delete(htmp);
        end
        
        xltmp = xlim;
        yltmp = ylim;
        g.vars.dimratio = (xltmp(2) - xltmp(1)) / (yltmp(2) - yltmp(1));
        
        axis(g.vars.hBrain(tmpcond),'off');
        if ~isempty(g.condtitle)
            h = title(g.vars.hBrain(tmpcond),g.condtitle(tmpcond,:),'color','w');
            if ~isempty(g.condtitleformat)
                set(h, g.condtitleformat{:} );
            end
        end
        
        % Create Footer Panel axis
        if ~isempty( g.footerPanelData )
            g.vars.hFooterPanel(tmpcond) = axes('position', [0/nbconditions+maxcoordx/nbconditions*(tmpcond-1), 0, ...
                maxcoordx/nbconditions-0.05/nbconditions, ordinate-0.1].*s+q,'visible', g.visible,'color','none','xcolor','w','ycolor','w','parent',g.figurehandle);
            
            minordinate = min(g.footerPanelData(:));
            maxordinate = max(g.footerPanelData(:));
            
            axis(g.vars.hFooterPanel(tmpcond),'on'); % set (g.figurehandle, 'visible', g.visible);
            
            hold(g.vars.hFooterPanel(tmpcond),'on')
            
            if any(strcmpi(g.footerPanelPlotMode,'all'))
                % plot individual footer panel data traces

                if any(strcmpi(g.footerPanelPlotMode,'envelope')) && size(g.footerPanelData(:,:,tmpcond),1)>1
                    % plot envelope of footer panel data
                    plot(g.vars.hFooterPanel(tmpcond),g.footerPanelTimes, env(g.footerPanelData(:,:,tmpcond)), 'color',g.envColor, 'linewidth', 3*g.resmult);
                    hold(g.vars.hFooterPanel(tmpcond),'on');
                end
                colors = distinguishable_colors(size(g.footerPanelData,1),[0 0 0; 1 1 1; g.envColor]);
                %colors = hsv(size(g.footerPanelData,1));
                for k=1:size(g.footerPanelData,1)
                    plot(g.vars.hFooterPanel(tmpcond),g.footerPanelTimes,g.footerPanelData(k,:),'linewidth',g.resmult,'color',colors(k,:));
                end
                hold(g.vars.hFooterPanel(tmpcond),'off');
            end
            
            set(g.vars.hFooterPanel(tmpcond), 'ylim', [minordinate maxordinate]);
            set(g.vars.hFooterPanel(tmpcond), 'xlim', [g.footerPanelTimes(1) g.footerPanelTimes(end)]);
            set(g.vars.hFooterPanel(tmpcond),'xcolor','w','ycolor','w');
            set(g.vars.hFooterPanel(tmpcond),'color','none');
            
            % draw event markers
            if ~isempty(g.events)
                for i=1:length(g.events)
                    events = g.events{i};
                    
                    % set defaults
                    if length(events) < 5
                        events{5} = '';     end
                    if length(events) < 4
                        events{4} = 2;      end
                    if length(events) < 3
                        events{3} = ':';    end
                    if length(events) < 2
                        events{2} = 'r';    end
                    
                    
                    lh = vline(events{1},events{2},events{5},0,g.vars.hFooterPanel(tmpcond));
                    set(lh,'linestyle',events{3},'linewidth',events{4}*g.resmult);
                end
            end
            
            xlabel(g.vars.hFooterPanel(tmpcond),'Time (sec)', 'fontweight', 'bold', 'fontsize', 12*g.resmult,'color','w');
            set(g.vars.hFooterPanel(tmpcond), 'box', 'off','fontsize', 10*g.resmult);
            if tmpcond == 1
                ylabel(g.vars.hFooterPanel(tmpcond),g.envylabel, 'fontweight', 'bold', 'fontsize', 12*g.resmult,'color','w');
            end
            title(g.vars.hFooterPanel(tmpcond),g.footerPanelTitle,'color','w','fontsize',12*g.resmult);
            
            hold(g.vars.hFooterPanel(tmpcond),'off');
            
        end
        
        % create a little title bar in right-top corner
        connstr = g.ConnMethod;
        if length(g.collapsedFreqs)>1
            connstr = [connstr ' (' num2str(g.collapsedFreqs(1)) '-' num2str(g.collapsedFreqs(end)) 'Hz) '];
        end
        
        text((maxcoordx+(1.1-maxcoordx)/2)*s(1)+q(1), 0.90, ...
            [fastif(iscell(g.title),g.title,{g.title}) {connstr}], ...
            'HorizontalAlignment','center','fontsize',12*g.resmult, ...
            'units','normalized','fontweight','bold', ...
            'parent',findobj(g.figurehandle,'tag','figureBackground'), ...
            'color','w');
        
        % draw a 'compass' indicating the directions
        if g.makeCompass
            % %- IN PREP
            g.vars.hCompass(tmpcond) = axes(...
                'position', [-0.05/nbconditions+maxcoordx/nbconditions*(tmpcond-1), 0.8631, ...
                0.2/nbconditions-0.05/nbconditions, 0.15].*s+q, ...
                'visible', g.visible,'color','none','xcolor','w', ...
                'ycolor','w','parent',g.figurehandle);
            
            % z-directon (superior)
            [hCompass(1).line,hCompass(1).head]=arrow3d([0 0 0],[0 0 0.1],30,'cylinder',[0.4,0.2],[20,10],g.vars.hCompass(tmpcond));
            set(hCompass(1).head,'FaceColor','r','EdgeColor','r');
            % y-direction (nose)
            [hCompass(2).line,hCompass(2).head]=arrow3d([0 0 0],[0 0.1 0],30,'cylinder',[0.4,0.2],[20,10],g.vars.hCompass(tmpcond));
            set(hCompass(2).head,'FaceColor','g','EdgeColor','g');
            % x-direction (right ear)
            [hCompass(3).line,hCompass(3).head]=arrow3d([0 0 0],[0.1 0 0],30,'cylinder',[0.4,0.2],[20,10],g.vars.hCompass(tmpcond));
            set(hCompass(3).head,'FaceColor',[11 131 222]/255,'EdgeColor',[11 131 222]/255); % blue
            
            % draw labels
            text(0,0,.15,g.compassLabels{1},'color',[0.9 0 0],'parent',g.vars.hCompass(tmpcond), ...
                'fontsize',12*g.resmult,'horizontalalignment','center');
            text(0,.15,-0.01,g.compassLabels{2},'color',[0 0.9 0],'parent',g.vars.hCompass(tmpcond),...
                'fontsize',12*g.resmult,'horizontalalignment','center');
            text(.15,0,-0.01,g.compassLabels{3},'color',[11 131 222]/255-0.04,'parent',g.vars.hCompass(tmpcond),...
                'fontsize',12*g.resmult,'horizontalalignment','center');
            
            axis(g.vars.hCompass(tmpcond),'equal','off');
            patches = [hCompass.head];
            set(patches(2,:),'edgecolor','k');
            
            %             lh=legend(g.vars.hCompass(tmpcond),patches(1,:),{'Sup','Ant','Rlat'},...
            %                 'textcolor','w','color','none', ...
            %                 'orientation','vertical','location','WestOutside');
            %             pos=get(lh,'position');
            %             set(lh,'position',[-0.02 pos(2:end)]);
            %             legend(g.vars.hCompass(tmpcond),'boxoff');
            
            % link rotation properties of the brain and compass axes
            CameraPosition = campos(g.vars.hBrain(tmpcond));
            CameraUpVector = camup(g.vars.hBrain(tmpcond));
            hlink=linkprop([g.vars.hCompass(tmpcond) g.vars.hBrain(tmpcond)] ,{'CameraPosition','CameraUpVector'});
            key = 'graphics_linkprop';
            setappdata(g.vars.hBrain(tmpcond),key,hlink);
            campos(g.vars.hBrain(tmpcond),CameraPosition);
            camup(g.vars.hBrain(tmpcond),CameraUpVector);
            camproj(g.vars.hCompass(tmpcond),'perspective');
            axis(g.vars.hCompass(tmpcond),'vis3d');
        end
        
    end
    
%     % adjust the colormap as needed
%     if strcmpi(g.nodeColorPolarity, 'pos')
%         g.nodeColormap = g.nodeColormap(length(g.nodeColormap)/2:end,:); end
%     if strcmpi(g.edgeColorPolarity, 'pos')
%         g.edgeColormap = g.edgeColormap(length(g.edgeColormap)/2:end,:); end
    
    
    % draw captions if necessary
    % --------------------------
    countl = 1;
    
    if g.caption
        
        xlimnorm = (1.1-maxcoordx)/(maxcoordx/nbconditions) * g.xlimaxes;
        ylimnorm = 0.45/(1-ordinate) * g.ylimaxes;
        
        % create 3 spheres to show NodeSize variation
        % --------------------------------------------
        if strcmpi(g.modulateNodeSize,'on')
            
            g.vars.hlgnd(countl) = axes('position', ...
                [maxcoordx+(1.1-maxcoordx)/2, -0.05, (1.1-maxcoordx)/2, 0.3 ].*s+q,  ...
                'xlim', xlimnorm, 'ylim', ylimnorm,'visible', g.visible, 'color', 'w', ...
                'parent',g.figurehandle);
            
            
            %                 small   = max(0.1,g.nodeSizeDataRange(1));
            %                 mid     = (g.nodeSizeDataRange(2)+g.nodeSizeDataRange(1))/2;
            %                 large   = g.nodeSizeDataRange(2);
            %
            %                 small = scalevalue(small, g.nodeSizeDataRange, g.nodeSizeLimits);
            %                 small = 0.05*small*10;
            %                 small = g.diskscale*small*10;
            
            %                 mid = scalevalue(mid, g.nodeSizeDataRange, g.nodeSizeLimits);
            %                 mid = 0.05*mid*10;
            %                 mid = g.diskscale*mid*10;
            
            %                 large = scalevalue(large, g.nodeSizeDataRange, g.nodeSizeLimits);
            %                 large = 0.05*large*10;
            %                 large = g.diskscale*large*10;
            
            small = 0.1;
            mid   = 0.5;
            large = 0.8;
            
            % draw 3 spheres
            [xstmp ystmp zs] = sphere(15);
            l=sqrt(xstmp.*xstmp+ystmp.*ystmp+zs.*zs);
            normals = reshape([xstmp./l ystmp./l zs./l],[16 16 3]);
            nodeSizeScaled = large;
            xs2 = nodeSizeScaled*xstmp; ys2 = nodeSizeScaled*ystmp; zs2 = nodeSizeScaled*zs + 2;
            nodeSizeScaled = mid;
            xs1 = nodeSizeScaled*xstmp; ys1 = nodeSizeScaled*ystmp; zs1 = nodeSizeScaled*zs;
            nodeSizeScaled = small;
            xs3 = nodeSizeScaled*xstmp; ys3 = nodeSizeScaled*ystmp; zs3 = nodeSizeScaled*zs - 1.5;
            colorarray = repmat(reshape([1 1 1],  1,1,3), [size(zs,1) size(zs,2) 1]);
            handles = surf(g.vars.hlgnd(countl),xs1+1, ys1, zs1, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                'backfacelighting', 'lit', 'facelighting', g.facelighting, 'facecolor', 'interp', 'ambientstrength', 0.3); hold on;
            handles = surf(g.vars.hlgnd(countl),xs2+1, ys2, zs2, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                'backfacelighting', 'lit', 'facelighting', g.facelighting, 'facecolor', 'interp', 'ambientstrength', 0.3);
            handles = surf(g.vars.hlgnd(countl),xs3+1, ys3, zs3, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                'backfacelighting', 'lit', 'facelighting', g.facelighting, 'facecolor', 'interp', 'ambientstrength', 0.3);
            axis(g.vars.hlgnd(countl),'off');
            %                 lh(1)=camlight('left');
            %                 lh(2)=camlight('right');
            %                 view(g.vars.hlgnd(countl),[1 0 0]);
            lh=lightangle(0,45);
            set(lh,'parent',g.vars.hlgnd(countl));
            %                 lighting(g.facelighting);
            %                 material shiny;
            axis(g.vars.hlgnd(countl),'equal');
            set(g.vars.hlgnd(countl), 'zlim', [-2 3]);
            set(g.vars.hlgnd(countl), 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'off');
            
            if g.displayLegendLimitText
                text(-0.5, 0.42, {'',[g.NodeSizeMapping '  ']}, 'units','normalized','horizontalalignment','center','rotation',90,'fontsize', 12, 'fontweight', 'normal','parent',g.vars.hlgnd(countl));
                text(0, 1, 2, num2str(g.nodeSizeDataRange(2),'%0.2f'), 'fontweight', 'bold','parent',g.vars.hlgnd(countl),'tag','NodeSizeMax');    % upper limit
                text(0, 1, 0, num2str((g.nodeSizeDataRange(2)+g.nodeSizeDataRange(1))/2,'%0.2f'), 'fontweight', 'bold','parent',g.vars.hlgnd(countl),'tag','NodeSizeMid');   % midrange
                text(0, 1, -1.5, num2str(g.nodeSizeDataRange(1),'%0.2f'), 'fontweight', 'bold','parent',g.vars.hlgnd(countl),'tag','NodeSizeMin'); % lower limit
            else
                text(-0.3, 0.42, {'',[g.NodeSizeMapping '  ']}, 'units','normalized','horizontalalignment','center','rotation',90,'fontsize', 12, 'fontweight', 'normal','parent',g.vars.hlgnd(countl));
            end
            
            
            %       renderNodeSizeLegend(mean(xlimnorm), min(ylimnorm)+0.2, g); % see function at the end
            %       axis off;
            
            %                 set(g.vars.hlgnd(countl),'ylim',[0 5]);
            
            countl = countl + 1;
            
        end
        
        % create ball colormap for NodeColor variation
        % --------------------------------------------
        if strcmpi(g.modulateNodeColor,'on');
            
            
            g.vars.hlgnd(countl) = axes('position', [maxcoordx+(1.1-maxcoordx)/2, 0.29 , (1.1-maxcoordx)/2, 0.14].*s+q, ...
                'visible', g.visible, 'color', 'none' );
            
            
            if strcmpi(g.nodeColorPolarity, 'posneg') % negative ITCs (difference only) ?
                cbar([-1 1], [-1 1], g.nodeColormap, 'vert', 'circle', g, g.vars.hlgnd(countl));
                %           ylabel(g.modulateNodeColor, 'fontweight', 'bold');
                %           set(gca, 'ytick', [-1 0 1], 'yticklabel', [-1 0 1], 'xticklabel', [], 'box', 'off');
            else
                cbar([0 1], [0 1], g.nodeColormap, 'vert', 'circle', g, g.vars.hlgnd(countl));
                %           cbar( [0 1], [0 1], g.nodeColormap(length(g.nodeColormap)/2:end,:), 'vert', 'circle', g);
                %           ylabel({'',g.modulateNodeColor},'rotation',90, 'fontweight', 'normal');
                %           set(gca, 'ytick', [0 1], 'yticklabel', [0 1], 'xticklabel', [], 'box', 'off');
                %           set(gca, 'ytick', [0], 'yticklabel', [0], 'xticklabel', [], 'box', 'off');
            end
            axis(g.vars.hlgnd(countl),'off');
            set(g.vars.hlgnd(countl), 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'off');
            % yticks
            
            if g.displayLegendLimitText
                text(-0.2, 1, num2str(g.nodeColorDataRange(2),'%0.2f'), 'fontsize', 10,'units','normalized','parent',g.vars.hlgnd(countl),'tag','NodeColorMax'); % upper limit
                text(-0.2, 0.5, num2str((g.nodeColorDataRange(2)+g.nodeColorDataRange(1))/2, '%0.2f'), 'fontsize', 10,'units','normalized','parent',g.vars.hlgnd(countl),'tag','NodeColorMid'); % midrange
                text(-0.2, 0, num2str(g.nodeColorDataRange(1),'%0.2f'), 'fontsize', 10,'units','normalized','parent',g.vars.hlgnd(countl),'tag','NodeColorMin'); % lower limit
                text(-0.5, 0.55, {'',[g.NodeColorMapping ' ']}, 'horizontalalignment','center','rotation',90,'fontsize', 12, 'fontweight', 'normal','units','normalized','parent',g.vars.hlgnd(countl));
            else
                text(-0.3, 0.55, {'',[g.NodeColorMapping ' ']}, 'horizontalalignment','center','rotation',90,'fontsize', 12, 'fontweight', 'normal','units','normalized','parent',g.vars.hlgnd(countl));
            end
                
            
            countl = countl + 1;
        end
        
        % create 3 cylinders for EdgeSize variation
        % --------------------------------------------
        if strcmpi(g.modulateEdgeSize,'on')
            
            g.vars.hlgnd(countl) = axes('position', [maxcoordx+(1.1-maxcoordx)/2, 0.72,(1.1-maxcoordx)/2, 0.20 ].*s+q, ...
                'visible', g.visible, 'parent', g.figurehandle, 'units', 'normalized');
            
            small   = 0.01;
            mid     = 0.05;
            large   = 0.09;
            
            %                 small   = max(0.001,g.edgeSizeDataRange(1));
            %                 mid     = (g.edgeSizeDataRange(2)+g.edgeSizeDataRange(1))/2;
            %                 large   = g.edgeSizeDataRange(2);
            
            %
            %                 small = scalevalue(small, g.nodeSizeDataRange, g.nodeSizeLimits);
            %                 small = 0.05*small*10;
            %                 small = small*10;
            
            %                 mid = scalevalue(mid, g.nodeSizeDataRange, g.nodeSizeLimits);
            %                 mid = 0.05*mid*10;
            %                 mid = mid*10;
            
            %                 large = scalevalue(large, g.nodeSizeDataRange, g.nodeSizeLimits);
            %                 large = 0.05*large*10;
            %                 large = large*10;
            
            [xc yc zc] = cylinder( small, 20 );   % create unit-length cylinder
            colorarray = repmat(reshape([1 1 1],  1,1,3), [size(zc,1) size(zc,2) 1]);
            handles = surf(g.vars.hlgnd(countl),xc, yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', ...
                'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', g.facelighting);
            [xc yc zc] = adjustcylinder2( handles, [0 0 0], [1 0 0] );  % stretch and rotate cylinder to match start-end pnts
            
            hold(g.vars.hlgnd(countl),'on');
            
            [xc yc zc] = cylinder( mid, 20 );   % create unit-length cylinder
            handles = surf(g.vars.hlgnd(countl),xc, yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', ...
                'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', g.facelighting);
            [xc yc zc] = adjustcylinder2( handles, [0 0.5 0], [1 0.5 0] );  % stretch and rotate cylinder to match start-end pnts
            
            [xc yc zc] = cylinder( large, 20 );   % create unit-length cylinder
            handles = surf(g.vars.hlgnd(countl),xc, yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', ...
                'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', g.facelighting);
            [xc yc zc] = adjustcylinder2( handles, [0 1 0], [1 1 0] );  % stretch and rotate cylinder to match start-end pnts
            
            hold(g.vars.hlgnd(countl),'off');
            
            axis(g.vars.hlgnd(countl),'off');
            %                 lh(1)=camlight('left');
            lh = light('parent',g.vars.hlgnd(countl));
            lh = camlight(lh,'right');
            view(g.vars.hlgnd(countl),[0 0 1]);
            camorbit(g.vars.hlgnd(countl),-25,0,[0 1 0]);
            camproj(g.vars.hlgnd(countl),'perspective');
            
            %                 lh(3)=lightangle(45,0);
%             set(lh,'parent',g.vars.hlgnd(countl));
            %                 lighting(g.facelighting);
            %                 material shiny;
            axis(g.vars.hlgnd(countl),'equal');
            
            %                 renderEdgeSizeLegend([0.02 1], [0.04 0.96], 5, g, g.vars.hlgnd(countl)); % see function at the end
            
            if g.displayLegendLimitText
                text(-0.3, 0.9, num2str(g.edgeSizeDataRange(2),'%0.2f'), 'fontsize', 10,'units','normalized','color','w','parent',g.vars.hlgnd(countl),'tag','EdgeSizeMax'); % upper limit
                text(-0.3, 0.5, num2str((g.edgeSizeDataRange(2)+g.edgeSizeDataRange(1))/2, '%0.2f'), 'fontsize', 10,'units','normalized','color','w','parent',g.vars.hlgnd(countl),'tag','EdgeSizeMid'); % midrange
                text(-0.3, 0.1, num2str(g.edgeSizeDataRange(1),'%0.2f'), 'fontsize', 10,'units','normalized','color','w','parent',g.vars.hlgnd(countl),'tag','EdgeSizeMin'); % lower limit
                text(-0.6, 0.55, {'',[g.EdgeSizeMapping ' ']}, 'horizontalalignment','center','rotation',90,'fontsize', 12, 'fontweight', 'normal','units','normalized','parent',g.vars.hlgnd(countl));
            else
                text(-0.3, 0.55, {'',[g.EdgeSizeMapping ' ']}, 'horizontalalignment','center','rotation',90,'fontsize', 12, 'fontweight', 'normal','units','normalized','parent',g.vars.hlgnd(countl));
            end
            
            set(g.vars.hlgnd(countl), 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'off');
            
            ylim(g.vars.hlgnd(countl),[-0.1 1.1]);
            xlim(g.vars.hlgnd(countl),[-0.4 1]);
            
            countl = countl + 1;
            
            % create rectangular colormap for EdgeColor variation
            % --------------------------------------------------
            if strcmpi(g.modulateEdgeColor,'on')
                
                g.vars.hlgnd(countl) = axes('position', [maxcoordx+(1.1-maxcoordx)/2, 0.57 , (1.1-maxcoordx)/2, 0.14].*s+q, ...
                    'visible', g.visible, 'color', 'none', 'parent',g.figurehandle);
                
                % NOTE: enable/disable lines below this if legend coloring is not as expected
%                 ylim(g.vars.hlgnd(countl),[0 1]);
%                 xlim(g.vars.hlgnd(countl),[0 1]);
                
                if strcmpi(g.edgeColorPolarity, 'posneg')
                    cbar( [-1 1], [-1 1], g.edgeColormap, 'vert', '', g, g.vars.hlgnd(countl));
                else
                    cbar( [0.1 1], [0 1], g.edgeColormap, 'vert', '', g, g.vars.hlgnd(countl));
                end
                
                if g.displayLegendLimitText
                    text(-0.3, 1, num2str(g.edgeColorDataRange(2),'%0.2f'), 'fontsize', 10,'units','normalized','parent',g.vars.hlgnd(countl),'tag','EdgeColorMax'); % upper limit
                    text(-0.3, 0.5, num2str((g.edgeColorDataRange(2)+g.edgeColorDataRange(1))/2, '%0.2f'), 'fontsize', 10,'units','normalized','parent',g.vars.hlgnd(countl),'tag','EdgeColorMid'); % midrange
                    text(-0.3, 0, num2str(g.edgeColorDataRange(1),'%0.2f'), 'fontsize', 10,'units','normalized','parent',g.vars.hlgnd(countl),'tag','EdgeColorMin'); % lower limit
                    text(-0.5, 0.55, {'',[g.EdgeColorMapping ' ']}, 'horizontalalignment','center','rotation',90,'fontsize', 12, 'fontweight', 'normal','units','normalized','parent',g.vars.hlgnd(countl));
                else
                    text(-0.3, 0.55, {'',[g.EdgeColorMapping ' ']}, 'horizontalalignment','center','rotation',90,'fontsize', 12, 'fontweight', 'normal','units','normalized','parent',g.vars.hlgnd(countl));
                end
                
                set(g.vars.hlgnd(countl), 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'off');
                
            end
            
        end
    else
        maxcoordx = 1;
    end
    
    % Determine number of dual symmetric dipoles
    g.vars.NumDupNodes  = nnz(cellfun(@(x)size(x,1),g.vars.coordinates{1})==2);

    % Initialize node text
    % --------------------
    g.vars.INIT_LABELS = true;
    g.vars.hNodeText = nan(length(SELECTED),1+(g.vars.NumDupNodes>0));
    
    
    % Set the renderer
    % -----------------
    if strcmpi(g.opengl,'on') && ~strcmpi(get(g.figurehandle,'renderer'),'opengl')
        set(g.figurehandle, 'renderer', 'opengl');
    end
    
end % initialization block

alltimepoints = g.vars.alltimepoints;

% return now if we are only initializing the brainmovie
if strcmpi(g.mode,'init')
    return;
end


% Main loop, draw frames
% ----------------------------
for indeximage = g.vars.alltimepoints
    
    switch g.verb
        case 1
            fprintf('Processing image %d\n', indeximage);
        case 2
            waitbar(indeximage/length(g.vars.alltimepoints),g.vars.hwaitbar,sprintf('Rendering timepoint (%d/%d)...',indeximage,length(g.vars.alltimepoints)));
    end
    
    
    % produce background event flash
    % ------------------------------
    if ~isempty(g.vars.flashTimes)
        %axes(g.vars.hFigureBg); set (g.figurehandle, 'visible', g.visible);
        if ~isempty(find(indeximage == g.vars.allFlashIndices, 1))
            posf = find(indeximage == g.vars.allFlashIndices);
            set(g.vars.hFlashPatch, 'facecolor', g.vars.flashColor{posf}, 'visible','on');
        elseif posf == 0 % allow the color to stay 2 images
            set(g.vars.hFlashPatch, 'facecolor', g.backcolor,'visible','off');
        else
            posf = 0;
        end
    end
    
    for tmpcond=1:nbconditions
        g.vars.hbrainax = g.vars.hBrain(tmpcond);
        
        if ~g.speedy
            set (g.figurehandle, 'visible', g.visible);
        end
        
        themeopts =  hlp_struct2varargin(g.theme.graph);
        
        % draw edges
        % -----------------
        if strcmpi(g.modulateEdgeSize,'on')
            
            % Build the cylinders (connections/edges)
            q = 1;
            NumPointsCircum = 11; % num points around circumference of cylinder
            NumNodes = length(SELECTED);
            % NOTE: we are pre-allocating memory for too many edges --
            % should only allocate for non-NaN edges
            NumEdges = ((NumNodes+g.vars.NumDupNodes)^2)-NumNodes;  %NumNodes^2-NumNodes+g.vars.NumDupNodes*(NumNodes-1); %((NumNodes+g.vars.NumDupNodes)^2)-NumNodes;
            [xc yc zc] = deal(NaN(3*NumEdges,NumPointsCircum));
            [colorarray normals] = deal(NaN(size(xc,1), NumPointsCircum, 3));
            cylwidth = NaN(1,NumEdges);
            
            for index1 = SELECTED
                for index2 = SELECTED
                    
                    if index1==index2 ...
                            || (~g.causality && index2 < index1) ... % if connectivity is symm
                            || all(isnan(ALLEDGESIZE{ index1, index2, tmpcond}(indeximage))) % or if connection is absent
                        continue;
                    end
                    
                    % collapse data across frequency
                    [tmpEdgeSize tmpEdgeSize_peakidx] = hlp_collapseFrequencies( ...
                        ALLEDGESIZE{ index1, index2, tmpcond}, ...
                        g.collapsefun,FREQS,indeximage);
                    
                    if g.EdgeColorMappedToDirectionality
                        tmpEdgeColor = sign(index1-index2);
                        if ~tmpEdgeColor, tmpEdgeColor=1; end
                    else
                        % collapse data across frequency
                        [tmpEdgeColor tmpEdgeColor_peakidx] = hlp_collapseFrequencies( ...
                            ALLEDGECOLOR{index1, index2, tmpcond}, ...
                            g.collapsefun,FREQS,indeximage);
                    end
                    
                    pos1 = g.vars.coordinates{tmpcond}{ index1 };
                    pos2 = g.vars.coordinates{tmpcond}{ index2 };
                    curvature = g.circfactor(index1, index2);
                    
                    % create edge
                    [xc(q:q+1,:) yc(q:q+1,:) zc(q:q+1,:) colorarray(q:q+1,:,:) cylwidth normals(q:q+1,:,:)] ...
                        = makeEdge( pos1(1,:), pos2(1,:), ...
                        tmpEdgeSize, tmpEdgeColor, curvature, NumPointsCircum,g);
                    
                    q=q+3;
                    
                    if size(pos1,1)==2
                        % create second edge for dual equiv dipole
                        [xc(q:q+1,:) yc(q:q+1,:) zc(q:q+1,:) colorarray(q:q+1,:,:) cylwidth normals(q:q+1,:,:)] ...
                            = makeEdge( pos1(2,:), pos2, ...
                            tmpEdgeSize, tmpEdgeColor, curvature, NumPointsCircum,g);

                        q=q+3;
                    end
                    
                    if size(pos2,1)==2
                        % create second edge for dual equiv dipole
                        [xc(q:q+1,:) yc(q:q+1,:) zc(q:q+1,:) colorarray(q:q+1,:,:) cylwidth normals(q:q+1,:,:)] ...
                            = makeEdge( pos1, pos2(2,:), ...
                            tmpEdgeSize, tmpEdgeColor, curvature, NumPointsCircum,g);

                        q=q+3;
                    end
                    
%                     %% DEBUG
%                     g.vars.hConnections = surf(g.vars.hbrainax, xc, yc, zc, colorarray, ...
%                     'tag', 'tmpmov_connections', 'edgecolor', 'none', ...
%                     'backfacelighting', 'lit', 'facecolor', 'interp', ...
%                     'facelighting', g.facelighting, 'vertexnormals',normals, ...
%                     themeopts{:});
%                     
%                     if index1 > 6 && index2 > 5
%                         disp(''); end
%                     
%                     %% DEBUG
                end
            end
            
            % Render all edges at once
            if ~isfield(g.vars,'hConnections') || isempty(g.vars.hConnections)
                g.vars.hConnections = surf(g.vars.hbrainax, xc, yc, zc, colorarray, ...
                    'tag', 'tmpmov_connections', 'edgecolor', 'none', ...
                    'backfacelighting', 'lit', 'facecolor', 'interp', ...
                    'facelighting', g.facelighting, 'vertexnormals',normals, ...
                    themeopts{:});
            else
                set(g.vars.hConnections,'Xdata',xc,'Ydata',yc,'Zdata',zc,'Cdata',colorarray,'VertexNormals',normals);
%                 set(g.vars.hConnections,'Xdata',xc(1:q,:),'Ydata',yc(1:q,:),'Zdata',zc(1:q,:),'Cdata',colorarray(1:q,:,:),'VertexNormals',normals(1:q,:,:));
            end
            
            
            
%             if strcmpi(g.project3d, 'on')
%                 if g.colorshadow
%                     colorarray=colorarray.*repmat(reshape(g.projcolor, 1,1,3), [size(zc,1) size(zc,2) 1]);
%                 else
%                     colorarray = repmat(reshape(g.projcolor, 1,1,3), [size(zc,1) size(zc,2) 1]);
%                 end
%                 surf(g.vars.hbrainax,xc, yc, g.factproj(3)*ones(size(zc)), colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
%                 surf(g.vars.hbrainax,xc, g.factproj(2)*ones(size(yc)), zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
%                 surf(g.vars.hbrainax,g.factproj(1)*ones(size(xc)), yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
%             end
        
        end
        
        % draw nodes
        % ------------
        
        themeopts =  hlp_struct2varargin(g.theme.graph);
        
        % Build the nodes
        q = 1;
%         NumPointsCircum = 11; % num points around circumference of cylinder
        NumNodes = length(SELECTED);
        nodeSizeScaled  = nan(1,NumNodes);
        nodeColorScaled = nan(NumNodes,3);
        [xs ys zs] = deal(NaN(17*NumNodes+g.vars.NumDupNodes,16));
        [colorarray normals] = deal(NaN(size(xs,1), size(xs,2), 3));
        
        
        for index1=1:length(SELECTED);
            
            % node size
            [tmpNodeSize tmpNodeSize_peakidx] = hlp_collapseFrequencies(ALLNODESIZE{ index1, tmpcond}, ...
                g.collapsefun,FREQS,indeximage);
            
            % node color
            [tmpNodeColor tmpNodeColor_peakidx] = hlp_collapseFrequencies(ALLNODECOLOR{ index1, tmpcond}, ...
                g.collapsefun,FREQS,indeximage);
            
            % render node
            pos = g.vars.coordinates{tmpcond}{ index1 };
            [xs(q:q+15,:) ys(q:q+15,:) zs(q:q+15,:) colorarray(q:q+15,:,:) ...
               nodeColorScaled(index1,:) nodeSizeScaled(index1) normals(q:q+15,:,:)] = makeNode(pos(1,:),tmpNodeSize,tmpNodeColor,g);
           
            q = q + 16;
            
            if size(pos,1)==2
                % render dual symmetric node
                pos = g.vars.coordinates{tmpcond}{ index1 };
                [xs(q:q+15,:) ys(q:q+15,:) zs(q:q+15,:) colorarray(q:q+15,:,:) ...
                   nodeColorScaled(index1,:) nodeSizeScaled(index1) normals(q:q+15,:,:)] = makeNode(pos(2,:),tmpNodeSize,tmpNodeColor,g);

                q = q + 16;
            end
            
        end
        
        if ~isfield(g.vars,'hNodes') || isempty(g.vars.hNodes)
            g.vars.hNodes = surf(g.vars.hbrainax,xs, ys, zs, colorarray, 'tag', 'tmpmov', ...
                'EdgeColor','none', 'vertexnormals',normals, ...
                'backfacelighting', 'lit', 'facelighting', g.facelighting, ...
                'facecolor', 'interp',themeopts{:});
            
            if g.plotNodeLabels && g.vars.INIT_LABELS
                
                g.vars.INIT_LABELS = false;
                
                for index1=1:length(SELECTED)
                    
                    nodeCenterPos = g.vars.coordinates{tmpcond}{index1};
                    
                    g.vars.hNodeText(index1,1) = text( ...
                        nodeCenterPos(1,1), ...
                        nodeCenterPos(1,2), ...
                        double(nodeCenterPos(1,3)+1.2*nodeSizeScaled(index1)), ...
                        g.nodelabels{index1},'color','w','parent',g.vars.hbrainax,'fontsize',12*g.resmult,'tag','nodelabel','interpreter','none');
                    
                    if size(nodeCenterPos,1)==2
                        % draw label for symmetric dipole
                        g.vars.hNodeText(index1,2) = text( ...
                            nodeCenterPos(2,1), ...
                            nodeCenterPos(2,2), ...
                            double(nodeCenterPos(2,3)+1.2*nodeSizeScaled(index1)), ...
                            g.nodelabels{index1},'color','w','parent',g.vars.hbrainax,'fontsize',12*g.resmult,'tag','nodelabel');
                    end
                end
            end
 
        else
            % update node surfaces
            set(g.vars.hNodes,'Xdata',xs,'Ydata',ys,'Zdata',zs,'Cdata',colorarray,'vertexnormals',normals);
            
            % update nodelabel position
            if g.plotNodeLabels && ~g.vars.INIT_LABELS
                for index1 = 1:length(SELECTED)
                    nodeCenterPos = g.vars.coordinates{tmpcond}{ index1 };
                    
                    pos = get(g.vars.hNodeText(index1,1),'Position');
                    set(g.vars.hNodeText(index1,1),'Position',[pos(1:2) nodeCenterPos(1,3)+1.2*nodeSizeScaled(index1)] );
                    % handle dual dipole, if present
                    if size(g.vars.hNodeText,2)>1 && ~isnan(g.vars.hNodeText(index1,2))
                        pos = get(g.vars.hNodeText(index1,2),'Position');
                        set(g.vars.hNodeText(index1,2),'Position',[pos(1:2) nodeCenterPos(2,3)+1.2*nodeSizeScaled(index1)] );
                    end
                end
            end
            
        end
        
%         if strcmpi(g.project3d, 'on')
%             if g.colorshadow
%                 colorarray=colorarray.*repmat(reshape(g.projcolor, 1,1,3), [size(zs,1) size(zs,2) 1]);
%             else
%                 colorarray = repmat(reshape(g.projcolor, 1,1,3), [size(zs,1) size(zs,2) 1]);
%             end
%             surf(g.vars.hbrainax,xs, ys, g.factproj(3)*ones(size(zs)), colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
%             surf(g.vars.hbrainax,xs, g.factproj(2)*ones(size(ys)), zs, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
%             surf(g.vars.hbrainax,g.factproj(1)*ones(size(xs)), ys, zs, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
%         end
    
        

        % draw a bar for time probability
        % -------------------------------
        %         if ~isempty(g.rt)
        %             if ~isempty(g.rt{tmpcond})
        %                 ll = line([g.rthistloc(1)-g.rthistloc(3)/2 g.rthistloc(1)+g.rthistloc(3)/2], [g.rthistloc(2) g.rthistloc(2)]);
        %                 set(ll, 'linewidth', 2*g.resmult, 'color', 'k');
        %                 barheight = RTdist(tmpcond, indeximage)*g.rthistloc(4);
        %                 x1 = g.rthistloc(1)-0.65*g.rthistloc(3)/2;
        %                 x2 = g.rthistloc(1)+0.65*g.rthistloc(3)/2;
        %                 y1 = g.rthistloc(2);
        %                 y2 = g.rthistloc(2)-barheight;
        %                 ll = patch([x1 x1 x2 x2], [y1 y2 y2 y1], g.rthistcolor, 'linewidth', 2*g.resmult);
        %             end
        %         end
        
        
        % Render Layers
        % ----------------------
        layers = fieldnames(g.Layers);
        
        for layeridx = 1 : length(layers)
            
            curlayer = layers{layeridx};
            
            if ~isstruct(g.Layers.(curlayer)), continue; end
            
            % check if we haven't yet rendered the mesh for this layer...
            if g.Layers.(curlayer).transparency < 1 ...
                    && isempty(findobj(g.figurehandle,'tag',sprintf('%smesh',curlayer))) ...
                    
                % ... and if no mesh present, render it.
                
                children = get(g.figurehandle,'children');
                
                hax = findobj(g.figurehandle,'tag',['brain' num2str(tmpcond)]);
                
                CameraPosition = get(hax,'CameraPosition');
                
                %         set(g.figurehandle,'visible','off')
                hold on;
                if strcmpi(curlayer,'cortex')
                    Handle = hlp_plotAtlas(g.Layers.(curlayer).mesh,hax,g.Layers.(curlayer).color.arg_selection,g.Layers.(curlayer).color.colormapping,hlp_struct2varargin(g.theme.(curlayer)));
                else
                    Handle = hlp_plotmesh(g.Layers.(curlayer).mesh.faces, g.Layers.(curlayer).mesh.vertices,[],false,hax,g.Layers.(curlayer).color,hlp_struct2varargin(g.theme.(curlayer)));
                end
                
                set(Handle, 'facealpha',1-g.Layers.(curlayer).transparency,'tag',sprintf('%smesh',curlayer));
                
                if ~strcmpi(get(g.figurehandle,'Renderer'),'opengl')
                    set(g.figurehandle,'Renderer' ,'opengl');
                end
                
                %         set(g.figurehandle,'visible',g.visible)
                
                set(g.figurehandle,'children',children);
                
                set(hax,'CameraPosition',CameraPosition);
                
            end
            
        end
        
%         if ~isempty(g.title) && tmpcond == 1
%             if iscell(g.title), g.title = g.title{1}; end
%             set(g.figurehandle,'Name',['BrainMovie3D ' g.title]);
%         end
        
        %
        % update view
        % ----------------------------
        if ~isempty(g.rotationpath3d)
            
            if isempty(g.rotationpath3d.FramesPerCycle)
                g.rotationpath3d.FramesPerCycle = length(g.vars.alltimepoints);
            end
            
            angle = (g.vars.path3d_frameindex-1)/g.rotationpath3d.FramesPerCycle*360;
            camorbit(g.vars.hbrainax, cos(angle/180*pi)*g.rotationpath3d.AngleFactor, sin(angle/180*pi)*g.rotationpath3d.PhaseFactor );
            
            if length(g.vars.alltimepoints)==1
                % there is just a single frame
                g.vars.path3d_frameindex = g.vars.path3d_frameindex + 1;
            else
                g.vars.path3d_frameindex = indeximage;
            end
        end
        
        
        %                 if ~isempty(g.title) & tmpcond == 1
        %                     t = textsc(g.title{1},'title');
        %                     set(t,'VerticalAlignment','top', 'fontsize', 15);
        %                 end
        
        % update the Footer Panel
        % -----------------------
        if ~isempty( g.footerPanelData )
            
            hold(g.vars.hFooterPanel(tmpcond),'on');
            
            % draw line for current time point
            [dummy curtime] = getindex(g.footerPanelTimes,TIMES(indeximage));
            if isfield(g.vars,'hCurrentTime')
                delete(g.vars.hCurrentTime);
            end
            g.vars.hCurrentTime = vline(curtime,'w','',0,g.vars.hFooterPanel(tmpcond));
            set(g.vars.hCurrentTime,'linewidth', 2*g.resmult);
            
            % draw sliding window rectangle
            if ~isempty(g.windowLength)
                if isfield(g.vars,'hSlidingWindow')
                    delete(g.vars.hSlidingWindow);
                end
                g.vars.hSlidingWindow = hlp_vrect([curtime-g.windowLength/2 curtime+g.windowLength/2], ...
                    'axesHandle',g.vars.hFooterPanel(tmpcond), ...
                    'patchProperties', ...
                    {'FaceColor',[0.7 0.7 1],'FaceAlpha',0.2,'EdgeColor',[0.2 0.2 0.2],'EdgeAlpha',0.5});
            end
            
            hold(g.vars.hFooterPanel(tmpcond),'off');
        end
        
        % update the Legend text
        % -----------------------
        if g.displayLegendLimitText && g.updateLegendText && g.caption
            
            countl = 1;
            
            if strcmpi(g.modulateNodeSize,'on')
                children = get(g.vars.hlgnd(countl),'children');
                
                % update max
                htxt=findobj(children,'tag','NodeSizeMax');
                set(htxt,'String',num2str(g.nodeSizeDataRange(2),'%0.2f'));
                % update mid
                htxt=findobj(children,'tag','NodeSizeMid');
                set(htxt,'String',num2str((g.nodeSizeDataRange(2)+g.nodeSizeDataRange(1))/2, '%0.2f'));
                % update min
                htxt=findobj(children,'tag','NodeSizeMin');
                set(htxt,'String',num2str(g.nodeSizeDataRange(1),'%0.2f'));
                
                countl = countl + 1;
            end
            
            if strcmpi(g.modulateNodeColor,'on')
                children = get(g.vars.hlgnd(countl),'children');
                
                % update max
                htxt=findobj(children,'tag','NodeColorMax');
                set(htxt,'String',num2str(g.nodeColorDataRange(2),'%0.2f'));
                % update mid
                htxt=findobj(children,'tag','NodeColorMid');
                set(htxt,'String',num2str((g.nodeColorDataRange(2)+g.nodeColorDataRange(1))/2, '%0.2f'));
                % update min
                htxt=findobj(children,'tag','NodeColorMin');
                set(htxt,'String',num2str(g.nodeColorDataRange(1),'%0.2f'));
                
                countl = countl + 1;
            end
            
            if strcmpi(g.modulateEdgeSize,'on')
                children = get(g.vars.hlgnd(countl),'children');
                
                % update max
                htxt=findobj(children,'tag','EdgeSizeMax');
                set(htxt,'String',num2str(g.edgeSizeDataRange(2),'%0.2f'));
                % update mid
                htxt=findobj(children,'tag','EdgeSizeMid');
                set(htxt,'String',num2str((g.edgeSizeDataRange(2)+g.edgeSizeDataRange(1))/2, '%0.2f'));
                % update min
                htxt=findobj(children,'tag','EdgeSizeMin');
                set(htxt,'String',num2str(g.edgeSizeDataRange(1),'%0.2f'));
                
                countl = countl + 1;
            end

            if strcmpi(g.modulateEdgeColor,'on')
                children = get(g.vars.hlgnd(countl),'children');
                
                % update max
                htxt=findobj(children,'tag','EdgeColorMax');
                set(htxt,'String',num2str(g.edgeColorDataRange(2),'%0.2f'));
                % update mid
                htxt=findobj(children,'tag','EdgeColorMid');
                set(htxt,'String',num2str((g.edgeColorDataRange(2)+g.edgeColorDataRange(1))/2, '%0.2f'));
                % update min
                htxt=findobj(children,'tag','EdgeColorMin');
                set(htxt,'String',num2str(g.edgeColorDataRange(1),'%0.2f'));
                
                countl = countl + 1;
            end
            
        end
            
        % Set the lighting options
        g.vars.lights = hlp_bm_updateLights(g);
        
        if g.caption && isfield(g.vars,'hlgnd')
            % reset text fontcolor to white
            set(findobj(g.vars.hlgnd,'type','text'),'color',[0.99 0.99 0.99]);
        end
        
    end  % LOOP OVER CONDITIONS
    
    
    % put the time in the right bottom corner
    % --------------------------------------
    if g.showLatency
        if ~isfield(g.vars,'hlatency'),
            g.vars.hlatency = text(0.92, 0.05, sprintf('%0.3g sec', TIMES(indeximage)), 'unit', 'normalized','parent',g.vars.hFigureBg);
            set(g.vars.hlatency, 'fontsize', 12*g.resmult, 'horizontalalignment', 'right', 'tag', 'tmpmov', 'color', 'w');
        else
            set(g.vars.hlatency,'String', sprintf('%0.3g sec', TIMES(indeximage)));
        end
        %         uistack(hlatency,'up',10);
    end
    
    
    drawnow;
    
    
    if ~isempty(g.moviename)
        % save the file for a movie
        % -------------------------
        movframes = getframe(g.figurehandle);
        mov = addframe(mov,movframes);
    end
    
    if ~isempty(g.framefolder)
        fname = sprintf('%simage%4.4d.%s', g.framefolder, indeximage,g.framesout);
        switch(g.framesout)
            case {'pdf', 'eps', 'png','tif','jpg','bmp'}
                export_fig(fname,g.figurehandle,'-nocrop');
            otherwise
                saveas(g.figurehandle,fname);
                pause(0.2);
        end
        
    end
end

if ~isempty(g.moviename)
    mov = close(mov);
end

if g.verb==2
    close(g.vars.hwaitbar);
end

% EOF






% function to construct a node
% ----------------------------
function [xs ys zs colorarray nodeColorScaled nodeSizeScaled normals] = makeNode( nodeCenterPos, nodeSize, nodeColor, g);
% nodeCenterPos     center coord of the sphere
% nodeSize          radius
% nodeColor         color
% g                 preferences

% deal with dual dipoles
% ----------------------
% handles = [];
% if size(nodeCenterPos,1) > 1
%     [nodeSizeScaled1, nodeColorScaled1, handles1] = drawcircle( nodeCenterPos(1,:), nodeSize, nodeColor, g);
%     [nodeSizeScaled2, nodeColorScaled2, handles2] = drawcircle( nodeCenterPos(2,:), nodeSize, nodeColor, g);
%     nodeSizeScaled  = [ nodeSizeScaled1 nodeSizeScaled2];
%     nodeColorScaled = [ nodeColorScaled1; nodeColorScaled2];
%     handles  = [ handles1 handles2];
%     return;
% end

[xs ys zs] = deal(nan(16));
colorarray = nan([16 16 3]);
nodeColorScaled = nodeColor;
normals= nan([16 16 3]);

% Node size
nodeSizeScaled = scalevalue(nodeSize, g.nodeSizeDataRange, g.nodeSizeLimits);
if strcmpi(g.modulateNodeSize, 'off') || nodeSize==eps, nodeSizeScaled = g.nodeSizeLimits(1); end
nodeSizeScaled = 0.05*nodeSizeScaled*10;
nodeSizeScaled = g.diskscale*nodeSizeScaled*100;

if isnan(nodeSize) || nodeSizeScaled == 0
    return;
end

% Node color
nodeColorScaled = scalevalue(nodeColor, g.nodeColorDataRange, g.nodeColorLimits);
indexcolor = floor(max(nodeColorScaled-0.001,0)*length(g.nodeColormap))+1;
nodeColorScaled   = g.nodeColormap( indexcolor,: );
if isnan(nodeColor) || strcmpi(g.modulateNodeColor, 'off'), nodeColorScaled = [1 1 1]; end %g.nodeColorLimits(1);


if nodeSizeScaled > 0
    [xstmp ystmp zs] = sphere(15);
    l=sqrt(xstmp.*xstmp+ystmp.*ystmp+zs.*zs);
    normals = reshape([xstmp./l ystmp./l zs./l],[16 16 3]);
    xs = nodeCenterPos(1) + nodeSizeScaled*ystmp*g.vars.dimratio;
    ys = nodeCenterPos(2) + nodeSizeScaled*xstmp;
    zs = nodeCenterPos(3) + nodeSizeScaled*zs;
    colorarray = repmat(reshape(nodeColorScaled, 1,1,3), [size(zs,1) size(zs,2) 1]);
end


% function to construct an edge
% -----------------------------
function [xc yc zc colorarray cylwidth normals] = makeEdge( pos1, pos2, edgeSize, edgeColor, edgeCurvature, numPointsCircum,g)
% pos1, pos2		position of the points
% edgeSize       coherence power for width of the line
% edgeColor       coherence angle for color and speed of the line
% cirfact           curvature of the line
% g                 preference
% arrow should point from pos2 to pos1

[xc yc zc]              = deal(nan(2,numPointsCircum));
[colorarray normals]    = deal(nan(2,numPointsCircum,3));
cylwidth = 0;

if edgeSize == 0, return; end

% % deal with dual dipoles
% % ----------------------
% if size(pos1,1) > 1
%     % draw first cylinder
%     [xc yc zc colorarray cylwidth normals] = ...
%         drawconnections( pos1(1,:), pos2, edgeSize, edgeColor, edgeCurvature, g);
%
%     % insert NaNs
%     [xc(end+1,:) yc(end+1,:) zc(end+1,:)] = deal(nan(1,size(xc,2)));
%
%     % draw second cylinder
%     [xc(end+1:end+3,:) yc(end+1:end+3,:) zc(end+1:end+3,:) colorarray cylwidth normals] = ...
%         drawconnections( pos1(2,:), pos2, edgeSize, edgeColor, edgeCurvature, g);
%
%     handles = [ handles1 handles2];
%     return;
% end
%
% if size(pos2,1) > 1
%     handles1 = drawconnections( pos1, pos2(1,:), edgeSize, edgeColor, edgeCurvature, g);
%     handles2 = drawconnections( pos1, pos2(2,:), edgeSize, edgeColor, edgeCurvature, g);
%     handles = [ handles1 handles2];
%     return;
% end

% if the two nodes are too close then do not draw the line
% --------------------------------------------------------
distance = sqrt(sum((pos1-pos2).^2));
if distance < 0.05*(g.ylimaxes(2) - g.ylimaxes(1))
    return;
end

if strcmpi(g.crossfcoh, 'off') || isnan(edgeSize)
    return; 
end

% Edge size
% ---------
tmpthick = scalevalue( edgeSize, g.edgeSizeDataRange, g.edgeSizeLimits );
tmpthick = 30*tmpthick;
if edgeSize == 0, tmpthick = 0; end

% Edge color
% ----------
nodeColorScaled = scalevalue( edgeColor, g.edgeColorDataRange, g.edgeColorLimits );
if strcmpi(g.modulateEdgeColor, 'off'), nodeColorScaled  = 0; end
indexcolor = floor(max(nodeColorScaled-0.001,0)*length(g.edgeColormap))+1;
nodeColorScaled   = g.edgeColormap( indexcolor,: );
if isnan(edgeColor), nodeColorScaled = [1 1 1]; end

if tmpthick > 0
    if g.causality
        cylwidth = [0.1 g.resmult*tmpthick/300*100];    % tapered cylinder pointing from pos2->pos1
    else
        cylwidth=g.resmult*tmpthick/300*100;
    end
    %[xc yc zc] = cylinder( cylwidth, 10 );   % create unit-length cylinder
    [xc yc zc] = cylinder2P(cylwidth, numPointsCircum, 2, [pos1(1) pos1(2) pos1(3)], [pos2(1) pos2(2) pos2(3)]);
    colorarray = repmat(reshape(nodeColorScaled, 1,1,3), [size(zc,1) size(zc,2) 1]);
    %     handles = surf(g.vars.hbrainax,xc, yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', ...
    %         'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', g.facelighting,themeopts{:});
    %     %[xc yc zc] = adjustcylinder2( handles, [pos1(1) pos1(2) pos1(3)], [pos2(1) pos2(2) pos2(3)] );  % stretch and rotate cylinder to match start-end pnts
    %
    
    if nargout > 5
        % compute cylinder normals (have to bias normal closer to sphere
        % to get a specular point
        % TM: this is the part that creates the "spotlight" effect on each
        % arc
        cx = mean(xc,2); cx = [(3*cx(1)+cx(2))/4; (cx(1)+3*cx(2))/4];
        cy = mean(yc,2); cy = [(3*cy(1)+cy(2))/4; (cy(1)+3*cy(2))/4];
        cz = mean(zc,2); cz = [(3*cz(1)+cz(2))/4; (cz(1)+3*cz(2))/4];
        tmpx = xc - repmat(cx, [1 11]);
        tmpy = yc - repmat(cy, [1 11]);
        tmpz = zc - repmat(cz, [1 11]);
        l=sqrt(tmpx.^2+tmpy.^2+tmpz.^2);
        normals = reshape([tmpx./l tmpy./l tmpz./l],[2 11 3]);
    end
    %     set( handles, 'vertexnormals', normals);
    
%         if edgeColor>0,k=2; else k=1; end
%         cx = mean(xc,2); cx = cx(k)-3;
%         cy = mean(yc,2); cy = cy(k); %[(cy(1)+cy(2))/4; (cy(1)+cy(2))/4];
%         cz = mean(zc,2); cz = cz(k); %[(cz(1)+cz(2))/4; (cz(1)+cz(2))/4];
%         tmpx = xc - repmat(cx, [2 11]);
%         tmpy = yc - repmat(cy, [2 11]);
%         tmpz = zc - repmat(cz, [2 11]);
%         l=sqrt(tmpx.^2+tmpy.^2+tmpz.^2);
%         normals = reshape([tmpx./l tmpy./l tmpz./l],[2 11 3]);
        
end





% ***************************************************************************************
%                              Caption and tests
% ***************************************************************************************

% function to draw Node spheres
% -------------------------------------
function renderNodeSizeLegend(posx, posy, g)

NBCIRCLE = 3;
coordy = posy;
nodesizeLimits = [ ceil( g.nodeSizeDataRange(1) ) 0 floor( g.nodeSizeDataRange(2) ) ];
xlim = get(gca, 'xlim');
ylim = get(gca, 'ylim');

for i=1:NBCIRCLE
    [nodeSizeScaled] = drawcircle( [posx coordy], nodesizeLimits(i), 0, g);
    if i == 1, nodeSizeScaledori = nodeSizeScaled; end
    
    if i == NBCIRCLE
        hlatency = text( 1.4*(xlim(2) - xlim(1))+xlim(1), coordy , sprintf('%2.1g dB', nodesizeLimits(i)));
    else hlatency = text( 1.4*(xlim(2) - xlim(1))+xlim(1), coordy , sprintf('%2.1g', nodesizeLimits(i)));
    end
    set(hlatency, 'fontsize', 10*g.resmult, 'horizontalalignment', 'left', 'fontweight', 'bold');
    coordy = coordy + nodeSizeScaled + 0.2*(ylim(2)-ylim(1));
    
    %command2 = sprintf('print -depsc -loose scale%d.eps', i);
    %eval(command2);
    %cla;
end
set(gca, 'xlim', xlim, 'ylim', ylim-nodeSizeScaledori, 'clipping', 'off', 'fontsize', 10*g.resmult);


% function to draw edge cylinders
% ---------------------------------------
function renderEdgeSizeLegend(posx, posy, thickness,g,axhandle)

if nargin<5
    axhandle = gca;
end

compter = -5;
for i=linspace( posy(1), posy(2), 11)
    % TODO: replace this with rendering of 3D cylinders
    
    superline( [ posx(1) posx(2) ], [ i i ], 'b', thickness*g.resmult, mod(compter/10, 1));
    compter = compter + 1;
end

axis(axhandle,'off');
set(axhandle, 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'on');


% colorbar special
% ----------------
function cbar( X, Y, colors, orientation, style, g, axhandle)
% colors = colors to plot
% orientation = 'vert' or 'horiz'
% style = shape of the colorbar, 'circle' = circle, bar otherwise

if nargin<7
    axhandle = gca;
end

NSEGMENTS = size(colors,1)-1;
compter = 0;
switch lower(orientation)
    case 'horiz'
        inc = (X(2)-X(1))/NSEGMENTS;
        for i=linspace(X(1),X(2)-inc,NSEGMENTS);
            compter = compter + 1;
            hold on;
            h = fill( [i i i+inc i+inc], [Y(1) Y(2) Y(2) Y(1)], colors(size(colors,1)+1-compter, :), 'parent',axhandle);
            set(h, 'edgecolor', 'none');
        end
    case 'vert'
        inc = (Y(2)-Y(1))/NSEGMENTS;
        for i=linspace(Y(1),Y(2)-(Y(2)-Y(1))/NSEGMENTS,NSEGMENTS);
            hold on;
            switch style
                case 'circle',
                    mid     = (Y(2)-Y(1))/2;
                    angle   = acos( compter/NSEGMENTS*2-1);
                    angle1  = acos( (compter+1)/NSEGMENTS*2-1);
                    coordx1 = mid - sin( angle )*mid;
                    coordx2 = mid + sin( angle )*mid;
                    coordx3 = mid + sin( angle1 )*mid;
                    coordx4 = mid - sin( angle1 )*mid;
                    coordx = real([coordx1 coordx2 coordx3 coordx4]);
                otherwise,	coordx = [X(1) X(2) X(2) X(1)];
            end
            compter = compter + 1;
            h = fill( coordx, [i i i+inc i+inc], colors(compter, :), 'parent',axhandle);
            set(h, 'edgecolor', 'none');
        end
    otherwise
        disp('Orientation has to be ''vert'' or ''horiz''');
end
set(axhandle, 'fontsize', 10*g.resmult);
if strcmp(style, 'circle'), axis(axhandle,'square'); end


% draw vertical lines
% -------------------
function drawvert(tmpev, tmpcond, coords,h)

if nargin<4
    h = gca;
end

if isstruct(tmpev) || isstruct(tmpev{1})
    
    % cooper envert
    %--------------
    if length(tmpev) > 1,
        verts = tmpev{ tmpcond };
    else   verts = tmpev{1};
    end
    
    for v=verts,
        if isstruct(v), ev = v;
        else,           ev.time = v;  ev.color = 'k'; ev.style = '-';
        end
        
        phandle = plot(h,[ev.time ev.time], coords, ev.style, 'linewidth', 1);
        set(phandle,'color',ev.color);
    end
else
    % standard envvert
    % ----------------
    for index = 1:length(tmpev)
        if ~iscell(tmpev{index}),
            plot(h,[tmpev{index} tmpev{index}], coords, 'k');
        else
            phandle = plot(h,[tmpev{index}{1} tmpev{index}{1}], coords, 'k');
            if length(tmpev{index}) > 2
                set(phandle,tmpev{index}{2:end});
            end
        end
    end
end

% % check the flux
% % --------------
% for indeximage = 1:nwin-7
%     index1 = 1;
%     index2 = 2;
%     % determine color = coherence phase
%     tmpcrossf = ALLEDGECOLOR     { index1, index2, 1 };
%     tmpvalue  = mean(tmpcrossf( 1:2, indeximage));
%     nodeColorScaled  = colormaphsv( ceil((tmpvalue+180)/360*63) + 1, : );    % index for color
%
%     % absolute value to 90 degree determine speed
%     speed = 1 - abs(90 - abs(tmpvalue))/90; % speed from 1 to 0
%     currentphase(index1, index2) = currentphase(index1, index2) + sign(tmpvalue)*speed/3; % 1 cycle in 5 images at max speed
%
%     superline( [ 2 1] , [ 1+indeximage 0.8+indeximage], 5, nodeColorScaled, mod(currentphase(index1, index2),1));
% end
% return;

% scale values for node, connector and color
% ------------------------------------------
function value = scalevalue( value, datarange, limits)

value = (value-datarange(1))/(datarange(2)-datarange(1));
value = max(min(value,1),0);
value = value * (limits(2)-limits(1)) + limits(1);










% circle() - draw a circle in the current Matlab axes
%
% Usage:
%   >> [linehandles fillhandle] = circle(X,Y,radius,colorfill,coloredge,oriangle,...
%                                         endangle,dashed,thickness,segments);
% Inputs:
%   X, Y       - circle center coordinates
%   radius     - circle radius. Can be a vector of 2 values, one
%                for the x dimension and one the y dimension.
%   colorfill  - circle fill color (default:0=none)
%   coloredge  - circle edge color (default:black; 0:no edge)
%   oriangle   - starting angle (default:0 degrees)
%   endangle   - ending angle (default:360 degrees)
%   dashed     - 0=no, 1=yes (default: no)
%   thickness  - thickness of boarder (default: line default)
%   segments   - number of line segments (default:50)
%
% Outputs:
%   linehandles - handle to the lines of the circle object
%   fillhandle  - handle to the interior of the circle object

% Author: Arnaud Delorme, CNL / Salk Institute, 2001

% This program is free software; you can redistribute it and/or modify it.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function [h, h2] = circle( X, Y, radius, colorfill, coloredge, oriangle, endangle, dashed, thickness, segments, h);

if nargin < 3
    error('Not enough arguments. Type help circle');
    return;
end
if nargin < 4
    colorfill = 0;
end
if nargin < 5
    coloredge = 'k';
end
if nargin < 6
    oriangle = 0;
end
if nargin < 7
    endangle = 360;
end
if nargin < 8
    dashed = 0;
end
if nargin < 9
    thickness = 0;
end
if nargin < 10
    segments = 50;
end
if nargin < 11
    h = gca;
end
if any(radius <= 0)
    return;
end

A = linspace(oriangle/180*pi, endangle/180*pi, segments);

% draw surface
% ------------
if any(colorfill)
    A = linspace(oriangle/180*pi, endangle/180*pi, segments);
    h2 = patch( X + cos(A)*radius(1), Y + sin(A)*radius(end), zeros(1,segments), colorfill, 'parent', h);
    set(h2, 'FaceColor', colorfill);
    set(h2, 'EdgeColor', 'none');
end

% draw lines
% ----------
if dashed
    compt=0;
    for i=1:3:segments-2
        compt = compt+1;
        h(compt) = line( X + cos(A(i:i+1))*radius(1), Y + sin(A(i:i+1))*radius(end));
    end
else
    h = line( X + cos(A)*radius(1), Y + sin(A)*radius(end), 'parent', h);
end
set(h, 'Color', coloredge);

if thickness ~= 0
    set(h, 'Linewidth', thickness);
end
return;
