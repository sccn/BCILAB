
function [bg bgh causality] = vis_biograph(varargin)
%
% Create an interactive 3D BrainMovie from a connectivity matrix. See [1]
% for more details on the Interactive BrainMovie3D.
%
% Inputs:
% 
%       ALLEEG:     Array of EEGLAB datasets
%       Conn:       SIFT Connectivity Structure
%
% Optional:
%
%     Stats:                        Name of variable in base containing statistics                                                                    
%                                   Input Data Type: string                                                                                           
% 
%     ConnectivityMethod:           Connectivity Measure to visualize                                                                                 
%                                   Possible values: {'Determined_by_data'}                                                                           
%                                   Default value  : 'Determined_by_data'                                                                             
%                                   Input Data Type: string                                                                                           
% 
%     MovieTimeRange:               Time Range for Movie [Min Max] (sec)                                                                              
%                                   Specified w.r.t to event time (e.g., [-1 2]). Leave blank for complete epoch.                                     
%                                   Input Data Type: real number (double)                                                                             
% 
%     FrequenciesToCollapse:        Frequencies over which to collapse (Hz)                                                                           
%                                   E.g., [1:50] for 1-50Hz. Leave blank for all frequencies                                                          
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%     FreqCollapseMethod:           Method for collapsing frequency dimension                                                                         
%                                   Possible values: {'mean','max','peak','integrate'}                                                                
%                                   Default value  : 'mean'                                                                                           
%                                   Input Data Type: string                                                                                           
% 
%     TimeResamplingFactor:         Time resampling factor                                                                                            
%                                   If 0, don't resample. If < 1, downsample timecourse by this factor. If > 1, upsample by this                      
%                                   factor. Uses resample() from Sigproc Toolbox                                                                      
%                                   Input Range  : [0  20]                                                                                            
%                                   Default value: 0                                                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%     SubtractConditions:           Subtract conditions                                                                                               
%                                   If true, then plot difference between conditions. If false, then render two brainmovies                           
%                                   side-by-side.                                                                                                     
%                                   Input Data Type: boolean                                                                                          
% 
%     NodeLabels:                   List of labels for each node. e.g., {'Node1','Node2',...}                                                         
%                                   Leave blank to omit labels.                                                                                       
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%     NodesToExclude:               Exclude these sources from Brainmovie                                                                             
%                                   Specify using the Name/ID of the source to exclude.                                                               
%                                   Input Data Type: boolean                                                                                          
% 
%     EdgeColorMapping:             Specify mapping for edge color                                                                                    
%                                   This determines how we index into the colormap. If 'None', edge color is not modulated. If                        
%                                   'Connectivity', use connectivity strength. If 'PeakFreq', use index of peak frequency                             
%                                   Possible values: {'None','Connectivity','PeakFreq','Directionality'}                                              
%                                   Default value  : 'Connectivity'                                                                                   
%                                   Input Data Type: string                                                                                           
% 
%     EdgeSizeMapping:              Specify mapping for edge size                                                                                     
%                                   If 'None', edges are not rendered. If 'Connectivity', use connectivity strength. If                               
%                                   'ConnMagnitude', use connectivity magnitude (absval). If 'PeakFreq', use index of peak frequency.                 
%                                   If 'Directionality', map directionality to the lower and upper extremes of the colormap (e.g.,                    
%                                   i->j: blue, j->i: red)                                                                                            
%                                   Possible values: {'None','ConnMagnitude','Connectivity'}                                                          
%                                   Default value  : 'ConnMagnitude'                                                                                  
%                                   Input Data Type: string                                                                                           
% 
%     NodeColorMapping:             Specify mapping for node color                                                                                    
%                                   This determines how we index into the colormap. Options are as follows. None: node color is not                   
%                                   modulated. Outflow: sum connectivity strengths over outgoing edges. Inflow: sum connectivity                      
%                                   strengths over incoming edges. CausalFlow: Outflow-Inflow. Asymmetry Ratio: node colors are defined               
%                                   by the equation C = 0.5*(1 + outflow-inflow/(outflow+inflow)). This is 0 for exclusive inflow, 1                  
%                                   for exclusive outflow, and 0.5 for balanced inflow/outflow                                                        
%                                   Possible values: {'None','Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'}  
%                                   Default value  : 'Outflow'                                                                                        
%                                   Input Data Type: string                                                                                           
% 
%     NodeSizeMapping:              Specify mapping for node size. Options are as follows:                                                            
%                                   None: node size is not modulated.                                                                                 
%                                   Outflow: sum connectivity strengths over outgoing edges.                                                          
%                                   Inflow: sum connectivity strengths over incoming edges.                                                           
%                                   CausalFlow: Outflow-Inflow.                                                                                       
%                                   Asymmetry Ratio: node size is defined by the equation C = 0.5*(1 +                                                
%                                   outflow-inflow/(outflow+inflow)). This is 0 for exclusive inflow, 1 for exclusive outflow, and 0.5                
%                                   for balanced inflow/outflow                                                                                       
%                                   Possible values: {'None','Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'}  
%                                   Default value  : 'Outflow'                                                                                        
%                                   Input Data Type: string                                                                                           
% 
%     Baseline:                     Time range of baseline [Min Max] (sec)                                                                            
%                                   Will subtract baseline from each point. Leave blank for no baseline.                                              
%                                   Input Data Type: real number (double)                                                                             
% 
%     NormalizeConn:                Normalize edge and node values to [0 1]                                                                           
%                                   Values mapped to edge/node width and color are devided by max to put in [0 1] range. Recommended!                 
%                                   Input Data Type: boolean                                                                                          
% 
%     UseStatistics:                Use Statistics                                                                                                    
%                                   Input Data Type: boolean                                                                                          
%     --------------                                                                                                                                  
% 
%         Thresholding:             Type of thresholding for stats                                                                                    
%                                   If 'both' then stats should have upper and lower thresholds                                                       
%                                   Possible values: {'single','both','lessthan'}                                                                     
%                                   Default value  : 'single'                                                                                         
%                                   Input Data Type: string                                                                                           
% 
%         AlphaSignificance:        Significance threshold. e.g., 0.05 for p<0.05                                                                     
%                                   Input Range  : [0  1]                                                                                             
%                                   Default value: 0.05                                                                                               
%                                   Input Data Type: real number (double)                                                                             
% 
%     PercentileThreshold:          Percentile threshold                                                                                              
%                                   Fraction of "strongest" connections to display. E.g: PercentileThreshold=0.05 will display only the               
%                                   top 5% of connections                                                                                             
%                                   Input Range  : [0  1]                                                                                             
%                                   Default value: n/a                                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%     AbsoluteThreshold:            Exact threshold                                                                                                   
%                                   If a single value, then render only connections with strength above this threshold. If [low hi]                   
%                                   then render only connections with strength between [lo hi]. Overrides PercentileThreshold                         
%                                   Input Data Type: real number (double)                                                                             
% 
%     FooterPanelDisplaySpec:       Configure footer panel displayed at the bottom of the figure                                                      
%                                   If 'off', don't render footer. If 'ICA_ERP_Envelope', then display the ERP envelope of                            
%                                   backprojected components. If 'Chan_ERP_Envelope' then display the ERP envelope of selected channels               
%                                   Possible values: {'off','ICA_ERPenvelope','Chan_ERPenvelope'}                                                     
%                                   Default value  : 'off'                                                                                            
%                                   Input Data Type: string                                                                                           
%     -----------------------                                                                                                                         
% 
%         icaenvelopevars:          Select components to use in the display                                                                           
%                                   Input Data Type: boolean                                                                                          
% 
%         backprojectedchans:       List of channels to use in the backprojection                                                                     
%                                   Input Data Type: boolean                                                                                          
% 
%         chanenvelopevars:         Select channels to use in the display                                                                             
%                                   Input Data Type: boolean                                                                                          
% 
%     BrainMovieOptions:            Additonal options for rendering the brainmovie                                                                    
%                                   Input Data Type: string                                                                                           
%     ------------------                                                                                                                              
% 
%         Visibility:               Figure visibility when rendering movie                                                                            
%                                   If 'on,' render frames on screen (slower). If 'off,' keep them hidden (faster).                                   
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'on'                                                                                             
%                                   Input Data Type: string                                                                                           
% 
%         LatenciesToRender:        Subset of latencies to render (sec)                                                                               
%                                   Must be in TimeRange. Can be a vector The time point closest to the latency given are plotted. If                 
%                                   empty, render all latencies in TimeRange.                                                                         
%                                   Input Data Type: real number (double)                                                                             
% 
%         FramesToRender:           Vector of frame indices to compute                                                                                
%                                   E.g. [1:2] only computes the first two frames. If empty, render all frames                                        
%                                   Input Data Type: real number (double)                                                                             
% 
%         FigureHandle:             Handle to a figure to render brainmovie in                                                                        
%                                   Input Data Type: real number (double)                                                                             
% 
%         RotationPath3D:           Specify the rotation path for the BrainMovie                                                                      
%                                   Possible values: {'none','automatic','manual'}                                                                    
%                                   Default value  : 'none'                                                                                           
%                                   Input Data Type: string                                                                                           
%         ---------------                                                                                                                             
% 
%             AngleFactor:          Angle multiplicative factor                                                                                       
%                                   Input Data Type: real number (double)                                                                             
% 
%             PhaseFactor:          Phase multiplicative factor                                                                                       
%                                   Input Data Type: real number (double)                                                                             
% 
%         InitialView:              3D static starting view for the movie                                                                             
%                                   See 'help view'. Default is [50 36].                                                                              
%                                   Input Data Type: real number (double)                                                                             
% 
%         ProjectGraphOnMRI:        Project the 3D graph onto the 2D MRI slices                                                                       
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'off'                                                                                            
%                                   Input Data Type: string                                                                                           
% 
%         RenderCorticalSurface:    Superimpose semi-transparent smoothed cortex on brainmovie                                                        
%                                   UseOpenGL must be set to "on"                                                                                     
%                                   Input Data Type: boolean                                                                                          
%         ----------------------                                                                                                                      
% 
%             VolumeMeshFile:       Filename or path to mesh volume file                                                                              
%                                   Input Data Type: string                                                                                           
% 
%             Transparency:         Transparency of the cortical surface                                                                              
%                                   Real number in [0 1], where 0=opaque, 1=transparent.                                                              
%                                   Input Range  : [0  1]                                                                                             
%                                   Default value: 0.7                                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%             VolumeMeshFile:       Filename or path to mesh volume file                                                                              
%                                   Input Data Type: string                                                                                           
% 
%             Transparency:         Transparency of the cortical surface                                                                              
%                                   Real number in [0 1], where 0=opaque, 1=transparent.                                                              
%                                   Input Range  : [0  1]                                                                                             
%                                   Default value: 0.7                                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%         UseOpenGL:                OpenGL usage                                                                                                      
%                                   OpenGL may cause rendering problems with MacOSX                                                                   
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'on'                                                                                             
%                                   Input Data Type: string                                                                                           
% 
%         EventFlashTimes:          Vector of time indices at which the background flashes                                                            
%                                   Specify the color of the flash with a cell array of [1,2] cell arrays. Ex. { { 200 'y' } { 1500 '5'               
%                                   }} will generate two flashes, yellow at 200 ms and red at 1500 ms                                                 
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%         Square:                   ?                                                                                                                 
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'on'                                                                                             
%                                   Input Data Type: string                                                                                           
% 
%         DisplayLegendPanel:       Display legends in BrainMovie                                                                                     
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'on'                                                                                             
%                                   Input Data Type: string                                                                                           
% 
%         ShowLatency:              Display latency of current frame                                                                                  
%                                   This will render in lower left corner.                                                                            
%                                   Input Data Type: boolean                                                                                          
% 
%         DisplayRTProbability:     Display reaction time probabilty (if RT available)                                                                
%                                   This will render a small bar the height of which will vary based on the probability of response.                  
%                                   Input Data Type: boolean                                                                                          
% 
%         BackgroundColor:          Background color                                                                                                  
%                                   Can use any allowable Matlab color specification (see 'help ColorSpec').                                          
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%         GraphColorAndScaling:     Graph and Color Scaling                                                                                           
%                                   Options for coloring and scaling components of the directed graph                                                 
%                                   Input Data Type: string                                                                                           
%         ---------------------                                                                                                                       
% 
%             NodeSizeLimits:       [Min Max] limits for node size (pixels)                                                                           
%                                   Input Data Type: real number (double)                                                                             
% 
%             NodeColorLimits:      [Min Max] limits for node color (colormap index)                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%             EdgeSizeLimits:       [Min Max] limits for edge size (pixels)                                                                           
%                                   Input Data Type: real number (double)                                                                             
% 
%             EdgeColorLimits:      [Min Max] limits for edge color (colormap index)                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%             NodeSizeDataRange:    [Min Max] range of node size data                                                                                 
%                                   Input Data Type: real number (double)                                                                             
% 
%             NodeColorDataRange:   [Min Max] range of node color data                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%             EdgeSizeDataRange:    [Min Max] range for edge size data                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%             EdgeColorDataRange:   [Min Max] limits for edge color data                                                                              
%                                   Input Data Type: real number (double)                                                                             
% 
%             CenterDataRange:      Make 0 in the center of the colormap/datarange                                                                    
%                                   Input Data Type: boolean                                                                                          
% 
%             EdgeColorMap:         Expression defining the colormap for edges                                                                        
%                                   E.g., jet(64). See 'help colormap'.                                                                               
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%             NodeColormap:         Expression defining the colormap for nodes                                                                        
%                                   E.g., jet(64). See 'help colormap'.                                                                               
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%             DiskScalingFactor:    Numeric value that scales the size of disks                                                                       
%                                   Input Range  : [0  Inf]                                                                                           
%                                   Default value: 0.3                                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%             MagnificationFactor:  Magnification factor for graphics                                                                                 
%                                   Input Range  : [0  Inf]                                                                                           
%                                   Default value: 1                                                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%         OutputFormat:             Options for saving the movie/figures                                                                              
%                                   Input Data Type: string                                                                                           
%         -------------                                                                                                                               
% 
%             ImageOutputDirectory: Output directory to save images                                                                                   
%                                   If 'prompt', then you will be prompted to select the folder from a dialog. If blank, don't save                   
%                                   images                                                                                                            
%                                   Input Data Type: string                                                                                           
% 
%             ImageOutputFormat:    Format for saving images                                                                                          
%                                   Possible values: {'jpg','eps','ppm','tif','ai','bmp','emf','pbm','pcx','pdf','pgm','png','fig'}                   
%                                   Default value  : 'jpg'                                                                                            
%                                   Input Data Type: string                                                                                           
% 
%             MovieOutputFilename:  Movie filename                                                                                                    
%                                   E.g, 'movie.avi'. If 'prompt', then you will be prompted to select the file from a dialog. If                     
%                                   blank, don't save movie.                                                                                          
%                                   Input Data Type: string                                                                                           
% 
%             MovieOpts:            Cell array of movie options for avifile function                                                                  
%                                   See "help avifile".                                                                                               
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%             ImageSize:            Image size (pixels)                                                                                               
%                                   Input should be [widthcond height]. If more than one condition is being plotted horizontally, then                
%                                   widthcond is the width of each condition subplot                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%         mri:                      Dipplot MRI structure                                                                                             
%                                   Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a                
%                                   path to a Matlab file containing MRI structure. Default uses MNI brain.                                           
%                                   Input Data Type: string                                                                                           
% 
%         DipoleCoordinateFormat:   Coordinate format for dipplot                                                                                     
%                                   Possible values: {'spherical','mni'}                                                                              
%                                   Default value  : 'spherical'                                                                                      
%                                   Input Data Type: string                                                                                           
% 
%         DipplotOptions:           Additional dipplot options                                                                                        
%                                   Cell array of <'name',value> pairs of additional options for dipplot (see 'doc dipplot')                          
%                                   Input Data Type: any evaluable Matlab expression.  
%        
% Outputs:
%
%       g:          Configuration structure. Any brainmovie can be replicated 
%                   via the command vis_causalBrainMovie3D(ALLEEG,Conn,g);
%       handles:    Handles to figure and other objects
%       BMout:      Reserved for future use
%
%
% See Also: pop_vis_causalBrainMovie3D(), brainmovie3d_causal(), 
%
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen and Arnaud Delorme, 2010, SCCN/INC, UCSD. 
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



[bg bgh causality] = deal([]);

subconddef = false;
haspower   = false;

% extract some stuff from inputs for arg defaults
Conn = arg_extract(varargin,'Conn',2);
if ~isempty(Conn)
    Conn = Conn(1);
    connnames   = hlp_getConnMethodNames(Conn);
    conndef     = connnames{1};
    freqrange   = [Conn.freqs(1) Conn.freqs(end)];
    if length(Conn.freqs)>1
        freqdef     = ['[' num2str(freqrange(1)) ':' num2str(freqrange(end)) ']'];
    else
        freqdef   = Conn.freqs;
    end
        
    timerange   = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
    if length(Conn.erWinCenterTimes)>1
        timedef     = [timerange(1) timerange(end)];
    else
        timedef = Conn.erWinCenterTimes;
    end
    
    if length(Conn)>1
        subconddef = true;
    else
        subconddef = false;
    end
    
    haspower = isfield(Conn,'S');
    
    clear Conn;
end

ALLEEG = arg_extract(varargin,{'EEG','ALLEEG'},1);
[MyComponentNames MyChannelNames defcoordformat] = deal({});
if ~isempty(ALLEEG)
    ALLEEG = ALLEEG(1);
    if isfield(ALLEEG.CAT,'curComponentNames') && ~isempty(ALLEEG.CAT.curComponentNames)
        MyComponentNames = ALLEEG.CAT.curComponentNames;
    else
        MyComponentNames = ALLEEG.CAT.curComps;
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    end
    
    if isfield(ALLEEG,'chanlocs') && ~isempty(ALLEEG.chanlocs)
        MyChannelNames = {ALLEEG.chanlocs.labels};
    else
        MyChannelNames = strtrim(cellstr(num2str((1:ALLEEG.nbchan)'))');
    end
    
    if isfield(ALLEEG,'dipfit')
        defcoordformat = lower(ALLEEG.dipfit.coordformat);
    end
    
    clear ALLEEG;
end

usestatsdef = [];  % false


% ensure we have row vectors
MyComponentNames = MyComponentNames(:)';
MyChannelNames = MyChannelNames(:)';

NodeColorMappingOpts = [{'None','Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'} fastif(haspower,'Power',[])];
NodeSizeMappingOpts  = [{'None','Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'} fastif(haspower,'Power',[])];
EdgeColorMappingOpts = {'None','Connectivity','PeakFreq','Directionality'};
EdgeSizeMappingOpts  = {'None','ConnMagnitude','Connectivity'};

g = arg_define([0 2],varargin, ...
    arg_norep({'ALLEEG','EEG'},mandatory),...
    arg_norep({'Conn'},mandatory),...
    arg_nogui({'stats','Stats'},'',[],'Name of variable in base containing statistics.','type','char','shape','row'), ...
    arg({'connmethod','ConnectivityMethod'},conndef,connnames,'Connectivity Measure to visualize','shape','row','cat','DataProcessing'), ...
    arg({'timeRangeToCollapse','TimeRangeToCollapse'},timedef,[],['Time Range(s) over which to collapse (sec).' ...
                                                                  'Each row is a [Min Max] time range in sec.' ...
                                                                  'A separate graph will be made for each time range.' ...
                                                                  'If timeRangeToCollapse is a vector, then a graph is made at each timepoint in the vector.' ...
                                                                  'If this is an event-locked dataset, then times should be specified w.r.t to event time (e.g., [-1 2]).' ...
                                                                  'Leave blank to collapse over complete epoch.'], ...
                                                                  'shape','matrix','type','denserealdouble','cat','DataProcessing'), ...
    arg({'freqsToCollapse','FrequenciesToCollapse'},freqdef,[],'Frequencies over which to collapse (Hz). E.g., [1:50] for 1-50Hz. Leave blank for all frequencies','type','expression','shape','row','cat','DataProcessing'), ...
    arg({'collapsefun','FreqCollapseMethod'},'mean',{'mean','max','peak','integrate'},'Method for collapsing frequency dimension.','cat','DataProcessing'),...
    arg({'resample','TimeResamplingFactor'},0,[0 20],'Time resampling factor. If 0, don''t resample. If < 1, downsample timecourse by this factor. If > 1, upsample by this factor. Uses resample() from Sigproc Toolbox','cat','DataProcessing'), ...
    arg({'subtractconds','SubtractConditions'},subconddef,[],'Subtract conditions. If true, then plot difference between conditions. If false, then render two brainmovies side-by-side.','cat','DataProcessing'), ...
    arg_subtoggle({'showNodeLabels','ShowNodeLabels'},{},...
        { ...
            arg({'nodelabels','NodeLabels'},MyComponentNames,[],'List of labels for each node. e.g., {''Node1'',''Node2'',...}. Leave blank to omit labels.','cat','DisplayProperties'), ...
        },'Show node labels','cat','DisplayProperties'), ...
    arg({'nodesToExclude','NodesToExclude'},false,MyComponentNames,'Exclude these sources from Brainmovie. Specify using the Name/ID of the source to exclude.','type','logical','cat','DisplayProperties'),...
    arg({'edgeColorMapping','EdgeColorMapping'},'Connectivity',EdgeColorMappingOpts,{'Specify mapping for edge color.', 'This determines how we index into the colormap.', ...
                    sprintf(['\n' ...
                             '-None: edges are not colored.\n' ...
                             '-Connectivity: use connectivity strength.\n' ...
                             '-PeakFreq: use index of peak frequency.\n'])}, ...
                             'cat','DisplayProperties'),...
    arg({'edgeSizeMapping','EdgeSizeMapping'},'ConnMagnitude',EdgeSizeMappingOpts,{'Specify mapping for edge size.', ...
                    sprintf(['\n' ...
                             '-None: edges are not rendered.\n' ...
                             '-Connectivity: use connectivity strength.\n' ...
                             '-ConnMagnitude: use connectivity magnitude (absval).\n' ...
                             '-PeakFreq: use index of peak frequency.\n' ...
                             '-Directionality: map directionality to the lower and upper extremes of the colormap (e.g., i->j: blue, j->i: red)'])}, ...
                             'cat','DisplayProperties'),...
    arg({'nodeColorMapping','NodeColorMapping'},'Outflow',NodeColorMappingOpts,{'Specify mapping for node color.', 'This determines how we index into the colormap.', ...
            sprintf(['\n-None: node color is not modulated.\n' ...
                     '-Outflow: sum connectivity strengths over outgoing edges.\n' ...
                     '-Inflow: sum connectivity strengths over incoming edges.\n' ...
                     '-CausalFlow: Outflow-Inflow.\n' ...
                     '-Asymmetry Ratio: node colors are defined by the equation C = 0.5*(1 + outflow-inflow/(outflow+inflow)). This is 0 for exclusive inflow, 1 for exclusive outflow, and 0.5 for balanced inflow/outflow'])}, ...
                     'cat','DisplayProperties'), ...
    arg({'nodeSizeMapping','NodeSizeMapping'},'Outflow',NodeSizeMappingOpts,{'Specify mapping for node size.', ...
            sprintf(['\n-None: node size is not modulated.\n' ...
                     '-Outflow: sum connectivity strengths over outgoing edges.\n' ...
                     '-Inflow: sum connectivity strengths over incoming edges.\n' ...
                     '-CausalFlow: Outflow-Inflow.\n' ...
                     '-Asymmetry Ratio: node colors are defined by the equation C = 0.5*(1 + outflow-inflow/(outflow+inflow)). This is 0 for exclusive inflow, 1 for exclusive outflow, and 0.5 for balanced inflow/outflow'])}, ...
                     'cat','DisplayProperties'), ...
    arg({'baseline','Baseline'},[],[],'Time range of baseline [Min Max] (sec). Will subtract baseline from each point. Leave blank for no baseline.','cat','DataProcessing'),...
    arg_nogui({'normalize','NormalizeConn'},true,[],'Normalize edge and node values to [0 1]. Values mapped to edge/node width and color are devided by max to put in [0 1] range. Recommended!','cat','DataProcessing'), ...
    arg_subtoggle({'useStats','UseStatistics'},usestatsdef, ...
            { ... 
            arg({'threshside','Thresholding'},'single',{'single','both','lessthan'},'Type of thresholding for stats. If ''both'' then stats should have upper and lower thresholds'), ...
            arg({'alpha','AlphaSignificance'},0.05,[0 1],'Significance threshold. e.g., 0.05 for p<0.05') ...
            }, 'Use Statistics','cat','Thresholding'), ...
    arg({'prcthresh','PercentileThreshold'},[],[0 1],'Percentile threshold. Fraction of "strongest" connections to display. E.g: PercentileThreshold=0.05 will display only the top 5% of connections','type','denserealdouble','cat','Thresholding'), ...
    arg({'absthresh','AbsoluteThreshold'},[],[],'Exact threshold. If a single value, then render only connections with strength above this threshold. If [low hi] then render only connections with strength between [lo hi]. Overrides PercentileThreshold','type','denserealdouble','shape','row','cat','Thresholding'), ...
    arg_sub({'BMopts','BrainMovieOptions'},[], ...
        { ...
            arg({'size','ImageSize'},[600 600],[],'Image size (pixels). Input should be [widthcond height]. If more than one condition is being plotted horizontally, then widthcond is the width of each condition subplot','cat','DisplayProperties') ...
            arg({'visible','Visibility'},'on',{'on','off'},'Figure visibility when rendering movie. If ''on,'' render frames on screen (slower). If ''off,'' keep them hidden (faster).','cat','DisplayProperties'), ...
            arg_nogui({'latency','LatenciesToRender'},[],[],'Subset of latencies to render (sec). Must be in TimeRange. Can be a vector The time point closest to the latency given are plotted. If empty, render all latencies in TimeRange.','cat','DataProcessing'), ...
            arg_nogui({'figurehandle','FigureHandle'},[],[],'Handle to a figure to render brainmovie in'), ...
            arg({'makeCompass','MakeCompass'},false,[],{'Label cardinal directions.','e.g. Anterior, Posterior, Right, Left','cat','DisplayProperties'}), ... 
            arg({'caption','DisplayLegendPanel'},true,[],'Display legends in BrainMovie','cat','DisplayProperties'), ...
            arg({'displayLegendLimitText','DisplayLegendLimitText'},true,[],'Display data limits in the legend.'), ...
            arg({'showLatency','ShowLatency'},true,[],'Display latency of current frame. This will render in lower left corner.','cat','DisplayProperties'), ...
            arg({'dispRT','DisplayRTProbability'},false,[],'Display reaction time probabilty (if RT available). This will render a small bar the height of which will vary based on the probability of response.','cat','DisplayProperties'), ...
            arg({'backcolor','BackgroundColor'},[0 0 0],[],'Background color.  Can use any allowable Matlab color specification (see ''help ColorSpec'').','shape','row','type','expression','cat','DisplayProperties'), ...
            arg_sub({'graphColorAndScaling','GraphColorAndScaling'},{}, ...
                { ...
                arg({'nodeSizeLimits','NodeSizeLimits'},[0.1 1],[],'[Min Max] limits for node size (pixels).','shape','row','cat','DisplayProperties'), ...
                arg({'nodeColorLimits','NodeColorLimits'},[0 1],[],'[Min Max] limits for node color (colormap index).','shape','row','cat','DisplayProperties'), ...
                arg({'edgeSizeLimits','EdgeSizeLimits'},[0.1 1],[],'[Min Max] limits for edge size (pixels)','shape','row','cat','DisplayProperties'), ...
                arg({'edgeColorLimits','EdgeColorLimits'},[0 1],[],'[Min Max] limits for edge color (colormap index).','shape','row','cat','DisplayProperties'), ...
                arg({'nodeSizeDataRange','NodeSizeDataRange'},[],[],'[Min Max] range of node size data.','shape','row','cat','DisplayProperties'), ...
                arg({'nodeColorDataRange','NodeColorDataRange'},[],[],'[Min Max] range of node color data.','shape','row','cat','DisplayProperties'), ...
                arg({'edgeSizeDataRange','EdgeSizeDataRange'},[],[],'[Min Max] range for edge size data.','shape','row','cat','DisplayProperties'), ...
                arg({'edgeColorDataRange','EdgeColorDataRange'},[],[],'[Min Max] limits for edge color data.','shape','row','cat','DisplayProperties'), ...
                arg({'centerDataRange','CenterDataRange'},false,[],'Make 0 in the center of the colormap/datarange','cat','DisplayProperties'), ...
                arg({'edgeColormap','EdgeColorMap'},'jet(64)',[],'Expression defining the colormap for edges. E.g., jet(64). See ''help colormap''.','type','expression','shape','row','cat','DisplayProperties'), ...
                arg({'nodeColormap','NodeColorMap'},'jet(64)',[],'Expression defining the colormap for nodes. E.g., jet(64). See ''help colormap''.','type','expression','shape','row','cat','DisplayProperties'), ...
                arg({'diskscale','DiskScalingFactor'},0.3,[0 Inf],'Numeric value that scales the size of disks.','cat','DisplayProperties'), ...
                arg({'magnify','MagnificationFactor'},1,[0 Inf],'Magnification factor for graphics','cat','DisplayProperties') ...
                },'Graph and Color Scaling. Options for coloring and scaling components of the directed graph'), ...
            arg_sub({'outputFormat','OutputFormat'},{}, ...   
                { ...
                arg({'framefolder','ImageOutputDirectory'},[],[],'Output directory to save images. If ''prompt'', then you will be prompted to select the folder from a dialog. If blank, don''t save images','shape','row','type','char'), ...
                arg({'framesout','ImageOutputFormat'},'jpg',{'jpg','eps','ppm','tif','ai','bmp','emf','pbm','pcx','pdf','pgm','png','fig'},'Format for saving images'), ...
                arg({'moviename','MovieOutputFilename'},[],[],'Movie filename. E.g, ''movie.avi''. If ''prompt'', then you will be prompted to select the file from a dialog. If blank, don''t save movie.','shape','row','type','char'),...
                arg({'movieopts','MovieOpts'},{'videoname',''},[],'Cell array of movie options for avifile function. See "help avifile".','type','expression','shape','row'), ...
                arg({'size','ImageSize'},[],[],'Image size (pixels). Input should be [widthcond height]. If more than one condition is being plotted horizontally, then widthcond is the width of each condition subplot','cat','DisplayProperties') ...
                },'Options for saving the movie/figures','cat','MovieOutput'), ...
            arg({'mri'},'standard_BEM_mri.mat',[],'Dipplot MRI structure. Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a path to a Matlab file containing MRI structure. Default uses MNI brain.','type','char','shape','row'), ...
            arg({'plotimgs','ShowMRISlices'},false,[],'Display axial/saggital/coronal MRI slices on the X,Y,Z axes'), ...
            arg({'coordformat','DipoleCoordinateFormat'},defcoordformat,{'spherical','mni'},'Coordinate format for dipplot','type','char','shape','row'), ...
            arg({'dipplotopt','DipplotOptions'},'{}','','Additional dipplot options. Cell array of <''name'',value> pairs of additional options for dipplot (see ''doc dipplot'')','type','expression','shape','row'), ...
            arg({'bmopts_suppl','BrainMovieSuppOptions'},'{}','','Additional options to brainmovie3d_causal.m. Cell array of <''name'',value> pairs of additional options (see ''brainmovie3d_causal.m'')','type','expression','shape','row'), ...
            arg_nogui('renderBrainMovie',true,[],'Special option. Determines whether to actually render the brainmovie (or just return values)'), ...
            arg_nogui('speedy',false,[],'Special option. Determines whether we are in "speedy" render mode with limited legend, font, etc elements'), ...
            arg_nogui({'mode','RenderMode'},'init_and_render',{'init','render','init_and_render'},'BrainMovie render mode. If set to ''init'', movie will only be initialized and internal state returned in BMout.BMopts.vars. If ''render'' then movie will only be rendered using {''BMopts'',BMout.BMopts{:}}. If ''init_and_render'', then both initialization and rendering is carried out'), ...
            arg_nogui({'vars','InternalStateVariables'},[],[],'Internal state variables set during initialization') ...
            }, ...
    'Additional options for rendering the brainmovie','cat','DisplayProperties') ...
    );

% copy the stats structure from base to current workspace
if ~isempty(g.stats)
    g.stats = evalin('base',g.stats);
end


% VALIDATE INPUTS
% ---------------------------------------------
if any(cellfun(@isempty,{g.ALLEEG.dipfit}))
    error('SIFT:vis_causalBrainMovie3D','In order to use BrainMovie3D, source locations must be stored in EEG.dipfit');
end

if any(ismember_bc({g.nodeSizeMapping, g.nodeColorMapping},'power')) && ~isfield(g.Conn,'S')
    error('To modulate node color/size by power you must have pre-computed the spectral density measure (Conn object must contain field ''S'')');
end

if ~ismember_bc(g.collapsefun,{'max','peak'}) && strcmpi(g.edgeColorMapping,'peakfreq')
    error('To use PeakFreq EdgeColorMapping, you must select ''max'' or ''peak'' as the FreqCollapseMethod');
end

if length(setdiff_bc(MyComponentNames,g.nodesToExclude)) < 2
    error('You must include at least two nodes in the BrainMovie');
end
    
% if no stats present, disable stats
if isempty(g.stats)
    g.useStats.arg_selection = false;  end
    
% check that frequencies are in valid range
if any(g.freqsToCollapse < freqrange(1)) || ...
   any(g.freqsToCollapse > freqrange(end))
    error('FreqsToCollapse contains elements outside valid range [%0.1f %0.1f]', ...
        freqrange(1),freqrange(end));
end
    
% reset default time range
if isempty(g.timeRange)
    g.timeRange   = [g.Conn.erWinCenterTimes(1) g.Conn.erWinCenterTimes(end)];
end

% check that times are in valid range
if g.timeRange(1)+1e-5 < timerange(1) || ...
   g.timeRange(2)-1e-5 > timerange(end)
    error('TimeRange contains elements outside valid range [%0.1f %0.1f]', ...
        timerange(1),timerange(end));
end

% check that latencies are in TimeRange
if any(g.BMopts.latency < g.timeRange(1)) || ...
   any(g.BMopts.latency > g.timeRange(end))
    error('Latency (%0.1f s) is outside the specified TimeRange [%0.1f %0.1f]. Please provide a valid latency', ...
            g.BMopts.latency,g.timeRange(1),g.timeRange(1));
end

% ---------------------------------------------

    

% HANDLE DEFAULTS
% ---------------------------------------------

% do some cleanup
g.nodeSizeMapping   = lower(g.nodeSizeMapping);
g.nodeColorMapping  = lower(g.nodeColorMapping);
g.edgeSizeMapping   = lower(g.edgeSizeMapping);
g.edgeColorMapping  = lower(g.edgeColorMapping);


if isempty(MyComponentNames)
    if isfield(g.ALLEEG.CAT,'curComponentNames') && ~isempty(g.ALLEEG.CAT.curComponentNames)
        MyComponentNames = g.ALLEEG.CAT.curComponentNames;
    else
        MyComponentNames = g.ALLEEG.CAT.curComps;
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    end
end
    
if isempty(g.connmethod)
    % if no connectivity methods specified, select first one
    g.connmethod = hlp_getConnMethodNames(g.Conn(1));
    g.connmethod = g.connmethod{1};                         end
if isempty(g.freqsToCollapse)
    g.freqsToCollapse = g.Conn(1).freqs;                    end
if isempty(g.timeRange)
    g.timeRange = g.Conn(1).erWinCenterTimes([1 end]);      end


% make a copy for convenience
erWinCenterTimes = g.Conn(1).erWinCenterTimes;

% obtain factors for resampling
if g.resample
    [p q] = rat(g.resample);
    g.resample = [p q];
end

if g.resample
    erWinCenterTimes = resample(erWinCenterTimes,g.resample(1),g.resample(2));
end


% number of sources
N=g.ALLEEG(1).CAT.nbchan-length(g.nodesToExclude);

% indices of nodes we will exclude
nodeIndicesToExclude = find(ismember_bc(MyComponentNames,g.nodesToExclude));

% get the indices of frequencies to collapse
freqIndicesToCollapse = getindex(g.Conn(1).freqs,g.freqsToCollapse);

% get indices of desired time windows ...
timeIndices =  getindex(erWinCenterTimes,g.timeRange(1)):getindex(erWinCenterTimes,g.timeRange(2));

% ... and construct the time vector for the Brainmovie
BrainMovieTimeRangeInSec = erWinCenterTimes(timeIndices);

% ---------------------------------------------


% TRANSFORM DATA INTO BRAINMOVIE FORMAT
% ---------------------------------------------
for cnd = 1:length(g.Conn)
    % select the connectivity measure
    Conn = g.Conn(cnd).(g.connmethod);
    
    % select the appropriate stats
    if ~isempty(g.stats)
        Stats = g.stats(cnd).(g.connmethod);
    end

    % remove the baseline
    if ~isempty(g.baseline)
        Conn = hlp_rmbaseline(Conn,g.baseline,g.Conn(cnd).erWinCenterTimes);
    end

    % Apply statistics and thresholding
    if g.useStats.arg_selection
        if isempty(Stats)
            error('SIFT:StatsThresholding:NoStats','You must supply a ''stats'' structure. See help');
            return;
        end
        % apply statistical threshold

        switch g.useStats.threshside
            case 'both'
                Conn((Conn > squeeze(Stats.thresh(1,:,:,:,:,:))) & (Conn < squeeze(Stats.thresh(2,:,:,:,:,:)))) = 0;
            case 'single'
                Conn(abs(Conn) < abs(Stats.thresh)) = 0;
            case 'lessthan'
                Conn(Conn < squeeze(Stats.thresh(2,:,:,:,:,:))) = 0;
        end
    end

    
    % Format edge data (size/color)
    % -----------------------------
    [EdgeSize EdgeColor] = deal(zeros(N,N,length(timeIndices),class(Conn(1))));
    for ch1=1:N
        for ch2=1:N

            if ch1==ch2 ...
               || ~isempty(intersect_bc(nodeIndicesToExclude,[ch1 ch2])) ...
               || all(ismember_bc({g.edgeSizeMapping, g.edgeColorMapping},'none'))
           
                continue;
            end
            
            % extract data
            causality = squeeze(Conn(ch1,ch2,:,:));
            
            % make row vector if necessary
            if length(g.freqsToCollapse)==1 && size(causality,1)>1 && size(causality,2)==1
                causality = causality';
            end
            
            % collapse matrix across frequencies for selected time range
            if length(g.freqsToCollapse)>1
                % collapse across frequency
                [causality peakidx] = hlp_collapseFrequencies( causality, ...
                  g.collapsefun,freqIndicesToCollapse,timeIndices,g.freqsToCollapse(2)-g.freqsToCollapse(1));

            else
                % multiple freqs, but only one frequency selected
                causality = causality(freqIndicesToCollapse,timeIndices);
                peakidx = freqIndicesToCollapse;
            end

            if g.resample
                causality   = resample(causality.',g.resample(1),g.resample(2)).';
                peakidx     = resample(peakidx.',g.resample(1),g.resample(2)).';
            end

            
            switch g.edgeSizeMapping
                case 'connectivity'
                    EdgeSize(ch1,ch2,:) = causality;
                case 'peakfreq'
                    EdgeSize(ch1,ch2,:) = peakidx;
                case 'connmagnitude'
                    EdgeSize(ch1,ch2,:) = abs(causality);    
%                 case 'none'
%                     EdgeSize(ch1,ch2,:) =
%                     zeros(1,length(timeIndices),class(Conn(1)));
            end
            
            if strcmpi(g.edgeSizeMapping,g.edgeColorMapping)
                EdgeColor(ch1,ch2,:) = EdgeSize(ch1,ch2,:);
            else
                switch g.edgeColorMapping
                    case 'connectivity'
                        EdgeColor(ch1,ch2,:) = causality;
                    case 'peakfreq'
                        EdgeColor(ch1,ch2,:) = peakidx;
%                     case 'none'
%                         EdgeColor(ch1,ch2,:) = zeros(1,length(timeIndices),class(Conn(1)));
                    case 'directionality'
                        EdgeColor(ch1,ch2,:) = causality;
                    case 'connmagnitude'
                        EdgeColor(ch1,ch2,:) = abs(causality);
                end
            end
 
            
        end % for ch2
    end % for ch1

    
    if ~isempty(g.absthresh)
        if isscalar(g.absthresh)
            edgeSizeThresh = [-Inf g.absthresh];
        else
            edgeSizeThresh = g.absthresh;
        end
        
        EdgeSize(EdgeSize   > edgeSizeThresh(1)     ...
               & EdgeSize   < edgeSizeThresh(2))    = 0;   
    end
    
    if ~isempty(g.prcthresh) && g.prcthresh < 1
        % apply percentile thresholding of the edges
        prcthresh = [0 1-g.prcthresh]*100;
        % get thresholds
        edgeSizeThresh  = prctile(abs(EdgeSize(:)) ,prcthresh);
        EdgeSize(abs(EdgeSize)   > edgeSizeThresh(1)     ...
               & abs(EdgeSize)   < edgeSizeThresh(2))    = 0;      
    end
    
    
    if ~all(ismember_bc({g.nodeSizeMapping, g.nodeColorMapping},'none'))
        % we will need the connectivity data in matrix format so we can
        % obtained graph statistics for mapping Node Size/Color
        
        if strcmpi(g.edgeSizeMapping,'connectivity')
            causality = squeeze(EdgeSize(:,:,:));
        elseif strcmpi(g.edgeColorMapping,'connectivity')
            causality = squeeze(EdgeColor(:,:,:));
        else % Conn matrices have not been collapsed ...
             %  ... we have to collapse matrix across frequencies
             causality = zeros(N,N,length(timeIndices));
             for ch1=1:N
                for ch2=1:N
                    [causality(ch1,ch2,:)] = hlp_collapseFrequencies( ...
                            squeeze(Conn(ch1,ch2,:,:)), g.collapsefun, ...
                            freqIndicesToCollapse,timeIndices,g.freqsToCollapse(2)-g.freqsToCollapse(1));
                end
             end
        end
    end


    % Format node data (size/color)
    [NodeSize NodeColor] = deal(zeros(N,length(timeIndices),class(Conn(1))));
    
    if ~all(ismember_bc({g.nodeSizeMapping, g.nodeColorMapping},'none'))
        
        for ch1=1:N
            if ismember_bc(ch1,nodeIndicesToExclude)
                continue;
            end

            if strcmpi(g.nodeSizeMapping,'power')
                % map power --> nodeSize
                NodeSize(ch1,:) = hlp_collapseFrequencies( ...
                            squeeze(g.Conn.S(ch1,ch1,:,:)), g.collapsefun, ...
                            freqIndicesToCollapse,timeIndices,g.freqsToCollapse(2)-g.freqsToCollapse(1));
            else
                % compute desired graph measure for this node and map to size
                % we add eps=10^-16 just to ensure that values are never exactly
                % zero (or they will get set to NaN later)
                NodeSize(ch1,:) = hlp_computeGraphMeasure(causality,ch1,g.BMopts.selected,g.nodeSizeMapping)+eps;
            end
        
            
            if strcmp(g.nodeSizeMapping,g.nodeColorMapping)
                % no need to re-compute anything for node color
                    NodeColor(ch1,:) = NodeSize(ch1,:);
            else
                
                if strcmpi(g.nodeColorMapping,'power')
                    % map power --> nodeColor
                    NodeColor(ch1,:) = hlp_collapseFrequencies( ...
                                squeeze(g.Conn.S(ch1,ch1,:,:)), g.collapsefun, ...
                                freqIndicesToCollapse,timeIndices,g.freqsToCollapse(2)-g.freqsToCollapse(1));
                else
                    % compute desired graph measure for this node and map to color
                    % we add eps=10^-16 just to ensure that values are never exactly
                    % zero (or they will get set to NaN later)
                    NodeColor(ch1,:) = hlp_computeGraphMeasure(causality,ch1,g.BMopts.selected,g.nodeColorMapping)+eps;
                end
            end

            
            if g.resample
                NodeSize(ch1,:)  = resample(NodeSize(ch1,:).',g.resample(1),g.resample(2)).';
                NodeColor(ch1,:) = resample(NodeColor(ch1,:).',g.resample(1),g.resample(2)).';
            end
        end % loop over channels
        
    end
    
    
    % check for negative values
%     HaveNegativeValues = (any(EdgeSize(:)<0) ||  any(NodeColor(:)<0));
    
    % set to NaN any nonsignificant values
    NodeSize(NodeSize==0)   = nan;
    NodeColor(NodeColor==0) = nan;
    EdgeSize(EdgeSize==0)   = nan;
    EdgeColor(EdgeColor==0) = nan;
    

end  % loop over conditions


% subtract multiple conditions
if g.subtractconds && size(EdgeSizeCell,3)==2  
    for ch1=1:N
        for ch2=1:N
            if ch1==ch2, continue; end
                tmpedgesize{ch1,ch2,1}   = EdgeSizeCell{ch1,ch2,1}  - EdgeSizeCell{ch1,ch2,2};
                tmpedgecolor{ch1,ch2,1}  = EdgeColorCell{ch1,ch2,1} - EdgeColorCell{ch1,ch2,2};
        end
    end
    EdgeSizeCell    = tmpedgesize;
    EdgeColorCell   = tmpedgecolor;

    for ch1=1:N
        tmpnodesize{ch1,1}     = NodeSizeCell{ch1,1}  - NodeSizeCell{ch1,2};
        tmpnodecolor{ch1,1}    = NodeColorCell{ch1,1} - NodeColorCell{ch1,2};
    end
    
    NodeSizeCell    = tmpnodesize;
    NodeColorCell   = tmpnodecolor;
    
    g.BMopts.condtitle = sprintf('%s - %s', g.ALLEEG(1).condition,g.ALLEEG(2).condition); 
else
    g.BMopts.condtitle    = {g.ALLEEG.condition};
end



% SETUP THE BRAINMOVIE OPTIONS
% ---------------------------------------------
BMopts              = g.BMopts;
BMopts.causality    = true;
if g.showNodeLabels.arg_selection
    BMopts.nodelabels   = g.showNodeLabels.nodelabels;
else
    BMopts.nodelabels = {''};
end
        
BMopts.title = BMopts.condtitle;
% extract the coordinates of dipoles
% for cnd=1:length(g.Conn)
%     coords{cnd} = {g.ALLEEG(cnd).dipfit.model.posxyz};
%     coords{cnd} = coords{cnd}(g.ALLEEG(cnd).CAT.curComps);
% end
coords = {g.ALLEEG(cnd).dipfit.model.posxyz};
coords = coords(g.ALLEEG(cnd).CAT.curComps);
BMopts.coordinates = coords;

BMopts.windowLength = g.ALLEEG(cnd).CAT.MODEL.winlen;

% determine whether to render and/or modulate properties of edges/nodes 
BMopts.modulateEdgeColor     = fastif(strcmpi('none',g.edgeColorMapping),'off','on');    % do/don't modulate edgecolor
BMopts.modulateEdgeSize      = fastif(strcmpi('none',g.edgeSizeMapping),'off','on');     % do/don't render edges
BMopts.modulateNodeColor     = fastif(strcmpi('none',g.nodeColorMapping),'off','on');    % do/don't modulate nodecolor
BMopts.modulateNodeSize      = fastif(strcmpi('none',g.nodeSizeMapping),'off','on');     % do/don't modulate nodesize

SLASH = fastif(isunix,'/','\');

% Select the folder to save frames
if strcmpi(BMopts.outputFormat.framefolder,'prompt')
    framefolder = uigetdir(pwd,'Select directory to save all movie frames');
    if framefolder~=0
        BMopts.outputFormat.framefolder = framefolder;
    else
        BMopts.outputFormat.framefolder = '';
    end
elseif ~isempty(BMopts.outputFormat.framefolder) && ~isdir(BMopts.outputFormat.framefolder)
    resp = questdlg2(sprintf('You asked to save figures here:\n ''%s''\nThis folder does not exist. Should I create it?',BMopts.outputFormat.framefolder),'Save BrainMovie Figures','Yes','No','No');
    if strcmpi(resp,'Yes')
        mkdir(BMopts.outputFormat.framefolder);
    else
        BMopts.outputFormat.framefolder = '';
    end
end
    

% select the movie filename
if strcmpi(BMopts.outputFormat.moviename,'prompt')
    moviename = uiputfile('*.mov, *.avi','Create output movie file');
    if moviename~=0
        BMopts.outputFormat.moviename = moviename;
    else
        BMopts.outputFormat.moviename = '';
    end
end


BMopts.footerPanelTimes     = envtimes;
BMopts.footerPanelData      = envdata;
if ~strcmpi(g.footerPanelSpec.arg_selection,'off')
    BMopts.footerPanelPlotMode = g.footerPanelSpec.plottingmode;
else
    BMopts.footerPanelPlotMode = [];
end

% BMopts.footerPanelData      = footerPanelData;

% get latencies for movie
if isempty(g.BMopts.latency) && isempty(g.BMopts.frames)
    % default: use all latencies in time range
    BMopts.latency = BrainMovieTimeRangeInSec;
else
    % convert to ms
    BMopts.latency = g.BMopts.latency;
end


if strcmpi(g.edgeColorMapping,'Directionality')
    BMopts.EdgeColorMappedToDirectionality = true;
end


% check if we have negative values (and need +- colormaps)
for i=1:length(NodeColorCell)
    if ~isempty(NodeColorCell{i})
        if any(NodeColorCell{i} < 0)
            BMopts.nodeColorPolarity = 'posneg';
        end
    end
end
for i=1:length(EdgeColorCell(:))
    if ~isempty(EdgeColorCell{i})
        if any(EdgeColorCell{i} < 0)
            BMopts.edgeColorPolarity = 'posneg';
        end
    end
end
% 
% if any(any(any(cell2mat(NodeColorCell)<0)))
%     BMopts.nodeColorPolarity = 'posneg'; end
% if any(any(any(cell2mat(EdgeColorCell)<0)))
%     BMopts.edgeColorPolarity = 'posneg'; end

BMopts.NodeColorMapping = g.nodeColorMapping;
BMopts.EdgeColorMapping = g.edgeColorMapping;
BMopts.NodeSizeMapping  = g.nodeSizeMapping;
BMopts.EdgeSizeMapping  = g.edgeSizeMapping;
BMopts.ConnMethod       = g.connmethod;
BMopts.collapsedFreqs   = g.freqsToCollapse;



%% Draw the biograph
% -------------------------------------------------------------------------

% Plot graphs for each of the time windows

for t=1:NumTimeWindows
    
    % get the [nchs_to x nchs_from] connectivity matrix for time window t
    C = causality(:,:,t);
    
    if issym(C)
        csym = true;
        C = triu(C);
    else
        csym = false;
    end
    
    if all(C==0)
        fprintf('No non-zero connections at time [%s]. Graph will not be constructed\n', ...
                regexprep(num2str(TimeRangeToCollapse(t,:)),'\s*','  '));
        bg{t}   = [];
        bgh{t}  = [];
        continue;
    else
    
        % create the biograph
        bg{t} = biograph(C',NodeLabels,'ShowWeights','off',...
                                       'ShowTextInNodes','Label', ...
                                       'EdgeType','curved', ...
                                       'NodeAutoSize','off', ...
                                       'LayoutType',fastif(ischar(LayoutType),LayoutType,'radial'));
    end
    
    if PlotFigures
        % plot the graph
        if ismatrix(LayoutType) || strcmpi(LayoutType,'radial')
            bh = bgCircleGraph(bg{t},NodeSize); % view(bg{t});
        else
            bh = view(bg{t});
        end
        
        if csym % symmetric connectivity
            set(bh,'showArrows','off');
%             dolayout(bh,'pathsOnly',true);
        end
        
        % bh=view(biograph(C',ComponentNames,'ShowWeights','on','LayoutType','radial','ShowTextInNodes','Label','Scale',1.5)); %'hierarchical', 'radial'
        
        % if colorclusters
        %     nodecolors = distinguishable_colors(numclusters,[1 1 1; 0 0 0; 0 0 1]);
        %     % make each cluster a different color
        %     for i=1:numclusters
        %         set(bh.Nodes(myclustering(myclustering(:,2)==i,1)),'Color',nodecolors(i,:));
        %     end
        % end
        
        for i=1:length(bh.Nodes)
            set(bh.Nodes(i),'label',NodeLabels{i},'FontSize',12);
        end
        
        % set the node locations
        if ismatrix(LayoutType)
            for i=1:length(bh.Nodes)
                set(bh.Nodes(i),'Position',LayoutType(i,:));
            end
            
            % redraw the edges of the graph
            bh.dolayout('pathsonly',true);
            bh.hgUpdate
            bh.hgReDraw
        end
        
        cmaplen = 128;           %length(find(C_cluster>0))
        edgecolors = mycolormap; %jet(cmaplen);
        
        cmin = EdgeSizeMin;     %0;  min(C(C>0));
        cmax = EdgeSizeMax;     %    max(C(:));
        lmax = EdgeSizeMax;     %    max(C(:));
        lmin = EdgeSizeMin;     %0;  min(C(C>0));
        
        for c1=1:size(C,1)
            for c2=1:size(C,1)
                cval=C(c1,c2);
                
                if cval>0
                    % determine line color
                    index = min(cmaplen,fix((cval-cmin)/(cmax-cmin)*(cmaplen-1))+1);
                    LineColor = edgecolors(index,:);
                    
                    % determine linewidth
                    
                    LineWidth = min(lMAX,1+(cval-lmin)*(lMAX-lMIN)/(lmax-lmin+eps));
                    c2name = NodeLabels{c2};
                    c1name = NodeLabels{c1};
                    set(findobj(bh.Edges,'ID',sprintf('%s -> %s',c2name,c1name)),'LineColor',LineColor,'LineWidth',LineWidth);
                end
                
            end
        end
        
        % Add a few more things to the graph
        fighandle = get(bh.hgAxes,'parent');
        
        if MakeLegend
            
            % create legend
            
            timestr = regexprep(num2str(TimeRangeToCollapse(t,:)),'\s*','  ');
            freqstr = regexprep(num2str(FreqRangeToCollapse),'\s*','  ');
            
            edgeleg     = sprintf('EdgeSize:\t%s  ',connmethod);
            nodeleg     = fastif(isempty(UnivariateGraphMeasure),{},{sprintf('NodeSize:\t%s  ',UnivariateGraphMeasure)});
            timesleg    = fastif(isnan(TimeRangeToCollapse),{}, ...
                                 sprintf('Times:       [%s] sec ',fastif(isempty(timestr),Conn.erWinCenterTimes,num2str(timestr))));
            freqleg     = fastif(isnan(FreqRangeToCollapse),{}, ...   
                                sprintf('Freqs:       [%s] Hz ',  fastif(isempty(freqstr),num2str(Conn.freqs),freqstr)));
            legend = [...
                      edgeleg, ...
                      nodeleg, ...
                      timesleg, ...
                      freqleg ...
                    ];
            
            % TODO: need to set up callbacks to allow the legend to be moved
            % NOTE: if we want the legend to be preserved in 'print to figure',
            %       then we must use the 'text' version (no FaceAlpha available)
            %         ha = annotation(fighandle,'textbox',[0.6895, 0.8784,0.1,0.2],'FaceAlpha',0.5,'FitBoxToText','on','String',legend,'units','normalized','EdgeColor','k','BackGroundColor',[1 0.7 0.7],'FontSize',10);
            ht=text(0.6895, 0.8784,legend,'parent',bh.hgAxes,'units','normalized','EdgeColor','k','BackGroundColor',[1 0.7 0.7],'FontSize',10); % 'BackGroundColor',[1 0.7 0.7]
            %         pos = get(ht,'Extent');
            %         axespos = get(bh.hgAxes,'Position');
            % transform pos to data units
            %         pos = pos.*repmat(axespos(3:4),1,2);
            %         hp = patch([pos(1) pos(1) pos(1)+pos(3) pos(1)+pos(3)],[pos(2) pos(2)+pos(4) pos(2)+pos(4) pos(2)],[1 0.7 0.7],'FaceAlpha',0.3,'parent',bh.hgAxes);
        end
        
        set(fighandle,'Name',sprintf('Connectivity Graph %d',t));
        
        %         set(fighandle,'Name',sprintf('Times: [%0.3g %0.3g] sec, Freqs: [%0.4g %0.4g] Hz',...
        %             TimeRangeToCollapse(t,1),TimeRangeToCollapse(t,2),FreqRangeToCollapse(1),FreqRangeToCollapse(2)));
        if NumTimeWindows>1
            set(fighandle,'WindowStyle','docked');
        end
        
        bgh{t} = bh;
        
    end

end

% EOF






        
        

