
function [g handles BMout] = vis_causalBrainMovie3D(varargin)
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
% Thanks to: Arnaud Delorme for contributing the original brainmovie3D
%            Alejandro Ojeda for helping with LONI Atlas meshes
%            Zeynep Akalin Acar for CSF/Scalp/Skull meshes
%            Nima Bigdely-Shamlo and Christian Kothe for helpful visualization input

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



[g handles BMout] = deal([]);

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

NodeColorMappingOpts = [{'None','Outflow','Mag_Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'} fastif(haspower,'Power',[])];
NodeSizeMappingOpts = [{'None','Outflow','Mag_Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'} fastif(haspower,'Power',[])];
EdgeColorMappingOpts = {'None','Connectivity','PeakFreq','Directionality'};
EdgeSizeMappingOpts = {'None','ConnMagnitude','Connectivity'};

g = arg_define([0 2],varargin, ...
    arg_norep({'ALLEEG','EEG'},mandatory),...
    arg_norep({'Conn'},mandatory),...
    arg_nogui({'stats','Stats'},'',[],'Name of variable in base containing statistics.','type','char','shape','row'), ...
    arg({'connmethod','ConnectivityMethod'},conndef,connnames,'Connectivity Measure to visualize','shape','row','cat','DataProcessing'), ...
    arg({'timeRange','MovieTimeRange'},timedef,[],'Time Range for Movie [Min Max] (sec). Specified w.r.t to event time (e.g., [-1 2]). Leave blank for complete epoch.','shape','row','type','denserealdouble','cat','DataProcessing'), ...
    arg({'freqsToCollapse','FrequenciesToCollapse'},freqdef,[],'Frequencies over which to collapse (Hz). E.g., [1:50] for 1-50Hz. Leave blank for all frequencies','type','expression','shape','row','cat','DataProcessing'), ...
    arg({'collapsefun','FreqCollapseMethod'},'mean',{'mean','max','peak','integrate','absmax','minmax'},'Method for collapsing frequency dimension.','cat','DataProcessing'),...
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
    arg_subswitch({'footerPanelSpec','FooterPanelDisplaySpec'},'off', ...
        {'off' { ...
            arg_norep({'dummy1'},[],[],'dummy') ...
            }, ...
         'ICA_ERPenvelope' {...
            arg({'icaenvelopevars','ICs'},true,MyComponentNames,'Select components to use in the display','type','logical'), ...
            arg({'backprojectedchans','BackProjectToChans'},true,MyChannelNames,'List of channels to use in the backprojection'), ...
            arg({'plottingmode','PlotMode'},{'all', 'envelope'},{'all','envelope'},{'Plotting Mode.' '''All'' plots all individual traces.' '''Envelope'' plots the envelope.'},'type','logical'), ...
            arg({'envColor','EnvelopeColor'},[1 0 0],[],'Envelope color [R G B]'), ...
            } ...
          'Chan_ERPenvelope' {...
            arg({'chanenvelopevars','Channels'},true,MyChannelNames,'Select channels to use in the display','type','logical'), ...
            arg({'plottingmode','PlotMode'},{'all', 'envelope'},{'all','envelope'},{'Plotting Mode.' '''All'' plots all individual traces.' '''Envelope'' plots the envelope.'},'type','logical'), ...
            arg({'envColor','EnvelopeColor'},[1 0 0],[],'Envelope color [R G B]'), ...
            } ...
           'GraphMetric' { ...
            arg({'metric','Metric'},'NodeColor',{'NodeColor','NodeSize'},'Select a causal metric to display.'), ...
            arg({'nodes','Nodes'},true,MyComponentNames,'Select components to use in the display','type','logical'), ...
            arg({'plottingmode','PlotMode'},{'all', 'envelope'},{'all','envelope'},{'Plotting Mode.' '''All'' plots all individual traces.' '''Envelope'' plots the envelope.'},'type','logical'), ...
            arg({'envColor','EnvelopeColor'},[1 0 0],[],'Envelope color [R G B]'), ...
            } ...
        },'Configure footer panel displayed at the bottom of the figure. If ''off'', don''t render footer. If ''ICA_ERP_Envelope'', then display the ERP envelope of backprojected components. If ''Chan_ERP_Envelope'' then display the ERP envelope of selected channels. If ''GraphMetric'' then display a graph theoretic metric (net metric value (integral) across all nodes if multiple nodes selected)','cat','DisplayProperties'), ...
    arg_sub({'BMopts','BrainMovieOptions'},[], ...
        { ...
            arg({'size','ImageSize'},[600 600],[],'Image size (pixels). Input should be [widthcond height]. If more than one condition is being plotted horizontally, then widthcond is the width of each condition subplot','cat','DisplayProperties') ...
            arg({'visible','Visibility'},'on',{'on','off'},'Figure visibility when rendering movie. If ''on,'' render frames on screen (slower). If ''off,'' keep them hidden (faster).','cat','DisplayProperties'), ...
            arg_nogui({'latency','LatenciesToRender'},[],[],'Subset of latencies to render (sec). Must be in TimeRange. Can be a vector The time point closest to the latency given are plotted. If empty, render all latencies in TimeRange.','cat','DataProcessing'), ...
            arg_nogui({'frames','FramesToRender'},[],[],'Vector of frame indices to compute. E.g. [1:2] only computes the first two frames. If empty, render all frames','cat','DataProcessing'), ...
            arg_nogui({'figurehandle','FigureHandle'},[],[],'Handle to a figure to render brainmovie in'), ...
            arg({'cameraMenu','ShowCameraMenu'},false,[],'Show Camera Menu.'), ...
            arg_subtoggle({'rotationpath3d','RotationPath3D'},[], ...
                {
                arg({'AngleFactor','VerticalRotationFactor'},1,[],'Vertical rotation multiplicative factor. Larger values produce more rotation around the x-axis (e.g. Anterior-Posterior axis).'), ...
                arg({'PhaseFactor','HorizontalRotationFactor'},0.75,[],'Horizontal rotation multiplicative factor. Larger values produce larger degrees of rotation around the z-axis (e.g. Dorsal-Rostral axis)'), ...
                arg({'FramesPerCycle'},[],[],'Number of image frames per rotation cycle. Default is all frames in sequence such that exactly one cycle will be completed in the movie. If set to 1, then continuous rotation ensues') ...
                },'Specify the rotation path for the BrainMovie.','cat','DisplayProperties'), ...
            arg({'view','InitialView'},[50 36],[],'3D static starting view for the movie.  See ''help view''. Default is [50 36].','cat','DisplayProperties'), ...
            arg({'makeCompass','MakeCompass'},false,[],{'Label cardinal directions.','e.g. Anterior, Posterior, Right, Left','cat','DisplayProperties'}), ... 
            arg({'compassLabels','CompassLabels'},{'Dorsal ','Anterior ','Right Lat '},[],'Labels for compass axes (Z, Y, X)'), ...
            arg({'project3d','ProjectGraphOnMRI'},'off',{'on','off'},'Project the 3D graph onto the 2D MRI slices','cat','DisplayProperties'), ...
            arg_sub({'theme','Theme'},{},@hlp_getBrainMovieTheme,'Brainmovie Color Theme. Each theme sets defaults for the graphic'), ...
            arg_sub('Layers', [], ...
            { ...
            arg_subtoggle({'scalp','Scalp','plotscalp'},{}, ...
                        { 
                        arg({'scalpres','Resolution'},'high',{'full','high','mid','low','custommesh'},'Resolution for the BEM scalp mesh. Full/High/Mid/Low corresponds to 100/50/10/1% of original mesh resolution.'), ...
                        arg({'volumefile','CustomMesh'},[],[],'Mesh structure. This only applies if ''Resolution'' is chosen to be ''custommesh''. This is a struct with fields ''vertices'' and ''faces'' defining the surface mesh. This can also be (1) a filename or path to custom mesh volume file or (2) the name of a variable in the base workspace containing the mesh or (3) an expression to evaluate in the workspace which returns the mesh structure.','type','expression','shape','row'), ...
                        arg({'scalptrans','Transparency','scalptransparency'},0.9,[0 1],'Transparency of the surface. Real number in [0 1], where 0=opaque, 1=transparent.'), ...
                        arg({'scalpcolor','Color'},[1,.75,.65],[],'Layer color. Can be a single vector, or a [num_vertices x 3] matrix of colors','shape','matrix','cat','Layers'), ...
                        },'Superimpose semi-transparent mesh of scalp on brainmovie. UseOpenGL must be set to "on"','cat','Layers'), ...
            arg_subtoggle({'skull','Skull','plotskull'},[], ...
                        { 
                        arg({'skullres','Resolution'},'high',{'full','high','mid','low','custommesh'},'Resolution for the BEM skull mesh. Full/High/Mid/Low corresponds to 100/50/10/1% of original mesh resolution.'), ...
                        arg({'volumefile','CustomMesh'},[],[],'Mesh structure. This only applies if ''Resolution'' is chosen to be ''custommesh''. This is a struct with fields ''vertices'' and ''faces'' defining the surface mesh. This can also be (1) a filename or path to custom mesh volume file or (2) the name of a variable in the base workspace containing the mesh or (3) an expression to evaluate in the workspace which returns the mesh structure.','type','expression','shape','row'), ...
                        arg({'skulltrans','Transparency','skulltransparency'},0.8,[0 1],'Transparency of the surface. Real number in [0 1], where 0=opaque, 1=transparent.'), ...
                        arg({'skullcolor','Color'},[1,.75,.65],[],'Layer color. Can be a single vector, or a [num_vertices x 3] matrix of colors','shape','matrix','cat','Layers'), ...
                        },'Superimpose semi-transparent mesh of skull on brainmovie. UseOpenGL must be set to "on"','cat','Layers'), ...
            arg_subtoggle({'csf','CSF','plotcsf'},[], ...
                        { 
                        arg({'csfres','Resolution'},'high',{'full','high','mid','low','custommesh'},'Resolution for the BEM csf mesh. Full/High/Mid/Low corresponds to 100/50/10/1% of original mesh resolution.'), ...
                        arg({'volumefile','CustomMesh'},[],[],'Mesh structure. This only applies if ''Resolution'' is chosen to be ''custommesh''. This is a struct with fields ''vertices'' and ''faces'' defining the surface mesh. This can also be (1) a filename or path to custom mesh volume file or (2) the name of a variable in the base workspace containing the mesh or (3) an expression to evaluate in the workspace which returns the mesh structure.','type','expression','shape','row'), ...
                        arg({'csftrans','Transparency','csftransparency'},0.8,[0 1],'Transparency of the surface. Real number in [0 1], where 0=opaque, 1=transparent.'), ...
                        arg({'csfcolor','Color'},[1,.75,.65],[],'Layer color. Can be a single vector, or a [num_vertices x 3] matrix of colors','shape','matrix','cat','Layers'), ...
                        },'Superimpose semi-transparent mesh of CSF layer on brainmovie. UseOpenGL must be set to "on"','cat','Layers'), ...
            arg_subtoggle({'cortex','Cortex','plotcortex'},{}, ...
                        {
                        arg({'cortexres','Resolution'},'mid',{'full','high','mid','low','smooth','custommesh'},'Resolution for the LONI mesh atlas. Full/High/Mid/Low corresponds to 100/50/10/1% of original mesh resolution.'), ...
                        arg({'volumefile','CustomMesh'},[],[],'Mesh structure. This only applies if ''Resolution'' is chosen to be ''custommesh''. This is a struct with fields ''vertices'' and ''faces'' defining the surface mesh. This can also be (1) a filename or path to custom mesh volume file or (2) the name of a variable in the base workspace containing the mesh or (3) an expression to evaluate in the workspace which returns the mesh structure.','type','expression','shape','row'), ...
                        arg({'cortextrans','Transparency'},0.8,[0 1],'Transparency of the surface. Real number in [0 1], where 0=opaque, 1=transparent.'), ...
                        arg_subswitch({'cortexcolor','Color'},'LONI_Atlas', ...
                            {'LONI_Atlas' { ...
                                arg({'colormapping','Color'},{'jet'},[],'Color according to LONI Atlas. If empty then use colortable from LONI database. Can also be one of the following: An [N x 3] colortable defining the color atlas (N = number of atlas labels). A string containing a function to evaluate to return the colortable (e.g. ''jet(112)''). A cell array containing a colormap function name which accepts as input the number of labels (automatically determined) (e.g. {''jet''})','type','expression','shape','row','cat','Layers'), ...
                                } ...
                             'Constant' { ...
                                arg({'colormapping','Color'},[1 1 1],[],'Layer color. Can be a single vector, or a [num_vertices x 3] matrix of colors','shape','matrix','cat','Layers','type','expression'), ...
                                } ...
                            },'Layer color','cat','Layers') ...
                         },'Superimpose semi-transparent cortex mesh on brainmovie. UseOpenGL must be set to "on"','cat','Layers'), ...
            arg_subtoggle({'custom','CustomMesh'},[], ...
                        { 
                        arg({'volumefile','VolumeMeshFile'},[],[],'Mesh structure. This is a struct with fields ''vertices'' and ''faces'' defining the surface mesh. This can also be (1) a filename or path to custom mesh volume file or (2) the name of a variable in the base workspace containing the mesh or (3) an expression to evaluate in the workspace which returns the mesh structure.','type','expression','shape','row'), ...
                        arg({'meshtrans','Transparency','transparency'},0.8,[0 1],'Transparency of the surface. Real number in [0 1], where 0=opaque, 1=transparent.'), ...
                        arg({'meshcolor','Color'},[1,.75,.65],[],'Layer color. Can be a single vector, or a [num_vertices x 3] matrix of colors','shape','matrix','cat','Layers'), ...
                        },'Superimpose custom semi-transparent mesh on brainmovie. UseOpenGL must be set to "on"','cat','Layers'), ...            
            },'Enable Visualization Layers','cat','DisplayProperties'), ...
            arg({'facelighting','FaceLighting'},'phong',{'flat','gouraud','phong','none'},{'Lighting algorithm.', '''flat'' produces uniform lighting across each of the faces of the object. Fastest to render, but poorest quality','''gouraud'' calculates the vertex normals and interpolates linearly across the faces. Good quality and faster to render than phong. The default choice for most applications.','''phong'' interpolates the vertex normals across each face and calculates the reflectance at each pixel. Phong lighting generally produces better results than Gouraud lighting, but it takes longer to render.'}), ...
            arg({'opengl','UseOpenGL'},'on',{'on','off'},'OpenGL usage. OpenGL may cause rendering problems with MacOSX','cat','DisplayProperties'), ...
            arg({'events','EventMarkers'},{{0 'r' ':' 2}},[],'Event marker time and style. Specify event markers with a cell array of {time, linecolor, linestyle, linewidth} cell arrays. Ex. { { 0.2 ''y'' '':'' 2} { 1.5 ''r'' '':'' 2}} will render two dotted-line event makers, yellow at 200 ms and red at 1500 ms','type','expression','shape','row','cat','DisplayProperties'), ...
            arg({'flashEvents','FlashEvents'},false,[],'Generate a colored flash at event latencies. The color used will match the specified colors of the event markers in Events field above.'), ...
            arg_nogui({'square','Square'},'on',{'on','off'},'?','cat','Miscellaneous'), ...
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
                arg({'size','ImageSize'},[600 600],[],'Image size (pixels). Input should be [widthcond height]. If more than one condition is being plotted horizontally, then widthcond is the width of each condition subplot','cat','DisplayProperties') ...
                },'Options for saving the movie/figures','cat','MovieOutput'), ...
            arg({'mri'},'standard_BEM_mri.mat',[],'Dipplot MRI structure. Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a path to a Matlab file containing MRI structure. Default uses MNI brain.','type','char','shape','row'), ...
            arg({'plotimgs','ShowMRISlices'},false,[],'Display axial/saggital/coronal MRI slices on the X,Y,Z axes'), ...
            arg({'coordformat','DipoleCoordinateFormat'},defcoordformat,{'manual','nose_x+','spherical','mni','ctf','auto'},'Coordinate format for dipplot. If ''manual'', no coordinate transformation is applied','type','char','shape','row'), ...
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

if isempty(g.BMopts.outputFormat.size)
    g.BMopts.outputFormat.size = g.BMopts.size;
end

if ~isfield(g.BMopts,'renderBrainMovie')
    g.BMopts.renderBrainMovie = true;
end

if ~ismember_bc(g.collapsefun,{'max','peak','absmax','minmax'}) && strcmpi(g.edgeColorMapping,'peakfreq')
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

if ~isempty(g.BMopts.mri) 
    if evalin('base',['exist(''' g.BMopts.mri ''',''var'')'])==1
        % MRI variable is in workspace, so copy it
        g.BMopts.mri = evalin('base',g.BMopts.mri);
        
        if ~isfield(g.BMopts.mri,'anatomy')
            error('MRI structure is invalid format (see ''dipplot'' for more info)');
        end
    elseif isdir(fileparts(g.BMopts.mri))
        % User specified path to MRI file, so load it up
        tmp = load(g.BMopts.mri);
        fn = fieldnames(tmp);
        g.BMopts.mri = tmp.(fn{1});
    elseif exist(g.BMopts.mri,'file')
        % mri file is on the path but full path is not provided
        % so reconstruct the full path...
        g.BMopts.mri = which(g.BMopts.mri);
        % ... and load it
        tmp = load(g.BMopts.mri);
        fn = fieldnames(tmp);
        g.BMopts.mri = tmp.(fn{1});
    else
        % User specified an invalid path to MRI file
        error('Invalid path to MRI matlab file');
    end
end

% get the theme structure
g.BMopts.theme = hlp_getBrainMovieTheme(g.BMopts.theme);

% HANDLE LAYERS
% ---------------------------------------------

% handle scalp layer
if g.BMopts.Layers.scalp.arg_selection
    % determine the file type and load it
    switch g.BMopts.Layers.scalp.scalpres
        case 'full'
            filename = 'scalp_bem_mesh_fullres.mat';
        case 'high'
            filename = 'scalp_bem_mesh_50res.mat';
        case 'mid'
            filename = 'scalp_bem_mesh_10res.mat';
        case 'low'
            filename = 'scalp_bem_mesh_1res.mat';
        case 'custommesh'
            filename = g.BMopts.Layers.scalp.volumefile;
    end

    g = hlp_loadMesh(g,filename,'scalp');
    
    g.BMopts.Layers.scalp.color        = g.BMopts.Layers.scalp.scalpcolor;
    g.BMopts.Layers.scalp.transparency = g.BMopts.Layers.scalp.scalptrans;
else
    g.BMopts.Layers.scalp.transparency = 1;
end

% handle skull layer
if g.BMopts.Layers.skull.arg_selection
    % determine the file type and load it
    switch g.BMopts.Layers.skull.skullres
        case 'full'
            filename = 'skull_bem_mesh_fullres.mat';
        case 'high'
            filename = 'skull_bem_mesh_50res.mat';
        case 'mid'
            filename = 'skull_bem_mesh_10res.mat';
        case 'low'
            filename = 'skull_bem_mesh_1res.mat';
        case 'custommesh'
            filename = g.BMopts.Layers.skull.volumefile;
    end

    g = hlp_loadMesh(g,filename,'skull');
    
%     % check if file is in the base workspace, and if so, load from there
%     if evalin('base',sprintf('exist(''%s'',''var'');',filename))
%         g.BMopts.Layers.skull.mesh = evalin('base',filename);
%     else
%         g.BMopts.Layers.skull.mesh = load(sprintf('%s.mat',filename));
%         fn = fieldnames(g.BMopts.Layers.skull.mesh);
%         if length(fn)==1 && isstruct(g.BMopts.Layers.skull.mesh.(fn{1}))
%             g.BMopts.Layers.skull.mesh = g.BMopts.Layers.skull.mesh.(fn{1});
%         end
%     end
    
    g.BMopts.Layers.skull.color        = g.BMopts.Layers.skull.skullcolor;
    g.BMopts.Layers.skull.transparency = g.BMopts.Layers.skull.skulltrans;
else
    g.BMopts.Layers.skull.transparency = 1;
end

% handle csf layer
if g.BMopts.Layers.csf.arg_selection
    % determine the file type and load it
    switch g.BMopts.Layers.csf.csfres
        case 'full'
            filename = 'csf_bem_mesh_fullres.mat';
        case 'high'
            filename = 'csf_bem_mesh_50res.mat';
        case 'mid'
            filename = 'csf_bem_mesh_10res.mat';
        case 'low'
            filename = 'csf_bem_mesh_1res.mat';
        case 'custommesh'
            filename = g.BMopts.Layers.csf.volumefile;
    end

    g = hlp_loadMesh(g,filename,'csf');
    
%     % check if file is in the base workspace, and if so, load from there
%     if evalin('base',sprintf('exist(''%s'',''var'');',filename))
%         g.BMopts.Layers.csf.mesh = evalin('base',filename);
%     else
%         g.BMopts.Layers.csf.mesh = load(sprintf('%s.mat',filename));
%         fn = fieldnames(g.BMopts.Layers.csf.mesh);
%         if length(fn)==1 && isstruct(g.BMopts.Layers.csf.mesh.(fn{1}))
%             g.BMopts.Layers.csf.mesh = g.BMopts.Layers.csf.mesh.(fn{1});
%         end
%     end
    
    g.BMopts.Layers.csf.color        = g.BMopts.Layers.csf.csfcolor;
    g.BMopts.Layers.csf.transparency = g.BMopts.Layers.csf.csftrans;
else
    g.BMopts.Layers.csf.transparency = 1;
end

% handle LONI layer
if g.BMopts.Layers.cortex.arg_selection
    % determine the file type and load it
    switch g.BMopts.Layers.cortex.cortexres
        case 'full'
            filename = 'LONImesh_fullres.mat';
        case 'high'
            filename = 'LONImesh_50res.mat';
        case 'mid'
            filename = 'LONImesh_10res.mat';
        case 'low'
            filename = 'LONImesh_1res.mat';
        case 'smooth'
            filename = 'LONImesh_smooth_no_coreg.mat';
        case 'custommesh'
            filename = g.BMopts.Layers.cortex.volumefile;
    end

    g = hlp_loadMesh(g,filename,'cortex');
    
%     % check if file is in the base workspace, and if so, load from there
%     if evalin('base',sprintf('exist(''%s'',''var'');',filename))
%         g.BMopts.Layers.cortex.mesh = evalin('base',filename);
%     else
%         g.BMopts.Layers.cortex.mesh = load(sprintf('%s.mat',filename));
%         fn = fieldnames(g.BMopts.Layers.cortex.mesh);
%         if length(fn)==1 && isstruct(g.BMopts.Layers.cortex.mesh.(fn{1}))
%             g.BMopts.Layers.cortex.mesh = g.BMopts.Layers.cortex.mesh.(fn{1});
%         end
%     end
    
    g.BMopts.Layers.cortex.color        = g.BMopts.Layers.cortex.cortexcolor;
    g.BMopts.Layers.cortex.transparency = g.BMopts.Layers.cortex.cortextrans;
else
    g.BMopts.Layers.cortex.transparency = 1;
end

% handle custom layer
if g.BMopts.Layers.custom.arg_selection
    % load the mesh
    
    filename = g.BMopts.Layers.custom.volumefile;
    g = hlp_loadMesh(g,filename,'custom');
    
    g.BMopts.Layers.custom.color        = g.BMopts.Layers.custom.meshcolor;
    g.BMopts.Layers.custom.transparency = g.BMopts.Layers.custom.meshtrans;
else
    g.BMopts.Layers.custom.transparency = 1;
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
    
    if isfield(g.ALLEEG,'chanlocs')
        MyChannelNames = {g.ALLEEG.chanlocs.labels};
    else
        MyChannelNames = strtrim(cellstr(num2str((1:g.ALLEEG.nbchan)'))');
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

% identify the indices of the components to keep
g.BMopts.selected = find(~ismember_bc(MyComponentNames,g.nodesToExclude));

% obtain factors for resampling
if g.resample
    [p q] = rat(g.resample);
    g.resample = [p q];
end

if g.resample
    erWinCenterTimes = resample(erWinCenterTimes,g.resample(1),g.resample(2));
end


% number of sources
N=g.ALLEEG(1).CAT.nbchan;

% indices of nodes we will exclude
nodeIndicesToExclude = find(~ismember_bc(1:N,g.BMopts.selected));

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
%     if g.absthresh > 0
%         % absolute threshold
%         Conn(Conn < g.absthresh) = 0;
%     elseif g.prcthresh<1
%         % simple percentile threshold
%     %     prcthresh = repmat(prctile(c.Conn,4),(1-g.prcthresh)*100,3),[1 1 size(Conn,3) size(Conn,4)]);
    if g.useStats.arg_selection
        if isempty(Stats)
            error('SIFT:StatsThresholding:NoStats','You must supply a ''stats'' structure. See help');
            return;
        end
        % apply statistical threshold
%         Conn(Stats.p > g.useStats.alpha) = 0;

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
           
                % don't compute data for this channel pair
%                 EdgeSize(ch1,ch2,:)=zeros(1,length(timeIndices),class(Conn(1)));
%                 EdgeColor(ch1,ch2,:) = EdgeSize{ch1,ch2};
                continue;
            end
            
            % extract data
            causality = squeeze(Conn(ch1,ch2,:,:));
            
            % make row vector if necessary
            if length(g.freqsToCollapse)==1 && size(causality,1)>1 && size(causality,2)==1
                causality = causality';
            end
%             
%             if any(size(causality)==1)
%                 causality = causality(:)';
%             end
            
            % collapse matrix across frequencies for selected time range
            if length(g.freqsToCollapse)>1
                % collapse across frequency
                [causality peakidx] = hlp_collapseFrequencies( causality, ...
                  g.collapsefun,freqIndicesToCollapse,timeIndices,g.freqsToCollapse(2)-g.freqsToCollapse(1));
%             elseif size(causality,2)==1
%                 % only one frequency
%                 causality = causality(:,freqIndicesToCollapse);
%                 peakidx = freqIndicesToCollapse;
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
%         if any(EdgeColor(:) < 0)
%             % two-sided thresholds
%             prcthresh = [g.prcthresh/2 1-g.prcthresh/2]*100;
%         else
%             prcthresh = [0 1-g.prcthresh]*100;
%         end
        
        prcthresh = [0 1-g.prcthresh]*100;

        % get thresholds
        edgeSizeThresh  = prctile(abs(EdgeSize(:)) ,prcthresh);
%         edgeColorThresh = prctile(abs(EdgeColor(:)),prcthresh);
        
        EdgeSize(abs(EdgeSize)   > edgeSizeThresh(1)     ...
               & abs(EdgeSize)   < edgeSizeThresh(2))    = 0;   
%         EdgeColor(abs(EdgeColor) > edgeColorThresh(1) 	...
%                 & abs(EdgeColor) < edgeColorThresh(2))   = 0;     
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

    %             NodeSize(ch1,:) = zeros(size(EdgeSize{1,2}),class(Conn(1)));
    %             NodeColor(ch1,:) = NodeSize{ch1,1};
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
    
    
    if nargout>2
        BMout.NodeSize  = NodeSize;
        BMout.NodeColor = NodeColor;
        BMout.EdgeSize  = EdgeSize;
        BMout.EdgeColor = EdgeColor;
    end
        

    % convert to Brainmovie format
    % ----------------------------
    for ch1=1:N
        for ch2=1:N
            if ch1==ch2, continue; end
                EdgeSizeCell{ch1,ch2,cnd}   = squeeze(EdgeSize(ch1,ch2,:))';
                EdgeColorCell{ch1,ch2,cnd}  = squeeze(EdgeColor(ch1,ch2,:))';
        end
    end
    
    for ch1=1:N
        NodeSizeCell{ch1,cnd}     = squeeze(NodeSize(ch1,:));
        NodeColorCell{ch1,cnd}    = squeeze(NodeColor(ch1,:));
    end
    
    
    % Create footer panel data
    % ---------------------------
    if strcmp(g.footerPanelSpec.arg_selection,'off')
        envdata     = [];
        envtimes    = [];
    else
        switch g.footerPanelSpec.arg_selection
            
            case 'Chan_ERPenvelope'
                
                % compute ERP of selected channels
                erpvars     = find(ismember_bc(MyChannelNames, ...                 % select chans
                                       g.footerPanelSpec.chanenvelopevars));
                if isempty(erpvars)
                    g.footerPanelSpec.arg_selection = 'off';
                    break; 
                end
                erpdata     = mean(g.ALLEEG(cnd).data(erpvars,:,:),3);         % compute ERP
                
                erptimes    = g.ALLEEG(cnd).times/1000;  % convert to sec
                erptlims    = getindex(erptimes,g.timeRange);
                erpdata     = erpdata(:,erptlims(1):erptlims(2));
                erptimes    = erptimes(erptlims(1):erptlims(2));
        
                % form the footer panel title string
                tmpstr = '[Chan ';
                for k=1:min(3,length(erpvars))
                    tmpstr = [tmpstr MyChannelNames{erpvars(k)} ','];
                end
                tmpstr(end) = [];
                if k<length(erpvars)
                    tmpstr = [tmpstr '...'];
                end
                tmpstr = [tmpstr '] '];
                g.BMopts.footerPanelTitle = sprintf('ERP%s: %s', ...
                        fastif(length(erpvars)>1,'s',''), ...
                        tmpstr);
                    
            case 'ICA_ERPenvelope'
                
                % compute ERP of selected backprojected components
                erpvars     = g.ALLEEG(cnd).CAT.curComps( ...
                                ismember_bc(MyComponentNames(g.BMopts.selected),...  % select comps
                                       g.footerPanelSpec.icaenvelopevars));
                if isempty(erpvars)
                    g.footerPanelSpec.arg_selection = 'off';
                    break; 
                end               
                erpchans    = find(ismember_bc(MyChannelNames,g.footerPanelSpec.backprojectedchans));
                erpdata     = mean(g.ALLEEG(cnd).icaact(erpvars,:,:),3);       % obtain ERP
                if ~isempty(erpchans)
                    erpdata     = g.ALLEEG(cnd).icawinv(erpchans,erpvars)*erpdata;        % backproject ERPs to scalp
                end
                
                erptimes    = g.ALLEEG(cnd).times/1000;  % convert to sec;
                erptlims    = getindex(erptimes,g.timeRange);
                if diff(erptlims)~=0
                    erpdata     = erpdata(:,erptlims(1):erptlims(2));
                    erptimes    = erptimes(erptlims(1):erptlims(2));
                end
                
                % form the footer panel title string
                tmpstr = '[IC ';
                for k=1:min(3,length(erpvars))
                    tmpstr = [tmpstr MyComponentNames{g.ALLEEG(cnd).CAT.curComps==erpvars(k)} ','];
                end
                tmpstr(end) = [];
                if k<length(erpvars)
                    tmpstr = [tmpstr '...'];
                end
                tmpstr = [tmpstr '] '];
                if ~isempty(erpchans)
                    tmpstr      = [tmpstr '-> [Chan '];
                    for k=1:min(3,length(erpchans))
                        tmpstr = [tmpstr MyChannelNames{erpchans(k)} ','];
                    end
                    tmpstr(end) = [];
                    if k<length(erpchans)
                        tmpstr = [tmpstr '...'];
                    end
                    tmpstr = [tmpstr '] '];
                end
                g.BMopts.footerPanelTitle = sprintf('%sERP%s: %s', ...
                        fastif(~isempty(erpchans),'Back-projected ',''), ...
                        fastif(length(erpvars)>1 || length(erpchans)>1,'s',''), ...
                        tmpstr);
                
            case 'GraphMetric'
%                 erpvars     = g.ALLEEG(cnd).CAT.curComps( ...
%                                 ismember_bc(MyComponentNames(g.BMopts.selected),...  % select comps
%                                        g.footerPanelSpec.nodes));
                erpvars     = find(ismember_bc(MyComponentNames(g.BMopts.selected),...  % select comps
                                            g.footerPanelSpec.nodes));
                if isempty(erpvars)
                    g.footerPanelSpec.arg_selection = 'off';
                    break; 
                end
                % compute graph metric over selected nodes
                eval(sprintf('tmp=%s;',g.footerPanelSpec.metric));
                tmp(isnan(tmp) | isinf(tmp))=0;
                erpdata = tmp(erpvars,:);
                
                % if all nan, set to zeros
                if all(isnan(erpdata(:)))
                    erpdata = zeros(size(erpdata));
                end
                
                erptimes    = erWinCenterTimes(timeIndices);
                erptlims    = getindex(erptimes,g.timeRange);
                
                
                erpdata     = erpdata(:,erptlims(1):erptlims(2));
                erptimes    = erptimes(erptlims(1):erptlims(2));
                
                % if all values are the same, create a slight difference so
                % we can determine ylimits.
                if all(erpdata==erpdata(1))
                    erpdata(1)=erpdata(1)*(1-1/1000);
                    erpdata(end) = erpdata(end)*(1+1/1000);
                end
                
                switch g.footerPanelSpec.metric
                    case 'NodeColor'
                        g.BMopts.envylabel = g.nodeColorMapping;
                    case 'NodeSize'
                        g.BMopts.envylabel = g.nodeSizeMapping;
                end 
        end

        if g.resample
            % resample the envelope timecourse to make consistent
            envdata(:,:,cnd)    = resample(erpdata' ,g.resample(1),g.resample(2))';
            envtimes(:,cnd)     = resample(erptimes,g.resample(1),g.resample(2));
        else
            envdata(:,:,cnd)    = erpdata;
            envtimes(:,cnd)     = erptimes;
        end
        
        % copy color spec
        g.BMopts.envColor = g.footerPanelSpec.envColor;
    end
    
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


if ~strcmpi(g.footerPanelSpec.arg_selection,'off')
    BMopts.footerPanelTimes     = envtimes;
    BMopts.footerPanelData      = envdata;
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


% setup the rotation path for the movie
if ~BMopts.rotationpath3d.arg_selection
    BMopts.rotationpath3d = [];
end
    
% if ~BMopts.plotCortex.arg_selection
%     % don't superimpose cortex
%     BMopts.cortexTransparency = 1;
% end

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

if strcmpi(BMopts.coordformat,'nose_x+')
    BMopts.dipplotopt = [BMopts.dipplotopt {'transform',[0 0 0 0 0 6*pi/4 1 1 1]}];
    BMopts.coordformat = 'mni';
elseif strcmpi(BMopts.coordformat,'manual')
    BMopts.dipplotopt = [BMopts.dipplotopt {'transform',[0 0 0 0 0 0 1 1 1]}];
    BMopts.coordformat = 'mni';
end

% convert options structure to ('name',value) arglist
bmargs = hlp_struct2varargin(hlp_flattenStruct(BMopts,Inf,{'mri','vars','Layers','rotationpath3d','theme'}),'suppress', ...
                     {'arg_selection','arg_direct','bmopts_suppl'});

bmargs = hlp_remDupNVPs([bmargs BMopts.bmopts_suppl], false);

if g.BMopts.renderBrainMovie
    % Call BrainMovie3D
    [dummy dummy g.BMopts] = brainmovie3d_causal( ...
        NodeSizeCell,NodeColorCell,EdgeSizeCell,EdgeColorCell, ...
        BrainMovieTimeRangeInSec,1,g.BMopts.selected,bmargs{:});
    
    handles.figurehandle = g.BMopts.figurehandle;
end

% EOF




function g = hlp_loadMesh(g,filename,layer)
% helper function to load a mesh
%
%
% g         - is the config structure
% layer     - is the fieldname of the layer we are working with 
%             (scalp, csf, skull, cortex, custom)
% filename  - determines where to find the mesh. filename is is either disk 
%             filename, base varname, or base-evaluable command string.
%             filename can also be an actual mesh structure.
%

% first check if we already have the mesh structure loaded
if isstruct(g.BMopts.Layers.(layer).volumefile)
        % volumefile contains the mesh
        g.BMopts.Layers.(layer).mesh = g.BMopts.Layers.(layer).volumefile;
else
    % the mesh structure is either:
    % (1) contained in a workspace variable or
    % (2) contained in a file on disk or
    % (3) returned as the result of evaluating volumefile as a command
    %     in the base workspace 

    % check if volumefile is a variable in the base workspace...
    if evalin('base',sprintf('exist(''%s'',''var'');',filename))
        % ...if so, copy it here
        g.BMopts.Layers.(layer).mesh = evalin('base',filename);
    elseif exist(filename,'file')==2
        % volumefile is a file on disk, load it
        g.BMopts.Layers.(layer).mesh = load(filename);
        fn = fieldnames(g.BMopts.Layers.(layer).mesh);
        if length(fn)==1 && isstruct(g.BMopts.Layers.(layer).mesh.(fn{1}))
            g.BMopts.Layers.(layer).mesh = g.BMopts.Layers.(layer).mesh.(fn{1});
        end
    else
        % try to evaluate the string in the workspace and get a mesh
        try
            g.BMopts.Layers.(layer).mesh = evalin('base',filename);
        catch err
            error('Unable to evaluate command %s in workspace\n. The following error was returned: %s\n',filename,err.message);
        end
    end
end



        
        

