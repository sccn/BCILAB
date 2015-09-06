
function p1 = hlp_plotmesh(faces, vertex, normal, newfig,hax,FaceColor,colortheme)
% plotmesh() - plot mesh defined by faces and vertex
%
% Usage: 
%     plotmesh(faces, vertex);
%
% Input:
%   faces   - array of N x 3. Each row defines a triangle. The 3 points
%             in each row are row indices in the matrix below.
%   vertex  - array of M x 3 points, (x = first colum; y=second colum
%             z=3rd column). Each row defines a point in 3-D.
%
% Optional input:
%   normal    - normal orientation for each face (for better lighting)
%   FaceColor - [1 x 3] color spec or [num_vertices x 3] matrix of colors
%               for each vertex
%   colortheme - a cell array with color theme options
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2003
% Modified: Tim Mullen, SCCN/INC/UCSD, 2011-2013
%
% Copyright (C) May 6, 2003 Arnaud Delorme, SCCN/INC/UCSD,
% arno@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
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


    if nargin < 2
        help plotmesh;
        return;
    end;
    if nargin < 3 
        normal = [];
    end;
    if nargin < 4
        figure; 
        axes;
    end;
    if nargin < 5
        hax = gca;
    end
    if nargin < 6
        FaceColor  = [.8 .55 .35]*1.1; % ~= ruddy Caucasian - pick your complexion!
    end
    if nargin < 7
        colortheme = {};
    end
    
    if any(any(faces == 0)), faces = faces+1; end;
 
    if isempty(normal)
        if size(FaceColor,1)==size(vertex,1)
            p1 = patch('vertices', vertex, 'faces', faces, ...
                       'parent',hax,'FaceVertexCdata',FaceColor,'facecolor','interp','edgecolor','none');
        else
            p1 = patch('vertices', vertex, 'faces', faces, ...
                       'facecolor', FaceColor,'parent',hax);
        end
    else
       if size(FaceColor,1)==size(vertex,1)
           p1 = patch('vertices', vertex, 'faces', faces, ...
                       'vertexnormals', normal,'parent',hax,'FaceVertexCdata',FaceColor,'facecolor','interp','edgecolor','none');
       else
            p1 = patch('vertices', vertex, 'faces', faces, ...
                       'facecolor', FaceColor, 'vertexnormals', normal,'parent',hax);
       end
    end;
    
    set(p1,'EdgeColor','none')
    
    % Lights
    %Lights = [-125  125  80; ...
    %          125  125  80; ...
    %          125 -125 125; ...
    %          -125 -125 125];    % default lights at four corners
    %for i = 1:size(Lights,1)
    %    hl(i) = light('Position',Lights(i,:),'Color',[1 1 1],...
    %                  'Style','infinite');
    %end
    %camlight left;
    
    
%     lightangle(45,30);
%     lightangle(45+180,30);

    
    if isempty(colortheme)
        set(p1, 'specularcolorreflectance', 0, 'specularexponent',50);

        set(p1,'DiffuseStrength',.6,'SpecularStrength',0,...
               'AmbientStrength',.4,'SpecularExponent',5);
    else
        set(p1,colortheme{:});
    end
    
    axis off;
