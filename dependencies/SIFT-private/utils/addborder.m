function img2 = addborder(img1, t, c, stroke)
% ADDBORDER draws a border around an image
%
%    NEWIMG = ADDBORDER(IMG, T, C, S) adds a border to the image IMG with
%    thickness T, in pixels. C specifies the color, and should be in the
%    same format as the image itself. STROKE is a string indicating the
%    position of the border:
%       'inner'  - border is added to the inside of the image. The dimensions
%                  of OUT will be the same as IMG.
%       'outer'  - the border sits completely outside of the image, and does
%                  not obscure any portion of it.
%       'center' - the border straddles the edges of the image.
%
% Example:
%     load mandrill
%     X2 = addborder(X, 20, 62, 'center');
%     image(X2);
%     colormap(map);
%     axis off          
%     axis image 

%    Eric C. Johnson, 7-Aug-2008
%
%
% Copyright (c) 2009, Eric Johnson
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The MathWorks, Inc. nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.



    % Input data validation
    if nargin < 4
        error('MATLAB:addborder','ADDBORDER requires four inputs.');
    end
    
    if numel(c) ~= size(img1,3)
        error('MATLAB:addborder','C does not match the color format of IMG.');
    end

    % Ensure that C has the same datatype as the image itself.
    % Also, ensure that C is a "depth" vector, rather than a row or column
    % vector, in the event that C is a 3 element RGB vector.
    c = cast(c, class(img1));
    c = reshape(c,1,1,numel(c));

    % Compute the pixel size of the new image, and allocate a matrix for it.
    switch lower(stroke(1))
        case 'i'
            img2 = addInnerStroke(img1, t, c);
        case 'o'
            img2 = addOuterStroke(img1, t, c);
        case 'c'
            img2 = addCenterStroke(img1, t, c);
        otherwise
            error('MATLAB:addborder','Invalid value for ''stroke''.');
    end
    
    
end


% Helper functions for each stroke type
function img2 = addInnerStroke(img1, t, c)

    [nr1 nc1 d] = size(img1);

    % Initially create a copy of IMG1
    img2 = img1;
    
    % Now fill the border elements of IMG2 with color C
    img2(:,1:t,:)           = repmat(c,[nr1 t 1]);
    img2(:,(nc1-t+1):nc1,:) = repmat(c,[nr1 t 1]);
    img2(1:t,:,:)           = repmat(c,[t nc1 1]);
    img2((nr1-t+1):nr1,:,:) = repmat(c,[t nc1 1]);

end

function img2 = addOuterStroke(img1, t, c)

    [nr1 nc1 d] = size(img1);

    % Add the border thicknesses to the total image size
    nr2 = nr1 + 2*t;
    nc2 = nc1 + 2*t;
    
    % Create an empty matrix, filled with the border color.
    img2 = repmat(c, [nr2 nc2 1]);
    
    % Copy IMG1 to the inner portion of the image.
    img2( (t+1):(nr2-t), (t+1):(nc2-t), : ) = img1;

end

function img2 = addCenterStroke(img1, t, c)

    % Add an inner and outer stroke of width T/2
    img2 = addInnerStroke(img1, floor(t/2), c);
    img2 = addOuterStroke(img2, ceil(t/2), c);    

end