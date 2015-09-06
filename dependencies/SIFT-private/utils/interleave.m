function r=interleave(x,y) 
% interleave two vectors, x and y

xc=num2cell(x); 
yc=num2cell(y); 
r=cell(1,numel(x)+numel(y)); 
r(1:2:2*numel(x))=xc; 
r(2:2:2*numel(y))=yc; 
r=[r{:}]; 