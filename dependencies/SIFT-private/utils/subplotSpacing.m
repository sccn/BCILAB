function h=subplotSpacing(h,numrows,numcols,space)
% subplotSpacing(h,numrows,numcols,space) shrinks or expands space between subplots
%
% h is a vector containing handles to all subplots
% numrows=number of numrows of subplots
% numcols=number of colums of subplots
% space=spacing between subplots (try 0.05 for a nice tight
% spacing)
% WARNING: you won't be able to undo this figure alteration
% NOTE: get(h,'position') returns four values: [left,bottom,width,height]
% NOTE: 'position' field of output h is modified to reflect new subplot
%       positions
%
% AUTHOR: Dan Pendleton, dep22@cornell.edu

[m n]=size(h);
if n~=numrows*numcols
    error('dimensions of subplot do not match number of elements of h')
end

for i=1:n
    B(i,:)=get(h(i),'position');
end

A=reshape((1:n)',numcols,numrows); 
A=A';

%pull columns of numrows together
ctr=1;
for j=1:numrows
   for k=2:numcols
        idx(ctr)=A(j,k);
        ctr=ctr+1;
    end
end

for i=1:length(idx)
    leftEdge=B(idx(i)-1,1)+B(idx(1)-1,3)+space;
    set(h(idx(i)),'position',[leftEdge B(idx(i),2) B(idx(i),3) B(idx(i),4)])
    B(idx(i),1)=leftEdge;
end

%pull rows together
clear idx
ctr=1;
for j=1:numrows-1
   for k=1:numcols
        idx(ctr)=A(j,k);
        ctr=ctr+1;
    end
end
for i=length(idx):-1:1
     bottom=B(idx(i)+numcols,2)+B(idx(i)+numcols,4)+space;
     set(h(idx(i)),'position',[B(idx(i),1) bottom B(idx(i),3) B(idx(i),4)]);
     B(idx(i),2)=bottom;
end