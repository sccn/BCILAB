function [d,A]=strdist(r,b,krk,cas)

%d=strdist(r,b,krk,cas) computes Levenshtein and editor distance 
%between strings r and b with use of Vagner-Fisher algorithm.
%   Levenshtein distance is the minimal quantity of character
%substitutions, deletions and insertions for transformation
%of string r into string b. An editor distance is computed as 
%Levenshtein distance with substitutions weight of 2.
%d=strdist(r) computes numel(r);
%d=strdist(r,b) computes Levenshtein distance between r and b.
%If b is empty string then d=numel(r);
%d=strdist(r,b,krk)computes both Levenshtein and an editor distance
%when krk=2. d=strdist(r,b,krk,cas) computes a distance accordingly 
%with krk and cas. If cas>0 then case is ignored.
%
%Example.
% disp(strdist('matlab'))
%    6
% disp(strdist('matlab','Mathworks'))
%    7
% disp(strdist('matlab','Mathworks',2))
%    7    11
% disp(strdist('matlab','Mathworks',2,1))
%    6     9

switch nargin
   case 1
      d=numel(r);
      return
   case 2
      krk=1;
      bb=b;
      rr=r;
   case 3
       bb=b;
       rr=r;
   case 4
      bb=b;
      rr=r;
      if cas>0
         bb=upper(b);
         rr=upper(r);
      end
end

if krk~=2
   krk=1;
end

d=[];
luma=numel(bb);	lima=numel(rr);
lu1=luma+1;       li1=lima+1;
dl=zeros([lu1,li1]);
dl(1,:)=0:lima;   dl(:,1)=0:luma;
%Distance
for krk1=1:krk
for i=2:lu1
   bbi=bb(i-1);
   for j=2:li1
      kr=krk1;
      if strcmp(rr(j-1),bbi)
         kr=0;
      end
   dl(i,j)=min([dl(i-1,j-1)+kr,dl(i-1,j)+1,dl(i,j-1)+1]);
   end
end
d=[d dl(end,end)];
end


