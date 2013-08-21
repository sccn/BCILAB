function data=padNaN(data)
% Creates a padded data matrix from input structural array of spike times
% pads with NaN
% Usage: data=padNaN(data)
% Input:
% data : structural array of spike times
% Output:
% data : data matrix (zero padded)
NC=length(data); 
fnames=fieldnames(data);
for c=1:NC;
   eval(['Nsp(c)=length(data(c).' fnames{1} ');']);
end;
dtmp(1:max(Nsp),1:NC)=NaN;
for c=1:NC;
    eval(['f=data(c).' fnames{1} ';'])
    dtmp(1:Nsp(c),c)=f;
end;
data=dtmp;
