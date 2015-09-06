function data = nanpad(data,morder)
% input: data is [pnts x chs x trials]
% output: data is [(pnts+morder+2)*trials x nchs] with morder+2 nans
%         between each trial

[pnts nch ntr] = size(data);

if ntr==1
    return;
end

nanc=nan;
data = cat(1,data,nanc(ones(morder+2,nch,ntr))); 
data = permute(data,[1 3 2]); % need it in [pnts,trials,chns]
data = reshape(data,(pnts+morder+2)*ntr,nch);