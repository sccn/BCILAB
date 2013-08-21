function EEG = smooth_prob(EEG, smoothlen)
if nargin<2
    smoothlen = 2;
end

llt_smooth = zeros(size(EEG.etc.amica.LLt));
smoothwnd = hann(smoothlen*EEG.srate); smoothwnd = smoothwnd/sum(smoothwnd);
for i = 1:size(llt_smooth,1)
    llt_smooth(i,:) = filtfilt(smoothwnd,1,EEG.etc.amica.LLt(i,:)')';
end
for i = 1:size(llt_smooth,1)
    EEG.etc.amica.LLt_smooth(i,:) = llt_smooth(i,:);
end

v_smooth = zeros(size(llt_smooth));
for m = 1:size(v_smooth,1)-1
    v_smooth(m,:) = 1./sum(exp(bsxfun(@minus,llt_smooth,llt_smooth(m,:))),1);
end
v_smooth(m+1,:) = ones(1,size(v_smooth,2)) - sum(v_smooth(1:m,:));
for i = 1:size(v_smooth,1)
    EEG.etc.amica.v_smooth(i,:) = v_smooth(i,:);
end

if ~isempty(EEG.etc.amica.v_cropped)
    EEG.etc.amica.v_smooth_cropped = cropvector(EEG.etc.amica.v_smooth,EEG.etc.amica.crop_start,EEG.etc.amica.crop_end);
end