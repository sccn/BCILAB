function [v_smooth LLt LLt_smooth] = flt_v_LLt(eeg, modelpath, smoothlen)
% attaches the per-model log-likelihood as new channels in the dataset
% should be run directly after an appropriate highpass

% if ~exp_begindef('filter') return; end
% if ~exist('smoothlen','var') smoothlen = 7.5; end

% subset, center and sphere the data using original model parameters
origm = loadmodout(modelpath);
%load([modelpath '/chanlocs.mat']);
data = double(eeg.data(:,:));
data = data - repmat(mean(data,2),1,size(data,2));
data = origm.S*data;

% then deactivate sphering etc in the mode
origm.data_mean = zeros(size(origm.data_mean));
origm.S = eye(size(origm.S));
% and compute smoothed likelihoods
smoothwnd = hann(round(smoothlen*eeg.srate)); smoothwnd = smoothwnd/sum(smoothwnd);
LLt = getmodLLt(data,origm);
LLt_smooth = filtfilt_fast(smoothwnd,1,LLt')';

 
% % write the data -- alternative way to compute the above results...
% floatwrite(data,[modelpath '/tempdata.fdt']);
% 
% % write control parameters
% copyfile([modelpath '/input.reference'], [modelpath '/input.param']);
% % append field_dim line
% try
%     f = fopen([modelpath '/input.param'],'a');
%     fprintf(f,'field_dim %d\n',size(eeg.data,2));
% catch,end
% fclose(f);
% 
% % clear contents of the out directory
% delete([modelpath '/out/*']);
% 
% % run amica to get LLt's
% system(['/data/common/amica/amica ' fullfile(modelpath,'input.param')]);
% 
% % load LLt's and filter them 
% m = loadmodout([modelpath '/out']);
% smoothwnd = hann(smoothlen*eeg.srate); smoothwnd = smoothwnd/sum(smoothwnd);
% LLt = filtfilt_fast(smoothwnd,1,m.LLt')';

%map to probabilities
v_smooth = zeros(size(LLt));
for k = 1:size(v_smooth,1)
    v_smooth(k,:) = 1./sum(exp(bsxfun(@minus,LLt_smooth,LLt_smooth(k,:))),1); end
