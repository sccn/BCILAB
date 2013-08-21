% finds the likelihood of data given mod structure (AMICA result.. It might contain multiple models)
%
% Inputs:    
%   EEG      - Input EEG dataset structure
%   mod      - AMICA structure
%
% Outputs:
%   LLt      - Log-likelihood of data. NxMxP matrix where N is the number of models in AMICA structure. 
%              M is the frames per epochs, P is the number of epochs. P can also be 1. (Continuous data) 
function LLt = findLLt(EEG,mod)

data = double(EEG.data(EEG.icachansind,:,:));
size2 = size(EEG.data,2);
size3 = size(EEG.data,3);
data = reshape(data,size(data,1),size(data,2)*size(data,3));
data = data - repmat(mean(data,2),1,size(data,2));
data = mod.S*data;

mod.S = eye(size(mod.S));
mod.data_mean = zeros(size(mod.data_mean));
LLt = getmodLLt(data,mod);
if size3 == 1
    LLt = squeeze(reshape(LLt,size(LLt,1),size2,size3));
else
    LLt = reshape(LLt,size(LLt,1),size2,size3);
end
