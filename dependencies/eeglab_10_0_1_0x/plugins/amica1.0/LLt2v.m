% Converts log-likelihood of data under a model (LLt) to model probability
% given data (v). If LLt is NxMxP then so is v. 
function v = LLt2v(LLt)
size2 = size(LLt,2);
size3 = size(LLt,3);
LLt = reshape(LLt,size(LLt,1),size2*size3);
v = zeros(size(LLt));
for m = 1:size(v,1)
    v(m,:) = 1./sum(exp(bsxfun(@minus,LLt,LLt(m,:))),1);
end
v = reshape(v,size(v,1),size2,size3);
end