function [V, COV2] = riemann_diag(COV,methodcov)
if nargin <2
    methodcov = 'cov';
end

C = mean_covariances(COV,'riemann');
T = Tangent_space(COV,C);

[V, ~] = eig(covariances(T,methodcov));

W = V(:,end)*inv(V(:,end)'*V(:,end))*V(:,end)';
T2 = W'*T;
% T2 = V(:,end)'*T;
% T2 = cat(1,zeros(size(T,1)-1,size(T,2)),T2);
% disp(size(T2));
% T = inv(V')*T2;

COV2 = UnTangent_space(T2,C);

[V ,~] = eig(COV2(:,:,1),COV2(:,:,2));

D = diag(V'*COV(:,:,1)*V);

V = V*diag(1./sqrt(D));