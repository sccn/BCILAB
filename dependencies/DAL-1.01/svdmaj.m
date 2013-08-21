function [U,S,V]=svdmaj(A, lambda,kinit,kstep)
% get some large SV's....

MM=min(size(A));
if max(size(A)) < 500 || kinit==MM
    [U,S,V]=svd(A);
else
    mm = inf;
    
    fprintf('[svdmaj]\n');
    
    kk=kinit-kstep;
    while mm>lambda && kk<MM
        kk=min(kk+kstep, MM);
        fprintf('kk=%d\n',kk);
        [U,S,V] = pca(A,kk,10); % [U,S,V]=lansvd(A, kk);
        mm=min(diag(S));
    end
end
