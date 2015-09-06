function [U,S,V]=svdmaj(A, lambda, varargin)
opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'kinit', 10, 'kstep', 2);

MM=min(size(A));

if max(size(A))<100 || opt.kinit==MM
  [U,S,V]=svd(A);
else

  mm = inf;

  fprintf('[svdmaj]\n');

  kk=opt.kinit/opt.kstep;
  while mm>lambda && kk<MM
    kk=min(kk*opt.kstep, MM);
    fprintf('kk=%d\n',kk);
    [U,S,V]=pca(A, kk, 10); % Using Mark Tygert's pca.m
    mm=min(diag(S));
  end
end