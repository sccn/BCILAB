%
% $Name:  $

n = 400; k = 7; d = 1; c = 3; e = 10;
fprintf('Sampling a %d-component %d-dimensional %d-separated mixture....\n',k,d,c);
[X,T] = gmmbvl_mixgen(n,n,k,d,c,e);
fprintf('------------------\n');

fprintf('Running 10 times normal EM with k-means initialization\n');
normal_em=[];
for i=1:10
  [W,M,R,Tlogl] = gmmbvl_em(X,T,k,0,0,0);
  normal_em = [normal_em Tlogl];
end
fprintf('Average log-likelihood %f with std. dev. %f best run: %f \n',mean(normal_em),std(normal_em), max(normal_em));

     max_k = 10;
fprintf('Running greedy EM\n');
[W,M,R,Tlogl] = gmmbvl_em(X,T,max_k,10,1,0);
title('Mixture model');


figure(2); clf;plot(Tlogl,'-o');hold on;
xlabel 'number of components'
ylabel 'log-likelihood of test set'
plot(repmat(mean(normal_em),max_k,1),'r')
plot(repmat(mean(normal_em)+std(normal_em),max_k,1),'r--')
plot(repmat(max(normal_em),max_k,1),'g')
plot(repmat(mean(normal_em)-std(normal_em),max_k,1),'r--')

[Tlogl,best_k] = max(Tlogl);
fprintf('Best number of components according to cross-validation: %d (yielding log-likelihood %f ) \n',best_k,Tlogl);


legend('Greedy EM','mean of normal EM', 'one standard deviation margins of normal EM','best result of normal EM',4);
title('Log-likelihood plots');
