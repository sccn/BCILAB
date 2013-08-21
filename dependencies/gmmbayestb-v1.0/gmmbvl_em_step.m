function [W,M,R] = gmmbvl_em_step(X,W,M,R,P,plo)
%gmmbvl_em_step - EM learning step for multivariate Gaussian mixtures
%
%[W,M,R] = gmmbvl_em_step(X,W,M,R,P,plo)
%  X - (n x d) matrix of input data
%  W - (k x 1) vector of mixing weights
%  M - (k x d) matrix of components means
%  R - (k x d^2) matrix of Cholesky submatrices of components covariances
%      in vector reshaped format. To get the covariance of component k:
%      Rk = reshape(R(k,:),d,d); S = Rk'*Rk;
%  P - (n x k) posterior probabilities of all components (from previous EM step)
%  plo - if 1 then plot ellipses for 2-d data
%returns
%  W - (k x 1) matrix of components priors
%  M - (k x d) matrix of components means
%  R - (k x d^2) matrix of Cholesky submatrices of components covariances

% Nikos Vlassis, 2000

%
% $Name:  $

[n,d] = size(X);


if plo 
	figure(1);
	if d == 1
		plot(X,zeros(n,1),'k+');
	else
		plot(X(:,1),X(:,2),'g+');
	end
	hold on;
end

Psum = sum(P,1);

for j = 1:length(W)
	if Psum(j) > eps
		% update mixing weight
		W(j) = Psum(j) / n;

		% update mean
		M(j,:) = P(:,j)' * X ./ Psum(j);
	
		% update covariance matrix
		Mj = repmat(M(j,:),n,1);
		Sj = ((X - Mj) .* repmat(P(:,j),1,d))' * (X - Mj) ./ ...
		   repmat(Psum(j),d,d);

		% check for singularities
		[U,L,V] = svd(Sj); 
		l = diag(L);
		if (min(l) > eps) & (max(l)/min(l) < 1e4)
			[Rj,p] = chol(Sj);
			if p == 0
				R(j,:) = Rj(:)';
			end
		end

		% plot ellipses
		if plo
			if d == 1
				x = linspace(min(X) - 3*max(R), ...
				   max(X) + 3*max(R), 500 )';
				Lx = gmmbvl_em_gauss(x,M,R);
				Fx = Lx*W;
				plot(x,Fx,'k-');
			else
				Rk = reshape(R(j,:),d,d); S = Rk'*Rk;l=svd(S);
				phi = acos(V(1,1));
				if V(2,1) < 0
					phi = 2*pi - phi;
				end
				plot(M(j,1),M(j,2),'k.',M(j,1),M(j,2),'k+');
				gmmbvl_ellipse( 2*sqrt(l(1)), 2*sqrt(l(2)), ...
				   phi, M(j,1), M(j,2),'k' );
			end
		end
	end
end

if plo
	if  d==2
		a = (max(X(:,1)) - min(X(:,1))) / 10;
		b = (max(X(:,2)) - min(X(:,2))) / 10;
		axis([min(X(:,1))-a max(X(:,1))+a min(X(:,2))-b max(X(:,2))+b]);
	end
	drawnow;
	hold off;
end


