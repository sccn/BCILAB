function [W,M,R] = gmmbvl_em_step_partial(X,W,M,R,P,n_all,plo)
%
% $Name:  $

[n,d] = size(X); n1=ones(n,1);d1=ones(1,d);
if plo
	figure(1), plot(X(:,1),X(:,2),'g+');
	hold on;
end

Psum = sum(P);

for j = 1:length(W)
	if Psum(j) > realmin
		W(j) = Psum(j) / n_all;
		M(j,:) = P(:,j)' * X ./ Psum(j);
		Mj = X-n1*M(j,:);
		Sj = (Mj .* (P(:,j)*d1))' * Mj / Psum(j);
		% check for singularities
		L = svd(Sj);  % get smallest eigenvalue
		if L(d) > realmin 
			[Rj,p] = chol(Sj);
			if p == 0
				R(j,:) = Rj(:)';
			end
		end
		% plot ellipses
		if plo
			[U,L,V] = svd(Sj); 
			phi = acos(V(1,1));
			if V(2,1) < 0
				phi = 2*pi - phi;
			end
			plot(M(j,1),M(j,2),'k.',M(j,1),M(j,2),'k+');
			
			% This code is commented out, because there is no
			% variable 'l' defined.
			%gmmbvl_ellipse( 2*sqrt(l(1)), 2*sqrt(l(2)), phi, ...
			%                M(j,1), M(j,2), 'k' );
		end

	end
end

if plo
	a = (max(X(:,1)) - min(X(:,1))) / 10;
	b = (max(X(:,2)) - min(X(:,2))) / 10;
	drawnow;
	hold off;
end

