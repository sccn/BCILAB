function [Mus, Covs, Ws]=gmmbvl_rand_split(P,X,M,R,sigma,F,W,nr_of_cand)
%
% $Name:  $

k       = size(R,1);
[n,d]   = size(X);

epsilon = 1e-2;      % threshold in relative loglikelihood improvement for convergence in local partial EM

[tmp,I] = max(P,[],2);

Mus = [];
Covs = [];
K = [];
Ws = [];
KL = [];



for i=1:k
		
	XI        = find(I==i);
	Xloc      = X(XI,:);
	start     = size(Mus,1);
	j=0;

	if length(XI) > 2*d  % generate candidates for this parent

		% number of candidates per parent component
		while j < nr_of_cand
			r  = randperm(length(XI));
			r  = r(1:2);
			if d==1
				cl = [Xloc-Xloc(r(1)) Xloc-Xloc(r(2))];
				[tmp,cl] = min(cl.^2,[],2);
			else			
				cl = gmmbvl_sqdist( Xloc', Xloc(r,:)' );
				[tmp,cl] = min(cl,[],2);
			end
			for guy = 1:2
				data = Xloc( find( cl==guy ), :);
				if size(data,1) > d
					j = j + 1;
					Mus  = [Mus; mean(data)];
					Rloc = cov(data) + eye(d)*eps;
					Rloc = chol(Rloc);
					Covs = [Covs; Rloc(:)'];
					Ws   = [Ws W(i)/2];
					Knew = zeros(n,1);
					Knew(XI) = gmmbvl_em_gauss( ...
					   Xloc,Mus(end,:),Covs(end,:) );
					K = [K Knew];
				end
			end
		end
	end


	last=size(Mus,1);
	if last > start
		% if candidates were added, do local partial EM
		alpha = Ws(start+1:last);
		K2 = K(XI,start+1:last);
		Mnew = Mus(start+1:last,:);
		Rnew = Covs(start+1:last,:);
		FF   = F(XI)*ones(1,last-start);
		PP   = FF.*(ones(length(XI),1)*(1-alpha)) + ...
		          K2.*(ones(length(XI),1)*alpha);
		Pnew = (K2.*(ones(length(XI),1)*alpha))./PP;
		OI   = ones(n,1);
		OI(XI) = 0;
		OI = find(OI==1);
		lpo   = sum(log(F(OI)));
		ll = sum(log(PP)) + length(OI)*log(1-alpha)+lpo;
		ll = ll/n;
		done = 0;
		iter = 1;
		
		while ~done
			[alpha,Mnew,Rnew] = gmmbvl_em_step_partial( ...
			   Xloc, alpha, Mnew, Rnew, Pnew, n, 0 );
			K2 = gmmbvl_em_gauss(Xloc,Mnew,Rnew);
			Fnew = FF.*(ones(length(XI),1)*(1-alpha)) + ...
			   K2.*(ones(length(XI),1)*alpha);
			old_ll = ll;
			ll = sum(log(Fnew))+length(OI)*log(1-alpha)+lpo;
			ll = ll/n;
			done = abs(max(ll/old_ll -1)) < epsilon;
			if iter > 20
				done=1;
			end;
			iter = iter+1;
			Pnew = (K2.*(ones(length(XI),1)*alpha))./Fnew;
		end   
		Pnew(find(Pnew<eps)) = eps;
		Pnew(find(Pnew==1)) = 1-eps;
		Ws(start+1:last) = alpha;
		Mus(start+1:last,:) = Mnew;
		Covs(start+1:last,:) = Rnew;
		KL = [KL n*log(1-alpha)-sum(log(1-Pnew))];
	end
end

I = [];
for i=1:length(Ws) % remove some candiates that are unwanted
	S = reshape(Covs(i,:),d,d);
	S = S'*S;
	S = min(eig(S));
	if (S<sigma/400 | Ws(i)<2*d/n  | Ws(i)>.99)
		I = [I i];
	end
end
Ws(I) = [];
KL(I) = [];
Mus(I,:) = [];
Covs(I,:) = [];


if isempty(Ws)
	Ws = 0;
else
	[logl sup] = max(KL);
	sup = sup(1);
	Mus = Mus(sup,:);
	Covs = Covs(sup,:);
	Ws = Ws(sup);
end


