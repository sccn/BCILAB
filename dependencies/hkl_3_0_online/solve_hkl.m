function output = solve_hkl(data,Y,lambda,loss,varargin);

n = size(Y,1);  % number of data points

% optional parameters
mingap = 1e-3;                  % duality gap
maxactive = 300;                % maximum number of active variables
display = 0;                    % 1: if display, 0 otherwise
hull =[];                       % for warm restarts
alpha =[];                      % for warm restarts
b =[];							% for warm restarts
eta =[];						% for warm restarts
max_difference_active = 80;		% max number of active variables needed to prove optimality
gapeta = 1e-3;                  % gap for eta
conjgrad = 0;                 

maxtoadd_nece = 10;              % number of variables added after nece
maxtoadd_suff = 40;              % number of variables added after suff

% READ OPTIONAL PARAMETERS
args = varargin;
nargs = length(args);
for i=1:2:nargs
	switch args{i},
		case 'display',        display = args{i+1};
		case 'mingap',        mingap = args{i+1};
		case 'hull',        hull = args{i+1};
		case 'alpha',        alpha = args{i+1};
		case 'b',               b = args{i+1};
		case 'eta',             eta = args{i+1};
		case 'maxactive',        maxactive = args{i+1};
		case 'max_difference_active',        max_difference_active = args{i+1};
		case 'maxtoadd_nece',        maxtoadd_nece = args{i+1};
		case 'maxtoadd_suff',        maxtoadd_suff = args{i+1};
		case 'gapeta',        gapeta = args{i+1};
		case 'conjgrad',        conjgrad = args{i+1};
	end
end

% INITIALIZATION OF HULL WHEN NONE IS GIVEN
if isempty(hull),
	switch data.dag_type
		case 'grid + input_space'
			hull = ones(1,data.p);
		case 'mkl + input_space'
			hull = ones(1,data.p);
		case 'bimkl + input_space'
			hull = ones(1,data.p);
	end
end

% SOLVE THE REDUCED HKL PROBLEM
% LOAD DATA FOR THE ENTIRE HULL

[data_hull,sources] = get_reduced_graph_weights_sources(data,hull);

if isempty(eta),
	eta_hull = ones(size(hull,1),1);
else
	eta_hull = eta;
end
% make sure normalization is correct
eta_hull = eta_hull / ( sum( (data_hull.weights).^2 .* eta_hull ));

% learn parameters
switch loss
	case 'square'
		if all(hull == 1 )
			% only root of the hull
			alpha_hull = 1/(n*lambda) * ( Y - mean(Y) );
			b_hull = mean(Y);
			eta_hull = 1 / data_hull.weights^2;
			zeta_hull = eta_hull;
			kappa_hull = 1;
		else
			switch data_hull.dag_type
				case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);
			end
		end

	case 'logistic'

		if all(hull == 1 )
			% only root of the hull
			b_hull = center_logistic(Y*0,Y);
			alpha_hull = - 1/(n*lambda) * ( 1/(1+exp(-b_hull)) - Y );
			eta_hull = 1 / data_hull.weights^2;
			zeta_hull = eta_hull;
			kappa_hull = 1;
		else
			switch data_hull.dag_type
				case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);
			end



		end


end


% check duality gap of the reduced problem (should be small!)
[ gap2,current_value] = compute_duality_gap(alpha_hull,data_hull,eta_hull,zeta_hull,kappa_hull,lambda);
active = find(eta_hull .* data_hull.weights.^2 - gapeta/  size(data_hull.X,2)>1e-12);





% START NECESSARY PHASE (only if active is full)
oldhull = [];
if ( length(eta_hull)>= maxactive ) || ( length(active) <= length(eta_hull) - max_difference_active ),
	necessary_phase = -Inf;
else
	if (length(eta_hull) == length(active) )
		necessary_phase = 1;
		if display
			fprintf('\nENTERING NECESSARY PHASE\n\n');
		end
	else
		necessary_phase = 0;
	end

end


while  necessary_phase > 0

	if display
		fprintf('number of active variables = %d - number of non zero = %d, gap =%f\n',size(data_hull.hull,1), length(active),gap2 );
	end

	if isempty(sources)
		% everything is selected: exit!
		necessary_phase = -Inf;
		break;
	end


	% compute necessary conditions with single added variables
	gapnec = zeros(size(sources,1),1);
	for i=1:size(sources,1)
		[Xloc,weightloc] = get_data_reduced(sources(i,:),data);
		gapnec(i) = norm( Xloc' * alpha_hull ).^2 * lambda/2;
		gapnec(i) = gapnec(i) / weightloc / weightloc;
	end


	if ~isempty(find(gapnec >= current_value + mingap ));
		% there are some single sources to add!
		% compute the ones to add
		oldhull = hull;
		[a,b]=sort(-gapnec);
		sources = sources(b,:);
		gapnec = gapnec(b,:);
        
		addhull_sources = find(gapnec >=   current_value + mingap );
		% do not include all of them
		addhull_sources = addhull_sources(1:min(maxtoadd_nece,length(addhull_sources)));


	else
		necessary_phase = 0;
		break;


	end

	% prepare new hull for next step
	if isempty(addhull_sources)
		necessary_phase = 0;
		break;
	end
	[data_hull,sources] = update_reduced_graph_weights_sources(data,data_hull,sources,addhull_sources);
	hull = data_hull.hull;
	eta_hull = [eta_hull; zeros(length(addhull_sources),1) ];


	% learn parameters
	switch loss
		case 'square'
			switch data_hull.dag_type
				case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);
			end


		case 'logistic'

			switch data_hull.dag_type
				case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
					[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);
			end

	end


	% check duality gap of the reduced problem (should be small!)
	[ gap2,current_value] = compute_duality_gap(alpha_hull,data_hull,eta_hull,zeta_hull,kappa_hull,lambda);


	if display
		fprintf('new eta = %f\n',eta_hull(end))
	end



	active = find(eta_hull .* data_hull.weights.^2 - gapeta/  size(data_hull.X,2)>1e-12);

	if ( length(eta_hull)>= maxactive ) || ( length(active) < length(eta_hull) ),
		if ( length(eta_hull)>= maxactive )
			necessary_phase = -Inf;
		else
			necessary_phase = 0;
		end
	end
end



if ~isinf(necessary_phase),
	if display
		fprintf('\nENTERING SUFFICIENT PHASE\n\n');
	end
	% START SUFFICIENT PHASE
	oldhull = [];
	sufficient_phase = 1;
	while  sufficient_phase > 0


		active = find(eta_hull .* data_hull.weights.^2 - gapeta/  size(data_hull.X,2)>1e-12);
		if display
			fprintf('number of active variables = %d - number of non zero = %d\n',size(data_hull.hull,1), length(active) );
		end

		if isempty(sources)
			% everything is selected: exit!
			sufficient_phase = Inf;
			break;
		end




		% compute sufficient condition
		switch data_hull.dag_type
			case {'mkl + input_space','mkl + kernels','bimkl + input_space','bimkl + kernels'}

				% compute necessary conditions  instead
				gapsuff = zeros(size(sources,1),1);
				for i=1:size(sources,1)
					[Xloc,weightloc] = get_data_reduced(sources(i,:),data);
					gapsuff(i) = norm( Xloc' * alpha_hull ).^2 * lambda/2;
					gapsuff(i) = gapsuff(i) / weightloc / weightloc;
				end

			case {'grid + input_space', 'grid + kernels'}
				gapsuff = check_sufficient_gaps_efficient(sources,alpha_hull,data,lambda,hull);
		end

		gap = max(gapsuff)-current_value;
		if display
			fprintf('sufficient gap = %e - number of violating sources = %d\n', gap, length(find(gapsuff >=   current_value + mingap)));
		end


		[a,b]=sort(-gapsuff);
		sources = sources(b,:);
		gapsuff = gapsuff(b);
		addhull_sources =  find(gapsuff >=   current_value + mingap );
		addhull_sources = addhull_sources(1:min(maxtoadd_suff,length(addhull_sources)));


		% prepare new hull for next step
		if isempty(addhull_sources)
			sufficient_phase = 0;
			break;
		end
		[data_hull,sources] = update_reduced_graph_weights_sources(data,data_hull,sources,addhull_sources);
		hull = data_hull.hull;
		eta_hull = [eta_hull; zeros(length(addhull_sources),1) ];



		% learn parameters
		switch loss
			case 'square'
				switch data_hull.dag_type
					case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}

 						[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
						[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
					case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
						[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
						[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);


				end

			case 'logistic'
				switch data_hull.dag_type
					case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}

						[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
						[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
					case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
						[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
						[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);
				end


		end



		% check duality gap of the reduced problem (should be small!)
		[ gap2,current_value] = compute_duality_gap(alpha_hull,data_hull,eta_hull,zeta_hull,kappa_hull,lambda);


		active = find(eta_hull .* data_hull.weights.^2 - gapeta/  size(data_hull.X,2)>1e-12);


		if size(hull,1)>maxactive, sufficient_phase=0;  end

		if ( length(eta_hull)>= maxactive ) || ( length(active) <= length(eta_hull) - max_difference_active ),

			sufficient_phase = 0;
		end


	end

end

% clean up the problem
 active = find(eta_hull .* data_hull.weights.^2 - gapeta/  size(data_hull.X,2)>1e-12);
 hull = hull(active,:);
 eta_hull = eta_hull(active);
 
 % take hull of new active set
 
 
data_hull = get_reduced_graph_weights_sources(data,hull);
eta_hull = eta_hull / ( sum( (data_hull.weights).^2 .* eta_hull ));

switch loss
	case 'square'
		switch data_hull.dag_type
			case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}

				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
			case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);

		end

	case 'logistic'


		switch data_hull.dag_type
			case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}

				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
			case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);
		end

end


% find active variables (ones with non zero eta)
active = find(eta_hull .* data_hull.weights.^2 - gapeta/  size(data_hull.X,2)>1e-12);


% compute hull of these active variables (to make sure we only keep hulls!)
hullactive = compute_hull_withinhull(hull,active,data_hull.dag_type);

hull = hull(hullactive,:);
eta_hull = eta_hull(hullactive);
data_hull= get_reduced_graph_weights_sources(data,hull);
eta_hull = eta_hull / ( sum( (data_hull.weights).^2 .* eta_hull ));
switch loss
	case 'square'
		switch data_hull.dag_type
			case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta/10,'conjgrad',conjgrad);
			case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_square_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta/10);

		end

	case 'logistic'


		switch data_hull.dag_type
			case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}

				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta,'conjgrad',conjgrad);
			case { 'grid + kernels', 'mkl + kernels', 'bimkl + kernels'}
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap,'display',display,'gapeta',gapeta);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);
				[alpha_hull,b_hull,eta_hull,zeta_hull,kappa_hull] = solve_hkl_logistic_dual_fast_kernels(eta_hull,data_hull,Y,lambda,'mingap',mingap/10,'display',display,'gapeta',gapeta);
		end

end





%output.W =W_hull;
output.b =b_hull;
output.eta =eta_hull;
output.zeta =zeta_hull;
output.kappa =kappa_hull;
output.hull=hull;
output.alpha=alpha_hull;
output.weights = data_hull.weights;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function active = compute_hull(hull,active);


% retake the hull, just to be sure we have not removed
toadd = [];
for i=1:size(active,1)
	toadd = [ toadd; ind2subv(hull(active(i),:),1:prod(hull(active(i),:)))];
end
toadd
toadd = unique([toadd; hull(active,:)],'rows');
% now maps
[tf,active]= ismember(toadd,hull,'rows');
active = sort(active);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function gapsuff = check_sufficient_gaps_efficient(sources,alpha,data,lambda,hull)



switch data.dag_type




	case { 'grid + input_space', 'grid + kernels'}
		p = data.p;
		n = data.n;
 		global K_cache;
		global max_cache;

		global Kdir_cache;
		global hull_cache;
		global Kdir_caches;
		global source_cache;
		global n_cache;
 
		for t=1:size(sources,1)

 			tt = int32(sources(t,:))';

			%                     temp = repmat(tt,1,size(source_cache,2)) - source_cache;
			%                     distance_to_existing = sum(temp~=0,1)';

			% distance_to_existing = find_existing_source_transpose(tt,source_cache);
			% [aold,bold] = min(distance_to_existing);
            [a,b] = find_existing_source_transpose_min(tt,source_cache);
            % [aold,a]
			%if (a>0) & (mod(t,10)==1), fprintf('%d ',length(find(distance_to_existing==a))); end

			if a==0
				% already exactly in the cache
				Knew = K_cache(:,b);
			else
				Knew = K_cache(:,b);
                % only changes the ones that need to be changed
				nonzero = find(tt-source_cache(:,b)~=0);
				for jj =1:length(nonzero);
					i = nonzero(jj);
					Knew = Knew .* Kdir_caches{i}{tt(i)} ./ Kdir_caches{i}{source_cache(i,b)};
				end
				n_cache = n_cache + 1;
				K_cache(:,mod(n_cache,max_cache)+1) =  Knew;
				source_cache(:,mod(n_cache,max_cache)+1) =  tt;
			end
			gloc =  vectorize_quad_single(Knew,alpha);
			gapsuff(t) = gloc * lambda / 2;
		end
end



function [ gap2,current_value] = compute_duality_gap(alpha,data_hull,eta_hull,zeta_hull,kappa_hull,lambda)
atKa = zeros(size(data_hull.hull,1),1);
for i=1:size(data_hull.hull,1)
	atKa(i) = norm( data_hull.X(:,data_hull.groups{i})' * alpha ).^2;
end
gap2 = 0;
for i=1:length(data_hull.groups)
	gap2 = gap2 + data_hull.weights(i) * sqrt( sum( zeta_hull(data_hull.affinity{i}).^2 .* atKa(data_hull.affinity{i}) ));
end
gap2 = lambda/2 * gap2 * gap2;
gap2 = gap2 ...
	- lambda * sum( zeta_hull .* atKa ) ...
	+ lambda/2 * max( data_hull.weights.^-2 .* ( kappa_hull.^2 * atKa ) );
current_value = lambda/2 * max( data_hull.weights.^-2 .* ( kappa_hull.^2 * atKa ) );


