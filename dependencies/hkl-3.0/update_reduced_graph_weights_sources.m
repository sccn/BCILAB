function  [data_hull,sources] = update_reduced_graph_weights_sources(data,data_hull,sources,addhull_sources);

initial_hull = data_hull.hull;
hull = data_hull.hull ;
addhull = sources(addhull_sources,:);
switch data.dag_type
	case 'grid + input_space'




		% update data by adding stuff
		ai1=size(data_hull.X,2)+1;
		for i1=1:size(addhull,1)
			Xloc = get_data_reduced(addhull(i1,:),data);
			data_hull.X = [ data_hull.X, Xloc ];
			data_hull.groups{i1+size(hull,1)} = ai1:ai1+size(Xloc,2)-1;
			ai1 = ai1 + size(Xloc,2);
		end
		hull = [ hull;  addhull];

		% get affinity structure
		affinity_hull = cell(1,size(hull,1));
		weights_hull = zeros(1,size(hull,1))';
		for i=1:size(hull,1)
			affinity_hull{i} = find(all(hull - repmat(hull(i,:),size(hull,1),1) >= 0,2));
			temp = 1;
			if all(hull(i,:)==1)

				weights_hull(i) = data.weight0;
			else
				for j=1:data.p, temp = temp * data.weightA^(hull(i,j)-1); end
				weights_hull(i) = temp;
			end
		end


		groupdescs_hull = cell(1,size(hull,1));
		for i1=1:size(hull,1)
			groupdescs_hull{i1}=[];
			for i2=1:length(affinity_hull{i1})
				groupdescs_hull{i1} = [ groupdescs_hull{i1},  data_hull.groups{affinity_hull{i1}(i2)} ];
			end
		end


		% sources = find_sources_complement_grid_fast_c(hull,data.qs);
		qs=data.qs;
		save temphull hull qs initial_hull sources addhull_sources;

		% update sources
		sources = update_sources(qs,sources,addhull_sources);
% 		% recompute from scratch
% 		sources_fast = compute_sources(hull,data.qs);
% 		sources_old =find_sources_complement_grid_fast_int_c( int32(hull'),int32(data.qs));


		data_hull.weights = weights_hull;
		data_hull.hull = hull;
		data_hull.affinity = affinity_hull;
		data_hull.groupdescs = groupdescs_hull;


	case 'mkl + input_space'


		% update data by adding stuff
		ai1=size(data_hull.X,2)+1;
		for i1=1:size(addhull,1)
			Xloc = get_data_reduced(addhull(i1,:),data);
			data_hull.X = [ data_hull.X, Xloc ];
			data_hull.groups{i1+size(hull,1)} = ai1:ai1+size(Xloc,2)-1;
			ai1 = ai1 + size(Xloc,2);
		end
		hull = [ hull;  addhull];

		% get affinity structure
		affinity_hull = cell(1,size(hull,1));
		weights_hull = zeros(1,size(hull,1))';
		for i=1:size(hull,1)
			affinity_hull{i} = find(all(hull - repmat(hull(i,:),size(hull,1),1) >= 0,2));
			temp = 1;
			if all(hull(i,:)==1)

				weights_hull(i) = data.weight0;
			else
				weights_hull(i) = data.weightA;
			end
		end


		groupdescs_hull = cell(1,size(hull,1));
		for i1=1:size(hull,1)
			groupdescs_hull{i1}=[];
			for i2=1:length(affinity_hull{i1})
				groupdescs_hull{i1} = [ groupdescs_hull{i1},  data_hull.groups{affinity_hull{i1}(i2)} ];
			end
		end







		%sources = find_sources_complement_mkl(hull,data.p,data.qs);
		included_kernels = find(max(hull,[],1)==2);
		temp = 1:data.p;
		temp(included_kernels) = [];
		sources = [];
		for i=1:length(temp)
			hloc = ones(1,data.p);

			hloc(temp(i))=2;
			sources =  [ sources; hloc];
		end


		data_hull.weights = weights_hull;
		data_hull.hull = hull;
		data_hull.affinity = affinity_hull;
		data_hull.groupdescs = groupdescs_hull;




	case 'bimkl + input_space'




		% update data by adding stuff
		ai1=size(data_hull.X,2)+1;
		for i1=1:size(addhull,1)
			Xloc = get_data_reduced(addhull(i1,:),data);
			data_hull.X = [ data_hull.X, Xloc ];
			data_hull.groups{i1+size(hull,1)} = ai1:ai1+size(Xloc,2)-1;
			ai1 = ai1 + size(Xloc,2);
		end
		hull = [ hull;  addhull];

		% get affinity structure
		affinity_hull = cell(1,size(hull,1));
		weights_hull = zeros(1,size(hull,1))';
		for i=1:size(hull,1)
			temp = 1;
			if all(hull(i,:)==1)
				affinity_hull{i}=1:size(hull,1);
				weights_hull(i) = data.weight0;
			else
				affinity_hull{i}=i;
				for j=1:data.p, temp = temp * data.weightA^(hull(i,j)-1); end
				weights_hull(i) = temp;
			end
		end


		groupdescs_hull = cell(1,size(hull,1));
		for i1=1:size(hull,1)
			groupdescs_hull{i1}=[];
			for i2=1:length(affinity_hull{i1})
				groupdescs_hull{i1} = [ groupdescs_hull{i1},  data_hull.groups{affinity_hull{i1}(i2)} ];
			end
		end







		p = data.p;

		included_kernels_one = zeros(p,1);
		included_kernels_two = zeros(p,p);
		for i=1:size(hull,1);
			tt = find(hull(i,:)==2);
			if length(tt)==1
				included_kernels_one(tt)=1;
			elseif length(tt)==2
				included_kernels_one(tt(1),tt(2)) = 1;
				included_kernels_one(tt(2),tt(1)) = 1;
			end
		end

		sources = [];
		for i=1:p
			if (included_kernels_one(i)==0),
				hloc = ones(1,p);
				hloc(i)=2;
				sources =  [ sources; hloc];
			end
		end
		for i=2:p
			for j=1:i-1
				if (included_kernels_two(i,j)==0),
					hloc = ones(1,p);
					hloc(i)=2; hloc(j)=2;
					sources =  [ sources; hloc];
				end
			end
		end



		data_hull.weights = weights_hull;
		data_hull.hull = hull;
		data_hull.affinity = affinity_hull;
		data_hull.groupdescs = groupdescs_hull;




end
data_hull.dag_type = data.dag_type;

% check that the total dimension of the input space is smaller than the
% number of observations. If not -> store kernel matrices
switch data.dag_type
	case {'grid + input_space', 'mkl + input_space'}

		if size(data_hull.X,2) > size(data_hull.X,1)
			n = size(data_hull.X,1);
			data_hull.kernels = zeros(n*(n+1)/2,size(hull,1),'single');
			for i=1:size(hull,1)
				data_hull.kernels(:,i) = symmetric_vectorize_single( single(data_hull.X(:,data_hull.groups{i})*data_hull.X(:,data_hull.groups{i})' ) );
			end
			switch data.dag_type
				case 'grid + input_space'
					data_hull.dag_type = 'grid + kernels';
				case 'mkl + input_space'
					data_hull.dag_type = 'mkl + kernels';
				case 'bimkl + input_space'
					data_hull.dag_type = 'bimkl + kernels';

			end

		end



end

























