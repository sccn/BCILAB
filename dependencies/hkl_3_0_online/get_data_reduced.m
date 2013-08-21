

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


function [Xloc,weightloc] = get_data_reduced(hull,data);

switch data.dag_type
	case 'grid + input_space'
		Xloc = get_X_cell_efficient(hull,data.Xs,data.ind_Xs,data.d_Xs,data.p,data.n);

		if nargout>1
			weightloc = 1;
			if all(hull ==1)

				weightloc= data.weight0;
			else
				for j=1:data.p, weightloc = weightloc * data.weightA^(hull(j)-1); end
			end
		end


	case 'mkl + input_space'
		if all(hull==1)
			% regular first one
			Xloc = get_X_cell_efficient(hull,data.Xs,data.ind_Xs,data.d_Xs,data.p,data.n);
			if nargout>1
				weightloc= data.weight0;
			end
		else
			Xloc = [];
			% include all row
			i2 = find(hull==2);
			for iii=0:data.q-1
				hloc = hull;
				hloc(i2) = hloc(i2) + iii;
				Xloc = [ Xloc, get_X_cell_efficient(hloc,data.Xs,data.ind_Xs,data.d_Xs,data.p,data.n) ];

			end
			if nargout>1
				weightloc= data.weightA;
			end
		end

	case 'bimkl + input_space'
		tt = find(hull==2);
		if length(tt)==0
			% regular first one
			Xloc = get_X_cell_efficient(hull,data.Xs,data.ind_Xs,data.d_Xs,data.p,data.n);
			if nargout>1
				weightloc= data.weight0;
			end
		elseif length(tt)==1
			Xloc = [];
			% include all row
			i2 = find(hull==2);
			for iii=0:data.q-1
				hloc = hull;
				hloc(i2) = hloc(i2) + iii;
				Xloc = [ Xloc, get_X_cell_efficient(hloc,data.Xs,data.ind_Xs,data.d_Xs,data.p,data.n) ];

			end
			if nargout>1
				weightloc= data.weightA;
			end
		else
			Xloc1 = [];
			% include all row
			i2 = tt(1);
			for iii=0:data.q-1
				hloc = hull;
				hloc(i2) = hloc(i2) + iii;
				Xloc1 = [ Xloc1, get_X_cell_efficient(hloc,data.Xs,data.ind_Xs,data.d_Xs,data.p,data.n) ];

			end


			Xloc2 = [];
			% include all row
			i2 = tt(2);
			for iii=0:data.q-1
				hloc = hull;
				hloc(i2) = hloc(i2) + iii;
				Xloc2 = [ Xloc2, get_X_cell_efficient(hloc,data.Xs,data.ind_Xs,data.d_Xs,data.p,data.n) ];

			end
			Xloc = [ Xloc1, Xloc2 ];

			if nargout>1
				weightloc= data.weightA^2;
			end
		end

end