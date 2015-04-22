function C = geodesic(C1,C2,t,metric)

if (nargin<4)||(isempty(metric))
    metric = 'euclid';
end

switch metric
    case 'riemann'
        C = riemann_geodesic(C1,C2,t);
    case 'logeuclid'
        C = logeuclidean_geodesic(C1,C2,t); 
    case 'opttransp'
        C = mintransp_geodesic(C1,C2,t);
    otherwise %euclidean metric
        C = euclidean_geodesic(C1,C2,t);
end