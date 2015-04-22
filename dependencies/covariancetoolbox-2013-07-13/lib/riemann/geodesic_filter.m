function filteredCovset = geodesic_filter(CovSet,C,W)

% passage dans le plan tangent
S = Tangent_space(CovSet,C);

%filtering
Out = (W*((W'*W)\W'))*S;

% passage dans la variété
filteredCovset = UnTangent_space(Out,C);
