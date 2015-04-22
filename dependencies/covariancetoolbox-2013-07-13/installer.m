% Installer for the Covariance toolbox


HOME = pwd;

if isunix
    path([HOME,'/lib'],path);
    path([HOME,'/lib/distance'],path);
    path([HOME,'/lib/geodesic'],path);
    path([HOME,'/lib/riemann'],path);
    path([HOME,'/lib/visu'],path);
    path([HOME,'/lib/estimation'],path);
    path([HOME,'/lib/mean'],path);
    path([HOME,'/lib/simulation'],path);
    path([HOME,'/lib/jointdiag'],path);
    path([HOME,'/lib/classification'],path);
    path([HOME,'/lib/potato'],path);
else
    path([HOME,'\lib'],path);
    path([HOME,'\lib\distance'],path);
    path([HOME,'\lib\geodesic'],path);
    path([HOME,'\lib\riemann'],path);
    path([HOME,'\lib\visu'],path);
    path([HOME,'\lib\estimation'],path);
    path([HOME,'\lib\mean'],path);
    path([HOME,'\lib\simulation'],path);
    path([HOME,'\lib\jointdiag'],path);
    path([HOME,'\lib\classification'],path);
    path([HOME,'\lib\potato'],path);
end    
disp('Covariance toolbox activated');
