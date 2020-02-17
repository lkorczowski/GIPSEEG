% Installer for the Covariance toolbox


HOMETMP = pwd;

if isunix
    path([HOMETMP,'/lib'],path);
    path([HOMETMP,'/lib/distance'],path);
    path([HOMETMP,'/lib/geodesic'],path);
    path([HOMETMP,'/lib/riemann'],path);
    path([HOMETMP,'/lib/visu'],path);
    path([HOMETMP,'/lib/estimation'],path);
    path([HOMETMP,'/lib/mean'],path);
    path([HOMETMP,'/lib/simulation'],path);
    path([HOMETMP,'/lib/jointdiag'],path);
    path([HOMETMP,'/lib/classification'],path);
    path([HOMETMP,'/lib/potato'],path);
else
    path([HOMETMP,'\lib'],path);
    path([HOMETMP,'\lib\distance'],path);
    path([HOMETMP,'\lib\geodesic'],path);
    path([HOMETMP,'\lib\riemann'],path);
    path([HOMETMP,'\lib\visu'],path);
    path([HOMETMP,'\lib\estimation'],path);
    path([HOMETMP,'\lib\mean'],path);
    path([HOMETMP,'\lib\simulation'],path);
    path([HOMETMP,'\lib\jointdiag'],path);
    path([HOMETMP,'\lib\classification'],path);
    path([HOMETMP,'\lib\potato'],path);
end    
clear HOMETMP
disp('Covariance toolbox activated');
