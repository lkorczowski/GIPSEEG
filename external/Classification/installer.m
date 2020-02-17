% ************************************************************************


p=mfilename('fullpath'); %find the path from this script
p=p(1:end-10);%remove 'installer'
addpath(genpath(p));
    
    disp('Classification toolbox successfully activated')
    clear p
    
