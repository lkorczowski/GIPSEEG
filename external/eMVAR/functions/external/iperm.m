%%%% EXTERNAL FUNCTION
% function taken from the LiNGAM package
% complete software may be downloaded from http://www.cs.helsinki.fi/group/neuroinf/lingam/


function q = iperm( p )
   
for i=1:length(p)
   q(i) = find(p==i); 
end

