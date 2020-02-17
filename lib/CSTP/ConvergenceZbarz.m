function Conv=ConvergenceZbarz(Zbarz1,Zbarz2)
% How are different two matrices of the same size
[Ne Nt]=size(Zbarz1);
Conv=norm(Zbarz2-Zbarz1,'fro')/(Ne*Nt);
end