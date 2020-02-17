function arrangeColor(A)
% this function prepare a colormap centred around zero according to maxima
% and minima of A

    m = max([abs(max(max(A))) abs(min(min(A)))]);
    caxis([-m m])
    colormap(redblue(128))
    colorbar
end