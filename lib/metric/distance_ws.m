function a = distance_ws(A,B)
% wasserstein distance according to eq.5 in Ning et al. (2013)
    a = trace(A)+trace(B)-2*(trace((A*B)^(1/2)));