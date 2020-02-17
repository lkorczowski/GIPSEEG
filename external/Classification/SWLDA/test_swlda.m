function y1 = test_swlda(Xtest,w,b)
% y1 = test_swlda(Xtest,w,b)

    y1 = zeros(size(Xtest,3),1);
    for j = 1:size(Xtest,3);
        y1(j) = sum(sum(w.*Xtest(:,:,j)'))+b;
    end