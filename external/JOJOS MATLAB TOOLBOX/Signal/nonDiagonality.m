function [NDiagFct]= nonDiagonality(A)
% This function computes the non-diagonality vector of a set of matrices A
% Input     -->     A ( N * N * K)
% Output    -->     NDiagFct (K * 1)
%
% HISTORY 
% *** 2011-02-24
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


    % Assertion on input matrices
    if size(A,1) ~= size(A,2)
        error('Attempt to compute non-diagonality function on non-square matrices')
    end
    
    % Calculation of cospectra non-diagonality:
    % (Sum square off_diagonal / sum square diagonal) / (N-1)
    [N N K] = size(A);
    SumDiag = zeros(K,1);
    SumOff  = zeros(K,1);
    for k = 1 : K
        for i = 1 : N
            for j = 1 : N
                if i == j
                    SumDiag(k) = SumDiag(k) + A(i,j,k)^2;
                else
                    SumOff(k) = SumOff(k) + A(i,j,k)^2;
                end
            end
        end
    end
    NDiagFct = SumOff./SumDiag / (N-1) ;
    NDiagFct = filter([1/3 1/3 1/3],1,NDiagFct); % moving average
    NDiagFct = filter([1/3 1/3 1/3],1,NDiagFct); % moving average

end