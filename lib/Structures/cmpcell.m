    function out=cmpcell(A,B)
    for indA=1:length(A)
        
        out(indA)=isequal(A{indA},B);
    end
    
    end
    