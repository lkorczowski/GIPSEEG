function drawSeparators(S)
% this function draw lines separators for submatrices whose size is 
% specified in vector S. It applies on current figure.

% Last version: 2011
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

    if(size(S,1) == 1) S=S'; end
    N = length(S);
    x = .5+reshape([zeros(N,1) sum(S)*ones(N,1)]',1,2*N);
    y = .5+reshape([cumsum(S) cumsum(S)]',1,2*N);
    for i=0:N-2
        line(x(2*i+1:2*i+2),y(2*i+1:2*i+2),'color','k', 'LineWidth', 2)
        line(y(2*i+1:2*i+2),x(2*i+1:2*i+2),'color','k', 'LineWidth', 2)
    end

end