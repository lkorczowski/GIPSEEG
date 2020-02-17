%% generation of vector FT surrogates (phase randomization)
% Note: the length of the series should be not be odd

%%% input:
% Y, M*N matrix of time series (each time series is in a row)

%%% output:
% Ys, M*N matrix of surrogate time series (each time series is in a row)

function Ys=surrVFT(Y)

Y=Y'; %transpose series (the following code works with column vectors)

[N,M]=size(Y);
Ys=NaN*ones(N,M);

for m=1:M
    Y(:,m)=Y(:,m)-mean(Y(:,m)); %tolgo la media (i surro sono a media nulla)
end

for cnt=1:M

    x=Y(:,cnt);
        
        % Procedures for FT surrogates
        fx=fft(x);
        szx=size(x);szx=szx(1);
        fasex=zeros(szx,1);
        mx=abs(fx);
        if szx/2==fix(szx/2)
            fasexa=2*pi*rand(szx/2-1,1)-pi;
            fasexb=-flipdim(fasexa,1);
            fasex(2:szx/2)=fasexa;
            fasex(szx/2+2:szx)=fasexb;
            fasex(1)=2*pi*rand(1)-pi;fasex(szx/2+1)=fasex(1);
        else
            szx2=fix(szx/2);
            fasexa=2*pi*rand(szx2,1)-pi;
            fasexb=-flipdim(fasexa,1);
            fasex(1:szx2)=fasexa;
            fasex(szx2+1)=2*pi*rand(1)-pi;
            fasex(szx2+2:szx)=fasexb;
        end
        fxs=mx.*(cos(fasex)+j*sin(fasex));
        xs=ifft(fxs);xs=real(xs);
%         xs=xs-mean(xs);

    Ys(:,cnt)=xs;
end

Ys=Ys'; % return to row vectors

