function Xs=powerspectra(x)
T=length(x);
    Xs=(abs(fft(x)).^2)/T;
