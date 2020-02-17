%% script PLVt (same freq, same phase)
t=0:0.01:10;
f1=1;f2=1.;
phi1=0;phi2=00;
a1=sin(2*pi*f1*t+phi1);
a2=sin(2*pi*f2*t+phi2);

dt=1;%sliding size
DT=10;%window length
N=length(a1);
k=1
while (dt*(k-1)+DT)<=N
samp=(1+dt*(k-1)):(dt*(k-1)+DT);
 out(k)=PLVt(a1(samp),a2(samp));
 k=k+1;
end

subplot(311);plot(smooth(out,5))

%% script PLVt (same freq, diff phase)
t=0:0.01:10;
f1=1;f2=1;
phi1=0;phi2=0.5;
a1=sin(2*pi*f1*t+phi1);
a2=sin(2*pi*f2*t+phi2);
dt=1;%sliding size
DT=10;%window length
N=length(a1);
k=1
while (dt*(k-1)+DT)<=N
samp=(1+dt*(k-1)):(dt*(k-1)+DT);
 out(k)=PLVt(a1(samp),a2(samp));
 k=k+1;
end

subplot(312);plot(smooth(out,5))

%% script PLVt (diff freq, same phase)
t=0:0.01:10;
f1=1;f2=5;
phi1=0;phi2=0.0;
a1=sin(2*pi*f1*t+phi1);
a2=sin(2*pi*f2*t+phi2);

dt=1;%sliding size
DT=10;%window length
N=length(a1);
k=1
while (dt*(k-1)+DT)<=N
samp=(1+dt*(k-1)):(dt*(k-1)+DT);
 out(k)=PLVt(a1(samp),a2(samp));
 k=k+1;
end

subplot(313);plot(smooth(out,5))