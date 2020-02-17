%% Test with Sine signals
ts=1e-3;
Ms=10;

ns=1e-1;
t=0:ts:Ms-ts;
S1=[];
for i=1:Ms;
    S1(i,:)=sin(i*t)+randn(1,Ms/ts)*ns;
end

%shuffle same cov matrix
C11=cov(S1')
C12=magic(Ms)/(Ms*Ms)*C11
C13=C11(:,randperm(Ms))

distance(C11,C13)
distance_riemann(C11,C13)

% different signal

S2=[];
for i=1:Ms;
    S2(i,:)=cos(i*t).^2+randn(1,Ms/ts)*ns;
end

C21=cov(S2')
C22=magic(Ms)/(Ms*Ms)*C21
C23=C21(:,randperm(Ms))


D=distance(C11,C21)
Dr=distance_riemann(C11,C21)

Sa=magic(Ms)*S1+magic(Ms)*S2;
Ca=cov(Sa')
Sb=magic(Ms)*S1+magic(Ms)*S2;
Cb=cov(Sb')
D=distance(Ca,Cb)
Dr=distance_riemann(Ca,Cb)
%% 