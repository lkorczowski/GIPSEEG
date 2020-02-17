% statistical checking
close all
r1 = 0 + 1.*randn(120,1);
r2 = 0.5 + 1.*randn(120,1);

hist([r1,r2])

mean(r1)
var(r1)
mean(r2)
var(r2)

[H,P] = ttest(r1,r2)