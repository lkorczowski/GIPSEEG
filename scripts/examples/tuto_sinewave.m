t=0:0.001:2*pi
s1=4*sin(t)/pi;s2=4*sin(3*t)/3/pi;s3=4*sin(5*t)/5/pi;s4=4*sin(7*t)/7/pi;
figure;
subplot(411);plot(t,[s1]')
subplot(412);plot(t,[s1+s2]')
subplot(413);plot(t,[s1+s2+s3]')
subplot(414);plot(t,[s1+s2+s3+s4]')