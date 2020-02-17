Dir='D:\GoogleDrive\Présentations & Communications\2016-12 Thesis EEG-Hyperscanning\thesis\images\'
%Dir='C:\Users\Oui-Oui\Mes Documents GIPSA\Présentations & Communications\2016-12 Thesis EEG-Hyperscanning\thesis\images\'
% Dir='.\test\'
%% Chap 5 Composite : Optimization Constraint
%This script illustrates a multivariate Gaussian distribution and its
%marginal distributions
close all
figure;
subplot(121)
x = [-pi:0.1:pi/1.5]; %a=-?, b=?
y = [-1.8:0.1:2.5]; %c=-2, d=2
[X,Y] = meshgrid(x,y); %notice X,Y have the same size
Z = sin(X) .* cos(Y);
surf(X,Y,-Z);
hold on;
%x = -2.5:0.2:2.5/2; %a=-?, b=?
%y = -2.5:0.2:2.5/2; %c=-2, d=2
[X2,Y2] = meshgrid(x,y); %notice X,Y have the same size
Z2 =(cos(X2-0.5)+sin(Y2-0.5))*0.3;
surf(X2,Y2,Z2);
hold off
subplot(122)

surf(X,Y,-Z+Z2);
figure
subplottight(2,3,1);surf(X,Y,Z);set(gca,'xtick',[],'ytick',[],'ztick',[]);title('f')
subplottight(2,3,4);contourf(Z,25);set(gca,'xtick',[],'ytick',[])
hold on;p=FastPeakFind(-Z);plot(p(1:2:end),p(2:2:end),'ro',...
'MarkerSize',6,...
 'MarkerFaceColor',[1 ,0.0,0.0]);hold off;
subplottight(2,3,2);surf(X,Y,Z2);set(gca,'xtick',[],'ytick',[],'ztick',[]);;title('f_b')
subplottight(2,3,5);contourf(Z2,25);set(gca,'xtick',[],'ytick',[])
hold on;p=FastPeakFind(-Z2);plot(p(1:2:end),p(2:2:end),'ro',...
'MarkerSize',6,...
 'MarkerFaceColor',[1 ,0.0,0.0]);hold off;
subplottight(2,3,3);surf(X,Y,Z+Z2);set(gca,'xtick',[],'ytick',[],'ztick',[]);title('f_c')
subplottight(2,3,6);contourf(Z+Z2,25);set(gca,'xtick',[],'ytick',[])
hold on;p=FastPeakFind(-(Z+Z2)*0.2);plot(p(1:2:end),p(2:2:end),'ro',...
'MarkerSize',6,...
 'MarkerFaceColor',[1 ,0.0,0.0]);hold off;
print(gcf, [ Dir '5_cost_fun_convergence'],'-dpng','-r450')

%%
figure; hold on; 
%Plot the samples on the "floor"
plot3(Samples(:,1),Samples(:,2),zeros(size(Samples,1),1),'k.','MarkerSize',2)
%Plot the 1,2, and 3-sigma ellipses slightly above the floor
%plot3(E1(1,:), E1(2,:), 1e-3+zeros(1,size(E1,2)),'Color','g','LineWidth',2);
%plot3(E2(1,:), E2(2,:), 1e-3+zeros(1,size(E2,2)),'Color','g','LineWidth',2);
plot3(E3(1,:), E3(2,:), 1e-3+zeros(1,size(E3,2)),'Color','g','LineWidth',2);

%Plot the histograms on the walls from the data in the middle
[n_x, xout] = hist(Samples(:,1),20);%Creates 20 bars
n_x = n_x ./ ( sum(n_x) *(xout(2)-xout(1)));%Normalizes to be a pdf
[~,~,~,x_Pos,x_Height] = makebars(xout,n_x);%Creates the bar points
plot3(x_Pos, Y(end)*ones(size(x_Pos)),x_Height,'-k')

%Now plot the other histograms on the wall
[n_y, yout] = hist(Samples(:,2),20);
n_y = n_y ./ ( sum(n_y) *(yout(2)-yout(1)));
[~,~,~,y_Pos,y_Height] = makebars(yout,n_y);
plot3(X(1)*ones(size(y_Pos)),y_Pos, y_Height,'-k')

%Now plot the 1-d pdfs over the histograms
plot3(X, ones(size(X))*Y(end), Z_x,'-b','LineWidth',2); 
plot3(ones(size(Y))*X(1), Y, Z_y,'-r','LineWidth',2);

%Make the figure look nice
grid on; view(45,55);
axis([X(1) X(end) Y(1) Y(end)])

print(gcf, [ Dir '2_joint_distribution'],'-dpng','-r450')
