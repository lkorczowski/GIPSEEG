clear all
Directory= 'D:\data\Hyperscanning\MARTI\Groups\results\BCI\'

% load([Directory 'MARTI_adaptive_20151128T232147.mat'])
load(['D:\data\Hyperscanning\EKATE\Groups\mat\' 'EKATE_adaptive_20151130T110608.mat'])

cell2mat(R.AUCall)

Index1(1,:)=(ParametersIndex(R.P,'classifier',R.Parameters.classifier{1}));
% Index1(2,:)=(ParametersIndex(R.P,'classifier',R.Parameters.classifier{2}));
% Index1(3,:)=(ParametersIndex(R.P,'classifier',R.Parameters.classifier{3}));

% Index2(1,:)=(ParametersIndex(R.P,'method_dist',R.Parameters.method_dist{1}));
% Index2(2,:)=(ParametersIndex(R.P,'method_dist',R.Parameters.method_dist{2}));
% 
% toPLOT=cell2mat(R.AUCall(:,Index1(1,:)&Index2(1,:)))
% plot(mean(toPLOT,1))
% hold all
% toPLOT=cell2mat(R.AUCall(:,Index1(1,:)&Index2(2,:)))
% plot(mean(toPLOT,1))
% hold off

disp('mean AUC hyper - 1 - 2')
mean(reshape([R.AUCall{1,:}],3,18)')
disp('mean single trial perf hyper - 1 - 2')

mean(reshape([R.R.perf],3,18)')

errC={R.R.errclass};
for indR=1:length(R.R)
    succes(indR)=count(errC{indR}>0)/length(errC{indR});
end
disp('mean repetition perf hyper - 1 - 2')
mean(reshape(succes,3,18)')
