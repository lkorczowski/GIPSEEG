function [weights, intercept] = train_swlda(Xtrain,Ytrain)
%  [weights, intercept] = train_swlda(Xtrain,Ytrain)
Responses = permute(Xtrain,[2 1 3]);
Regressor = Ytrain + 1;
%States = BCI2k.FeatExtraction.States;
maxFeature=60;

nchannels=size(Responses,2);
if size(Xtrain,3)>1
    NumTrials = size(Responses,3);
    ResponsesOutSizeX = size(Responses,1)*size(Responses,2);
    ResponsesReshaped = reshape(Responses, [ResponsesOutSizeX NumTrials])';
else
    NumTrials = size(Responses,1);
    ResponsesReshaped=Xtrain;
end
Label=2*(Regressor-1.5);

[b,se,pval,inmodel,sta] = stepwisefit( ...
    ResponsesReshaped,Label, ...
    'maxiter',maxFeature,'display','off','penter',0.1,'premove',0.15);%,'scale','on');
intercept= sta.intercept;

index = find(inmodel~=0);
Feature1 = b(index);


Feature = Feature1;%/sqrt(mean(abs(Feature1)));
chinmodel = reshape(inmodel,length(inmodel)/nchannels,nchannels);

[samp,ch1] = find(chinmodel==1);

MUD = zeros(size(Feature,1), 3);

MUD= [ch1 samp Feature];
weights = zeros (size( Responses,1),nchannels);

for i=1:size(MUD,1),
    weights(MUD(i,2),MUD(i,1))=MUD(i,3);
end