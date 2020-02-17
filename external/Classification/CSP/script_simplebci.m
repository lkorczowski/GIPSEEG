
Training=test{1}.s(1:end/2,:);
Tags_tr=(test{1}.Flash(1:end/2,:)+test{1}.Target(1:end/2,:))';
Fs=test{1}.Fs;
Test=test{1}.s(end/2:end,:);
Tags_te=(test{1}.Flash(end/2:end,:))';
Tags_teCHECK=(test{1}.Flash(end/2:end,:)+test{1}.Target(end/2:end,:))';
Yte=(Tags_teCHECK(Tags_teCHECK>0)-1)*2-1;
flt = @(f)(f>7&f<30).*(1-cos((f-(7+30)/2)/(7-30)*pi*4));
[S,T,w,b] = train_bci(Training, Fs, ...
    Tags_tr,[0 1],flt,3,128);

Testepoch= epoch_p300(Test',Tags_te,128);

for x=1:size(Testepoch,3)
y(x) = log(var(T*(Testepoch(:,:,x)'*S)))*w - b;
end
                % matrice de confusion
                ConfM=confusionmat(Yte',sign(y)');
                % courbe roc
                                %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
                [TPR,FPR,TH] =roc(Yte',y')
                plotroc(y,Yte)
                % Area under curve (higer is better, max=1);
                                %%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%
                [PerfX,PerfY,~,AUC{indNB},OPTROCPT] = perfcurve(Ytest',scores{indNB}',1);
                
plot((1:size(Testepoch,3))/Fs,[sign(y/sqrt(mean(y.*y)))- Yte]);
xlabel('time (seconds)'); ylabel('class');

axis([0 2 -2 2])

confusion(Yte,y)
count(sign(y)==Yte)/length(y)