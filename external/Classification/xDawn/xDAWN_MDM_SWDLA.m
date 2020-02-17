%% MDM 



clear all

disp('install Toolbox');

cd CovarianceToolbox/
installer
cd ../P300/
installer
cd ..


disp('Start Processing');
channels = [6 8 10:17];
users = [2 5 23 24 26];
method_mean = 'ld';
Nfilter = 2;
for i=users
    if i < 10
        Nuser = ['0' num2str(i)];
    else
        Nuser = num2str(i);
    end
    disp(Nuser);
    data = load_gijs_data(Nuser);
    data.training.s = data.training.s(:,channels);
    data.test.s = data.test.s(:,channels);
    % parameter

    Fs = data.Fs;
    window = fix(1*Fs);
    off = fix(0*Fs);

    % Compute feature
    Xtr = epoch_p300((data.training.s')',data.training.Flash,window,off);
    Ytr = data.training.Y;
    Xte = epoch_p300((data.test.s')',data.test.Flash,window,off);
    Yte = data.test.Y;
    clear Xtrdec
    for j=1:size(Xtr,3)
        Xtrdec(:,:,j) = resample(Xtr(:,:,j)',2,8)';
    end
    clear Xtedec
    for j=1:size(Xte,3)
        Xtedec(:,:,j) = resample(Xte(:,:,j)',2,8)';
    end
    
    % mdm
    y1 = mdm(Xte,Xtr,Ytr,method_mean);   
    aucmdm(i) = area_under_curve(y1',data.test.Y');
    disp(['MD AUC : ' num2str(aucmdm(i),'%10.2f')])
    
    %stepwisereg
    [~,~,w] = p300SWLDA(permute(Xtrdec,[2 1 3]),Ytr+1);
    y1 = zeros(size(Xtedec,3),1);
    for j = 1:size(Xtedec,3);
        y1(j) = sum(sum(w.*Xtedec(:,:,j)'));
    end        
    aucsw(i) = area_under_curve(y1',data.test.Y');        
    disp(['SW AUC : ' num2str(aucsw(i),'%10.2f')])
    
    %xdawn
    V = xdawn(data.training.s,data.training.Flash,data.training.Target,window);
    Ftrain  = zeros(Nfilter,size(Xtrdec,2),size(Xtrdec,3));
    for k=1:size(Ftrain,3)
        Ftrain(:,:,k) = V(:,1:Nfilter)'*Xtrdec(:,:,k);
    end
    Ftrain = reshape(Ftrain,size(Xtrdec,2)*Nfilter,[]);
    Ftest  = zeros(Nfilter,size(Xtedec,2),size(Xtedec,3));
    for k=1:size(Ftest,3)
        Ftest(:,:,k) = V(:,1:Nfilter)'*Xtedec(:,:,k);
    end
    Ftest = reshape(Ftest,size(Xtedec,2)*Nfilter,[]);
    
    % LDA classification on the xdawn result
    [~,y1] = reglda(Ftest,Ftrain,Ytr,'shcovft',{});
    aucxd(i) = area_under_curve(y1,data.test.Y');        
    disp(['XD AUC : ' num2str(aucxd(i),'%10.2f')])
    
    % BPM
    
    % classify with BPM using : 
    %   - Ftest : test features :  size(Ftest) = Nfeat x Nepochtest 
    %   - Ftrain : train features : size(Ftrain) = Nfeat x Nepochtrain
    %   - Ytr : Labels of the traning : size(Ytr) = Nepochtrain x 1
    [V D] = eig(cov(Ftrain')); % PCA
	
	%anova feature selection
    for k=1:size(Ftrain,1)
        p(k) = anova1(V(:,k)'*Ftrain,Ytr,'off');
    end
    [~,ix] = sort(p,'ascend');
	
    y1 = bpmClassify(V(:,ix(1:10))'*Ftest,V(:,ix(1:10))'*Ftrain,Ytr);
    
	aucbpm(i) = area_under_curve(1-y1,data.test.Y');        
    disp(['BPM AUC : ' num2str(aucbpm(i),'%10.2f')])
    
    
end
%save('results/Training-test.mat','aucmdm','aucsw','aucxd','users')
boxplot([aucmdm(users)',aucsw(users)',aucxd(users)' ,aucbpm(users)']);