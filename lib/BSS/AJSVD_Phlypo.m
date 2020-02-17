function [U,V,crit,critOFF,rhoU,rhoV] = AJSVDnew(CC,epsilon,U,V,t,Ur,Vr,orth,pl,noc,max_iter,Binitstep)
%,dUc,dVc,dU,dV

% function [U,V,crit,critoff,dU,dV,dUeucl,dVeucl] = AJSVDnew(C,epsilon,Uo,Vo,type,Uorig,Vorig,orthtype,plots,noc)
%
% Joint Approximate Bilinear Diagonalisation
%
% Description:
% ------------
%
% For all matrices C_k, we search matrices U and V, such that U' C_k V is
% approximately diagonal. All C_k need to have the same dimensions, not
% necessarily square.
%
% Inputs:
% ------- 
%
%   C       : a collection of matrices C 3D array, matrices are stacked as
%               C(:,:,k)
%   epsilon : the convergence criterion [1e-6]
%   Uo      : initial estimate for U [identity]
%   Vo      : initial estimate for V [identity]
%   type    : type of iterations
%               (0) Power Iterations
%               (1) Optimal Step Size
%               (2) Optimal Step Size of Normalised Objective
%               [3] Given's rotations (Pham), default
%               (4) Initialization step only
%   Uorig   : the original U used in the simulation for comparison only [Uo]
%   Vorig   : the original V used in the simulation for comparison only [Vo]
%   orth    : Orthogonalisation type
%               'none'  no orthogonalisation
%               ['lw']  Lödwin "symmetric orthogonalisation"
%               'gs'    Gram-Schmidt orthogonalisation
%   plots   : [0] no plots, nor written feedback, (1) plots of results and
%               user feedback
%   noc     : dimension reduction, reduce the search space to noc
%               components
%
% Outputs:
% --------
%
%   U       : final estimate of U
%   V       : final estimate of V
%   crit    : the criterion sum of diagonals (which is maximised)
%   critoff : the off-diagonal criterion proportional to
%               SUM Sum(off-diagonal squared)/Sum(diagonal squared)
%   dU      : the Riemannian corrected gradient of U
%   dV      : the Riemannian corrected gradient of V
%   dUeucl  : the Euclidean gradient of U
%   dVeucl  : the Euclidean gradient of V
%
% History:
% --------
% Algorithm developed by M. Congedo
%
% *** 2010-01-12
% Matlab code written by R. Phlypo

%% Initialisation
% NEW: 3D Array for speeding up the algorithms
[P,Q,K] = size(CC);
N = min(P,Q);

if nargin < 2 || isempty( epsilon), epsilon = 1e-6; end%epsilon = 1e-6; end
if nargin < 3 || isempty(U), U = eye(P); end
if nargin < 4 || isempty(V), V = eye(Q); end
if nargin < 5 || isempty(t), t = 3; end
if nargin < 6 || isempty(Ur), Ur = U; end
if nargin < 7 || isempty(Vr), Vr = V; end
if nargin < 8 || isempty(orth), orth = 'lw'; end
if nargin < 11 || isempty(max_iter), max_iter = 10; end
if nargin < 12 || isempty(Binitstep), Binitstep = 0; end
if ~isempty(orth) && strcmpi(orth,'gs')
    orth = 2;
elseif ~isempty(orth) && strcmpi(orth,'none')
    orth = 0;
else
    orth = 1;
end
if nargin < 9 || isempty(pl)
    pl = 0;
end
if nargin < 10 || isempty( noc), noc = N; end
% max_iter = 1e2;
N= noc;
%%%%%%%%%%%%%%%%%%%%% Start Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation U and V

% C = reshape(CC,P,Q*K);
% CT  = reshape(shiftdim(C,1),Q,P*K);

H = CC;
critOFF(1) =  ( sum(sum(sum(abs(H).^2)))/sum(diag(sum(abs(H).^2,3))) - 1)/(max(P,Q)-1)/K;

GU = U(:,1:N).'*Ur(:,1:N);
rhoU(1) = (sum(sum(abs(GU))./max(abs(GU)))-N + sum(sum(abs(GU'))./max(abs(GU')))-N)/2/(N-1)/N;
GV = V(:,1:N).'*Vr(:,1:N);
rhoV(1) = (sum(sum(abs(GV))./max(abs(GV)))-N + sum(sum(abs(GV'))./max(abs(GV')))-N)/2/(N-1)/N;

if Binitstep%t==4 || t==3 || t==5
    [U,D,nouse] = svd(reshape(CC,P,Q*K));%eig(reshape(CC,P,Q*K)*reshape(CC,P,Q*K)');
    [V,D,nouse] = svd(reshape(shiftdim(CC,1),Q,P*K));%eig(reshape(shiftdim(CC,1),Q,P*K)*reshape(shiftdim(CC,1),Q,P*K)');
else
    U = eye(P);
    V = eye(Q);
end
if pl
    h1 = figure; subplot(223)
    theta = (0:360)/180*pi;
    plot(cos(theta),sin(theta),'k'); hold on, axis equal
    quiver(zeros(1,P),zeros(1,P),Ur(1,:),Ur(2,:),0)
    %     h2 = figure;
    %     h3 = figure;
    figure(h1)
end

H = shiftdim(reshape(V.'*reshape(shiftdim(reshape(U'*reshape(CC,P,[]),P,Q,K),1),Q,[]),Q,K,P),2);
if t == 3 || t == 5
    orth = 0;
    CC = shiftdim(reshape((reshape(shiftdim(reshape(U'*reshape(CC,P,[]),P,Q,K),1),Q,[])'*V)',Q,K,P),2);
end
if Binitstep
    critOFF(2) = ( sum(sum(sum(abs(H).^2)))/sum(diag(sum(abs(H).^2,3))) - 1)/(max(P,Q)-1)/K;
    GU = U(:,1:N).'*Ur(:,1:N);
    rhoU(2) = (sum(sum(abs(GU))./max(abs(GU)))-N + sum(sum(abs(GU'))./max(abs(GU')))-N)/2/(N-1)/N;
    GV = V(:,1:N).'*Vr(:,1:N);
    rhoV(2) = (sum(sum(abs(GV))./max(abs(GV)))-N + sum(sum(abs(GV'))./max(abs(GV')))-N)/2/(N-1)/N;
end

%% Initial Objective Function Value
% h100 = figure;
crit([1 2],1) = 0;
if t==0 || t==1 || t==2
    for n = 1:N
        MVc{n} = V(:,n)*V(:,n)';
        MV{n} = zeros(P,P);
        for k = 1:K
            MV{n} = MV{n} + squeeze(CC(:,:,k))*MVc{n}*squeeze(CC(:,:,k))';
        end
        crit(1,1) = crit(1,1) + U(:,n)'*MV{n}*U(:,n);

        MUc{n} = U(:,n)*U(:,n)';
        MU{n} = zeros(Q,Q);
        for k = 1:K
            MU{n} = MU{n} + squeeze(CC(:,:,k))'*MUc{n}*squeeze(CC(:,:,k));
        end
        crit(2,1) = crit(2,1) + V(:,n)'*MU{n}*V(:,n);
    end
end


if pl
    subplot(221)
    plot(0,critOFF(1),'b+'), hold on
end


%% Actual Iterations
% h100 = figure;
% for i = 1:epsilon
msin = 1;
i = 0;

while i < max_iter-Binitstep %&& msin > epsilon %msin > 1e-3 &&
    i = i+1;
    %     fprintf('Iteration %i:\tepsilon = %.3f\n',i,msin)
    if pl
        st = ones(10,1)*'-';
        fprintf([st' 'iteration %i' st' '\n'],i)
    end
       
    if nargout > 2
        NablaU(:,:,i) = zeros(P);
        NablaV(:,:,i) = zeros(Q);

        U_old = U;
        V_old = V;
        crit([1 2],i+1) = 0;
    end

    for n = 1:N
        if t==0 || t==1 || t==2
            MVc{n} = V(:,n)*V(:,n)';
            MV{n} = zeros(P,P);
            for k = 1:K
                MV{n} = MV{n} + squeeze(CC(:,:,k))*MVc{n}*squeeze(CC(:,:,k))';
            end
        end
                
        switch t
            case 0 % Power Iterations
                U(:,n) = MV{n}*U(:,n);%/(U(:,n)'*MV{n}'*MV{n}*U(:,n))^.5;
                NablaU(:,:,i) = NablaU(:,:,i) + (MV{n}*U_old(:,n)-U_old(:,n))*[zeros(1,n-1) 1 zeros(1,P-n)];
                
            case 1 % Derivative Based
                alpha = -U(:,n)'*MV{n}*MV{n}*U(:,n)/( U(:,n)'*MV{n}'*MV{n}*MV{n}*U(:,n) )*sum(abs(MV{n}*U(:,n)).^2).^.5;
                d2Jdalpha2 = U(:,n)'*MV{n}*MV{n}*MV{n}*U(:,n)/(U(:,n)'*MV{n}'*MV{n}*U(:,n));
                if sign( d2Jdalpha2)<0
                    U(:,n) = U_old(:,n) + alpha*MV{n}*U(:,n)*sum((MV{n}*U(:,n)).^2).^-.5;
                end
                NablaU(:,:,i) = NablaU(:,:,i) + alpha*MV{n}*U(:,n)*sum((MV{n}*U(:,n)).^2).^-.5*[zeros(1,n-1) 1 zeros(1,P-n)];

            case 2 % Derivative Based Rayleigh Quotient
                beta = 2*( MV{n}*U(:,n)/(U(:,n)'*U(:,n)) - U(:,n)'*MV{n}*U(:,n)/(U(:,n)'*U(:,n))^2*U(:,n) );
                MVtilde = MV{n} - eye(size(MV{n}));
                c = [-beta'*beta*U(:,n)'*MV{n}*beta+U(:,n)'*beta*beta'*MV{n}*beta ...
                    beta'*MV{n}*beta*U(:,n)'*U(:,n) - beta'*beta*U(:,n)'*MV{n}*U(:,n) ...
                    U(:,n)'*U(:,n)*beta'*MV{n}*U(:,n) - beta'*U(:,n)*U(:,n)'*MV{n}*U(:,n)];
                alpha = (-c(2) + [1; -1]*sqrt(c(2)^2-4*c(1)*c(3)))/2/c(1);

                alpha = real(alpha);
                UU = U_old(:,n)*ones(1,2) + beta*alpha';
                PHI = diag(UU'*MV{n}*UU) ./ diag(UU'*UU);
                [nu,ix] = max(PHI);
                alpha = alpha(ix);

                U(:,n) = UU(:,ix);

                NablaU(:,:,i) = NablaU(:,:,i) + beta*alpha*[zeros(1,n-1) 1 zeros(1,P-n)];
            case 3 % Given's Rotations (Pham)
                if n>1, break; end
                msin = 0;
                for q = 1:N
                    for p = q+1:P
                        %                         fprintf('Updating U:\tTreating pair (%i,%i)\n',p,q)
                        if p<=N
                            temp = [sum(abs(CC(p,p,:)).^2)+sum(abs(CC(q,q,:)).^2), -sum(CC(p,p,:).*CC(q,p,:))+sum(CC(q,q,:).*CC(p,q,:)); ...
                                -sum(CC(p,p,:).*CC(q,p,:))+sum(CC(q,q,:).*CC(p,q,:)), sum(CC(p,q,:).^2)+sum(CC(q,p,:).^2)]/K;
                            [ev,di] = eig(temp);
                            [nu,ix] = max(diag(di));
                            ev = sign(ev(1,ix))*ev(:,ix)/sum(ev(:,ix).^2)^.5;
                            %                             CC_old = CC;
                            CC([p q],:,:) = reshape([ev(1) -ev(2); ev(2) ev(1)]*reshape(CC([p q],:,:),2,Q*K),2,Q,K);
                            %                             CC(p,:,:) = CC_old(p,:,:)*ev(1) - CC_old(q,:,:)*ev(2);
                            %                             CC(q,:,:) = CC_old(p,:,:)*ev(2) + CC_old(q,:,:)*ev(1);
                            U(:,[p q]) = U(:,[p q])*[ev(1) ev(2); -ev(2) ev(1)];
                        else
                            temp = [sum(CC(q,q,:).^2), -sum(CC(q,q,:).*CC(p,q,:)); ...
                                -sum(CC(q,q,:).*CC(p,q,:)), sum(CC(p,q,:).^2)]/K;
                            [ev,di] = eig(temp);
                            [nu,ix] = max(diag(di));
                            ev = sign(ev(1,ix))*ev(:,ix)/sum(ev(:,ix).^2)^.5;
                            %                             CC_old = CC;
                            CC([p q],:,:) = reshape([ev(1) -ev(2); ev(2) ev(1)]*reshape(CC([p q],:,:),2,Q*K),2,Q,K);
                            %                             CC(p,:,:) = ev(1)*CC_old(p,:,:) - ev(2)*CC_old(q,:,:);
                            %                             CC(q,:,:) = ev(2)*CC_old(p,:,:) + ev(1)*CC_old(q,:,:);
                            U(:,[p q]) = U(:,[p q])*[ev(1) ev(2); -ev(2) ev(1)];%[ev(1)*U(:,p) - ev(2)*U(:,q), ev(2)*U(:,p) + ev(1)*U(:,q)];
                        end
                        msin = max(msin, abs(ev(2)));
                    end
                end
                %                 figure(h100); subplot(1,4,[2:4]), cla, eegplot(U(:,1:N),0,[],256/1000); hold off;
                if nargout > 2
                    NablaU(:,:,i) = U - U_old;
                end
            case 4 % Initialization only
                if n>1, break; end
                NablaU(:,:,i) = U - U_old;
        end
    end

    if orth == 1
        [R,S,T] = svd(U);
        U = R*T';
    elseif orth == 2
        [U,R] = qr(U);
    end

    % Check the Criterion
    if t==1 || t==2
        for n = 1:N
            crit(1,i+1) = crit(1,i+1) + U(:,n)'*MV{n}*U(:,n);
        end
    end

    if pl % 2D Plotting
        subplot(223)
        quiver(U_old(1,:),U_old(2,:),U(1,:)-U_old(1,:),U(2,:)-U_old(2,:),0); hold on

        figure;

        theta = (0:360)/180*pi;
        plot(cos(theta),sin(theta),'color',[.7 .7 .7],'linewidth',2); hold on, axis equal
        plot(U_old(1,:),U_old(2,:),'ko','linewidth',2,'markersize',8)
        plot(U(1,:),U(2,:),'k*','linewidth',2,'markersize',8)
        if t==0
            quiver(U_old(1,:),U_old(2,:),squeeze(NablaU(1,:,i)),squeeze(NablaU(2,:,i)),0,'linewidth',2,'color',[0 0 0]); hold on
            line([U_old(1,:)+squeeze(NablaU(1,:,i)); U(1,:)], [U_old(2,:)+squeeze(NablaU(2,:,i)); U(2,:)],'linestyle',':','color',[0 0 0],'linewidth',1)
        end
        alpha = 0:.001:1;
        temp1 = angle(U_old(1,:)+sqrt(-1)*U_old(2,:));
        temp2 = angle(U(1,:)+sqrt(-1)*U(2,:));
        [temp,ix] = min([mod(2*pi-(temp2-temp1),2*pi); mod(temp2-temp1,2*pi)]);
        temp = sign(ix-1.5).*temp;

        plot(cos(ones(length(alpha),1)*temp1+alpha'*temp),sin(ones(length(alpha),1)*temp1+alpha'*temp),'k','linewidth',2)
        axis equal
        Mxy = 3;%max(abs([xlim ylim]));
        xlim([-Mxy Mxy])
        ylim([-Mxy Mxy])
        plot((-10:.1:10)'*Ur(1,:),(-10:.1:10)'*Ur(2,:),'linewidth',2,'linestyle','--','color',[.7 .7 .7])

        %         quiver(U_old(1,:),U_old(2,:),U(1,:)-U_old(1,:),U(2,:)-U_old(2,:),0);
        %         hold off
        %pause(.5)
    end

    for n = 1:N
        if ismember(t,[0 1 2])
            MUc{n} = U(:,n)*U(:,n)';
            MU{n} = zeros(Q,Q);
            for k = 1:K
                MU{n} = MU{n} + squeeze(CC(:,:,k))'*MUc{n}*squeeze(CC(:,:,k));
            end
        end
        switch t
            case 0
                V(:,n) = MU{n}*V(:,n);%/(V(:,n)'*MU{n}'*MU{n}*V(:,n))^.5;
                NablaV(:,:,i) = NablaV(:,:,i) + (MU{n}*V_old(:,n)-V_old(:,n))*[zeros(1,n-1) 1 zeros(1,Q-n)];

            case 1
                alpha = -V(:,n)'*MU{n}*MU{n}*V(:,n)/( V(:,n)'*MU{n}*MU{n}*MU{n}*V(:,n) )*sum((MU{n}*V(:,n)).^2).^.5;
                d2Jdalpha2 = V(:,n)'*MU{n}*MU{n}*MU{n}*V(:,n)/(V(:,n)'*MU{n}*MU{n}*V(:,n));

                if sign( d2Jdalpha2)<0
                    V(:,n) = V_old(:,n) + alpha*MU{n}*V(:,n)*sum((MU{n}*V(:,n)).^2).^-.5;
                end

                NablaV(:,:,i) = NablaV(:,:,i) + alpha*MU{n}*V(:,n)*sum((MU{n}*V(:,n)).^2).^-.5*[zeros(1,n-1) 1 zeros(1,Q-n)];
            case 2
                beta = 2*( MU{n}*V(:,n)/(V(:,n)'*V(:,n)) - V(:,n)'*MU{n}*V(:,n)/(V(:,n)'*V(:,n))^2*V(:,n) );
                c = [-beta'*beta*V(:,n)'*MU{n}*beta + V(:,n)'*beta*beta'*MU{n}*beta ...
                    beta'*MU{n}*beta*V(:,n)'*V(:,n) - beta'*beta*V(:,n)'*MU{n}*V(:,n) ...
                    V(:,n)'*V(:,n)*beta'*MU{n}*V(:,n) - beta'*V(:,n)*V(:,n)'*MU{n}*V(:,n)];
                alpha = (-c(2) + [1; -1]*sqrt(c(2)^2-4*c(1)*c(3)))/2/c(1);

                alpha = real(alpha);
                VV = V_old(:,n)*ones(1,2) + beta*alpha';
                PHI = diag(VV'*MU{n}*VV) ./ diag(VV'*VV);
                [nu,ix] = max(PHI);
                alpha = alpha(ix);

                V(:,n) = VV(:,ix);

                NablaV(:,:,i) = NablaV(:,:,i) + beta*alpha*[zeros(1,n-1) 1 zeros(1,Q-n)];
            case 3 % Given's rotations (Pham)
                if n>1, break; end
                %                 fprintf('\nIteration %i\n\n',i)
                for p = 1:N
                    for q = p+1:Q
                        %                         fprintf('Updating U:\tTreating pair (%i,%i)\n',p,q)
                        if q<=N
                            temp = [sum(CC(p,p,:).^2)+sum(CC(q,q,:).^2), sum(CC(q,q,:).*CC(q,p,:))-sum(CC(p,p,:).*CC(p,q,:)); ...
                                sum(CC(q,q,:).*CC(q,p,:))-sum(CC(p,p,:).*CC(p,q,:)), sum(CC(p,q,:).^2)+sum(CC(q,p,:).^2)]/K;
                            [ev,di] = eig(temp);
                            [nu,ix] = max(diag(di));
                            ev = sign(ev(1,ix))*ev(:,ix)/sum(ev(:,ix).^2)^.5;
                            %                             CC_old = CC;
                            CC(:,[p q],:) = shiftdim(reshape([ev(1) -ev(2); ev(2) ev(1)]*reshape(shiftdim(CC(:,[p q],:),1),2,P*K),2,K,P),2);
                            %                             CC(:,p,:) = ev(1)*CC_old(:,p,:) - ev(2)*CC_old(:,q,:);
                            %                             CC(:,q,:) = ev(2)*CC_old(:,p,:) + ev(1)*CC_old(:,q,:);
                            V(:,[p q]) = V(:,[p q])*[ev(1) ev(2); -ev(2) ev(1)];
                        else
                            temp =  [sum(CC(p,p,:).^2), -sum(CC(p,p,:).*CC(p,q,:)); ...
                                -sum(CC(p,p,:).*CC(p,q,:)), sum(CC(p,q,:).^2)]/K;
                            [ev,di] = eig(temp);
                            [nu,ix] = max(diag(di));
                            ev = sign(ev(1,ix))*ev(:,ix)/sum(ev(:,ix).^2)^.5;
                            %                             CC_old = CC;
                            CC(:,[p q],:) = shiftdim(reshape([ev(1) -ev(2); ev(2) ev(1)]*reshape(shiftdim(CC(:,[p q],:),1),2,P*K),2,K,P),2);
                            %                             CC(:,p,:) = ev(1)*CC_old(:,p,:) - ev(2)*CC_old(:,q,:);
                            %                             CC(:,q,:) = ev(2)*CC_old(:,p,:) + ev(1)*CC_old(:,q,:);
                            V(:,[p q]) = V(:,[p q])*[ev(1) ev(2); -ev(2) ev(1)];
                        end
                        msin = max(msin, abs(ev(2)));
                    end
                end
                %                 figure(h100); subplot(4,4,1), imagesc(V);
                if nargout > 2
                    NablaV(:,:,i) = V- V_old;
                end
            case 4
                if n>1, break; end
                NablaV(:,:,i) = V- V_old;
        end
    end
    if orth == 1
        [R,S,T] = svd(V);
        V = R*T';
    elseif orth == 2
        [V,R] = qr(V);
    end

    
    for pp = 1:size(U,2)
        U(:,pp) = abs(U(:,pp)).*( cos(atan2(imag(U(:,pp)),real(U(:,pp))) - atan2(imag(U(:,1)),real(U(:,1)))) + sqrt(-1)*sin(atan2(imag(U(:,pp)),real(U(:,pp))) - atan2(imag(U(:,1)),real(U(:,1)))) );
    end
    for pp = 1:size(V,2)
        V(:,pp) = abs(V(:,pp)).*( cos(atan2(imag(V(:,pp)),real(V(:,pp))) - atan2(imag(V(:,1)),real(V(:,1)))) + sqrt(-1)*sin(atan2(imag(V(:,pp)),real(V(:,pp))) - atan2(imag(V(:,1)),real(V(:,1)))) );
    end
    
    if nargout > 2

        if t==1 || t==2
            for n = 1:N
                crit(2,i+1) = crit(2,i+1) + V(:,n)'*MU{n}*V(:,n);
            end
        end

        if t~=3 && t~=5
            H = shiftdim(reshape((reshape(shiftdim(reshape(U'*reshape(CC,P,[]),P,Q,K),1),Q,[])'*V)',Q,K,P),2);
        else
            H = CC;
        end
        
        % This line can be commented out
%         ( sum(sum(sum(H.^2)))/sum(diag(sum(H.^2,3))) - 1)/(max(P,Q)-1)/K
        critOFF(i+1+Binitstep) =  ( sum(sum(sum(abs(H).^2)))/sum(diag(sum(abs(H).^2,3))) - 1)/(max(P,Q)-1)/K;

        derivU = squeeze(NablaU(:,:,i));
        derivV = squeeze(NablaV(:,:,i));

        derivUcorr = U-U_old;
        derivVcorr = V-V_old;
    end
    if pl
        figure(h1)
        subplot(221), hold on

        %     plot(i,sum(sum(NablaU(:,:,i).^2)),'r+');
        %     plot(i,sum(sum(NablaV(:,:,i).^2)),'g.');
        %     plot(i,crit(1,i),'ro')
        %     plot(i,crit(2,i),'go')
        plot(i,critOFF(i+1),'b+')

        subplot(222)
        imagesc(U'*Ur), colorbar

        subplot(224)
        imagesc(V'*Vr), colorbar

        %         figure(h3)
        %         plot(i,sum(sum(derivUcorr.^2))^.5,'*r'), hold on
        %         plot(i,sum(sum(derivVcorr.^2))^.5,'+g')
        %         figure(h1)
        %         pause%(.5)
    end

    if nargout > 2
        GU = U(:,1:N)'*Ur(:,1:N);

        rhoU(i+1+Binitstep) = (sum(sum(abs(GU))./max(abs(GU)))-N + sum(sum(abs(GU'))./max(abs(GU')))-N)/2/(N-1)/N;


        GV = V(:,1:N)'*Vr(:,1:N);

        rhoV(i+1+Binitstep) = (sum(sum(abs(GV))./max(abs(GV)))-N + sum(sum(abs(GV'))./max(abs(GV')))-N)/2/(N-1)/N;

        %     rhoU(i+1) = sum(acos(abs(diag(U(:,1:N)'*Ur(:,1:N)))))/N;
        %     rhoV(i+1) = sum(acos(abs(diag(V(:,1:N)'*Vr(:,1:N)))))/N;

        dU(i) = sum(sum(derivU.^2))^.5;
        dV(i) = sum(sum(derivV.^2))^.5;

        dUc(i) = sum(acos(diag(U'*U_old)));
        dVc(i) = sum(acos(diag(V'*V_old)));
    end
    % 	iter = iter+1;
%     sum(sum(sum(CC.^2))) - sum(diag(sum(CC.^2,3)))
%     if sum(sum(sum(abs(CC).^2))) - sum(diag(sum(abs(CC).^2,3))) < 1e-12
%         return
%     end
end

% if pl
%     subplot(223)
%     quiver(zeros(1,P),zeros(1,P),U(1,:),U(2,:),0,'r')
% end
