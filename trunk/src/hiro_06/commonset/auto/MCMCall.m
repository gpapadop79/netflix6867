%% Calculate mean

global q

q = 3;  

meanV = zeros(D,2);
meanV(1:d,1) = Myukernel.*n;
meanV(1:d,2) = n;
for k=1:dataSize
    meanV(data(k,1),1)=meanV(data(k,1),1)+data(k,3);
    meanV(data(k,1),2)=meanV(data(k,1),2)+1;
end
meanV(:,1) = meanV(:,1)./meanV(:,2);

%% MCMC

tic;

for kkk=1:nMCMC

    if(mod(kkk,10)==1)
        kkk
        
        Rmean = Rsum./(kkk-nburnin);
        toc;
        SE = (R-Rmean).^2;
        RMSE = sqrt(mean(mean(SE)))
    end

    
    %iterate over V
    Vvar = getD(X,sigmaML,sigmaPriorV);
    S = chol(Vvar);
    for i = 1:D
        myu = getMyuV(Vvar,Rmcmc,X,i,meanV,sigmaML);
        V(:,i) = S'*randn(3,1) + myu;
    end
    
    %iterate over X
    Xvar = getD(V,sigmaML,sigmaPriorX);
    S = chol(Xvar);
    for i = 1:N
        myu = getMyuX(Xvar,Rmcmc,V,i,meanV,sigmaML);
        X(:,i) = S'*randn(3,1) + myu;
    end
    
    %iterate over R
    for k = 1:sizeR
        i=freeData(k,1);
        j=freeData(k,2);
        Rmcmc(i,j) = sigmaML*randn() + meanV(i,1) + V(:,i)'*X(:,j);
    end
    
    if(kkk>nburnin)
        %Xsum = Xmean + X./(nMCMC-nburnin);
        %Vmean = Vmean + V./(nMCMC-nburnin);
        Rsum = Rsum + Rmcmc;
    end    
end




%% RMSE
% Rest=Vmean'*Xmean;
% for i=1:N
%     Rest(:,i) = Rest(:,i)+meanV(:,1);
% end
% 
% for i=1:numel(Rest)
%     if Rest(i)>5
%         Rest(i)=5;
%     elseif Rest(i)<1
%         Rest(i)=1;
%     end
%     Rest(i) = round(Rest(i));
% end

SE = (R-Rmean).^2;
RMSE = sqrt(mean(mean(SE)))


toc;
%clear Vps;
%clear Xps;
    