%% Initialize

temp=size(R);
D = temp(1);
N = temp(2);
q = 3;

sigmaPriorV = 0.5*eye(3);
sigmaPriorX = 0.5*eye(3);

nP=100;

alpha = -0.5/sigmaML^2;

Vps=zeros(nP,q+1,D);
Xps=zeros(nP,q+1,N);

Vps(1:nP,1:q,1:d) = Vmem(1:nP,1:q,1:d);
Xps(1:nP,1:q,1:n) = Xmem(1:nP,1:q,1:n);


%initialize
for k=1:nP
    Vps(k,1:q,d+1:D) = sigmaPriorV*randn(q,D-d);
    Xps(k,1:q,n+1:N) = sigmaPriorX*randn(q,N-n);
end

%initialize weights.
Vps(:,q+1,:)=1;
Xps(:,q+1,:)=1;

'initialized'

%% Calculate mean

temp = size(data);
dataSize=temp(1);

%calculate mean
meanV = zeros(D,2);
meanV(1:d,1) = Myukernel.*n;
meanV(1:d,2) = n;
for k=1:dataSize
    meanV(data(k,1),1)=meanV(data(k,1),1)+data(k,q);
    meanV(data(k,1),2)=meanV(data(k,1),2)+1;
end
meanV(:,1) = meanV(:,1)./meanV(:,2);

%% Filter

VW = zeros(1,nP);
XW = zeros(1,nP);
VWacc = zeros(1,nP);
XWacc = zeros(1,nP);

VpsOld = zeros(nP,q);
XpsOld = zeros(nP,q);

rej_b = 0.5;
rej_a = sqrt(1-rej_b^2);
rej_1a = 1-rej_a;

tic;

for k=1:dataSize
     logW=(data(k,3)-meanV(data(k,1),1)-Vps(1:nP,1:q,data(k,1))*Xps(1:nP,1:q,data(k,2))').^2;
     %Vps(1:nP,q+1,data(k,1)) = Vps(1:nP,q+1,data(k,1))+sum(logW,2);
     %Xps(1:nP,q+1,data(k,2)) = Xps(1:nP,q+1,data(k,2))+(sum(logW,1))';
     
     logVW = alpha*sum(logW,2);
     maxlogVW = max(logVW);
     VW = exp(logVW-maxlogVW);
     
     logXW = alpha*sum(logW,1);
     maxlogXW = max(logXW);
     XW = exp(logXW-maxlogXW);
          
     [Vchilda Vcov] = mymean(Vps(:,1:q,data(k,1)),VW);
     [Xchilda Xcov] = mymean(Xps(:,1:q,data(k,2)),XW');
%      Vsigma = rej_b*chol(Vcov);
%      Xsigma = rej_b*chol(Xcov);
      Vsigma = rej_b*diag(sqrt(var(Vps(:,1:q,data(k,1)),VW)));
      Xsigma = rej_b*diag(sqrt(var(Xps(:,1:q,data(k,2)),XW')));
         
     %resample
     VpsOld = Vps(:,1:q,data(k,1));
     XpsOld = Xps(:,1:q,data(k,2));      
     [VW VWi]= sort(VW,1,'descend');
     [XW XWi]= sort(XW,1,'descend');
     VWsum = sum(VW);
     XWsum = sum(XW);
     
     VW=VW./VWsum.*nP;
     XW=XW./XWsum.*nP;
     
     resV=VW(1);
     resX=XW(1);
     jV=1; 
     jX=1;
     for i=1:nP
         if(resV>0.5)
             resV = resV-1;
         else
             while(resV<=0.5)
                jV=jV+1;
                resV = resV+VW(jV);
             end
         end
         Vps(i,1:q,data(k,1)) = VpsOld(jV,:);
        
         if(resX>0.5)
             resX = resX-1;
         else
            while(resX<=0.5)
                jX=jX+1;
                resX = resX+XW(jX);
            end
         end
         Xps(i,1:q,data(k,2)) = XpsOld(jV,:);
     end     
    
     Vps(:,1:q,data(k,1));
     
%rejuvenation using 1PFS
    
    Vps(:,1:q,data(k,1)) = (rej_a.*Vps(:,1:q,data(k,1)) + rej_1a.*Vchilda) + randn(nP,q)*Vsigma';
    Xps(:,1:q,data(k,2)) = (rej_a.*Xps(:,1:q,data(k,2)) + rej_1a.*Xchilda) + randn(nP,q)*Xsigma';
    
%     Vps(:,1:q,data(k,1))
%     Xps(:,1:q,data(k,2))
         
%      VW = exp(alpha*sum(logW,2));
%      XW = exp(alpha*sum(logW,1));
%      VWacc(1) = VW(1);
%      XWacc(1) = XW(1);
%      for i=2:nP
%          VWacc(i)=VWacc(i-1)+VW(i);
%          XWacc(i)=XWacc(i-1)+XW(i);
%      end
%      
%     
%      for i=1:nP
%          Vps(i,1:q,data(k,1)) = VpsOld(resampleIndex(VWacc, rand()*VWacc(nP)),:);
%          Xps(i,1:q,data(k,2)) = XpsOld(resampleIndex(XWacc, rand()*XWacc(nP)),:);
%      end
%          
 
     %caluculate ESS
%      maxVW = min(Vps(1:nP,q+1,data(k,1)));
%      W_V = exp(-Vps(1:nP,q+1,data(k,1))+maxVW);
%      wSum=sum(W_V);
%      wSqSum=sum(W_V.^2);
%      ESS_V = wSum^2/wSqSum;
%      
%      
%      maxVX = min(Xps(1:nP,q+1,data(k,1)));
%      W_X = exp(-Xps(1:nP,q+1,data(k,2))+maxVX);
%      wSum=sum(W_X);
%      wSqSum=sum(W_X.^2);
%      ESS_X = wSum^2/wSqSum;
     
%      if k>1720
%          k
%          ESS_V
%          ESS_X
%      end
     
%     for i=1:nP
%         for j=1:nP
%             %w = exp(alpha*(data(k,3)-meanV(data(k,1),1)-Xps(j,1:q,data(k,2))*Vps(i,1:q,data(k,1))')^2);
%             %Vps(i,q+1,data(k,1)) = w.*Vps(i,q+1,data(k,1));
%             %Xps(j,q+1,data(k,2)) = w.*Xps(j,q+1,data(k,2));
%             w = alpha*(data(k,3)-meanV(data(k,1),1)-Xps(j,1:q,data(k,2))*Vps(i,1:q,data(k,1))')^2;
%             Vps(i,q+1,data(k,1)) = w + Vps(i,q+1,data(k,1));
%             Xps(j,q+1,data(k,2)) = w + Xps(j,q+1,data(k,2));
%         end
%     end
%         if(mod(k,100)==0)
%             k
%         end
end

toc;

% for k=1:dataSize
%      logW=(data(k,3)-meanV(data(k,1),1)-Vps(1:nP,1:q,data(k,1))*Xps(1:nP,1:q,data(k,2))').^2;
%      Vps(1:nP,q+1,data(k,1)) = Vps(1:nP,q+1,data(k,1))+sum(logW,2);
%      Xps(1:nP,q+1,data(k,2)) = Xps(1:nP,q+1,data(k,2))+(sum(logW,1))';
% end

% Vps(:,q+1,:) = Vps(:,q+1,:).*alpha;
% Xps(:,q+1,:) = Xps(:,q+1,:).*alpha;

%% mean

Vmean=zeros(q,D);
Xmean=zeros(q,N);
% Vw = zeros(1,D);
% Xw = zeros(1,N);
% 
% maxWV=zeros(1,D);
% maxWX=zeros(1,N);
% 
% maxWV(1,1:D) = max(Vps(:,4,:),[],1);
% maxWX(1,1:N) = max(Xps(:,4,:),[],1);
% 
% 
% for k=1:nP
%     for i=1:D
%         Vmean(:,i) = Vmean(:,i)+Vps(k,4,i).*Vps(k,1:q,i)';
%         Vw(i) = Vw(i) + Vps(k,4,i);
%     end
%     for i=1:N
%         Xmean(:,i) = Xmean(:,i)+Xps(k,4,i).*Xps(k,1:q,i)';
%         Xw(i) = Xw(i) + Xps(k,4,i);
%     end
% end
% 
% for k=1:q
%     Vmean(k,:) = Vmean(k,:)./Vw;
%     Xmean(k,:) = Xmean(k,:)./Xw;
% end

for i=1:D
    Vmean(:,i) = mean(Vps(:,1:q,i),1)';
end
for i=1:N
    Xmean(:,i) = mean(Xps(:,1:q,i),1)';
end


% Vmax=zeros(q,D);
% Xmax=zeros(q,N);
% 
% for i=1:D
%     [temp, index] = max(Vps(:,4,i));
%     Vmax(:,i) = Vps(index,1:q,i)';
% end
% for i=1:N
%     [temp, index] = max(Xps(:,4,i));
%     Xmax(:,i) = Xps(index,1:q,i)';
% end

%% RMSE
Rest=Vmean'*Xmean;
%Rest=Vmax'*Xmax;
for i=1:N
    Rest(:,i) = Rest(:,i)+meanV(:,1);
end

for i=1:numel(Rest)
    if Rest(i)>5
        Rest(i)=5;
    elseif Rest(i)<1
        Rest(i)=1;
    end
    %Rest(i) = round(Rest(i));
end

SE = (R-Rest).^2;
RMSE = sqrt(mean(mean(SE)))

col=86;
plot(1:D,R(:,col),1:D,Rest(:,col));

%clear Vps;
%clear Xps;
    