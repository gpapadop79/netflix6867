%% Initialize

temp=size(R);
D = temp(1);
N = temp(2);
q = 3;

sigmaPriorV = 0.5*eye(3);
sigmaPriorX = 0.5*eye(3);

nP=10;

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
     Vps(1:nP,q+1,data(k,1)) = Vps(1:nP,q+1,data(k,1))+sum(logW,2);
     Xps(1:nP,q+1,data(k,2)) = Xps(1:nP,q+1,data(k,2))+(sum(logW,1))';

%      if(mod(k,100)==0)
%           k
%      end
end

toc;

 Vps(:,q+1,:) = Vps(:,q+1,:).*alpha;
 Xps(:,q+1,:) = Xps(:,q+1,:).*alpha;

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


Vmax=zeros(q,D);
Xmax=zeros(q,N);

for i=1:D
    [temp, index] = max(Vps(:,4,i));
    Vmax(:,i) = Vps(index,1:q,i)';
end
for i=1:N
    [temp, index] = max(Xps(:,4,i));
    Xmax(:,i) = Xps(index,1:q,i)';
end

%% RMSE
%Rest=Vmean'*Xmean;
Rest=Vmax'*Xmax;
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
    