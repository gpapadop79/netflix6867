%% Initialize

load result;
sigmaML;
makeData;

d = 10;
n = 30;
q = 3;
temp=size(R);
D = temp(1);
N = temp(2);

nMCMC = 1000;
nburnin = 0;

sigmaPriorX = eye(q);
sigmaPriorV = eye(q);


%initialization
X = sigmaPriorX*randn(q,N);
V = sigmaPriorV*randn(q,D);
Xmean = zeros(q,N);
Vmean = zeros(q,D);
Rmean = zeros(D,N);
Xsum = zeros(q,N);
Vsum = zeros(q,D);
Rsum = zeros(D,N);

dataSize=size(data,1);

clear Rmcmc;
Rmcmc = zeros(D, N);
Rmcmc(1:d,1:n) = Rkernel;
for i=1:dataSize
    Rmcmc(data(i,1),data(i,2)) = data(i,3);
end

clear freeData;
count=1;
for i=1:D
    for j=1:N
        if (Rmcmc(i,j)==0)
            freeData(count,:)=[i j];
            Rmcmc(i,j)=randn()+3;
            count = count+1;
        end
    end
    i
end
sizeR = count-1;

'initialized'

MCMCall;
