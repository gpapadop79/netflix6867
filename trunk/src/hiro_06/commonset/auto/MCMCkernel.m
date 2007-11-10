global q n d

d = 10;
n = 30;
q = 3;

nMCMC = 110;
nburnin = 10;

sigmaPriorX = eye(q);
sigmaPriorV = eye(q);

%initialization
Xmem = zeros(nMCMC-nburnin,q,n);
Vmem = zeros(nMCMC-nburnin,q,d);
X = sigmaPriorX*randn(q,n);
V = sigmaPriorV*randn(q,d);
Xmean = zeros(q,n);
Vmean = zeros(q,d);

tic;

for k=1:nMCMC
    
    %iterate over V
    D = getD(X,sigmaML,sigmaPriorV);
    S = chol(D);
    for i = 1:d
        myu = getMyuV(D,Rkernel,X,i,Myukernel,sigmaML);
        V(:,i) = S'*randn(3,1) + myu;
    end
    
    %iterate over X
    D = getD(V,sigmaML,sigmaPriorX);
    S = chol(D);
    for i = 1:n
        myu = getMyuX(D,Rkernel,V,i,Myukernel,sigmaML);
        X(:,i) = S'*randn(3,1) + myu;
    end
    
    if(k>nburnin)
        Xmem(k-nburnin,:,:) = X;
        Vmem(k-nburnin,:,:) = V;
        Xmean = Xmean + X./(nMCMC-nburnin);
        Vmean = Vmean + V./(nMCMC-nburnin);
    end    
end

toc;

arg=1;
reconstruct;