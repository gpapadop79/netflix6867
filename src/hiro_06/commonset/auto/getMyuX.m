function myu = getMyuX(D, R, V, j,Myukernel,sigmaLik)

global q d

S = zeros(q,1);

for i=1:d
    S = S+ (R(i,j)-Myukernel(i)).*V(:,i);
end

myu = D*S./(sigmaLik^2);