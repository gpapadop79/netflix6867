function myu = getMyuV(D, R, X, i,Myukernel,sigmaLik)

global q n

S = zeros(q,1);

for j=1:n
    S = S+ (R(i,j)-Myukernel(i)).*X(:,j);
end

myu = D*S/(sigmaLik^2);