function D = getD(X, sigmaLik, sigmaPrior)

global q

sizeX = size(X);

Dinv = zeros(q);
for i=1:sizeX(2)
    Dinv = Dinv + X(:,i)*X(:,i)';
end

Dinv=Dinv/sigmaLik^2 + sigmaPrior;
%Dinv=Dinv/sigmaLik^2;

D = inv(Dinv);