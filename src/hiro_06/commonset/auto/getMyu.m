function myu = getMyuV(D, R, X, i)

global q n d

S = zeros(q);

for j=1:n
    S = S+ R(i,j)*X(:,j);
end

myu = D\S;