function [M condnum] = lrem_update(R, U)

[nm nu] = size(R);
dim = size(U,1);
M = zeros(nm, dim);

condn = zeros(nm,1);

for i = 1:nm
    r = R(i,:);
    X = U(:,~isnan(r));
    y = r(~isnan(r));
    if nargout > 1
       % condn(i) = cond(X*X');
    end
    M(i,:) = y/X;
end

condnum = max(condn);