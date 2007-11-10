function evaluateRMSE(R)

[U S V] = svd(R);
n=min(size(S));

for i=1:n
    Ui = U(:,1:i);
    Si = S(1:i,1:i);
    Vi = V(:,1:i);
    Ri = Ui*Si*Vi';
    
    SE = (R-Ri).^2;
    RMSE(i) = sqrt(mean(mean(SE)));
    
    disp(sprintf('%dth order approximation: RMSE = %f\n', i, RMSE(i)));
end

plot(1:51,RMSE,'.-');
    