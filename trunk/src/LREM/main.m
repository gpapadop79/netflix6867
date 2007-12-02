function [trainRMSE testRMSE] = main(R, sparsity, dim, plotresult)
% Netflix prize rating estimation using Linear Regression + EM-like
% algorithm.
%
% Written by Hiro Ono, 11/29/2007


N = 100;

[nm nu] = size(R);

[R_train_original R_test]=devideData(R, sparsity);

[R_train average]=subtract_mean(R_train_original);
averageR = repmat(average, 1, nu);

%randamly initialize 
rand('twister',sum(100*clock));
U = randn(dim, nu);

rmse_average = rmse(averageR,R_test);

for i = 1:N
    [M condm(i)] = lrem_update(R_train,U);
    [U condu(i)] = lrem_update(R_train',M');
    U = U';
    
    %rmse_hist_r2r(i) = rmse(real2rate(M*U+averageR),R);
    trainRMSE(i) = rmse(M*U+averageR,R_train_original);
    rmse_test_hist(i) = rmse(M*U+averageR,R_test);
    testRMSE(i) = rmse(real2rate(M*U+averageR),R_test);
    max_hist(i) = max(max(M*U));
    min_hist(i) = min(min(M*U));

end

if plotresult==1
    figure;
    hold on;
    plot(1:N, trainRMSE,'b',1:N,testRMSE,'r-.',1:N,repmat(rmse_average,1,N),':');
    legend('Training RMSE', 'Test RMSE', 'Average')

    figure;
    plot(condm);
    title('Condition number M');

    figure;
    plot(condu);
    title('Condition number U');
end