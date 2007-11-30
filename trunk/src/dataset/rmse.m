function rmse_out = rmse(R_est, R_test)

% Calculate RMSE (root mean square error)
%
% Input
%  R_est - estimated rating matrix
%  R_test - test data set (as rating matrix)
%
% The function only uses the entries (i,j) where R_test(i,J) is not zero
% nor NaN.
%
% Hiro Ono, 11/29/2007

[n m]=size(R_est);
SE = 0;
count = 0;
for i=1:n
    for j=1:m
        if (R_test(i,j) ~= 0 && ~isnan(R_test(i,j)))
            SE = SE + (R_est(i,j)-R_test(i,j)).^2;
            count=count+1;
        end
    end
end
rmse_out = sqrt(SE/count);
