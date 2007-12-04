clear;
load ../dataset/Netflix_subset.mat

sparsity = 0.5;

[R_train R_test]=devideData(R, sparsity);

[prum pui pmj like pumrij] = em(R_train, 10, 4, 100);

%%
[err errmean Rhat Rhatmean] = rmse(R_train, prum, pui, pmj);

