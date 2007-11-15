% function [trainingSet testSet] = devideData(R, sparsity)
% This function devides the data matrix into trainingSet and testSet
% Arguments
%  R - data matrix
%  sparsity - ratio of the number data in trainigSet and in R
% Outputs
%  trainingSet - training set data matrix. Missing elements are filled with
%  zero.
%  testSet - test set data. Missing elements are filled with
%  zero.
%
% Example:
% [train test] = devideData(R, 0.1);

function [trainingSet testSet] = devideData(R, sparsity)

trainingSet = zeros(size(R));
testSet = zeros(size(R));
rand('twister',sum(100*clock));
randM = rand(size(R));
indexTraining = find(randM <= sparsity);
indexTest = find(randM > sparsity);

trainingSet(indexTraining) = R(indexTraining);
testSet(indexTest) = R(indexTest);
