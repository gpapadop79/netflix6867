function [mean var] = mymean(x,W)

%weighted mean
S = sum(W);
mean = zeros(size(x));
WW = zeros(size(x));

for i=1:size(x,2)
    WW(:,i) = W;
    mean(:,i) = x(:,i)'*W/S;
end

xc = x - mean;  % Remove mean
var = (xc' * (xc.*WW)) / S;


% function mean = mymean(x,W)
% 
% %weighted mean
% S = sum(W);
% mean = zeros(size(x));
% 
% for i=1:size(x,2)
%     mean(:,i) = x(:,i)'*W/S;
% end


