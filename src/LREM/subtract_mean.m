function [R_out average]=subtract_mean(R_in)

[n m]=size(R_in);

average = zeros(n,1);
R_out = zeros(n,m);

for i=1:n
    r = R_in(i,:);
    average(i) = mean(r(r~=0));
    for j=1:m
        if R_in(i,j) == 0
            R_out(i,j) = NaN;
        else
            R_out(i,j) = R_in(i,j) - average(i);
        end
    end
end
