Rest=Vmean'*Xmean;

plot(1:d,Rkernel(1:d,arg),1:d,(Rest(1:d,arg)+Myukernel));

for i=1:n
    Rest(1:d,i) = Rest(1:d,i)+Myukernel;
end

SE = (Rkernel-Rest(1:d,1:n)).^2;
RMSE = sqrt(mean(mean(SE)))