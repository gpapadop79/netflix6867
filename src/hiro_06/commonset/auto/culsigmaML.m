d = 10;
n = 30;
q = 3;

temp = size(A);
D = temp(1)-1;
N = temp(2);
%d<n

Rkernel=A(2:d+1,1:n);
XkernelLabel=A(1,1:n);
R = A(2:D+1,1:N);

Myukernel(:,1)=mean(Rkernel,2);

S=zeros(d);
for i=1:d
    S=S+(Rkernel(:,i)-Myukernel)*(Rkernel(:,i)-Myukernel)';
end
S = S/d;

[U, Sigmas]=eig(S);
sigma=diag(Sigmas);

sigmaML = 0;
for i=1:d-q
    sigmaML = sigmaML + sigma(i);
end
sigmaML=sqrt(sigmaML/(d-q))
