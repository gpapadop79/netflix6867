function R_out = real2rate(R)

R_out=R;

R_out(R>5) = 5*ones(size(R(R>5)));
R_out(R<1) = ones(size(R(R<1)));