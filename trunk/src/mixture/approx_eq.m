function out = approx_eq(A, B, tol)

ineq1 = A < B+tol;
ineq2 = B-tol < A;

out = ineq1 & ineq2;