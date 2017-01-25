% This script runs illustrates some simple scaling and rank issues. 
% Stefan Wild. Argonne National Laboratory, January 2016.
clear all


% Try n = 3, 30, 150
n = 3;
vec  = 2*10.^(2*(1-(1:n)));
A = diag(vec); % This is the Hessian for probnum=3
matlab_computed_rank = rank(A)
actual_rank = sum(diag(A)>0)

min_diag = min(diag(A));
recomputed_matlab_rank = rank(A,min_diag/10)

if n==150
    B = spdiags(vec', 0, n,n); % Make B  a sparse diagonal matrix
    whos % Look at difference between A and B
    matlab_sparse_rank = sprank(A)
end


if 1==0
    probnum=3;
    func=@examplefun; 
    
    n=2;
    B = spdiags(vec(1:n)', 0, n,n); % Make B  a sparse diagonal matrix
    x=10*randn(n,1);
    g = 2*( x'.*10.^(2*(1-(1:n))) )';
    y=x-(B\g)
    f = func(y);
end