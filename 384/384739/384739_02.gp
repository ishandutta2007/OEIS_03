\\ E.g.f. A(x) satisfies A(x) = exp( x * A(x*A(x)^2) ).
a(n, k) = my(A=1); for(i=1, n, A = exp( x * subst(A, x, x*A^2) )); n! * polcoef(A^k, n);
for(n=0, 15, print1(a(n, 1),", "));

