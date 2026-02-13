M=30;

\\ G.f. A(x) satisfies A(x) = 1/(2/(1 + x*A(x))^3 - 1).
seq(n) = my(A=1); for(i=1, n, A=1/(2/(1 + x*A)^3 - 1) + x*O(x^n)); Vec(A);
seq(M) 