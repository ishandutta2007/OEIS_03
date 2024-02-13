\\ Coefficient of x^n in the expansion of ( 1/(1-x) * (1+x^2)^3 )^n.
a(n) = polcoef( ( 1/(1-x + x*O(x^n)) * (1+x^2)^3 )^n, n);
for(n=0, 23, print1(a(n), ", "));

\\ See A369262.
my(N=30, x='x+O('x^N)); Vec(exp(sum(k=1, N, a(k)*x^k/k))) 