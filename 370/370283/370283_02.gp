\\ Coefficient of x^n in the expansion of 1/( (1-x)^3 - x^2 )^n.
a(n) = polcoef( 1/( (1-x)^3 - x^2 + x*O(x^n) )^n, n);
for(n=0, 22, print1(a(n), ", "));

\\ See A369694.
my(N=30, x='x+O('x^N)); Vec(exp(sum(k=1, N, a(k)*x^k/k))) 
