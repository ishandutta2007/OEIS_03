M=20;

\\ E.g.f. A(x) satisfies A(x) = exp(x) * B(x*A(x)^2), where B(x) = 1 + x*B(x)^3 is the g.f. of A001764.
a(n, s, t) = {
    my(A=1, B=sum(k=0, n, binomial(s*k,k)/((s-1)*k+1) * x^k) + x*O(x^(n+1))); 
    for(i=1, n, A=exp(x + x*O(x^(n+1))) * subst(B, x, x*A^t) ); 
    Vec(serlaplace(A)); 
}
a(M, 3, 2)

\\ Let F(x) be the e.g.f. of A382000. F(x) = B(x*A(x)^2) = exp( 1/3 * Sum_{k>=1} binomial(3*k,k) * (x*A(x)^2)^k/k ).
b(n, s, t) = {
    my(A=1, B=sum(k=0, n, binomial(s*k,k)/((s-1)*k+1) * x^k) + x*O(x^(n+1))); 
    for(i=1, n, A=exp(x + x*O(x^(n+1))) * subst(B, x, x*A^t) ); 
    A; 
}

g = b(M, 3, 2);
f1 = g * exp(-x + x * O(x^(M+1)));
f2 = exp( 1/3 * sum(k=1, M, binomial(3*k,k) * (x*g^2)^k/k) );
Vec(serlaplace(f1))
Vec(serlaplace(f2))
serlaplace(f1 - f2)
