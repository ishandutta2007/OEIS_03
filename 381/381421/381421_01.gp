M=28;

\\ a(n) = Sum_{k=0..n} (k+1) * binomial(2*k,2*n-2*k).
b(n) = sum(k=0, n, (k+1)*binomial(2*k, 2*n-2*k));
for(n=0, M, print1(b(n),", "))     

\\ signature (4,-2,0,-11,0,-2,4,-1).
a(n) = if(n<8, b(n), 4*a(n-1)-2*a(n-2)-11*a(n-4)-2*a(n-6)+4*a(n-7)-a(n-8));
for(n=0, 25, print1(a(n)-b(n),", "))

