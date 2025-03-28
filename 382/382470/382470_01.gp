M=25;

\\ a(n) = Sum_{k=0..n} binomial(k+3,3) * binomial(2*k,2*n-2*k).
b(n) = sum(k=0, n, binomial(k+3,3)*binomial(2*k,2*n-2*k));
for(n=0, M, print1(b(n),", "))     

\\ signature (8,-20,16,-26,88,-48,24,-163,24,-48,88,-26,16,-20,8,-1).
a(n) = if(n<16, b(n), 8*a(n-1)-20*a(n-2)+16*a(n-3)-26*a(n-4)+88*a(n-5)-48*a(n-6)+24*a(n-7)-163*a(n-8)+24*a(n-9)-48*a(n-10)+88*a(n-11)-26*a(n-12)+16*a(n-13)-20*a(n-14)+8*a(n-15)-a(n-16));
for(n=0, 30, print1(a(n)-b(n),", "))

