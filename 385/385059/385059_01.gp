M=16;

a(n, k=-1) = if(n*k==0, 0^n, (-1)^n*k*sum(j=1, n, (-n+k)^(j-1)*binomial(n, j)*a(n-j, 2*j)));
for(n=0, M, print1(a(n),", "));

