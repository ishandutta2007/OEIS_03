M=15;

a(n) = n! * sum(k=0, n\3, (1/2)^k * (k+1)^(n-2*k-1) / (k! * (n-3*k)!));
for(n=0, M, print1(a(n), ", "))
