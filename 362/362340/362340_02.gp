M=21;

a(n) = n! * sum(k=0, n\2, (-k/2)^k / (k! * (n-2*k)!));
for(n=0, M, print1(a(n), ", ")); 
