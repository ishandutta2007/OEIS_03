\\ a(n) = n! * Sum_{k=1..floor(n/3)} k^(n-3*k-1)/(6^k * (n-3*k)!).
a(n) = n! * sum(k=1, n\3, k^(n-3*k-1)/(6^k * (n-3*k)!));
for(n=0, 28, print1(a(n), ", ")); 
