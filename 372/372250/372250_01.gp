\\ a(n) = n! * Sum_{k=0..floor(n/3)} (n-k+1)^(n-2*k-1) / (k! * (n-3*k)!).
a(n) = n! * sum(k=0, n\3, (n-k+1)^(n-2*k-1) / (k! * (n-3*k)!));
for(n=0, 18, print1(a(n), ", "))

