M=20;

\\ a(n) = 2 * (n+1)! * Sum_{k=0..n} (-1)^k * Stirling2(n,k)/(n-k+2)!.
a(n) = 2 * (n+1)! * sum(k=0, n, (-1)^k * stirling(n, k, 2) / (n-k+2)!);
for(n=0, M, print1(a(n), ", "));
