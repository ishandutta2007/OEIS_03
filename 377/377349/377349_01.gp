M=19;

\\ a(n) = Sum_{k=0..floor((2*n+1)/3)} (2*n-2*k)!/(2*n-3*k+1)! * |Stirling1(n,k)|.
a(n) = sum(k=0, (2*n+1)\3, (2*n-2*k)!/(2*n-3*k+1)! * abs(stirling(n,k,1)));
for(n=0, M, print1(a(n), ", ")) 